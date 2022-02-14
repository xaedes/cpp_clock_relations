#include "clock_relations/clock_relations.h"

#include <array>
#include <vector>
#include <iostream>
#include <cstdint>
#include <unordered_map>
#include <functional>

using ClockId = uint32_t;

struct ClockTime
{
    ClockId clock_id = 0;
    double time = 0;
};

struct LinearRelation
{
    // y = f(x)
    // y = x * scale + offset
    // x = (y - offset) / scale

    double offset = 0;
    double scale = 1;

    LinearRelation inverse() const
    {
        // x = (y - offset) / scale
        // x = y/scale - offset/scale
        LinearRelation inv;
        inv.scale = 1/scale;
        inv.offset = -offset/scale;
        return inv;
    }

    // as_time is like homogeneous vector entry 1.0.
    // setting it to false, treats x as duration, i.e.
    // as one dimensional direction vector instead of 
    // one dimensional position vector. only scale is 
    // applied, similar to how direction vectors are not
    // translated by transformation matrices in linear algebra,
    // but only scaled and rotated.
    double applyForward(double x, bool as_time = true) const
    {
        return x * scale + (as_time ? offset : 0);
    }
    double applyBackward(double y, bool as_time = true) const
    {
        return (as_time ? y - offset : y) / scale;
    }


};



struct MockClock
{
    ClockId id = 0;
    double current = 0;

    LinearRelation world_relation;

    ClockTime now() const {
        return ClockTime{id, nowRaw()};
    }
    double nowRaw() const {
        return current;
    }
    void advance(double world_duration) {
        const bool as_time = false;
        current += world_relation.applyForward(world_duration, as_time) ;
    }
    void reset(double world_time = 0) {
        const bool as_time = true;
        current = world_relation.applyForward(world_time, as_time) ;
    }
};

struct TimeCorrespondence
{
    double a;
    double b;
};

struct IdRelation
{
    ClockId reference = 0;
    ClockId id = 0;

    bool operator==(const IdRelation& rhs) const
    {
        return (
            (reference == rhs.reference)
         && (id == rhs.id)
        );
    }

    struct Hash
    {
        // see https://en.cppreference.com/w/cpp/utility/hash
        std::size_t operator()(IdRelation const& x) const noexcept
        {
            std::size_t h1 = std::hash<ClockId>{}(x.reference);
            std::size_t h2 = std::hash<ClockId>{}(x.id);
            return h1 ^ (h2 << 1); // or use boost::hash_combine
        }
    };

};

struct ClockTimeCorrespondence
{
    ClockTime a;
    ClockTime b;

    operator TimeCorrespondence() const {
        return TimeCorrespondence{
            a.time,
            b.time
        };
    }

    operator IdRelation() const {
        return IdRelation{
            a.clock_id,
            b.clock_id
        };
    }
};


struct ClockRelation
{
    IdRelation ids;
    LinearRelation linear{};
};

struct EstimateClockRelation
{
    void setup(IdRelation ids, LinearRelation guess = LinearRelation{})
    {
        m_ids = ids;
        m_estimate = guess;
    }
    IdRelation getIds() const {
        return m_ids;
    }

    LinearRelation getLinearRelation() const {
        return m_estimate;
    }

    ClockRelation getClockRelation() const {
        return ClockRelation{m_ids, m_estimate};
    }

    void reportCorrespondence(const TimeCorrespondence x, bool update=true)
    {
        bool didObserve = false;
        if (correspondences.size() > 0)
        {
            Observation observation{{
                correspondences.back(),
                x
            }};
            observations.push_back(observation);
            didObserve = true;
        }
        correspondences.push_back(x);

        if (didObserve && update)
        {
            updateEstimate();
        }
    }

    void updateEstimate()
    {
        // static thread_local std::vector<LinearRelation> observedRelations;
        // observedRelations.clear();
        LinearRelation sum{0,0};
        double infoSum = 0;
        // estimate relation to convert from clock_a to clock_b
        // estimate offset and scale for is_time in [0,1] and:
        // time_b = scale * time_a + is_time * offset
        for (int i=0; i<observations.size(); i++)
        {

            auto& observation = observations[i];
            // time from clock a and b at (observation local) steps 0 and 1
            auto a0 = observation.items[0].a;
            auto b0 = observation.items[0].b;
            auto a1 = observation.items[1].a;
            auto b1 = observation.items[1].b;

            auto da = a1 - a0;
            auto db = b1 - b0;
            if (abs(da) > 0)
            {
                auto scale = db / da;
                
                // b = scale * a + offset
                // offset = b - scale * a
                auto offset0 = b0 - scale * a0;
                auto offset1 = b1 - scale * a1;

                auto offset = (offset0 + offset1) * 0.5;

                // observedRelations.push_back(LinearRelation{offset, scale});
                const double weight = 1.0;
                sum.offset += weight * offset;
                sum.scale += weight * scale;
                infoSum += weight;
            }
            //observation.b.a
            //observation.b.b
        }
        if (infoSum > 0)
        {
            sum.offset /= infoSum;
            sum.scale /= infoSum;
            m_estimate = sum;
        }
    }

    std::vector<TimeCorrespondence> correspondences;

protected:
    struct Observation
    {
        std::array<TimeCorrespondence, 2> items;
    };
    std::vector<Observation> observations;

    IdRelation m_ids;
    LinearRelation m_estimate;

};


struct ClockRelationGraph
{
    std::vector<ClockId> nodes;
    std::unordered_map<ClockId, ClockId> edges;

    void insert(IdRelation edge, EstimateClockRelation& estimator)
    {

    }
};

struct EstimateClockRelations
{
    void reportCorrespondence(const ClockTimeCorrespondence x)
    {
        IdRelation key = static_cast<IdRelation>(x);
        if (estimates.count(key) == 0)
        {
            estimates[key].setup(key);
            graph.insert(key, estimates[key]);
        }
        estimates[key].reportCorrespondence(static_cast<TimeCorrespondence>(x));
    }
    ClockRelationGraph graph;
    std::unordered_map<IdRelation, EstimateClockRelation, IdRelation::Hash> estimates;
};

int main(int argc, char **argv)
{
    // track relation between different clocks
    // clock_relations::ClockRelations relations;
    EstimateClockRelations relations;
    
    // this example assumes a local system clock and multiple
    // other clocks. ticks from these other clocks are regularly 
    // reported to the local system. in real world this could be 
    // timestamps from external sensors each having its own clock.
    // the relations between external clocks and system clock are 
    // estimated using time correspondences of time points expressed
    // in different clocks.
    // in our example we use the correspondence between system clock
    // and external clock. this applies well to the real world use case
    // of external sensors reporting their own timestamps.

    // the estimated relations between those clocks can be used to 
    // transform timepoints from one clock into another.
    
    // each clock is identified by an unique identifier.

    // a clock that represents a system clock, which will be used as reference 
    // in TimeCorrespondence
    ClockId next_clock_id = 0;
    MockClock systemClock{next_clock_id++, 0.0, LinearRelation{0.0, 1.0}};

    // external clocks
    std::vector<MockClock> clocks;
    clocks.push_back(MockClock{next_clock_id++, 0.0, LinearRelation{20.0, 1 + 1e-6}});
    clocks.push_back(MockClock{next_clock_id++, 0.0, LinearRelation{0.0, 0.001}});
    clocks.push_back(MockClock{next_clock_id++, 0.0, LinearRelation{10.0, 0.001 - 1e-6}}); 
    clocks.push_back(MockClock{next_clock_id++, 0.0, LinearRelation{20.0, 0.001 + 1e-6}});

    std::vector<double> clockIntervals{
        1.0 / 30.0,
        1.0 / 60.0,
        1.0 / 10.0,
        1.0 / 1.0
    };

    // functions to reset and advance all mock clocks in sync
    // auto tick = [&systemClock, &clocks](){
    //     systemClock.tick();
    //     for (int i=0; i<clocks.size(); i++)
    //     {
    //         clocks[i].tick();
    //     }
    // };
    // auto reset = [&systemClock, &clocks](){
    //     systemClock.reset();
    //     for (int i=0; i<clocks.size(); i++)
    //     {
    //         clocks[i].reset();
    //     }
    // };

    // collect traces of time correspondences for each external clock
    double maxWorldTime = 5.0;
    std::vector<std::vector<ClockTimeCorrespondence>> clockTraces;
    clockTraces.resize(clocks.size());
    for (int i=0; i<clocks.size(); i++)
    {
        auto& trace = clockTraces[i];
        auto& clock = clocks[i];
        double interval = clockIntervals[i];

        systemClock.reset();
        clock.reset();
        int count = static_cast<int>(floor(maxWorldTime / interval));
        std::cout << "clock #" << i << "\n";
        for (int k=0; k<count; k++)
        {
            ClockTimeCorrespondence correspondence{
                systemClock.now(),
                clock.now()
            };
            trace.push_back(correspondence);
            systemClock.advance(interval);
            clock.advance(interval);
            std::cout << correspondence.a.time << " " << correspondence.b.time << "\n";
        }
    }

    for (int i = 0; i < clocks.size(); i++)
    {
        auto& trace = clockTraces[i];
        for (int k = 0; k < trace.size(); k++)
        {
            relations.reportCorrespondence(
                trace[k]
            );
        }
    }

    for (int i = 0; i < clocks.size(); i++)
    {
        auto& clock = clocks[i];
        IdRelation key{systemClock.id, clock.id};
        auto& estimator = relations.estimates[key];
        auto& relation = estimator.getClockRelation().linear;
        std::cout << "clock #" << i  << "\n";
        std::cout << "true      (offset, scale):" 
            << clock.world_relation.offset << " " 
            << clock.world_relation.scale << "\n";
        std::cout << "estimated (offset, scale):" 
            << relation.offset << " " 
            << relation.scale << "\n";
        std::cout << "\n";
    }
    
    /*    
    for (int k=0; k<count; k++)
    {
        std::cout << "#" << k;
        for (int i=0; i<times.size(); i++)
        {
            if (k == 0)
            {
                // output system clock time once per line, 
                // since it is the same for all times[...][k]
                std::cout << " " << times[i][k].a.time;
            }
            std::cout << " " << times[i][k].b.time;
        }
        std::cout << "\n";
    }
    */
    return 0;
}
