cmake_minimum_required(VERSION 3.8)

add_library(
    ${PROJECT_NAME} 
    STATIC
    src/clock_relations.cpp
)

target_include_directories(
    ${PROJECT_NAME}
    PUBLIC
        $<INSTALL_INTERFACE:include>    
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
)

add_executable(
    ${PROJECT_NAME}_example
    src/example.cpp
)

target_link_libraries(${PROJECT_NAME}_example PUBLIC ${PROJECT_NAME})
