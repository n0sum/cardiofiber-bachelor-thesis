cmake_minimum_required(VERSION 3.20)
project(Estimate)

set(CMAKE_CXX_STANDARD 20)

find_package(glm REQUIRED)

# estimation package
add_subdirectory(../cardiofiber cardiofiber)

add_executable(Estimate
        src/main.cpp
        src/ExperimentConfig.cpp
        src/io.cpp
        src/GroundTruth.cpp
)
target_include_directories(Estimate PRIVATE src include)
target_link_libraries(Estimate PRIVATE cardiofiber)

add_executable(ErrorMaps
        src/ErrorMaps.cpp
        src/io.cpp
)
target_include_directories(ErrorMaps PRIVATE src include)
target_link_libraries(ErrorMaps PRIVATE cardiofiber)

add_executable(Simulate
        src/Simulate.cpp
)
target_include_directories(Simulate PRIVATE src include)
target_link_libraries(Simulate PRIVATE cardiofiber)

add_executable(FiberView
        src/FiberView.cpp
        src/io.cpp
)
target_include_directories(FiberView PRIVATE src include)
target_link_libraries(FiberView PRIVATE cardiofiber)
