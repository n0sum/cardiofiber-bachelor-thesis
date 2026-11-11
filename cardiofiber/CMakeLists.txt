cmake_minimum_required(VERSION 3.20)
project(cardiofiber LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)

find_package(glm REQUIRED)

add_library(cardiofiber STATIC src/Mesh.cpp)
target_compile_features(cardiofiber PUBLIC cxx_std_20)
target_include_directories(cardiofiber PUBLIC include)
target_link_libraries(cardiofiber PUBLIC glm::glm)
