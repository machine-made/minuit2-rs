cmake_minimum_required(VERSION 3.16)
project(minuit2_ref_runner LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

find_package(Minuit2 CONFIG REQUIRED)

add_executable(ref_runner main.cpp)
target_link_libraries(ref_runner PRIVATE Minuit2::Minuit2)
