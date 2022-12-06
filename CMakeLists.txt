cmake_minimum_required(VERSION 3.22)
project(BytePairRam)

set(CMAKE_CXX_STANDARD 11)

find_package(ZLIB 1.2.8 REQUIRED)

add_executable(BytePairRam src/main.cpp include/parser.h src/parser.cpp include/nucleic_acid.h)

target_link_libraries(BytePairRam ZLIB::ZLIB)
