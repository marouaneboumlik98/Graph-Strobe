# Almost all CMake files should start with this
# You should always specify a range with the newest
# and oldest tested versions of CMake. This will ensure
# you pick up the best policies.
cmake_minimum_required(VERSION 3.2)

# This is your project statement. You should always list languages;
project(GraphAligner
        VERSION 1.0
        LANGUAGES CXX C)

# compile options
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -O3 -Wno-unused-parameter -DVERSION=1.0")

# If you set any CMAKE_ variables, that can go here.
# (But usually don't do this, except maybe for C++ standard)
# Find packages go here.

# information
message(PROJECT_SOURCE_DIR="${PROJECT_SOURCE_DIR}")
message(CMAKE_CURRENT_SOURCE_DIR="${CMAKE_CURRENT_SOURCE_DIR}/")
message(CMAKE_CURRENT_BINARY_DIR="${CMAKE_CURRENT_BINARY_DIR}")


# libraries retrieved via git
##########################################################################
include(FetchContent)

# sdsl
FetchContent_Declare(
        sdsl-lite
        GIT_REPOSITORY https://github.com/simongog/sdsl-lite.git
        GIT_TAG        v2.1.1 # adjust tag/branch/commit as needed
)

# marouane-modified strobemer
# requires that ssh access is configured properly
FetchContent_Declare(
        strobemer
        GIT_REPOSITORY https://github.com/blinard-BIOINFO/strobemers.git
        GIT_TAG        e22fa4691bbb071e31f812b1c3caea307b756d0e
        SOURCE_SUBDIR  strobemers_cpp
        OVERRIDE_FIND_PACKAGE
)

FetchContent_MakeAvailable(sdsl-lite strobemer)


# libraries directly loaded via their cmakefiles
##########################################################################

add_subdirectory(zstr)
add_subdirectory(concurrentqueue)
add_subdirectory(parallel-hashmap)


# packages searched on system
##########################################################################

set(Boost_USE_STATIC_LIBS  ON)
find_package(Boost 1.70 COMPONENTS program_options REQUIRED)

set(ZLIB_USE_STATIC_LIBS ON)
find_package(ZLIB REQUIRED)

set(PkgConfig_USE_STATIC_LIBS ON)
find_package(PkgConfig REQUIRED)

set(strobemer_USE_STATIC_LIBS ON)
find_package(strobemer REQUIRED)

# modules
##########################################################################
include(FindPkgConfig)
pkg_check_modules(JEMALLOC REQUIRED jemalloc)
pkg_check_modules(SPARSEHASH REQUIRED libsparsehash)

#target_sources(
#        bbhash
#        PRIVATE
#        ${CMAKE_CURRENT_SOURCE_DIR}/BBHash/BooPHF.h
#        ${CMAKE_CURRENT_SOURCE_DIR}/BBHash/MinimalHashmap.h
#        ${CMAKE_CURRENT_SOURCE_DIR}/BBHash/ParallelBB.h
#)
#set_target_properties(bbhash PROPERTIES LINKER_LANGUAGE CXX)

##########################################################################
##########################################################################

# library targets
##########################################################################


# create vg protobuf headers locally: e.g. vg.pb.cc and vg.pb.h
# this uses "protoc" compiler
##########################################################################

# search for protobug install
set(Protobuf_USE_STATIC_LIBS ON)
find_package(Protobuf REQUIRED)
#include_directories(${Protobuf_INCLUDE_DIRS})
#include_directories(${CMAKE_CURRENT_BINARY_DIR})
# create library from vg.proto
add_library(
        proto-objects
        OBJECT
        "${CMAKE_CURRENT_LIST_DIR}/src/vg.proto"
)
target_link_libraries(
        proto-objects
        PUBLIC
        protobuf::libprotobuf
)
set(PROTO_BINARY_DIR "${CMAKE_CURRENT_BINARY_DIR}/protobuf_generated")
file(MAKE_DIRECTORY ${PROTO_BINARY_DIR})  #for some reason this is necessary in some setups
target_include_directories(proto-objects PUBLIC "$<BUILD_INTERFACE:${PROTO_BINARY_DIR}>")
# protoc compilation
protobuf_generate(
        TARGET proto-objects
        IMPORT_DIRS "${CMAKE_CURRENT_LIST_DIR}/src"
        PROTOC_OUT_DIR "${PROTO_BINARY_DIR}"
)

# C library lisais, used by MEMfinder
add_library(
        sais
        STATIC
        MEMfinder/libsais/src/libsais.c
        MEMfinder/libsais/src/libsais64.c
)
set_target_properties(sais PROPERTIES LINKER_LANGUAGE C)  # !!! C not C++

# C++ library PartSortBWT, used by MEMfinder
add_library(
        PartSortBWT
        STATIC
        MEMfinder/PartSortBWT/src/PartSortBWT.cpp

)
set_target_properties(PartSortBWT PROPERTIES LINKER_LANGUAGE CXX)


# MEMfinder library
add_library(
        MEMfinder
        STATIC
        MEMfinder/src/FlatRanks.cpp
        MEMfinder/src/FMIndex.cpp
        MEMfinder/src/MEMfinder.cpp
        MEMfinder/src/RankBitvector.cpp
        MEMfinder/src/ReverseComplementView.cpp
        MEMfinder/src/Serialize.cpp
        MEMfinder/src/WaveletTree.cpp
        MEMfinder/src/RankBitvector.cpp
)
target_link_libraries(
        MEMfinder
        PRIVATE
        sais
        PartSortBWT
)
target_include_directories(
        MEMfinder
        PUBLIC
        #sais
        ${CMAKE_CURRENT_SOURCE_DIR}/MEMfinder/libsais/src
        #PartSortBWT
        ${CMAKE_CURRENT_SOURCE_DIR}/MEMfinder/PartSortBWT/src
)
set_target_properties(MEMfinder PROPERTIES LINKER_LANGUAGE CXX)

##########################################################################
##########################################################################

# the executable graphaligner
add_executable(
        graphaligner
        src/Aligner.cpp
        src/AlignerMain.cpp
        src/AlignmentGraph.cpp
        src/AlignmentSelection.cpp
        src/BigraphToDigraph.cpp
        src/CommonUtils.cpp
        src/DiploidHeuristic.cpp
        src/DNAString.cpp
        src/EValue.cpp
        src/fastqloader.cpp
        src/GfaGraph.cpp
        src/GraphAlignerWrapper.cpp
        src/MEMSeeder.cpp
        src/MinimizerSeeder.cpp
        src/ReadCorrection.cpp
        src/ThreadReadAssertion.cpp
)

target_include_directories(
        graphaligner
        PUBLIC
        #bbhash
        ${CMAKE_CURRENT_SOURCE_DIR}/BBHash
        #MEMfinder
        ${CMAKE_CURRENT_SOURCE_DIR}/MEMfinder/src
        #proto-generated headers
        ${PROTO_BINARY_DIR}/src

)

target_link_libraries(
        graphaligner
        PUBLIC
        zstr
        phmap
        proto-objects
        sdsl
        MEMfinder
        concurrentqueue
        Boost::program_options
        ZLIB::ZLIB
        ${JEMALLOC_LIBRARIES}
        ${SPARSEHASH_LIBRARIES}

)


# the executable graphaligner_strobe
add_executable(
        graphaligner_strobe
        src/Aligner.cpp
        src/AlignerMain.cpp
        src/AlignmentGraph.cpp
        src/AlignmentSelection.cpp
        src/BigraphToDigraph.cpp
        src/CommonUtils.cpp
        src/DiploidHeuristic.cpp
        src/DNAString.cpp
        src/EValue.cpp
        src/fastqloader.cpp
        src/GfaGraph.cpp
        src/GraphAlignerWrapper.cpp
        src/MEMSeeder.cpp
        src/MinimizerSeeder.cpp
        src/ReadCorrection.cpp
        src/ThreadReadAssertion.cpp
)

target_include_directories(
        graphaligner_strobe
        PUBLIC
        #bbhash
        ${CMAKE_CURRENT_SOURCE_DIR}/BBHash
        #MEMfinder
        ${CMAKE_CURRENT_SOURCE_DIR}/MEMfinder/src
        #proto-generated headers
        ${PROTO_BINARY_DIR}/src
        #
)

target_link_libraries(
        graphaligner_strobe
        PUBLIC
        strobemer
        zstr
        phmap
        proto-objects
        sdsl
        MEMfinder
        concurrentqueue
        Boost::program_options
        ZLIB::ZLIB
        ${JEMALLOC_LIBRARIES}
        ${SPARSEHASH_LIBRARIES}

)



