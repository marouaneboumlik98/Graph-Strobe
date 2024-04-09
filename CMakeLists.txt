# Almost all CMake files should start with this
# You should always specify a range with the newest
# and oldest tested versions of CMake. This will ensure
# you pick up the best policies.
cmake_minimum_required(VERSION 3.2)

# This is your project statement. You should always list languages;
# Listing the version is nice here since it sets lots of useful variables
project(GraphAligner
        VERSION 1.0
        LANGUAGES CXX)

#target_compile_features(GraphAligner PRIVATE cxx_std_17)
#add_compile_options(-Wall -Wextra -O3 -g -Wno-unused-parameter)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -O3 -g -Wno-unused-parameter")

# If you set any CMAKE_ variables, that can go here.
# (But usually don't do this, except maybe for C++ standard)
# Find packages go here.

# You should usually split this into folders, but this is a simple example

# This is a "default" library, and will match the *** variable setting.
# Other common choices are STATIC, SHARED, and MODULE
# Including header files here helps IDEs but is not required.
# Output libname matches target name, with the usual extensions on your system


include(FetchContent)
FetchContent_Declare(
        sdsl-lite
        GIT_REPOSITORY https://github.com/simongog/sdsl-lite.git
        GIT_TAG        v2.1.1 # adjust tag/branch/commit as needed
)
FetchContent_MakeAvailable(sdsl-lite)

#libraries with cmake config
add_subdirectory(zstr)
#add_subdirectory(libvgio)
add_subdirectory(parallel-hashmap)

#search packages sur system
find_package(Boost COMPONENTS program_options REQUIRED)
find_package(ZLIB REQUIRED)
find_package(PkgConfig REQUIRED)

include(FindPkgConfig)
pkg_check_modules(JEMALLOC REQUIRED jemalloc)
pkg_check_modules(SPARSEHASH REQUIRED libsparsehash)
#pkg_check_modules(SDSL REQUIRED sdsl)

# protobuf
find_package(Protobuf REQUIRED)
include_directories(${Protobuf_INCLUDE_DIRS})
include_directories(${CMAKE_CURRENT_BINARY_DIR})
#add_library(proto-objects OBJECT "${CMAKE_CURRENT_LIST_DIR}/src/vg.proto")
#protobuf_generate_cpp(PROTO_SRCS PROTO_HDRS src/vg.proto) # compilation

#target_sources(
#        bbhash
#        PRIVATE
#        ${CMAKE_CURRENT_SOURCE_DIR}/BBHash/BooPHF.h
#        ${CMAKE_CURRENT_SOURCE_DIR}/BBHash/MinimalHashmap.h
#        ${CMAKE_CURRENT_SOURCE_DIR}/BBHash/ParallelBB.h
#)
#set_target_properties(bbhash PROPERTIES LINKER_LANGUAGE CXX)

#target_sources(
#        concurrentqueue
#        PRIVATE
#        ${CMAKE_CURRENT_SOURCE_DIR}/concurrentqueue/blockingconcurrentqueue.h
#        ${CMAKE_CURRENT_SOURCE_DIR}/concurrentqueue/concurrentqueue.h
#        #${CMAKE_CURRENT_SOURCE_DIR}/concurrentqueue/internal/concurrentqueue_internal_debug.h
#)
#set_target_properties(concurrentqueue PROPERTIES LINKER_LANGUAGE CXX)


#

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
        MEMfinder/PartSortBWT/src/PartSortBWT.cpp
#        MEMfinder/libsais/src/libsais.c
#        MEMfinder/libsais/src/libsais16.c
#        MEMfinder/libsais/src/libsais64.c
        MEMfinder/src/RankBitvector.cpp
)
#target_include_directories(${PROJECT_NAME}
#        PRIVATE
#        # where the library itself will look for its internal headers
#        ${CMAKE_CURRENT_SOURCE_DIR}/MEMfinder/src
#        PUBLIC
#        # where top-level project will look for the library's public headers
#        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/MEMfinder/include>
#        # where external projects will look for the library's public headers
#        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
#)
set_target_properties(MEMfinder PROPERTIES LINKER_LANGUAGE CXX)


# create library protobuf locally using vg.proto
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
target_include_directories(proto-objects PUBLIC "$<BUILD_INTERFACE:${PROTO_BINARY_DIR}>")
protobuf_generate(
        TARGET proto-objects
        IMPORT_DIRS "${CMAKE_CURRENT_LIST_DIR}/src"
        PROTOC_OUT_DIR "${PROTO_BINARY_DIR}"
)

# Adding something we can run - the executable graphaligner
add_executable(
        graphaligner
        src/Aligner.cpp
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
        #src/vg.pb.cc
)

#include_directories(${parallel-hashmap_SOURCE_DIR})

target_include_directories(
        graphaligner
        PUBLIC
        #bbhash
        ${CMAKE_CURRENT_SOURCE_DIR}/BBHash
        #concurrentqueue
        ${CMAKE_CURRENT_SOURCE_DIR}/concurrentqueue
        #MEMfinder
        ${CMAKE_CURRENT_SOURCE_DIR}/MEMfinder/src
        ${CMAKE_CURRENT_SOURCE_DIR}/MEMfinder
)


# Make sure you link your targets with this command. It can also link libraries and
# even flags, so linking a target that does not exist will not give a configure-time error.

#link my local static librarires
target_link_libraries(
        graphaligner
        PRIVATE
        zstr
        phmap
        proto-objects
        sdsl
        MEMfinder
        #concurrentqueue
        #bbhash
        PUBLIC
        Boost::program_options
        ZLIB::ZLIB
        #PkgConfig::PkgConfig
        ${JEMALLOC_LIBRARIES}
        ${SPARSEHASH_LIBRARIES}

)




