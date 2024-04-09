# Almost all CMake files should start with this
# You should always specify a range with the newest
# and oldest tested versions of CMake. This will ensure
# you pick up the best policies.
cmake_minimum_required(VERSION 3.22)

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

#libraries with cmake config
add_subdirectory(zstr)
#add_subdirectory(libvgio)
add_subdirectory(parallel-hashmap)

#search pacjkages sur system
find_package(Boost COMPONENTS program_options REQUIRED)
find_package(ZLIB REQUIRED)
find_package(PkgConfig REQUIRED)

pkg_check_modules(JEMALLOC REQUIRED jemalloc)
pkg_check_modules(SPARSEHASH REQUIRED libsparsehash)

# protobuf
find_package(Protobuf REQUIRED)
include_directories(${Protobuf_INCLUDE_DIRS})
include_directories(${CMAKE_CURRENT_BINARY_DIR})
#protobuf_generate_cpp(PROTO_SRCS PROTO_HDRS foo.proto) # compilation

add_library(
        BBHash
        BBHash/BooPHF.h
        BBHash/MinimalHashmap.h
        BBHash/ParallelBB.h
)
set_target_properties(BBHash PROPERTIES LINKER_LANGUAGE CXX)

add_library(
        concurrentqueue
        concurrentqueue/blockingconcurrentqueue.h
        concurrentqueue/concurrentqueue.h
        concurrentqueue/internal/concurrentqueue_internal_debug.h
)
set_target_properties(concurrentqueue PROPERTIES LINKER_LANGUAGE CXX)

add_library(
        MEMfinder
        MEMfinder/src/FlatRanks.h
        MEMfinder/src/FMIndex.h
        MEMfinder/src/MEMfinder.h
        MEMfinder/src/RankBitvector.h
        MEMfinder/src/ReverseComplementView.h
        MEMfinder/src/Serialize.h
        MEMfinder/src/WaveletTree.h
        MEMfinder/PartSortBWT/src/PartSortBWT.h
        MEMfinder/libsais/src/libsais.h
        MEMfinder/libsais/src/libsais16.h
        MEMfinder/libsais/src/libsais64.h
)
set_target_properties(MEMfinder PROPERTIES LINKER_LANGUAGE CXX)

#add_library(
#        parallel-hashmap
#        parallel-hashmap/parallel_hashmap/btree.h
#        parallel-hashmap/parallel_hashmap/meminfo.h
#        parallel-hashmap/parallel_hashmap/phmap.h
#        parallel-hashmap/parallel_hashmap/phmap_base.h
#        parallel-hashmap/parallel_hashmap/phmap_bits.h
#        parallel-hashmap/parallel_hashmap/phmap_config.h
#        parallel-hashmap/parallel_hashmap/phmap_dump.h
#        parallel-hashmap/parallel_hashmap/phmap_fwd_decl.h
#        parallel-hashmap/parallel_hashmap/phmap_utils.h
#)
#set_target_properties(parallel-hashmap PROPERTIES LINKER_LANGUAGE CXX)

#add_library(
#        zstr
#        ./zstr/src/zstr_make_unique_polyfill.h
#        ./zstr/src/zstr.hpp
#        ./zstr/src/strict_fstream.hpp
#)
#set_target_properties(zstr PROPERTIES LINKER_LANGUAGE CXX)

# Link each target with other targets or add options, etc.

# Adding something we can run - Output name matches target name
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

# Make sure you link your targets with this command. It can also link libraries and
# even flags, so linking a target that does not exist will not give a configure-time error.

#link my local static librarires
target_link_libraries(
        graphaligner
        PRIVATE
        zstr
        parallel-hashmap
        MEMfinder
        concurrentqueue
        BBHash
        PUBLIC
        #${PROTO_SRCS} ${PROTO_HDRS}
        Boost::program_options
        ZLIB::ZLIB
        #PkgConfig::PkgConfig
        #${JEMALLOC_INCLUDE_DIRS}
        #${SPARSEHASH_INCLUDE_DIRS}

)



