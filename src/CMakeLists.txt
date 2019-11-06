add_library(bvOperators_core
        rank_support.cpp
        select_support.cpp
        waveletTree.cpp)

set(BVOPERATORS_DEBUG_CXXFLAGS "${BVOPERATORS_CXX_FLAGS};-g")
set(BVOPERATORS_RELEASE_CXXFLAGS "${BVOPERATORS_CXX_FLAGS};-O3")

target_include_directories(bvOperators_core PUBLIC $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/include>)
target_compile_options(bvOperators_core PUBLIC "$<$<AND:$<CONFIG:DEBUG>,$<COMPILE_LANGUAGE:CXX>>:${BVOPERATORS_DEBUG_CXXFLAGS}>")
target_compile_options(bvOperators_core PUBLIC "$<$<AND:$<CONFIG:RELEASE>,$<COMPILE_LANGUAGE:CXX>>:${BVOPERATORS_RELEASE_CXXFLAGS}>")

#[[
target_link_libraries(bvOperators_core
        Threads::Threads
        m
        z
        rt)
]]

add_executable(bvOperators BVOperator.cpp)
target_include_directories(bvOperators PUBLIC $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/include>)
target_link_libraries(bvOperators bvOperators_core)
target_compile_options(bvOperators PUBLIC "$<$<AND:$<CONFIG:DEBUG>,$<COMPILE_LANGUAGE:CXX>>:${BVOPERATORS_DEBUG_CXXFLAGS}>")
target_compile_options(bvOperators PUBLIC "$<$<AND:$<CONFIG:RELEASE>,$<COMPILE_LANGUAGE:CXX>>:${BVOPERATORS_RELEASE_CXXFLAGS}>")
