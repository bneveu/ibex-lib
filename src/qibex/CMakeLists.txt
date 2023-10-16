# see arithmetic/CMakeLists.txt for comments

target_sources (ibex PRIVATE
  
  ${CMAKE_CURRENT_SOURCE_DIR}/ibex_QibexOptimizer.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/ibex_QibexOptimizer.h
  ${CMAKE_CURRENT_SOURCE_DIR}/ibex_QibexLargestFirst.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/ibex_QibexLargestFirst.h
  ${CMAKE_CURRENT_SOURCE_DIR}/ibex_QibexSmearsum.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/ibex_QibexSmearsum.h
  ${CMAKE_CURRENT_SOURCE_DIR}/ibex_QibexSmearsumRelative.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/ibex_QibexSmearsumRelative.h
  ${CMAKE_CURRENT_SOURCE_DIR}/ibex_QibexRoundRobin.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/ibex_QibexRoundRobin.h
)

target_include_directories (ibex PUBLIC
                              "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>")
