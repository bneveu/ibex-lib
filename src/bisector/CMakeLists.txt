# see arithmetic/CMakeLists.txt for comments

target_sources (ibex PRIVATE
  ${CMAKE_CURRENT_SOURCE_DIR}/ibex_Bisection.h
  ${CMAKE_CURRENT_SOURCE_DIR}/ibex_BisectionPoint.h
  ${CMAKE_CURRENT_SOURCE_DIR}/ibex_Bsc.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/ibex_Bsc.h
  ${CMAKE_CURRENT_SOURCE_DIR}/ibex_LargestFirst.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/ibex_LargestFirst.h
  ${CMAKE_CURRENT_SOURCE_DIR}/ibex_LSmear.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/ibex_LSmear.h
  ${CMAKE_CURRENT_SOURCE_DIR}/ibex_NoBisectableVariableException.h
  ${CMAKE_CURRENT_SOURCE_DIR}/ibex_OptimLargestFirst.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/ibex_OptimLargestFirst.h
  ${CMAKE_CURRENT_SOURCE_DIR}/ibex_RoundRobin.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/ibex_RoundRobin.h
  ${CMAKE_CURRENT_SOURCE_DIR}/ibex_SmearFunction.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/ibex_SmearFunction.h
  ${CMAKE_CURRENT_SOURCE_DIR}/ibex_MinlpLargestFirst.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/ibex_MinlpLargestFirst.h
  ${CMAKE_CURRENT_SOURCE_DIR}/ibex_MinlpSmearSum.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/ibex_MinlpSmearSum.h
  ${CMAKE_CURRENT_SOURCE_DIR}/ibex_MinlpSmearSumRelative.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/ibex_MinlpSmearSumRelative.h
  ${CMAKE_CURRENT_SOURCE_DIR}/ibex_QibexLargestFirst.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/ibex_QibexLargestFirst.h
   ${CMAKE_CURRENT_SOURCE_DIR}/ibex_QibexSmearsumRelative.cpp
   ${CMAKE_CURRENT_SOURCE_DIR}/ibex_QibexSmearsumRelative.h
   ${CMAKE_CURRENT_SOURCE_DIR}/ibex_QibexRoundRobin.cpp
   ${CMAKE_CURRENT_SOURCE_DIR}/ibex_QibexRoundRobin.h
  )

target_include_directories (ibex PUBLIC
                              "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>")
