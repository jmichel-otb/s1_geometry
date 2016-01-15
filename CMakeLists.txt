cmake_minimum_required(VERSION 2.6)

# Orfeo ToolBox
find_package(OTB REQUIRED
  COMPONENTS OTBImageIO OTBVectorDataIO OTBProjection OTBStatistics)
if(OTB_FOUND)
  MESSAGE(STATUS "OTB found")
  include(${OTB_USE_FILE})
else(OTB_FOUND)
  message(FATAL_ERROR "OTB not found. Please set OTB_DIR")
endif(OTB_FOUND)


find_package ( Boost
  1.35.0
  REQUIRED
  COMPONENTS date_time
  )

MESSAGE("${Boost_LIBRARIES}")

add_executable(s1_inversemodel s1_inversemodel.cpp)
target_link_libraries(s1_inversemodel ${OTB_LIBRARIES} ${Boost_LIBRARIES})
