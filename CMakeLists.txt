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

include_directories(${CMAKE_CURRENT_SOURCE_DIR})

add_executable(s1_inversemodel s1_inversemodel.cpp ossimSarSensorModel.cpp ossimSentinel1SarSensorModel.cpp)
target_link_libraries(s1_inversemodel ${OTB_LIBRARIES} ${Boost_LIBRARIES})


#Test executables
add_executable(ossimSarSensorModelTest ossimSarSensorModelTest.cpp ossimSarSensorModel.cpp)
target_link_libraries(ossimSarSensorModelTest ${OTB_LIBRARIES} ${Boost_LIBRARIES})

add_executable(ossimSentinel1SarSensorModelTest ossimSentinel1SarSensorModelTest.cpp ossimSarSensorModel.cpp ossimSentinel1SarSensorModel.cpp)
target_link_libraries(ossimSentinel1SarSensorModelTest ${OTB_LIBRARIES} ${Boost_LIBRARIES})

add_executable(ossimTerraSarXSarSensorModelTest ossimTerraSarXSarSensorModelTest.cpp ossimSarSensorModel.cpp ossimTerraSarXSarSensorModel.cpp)
target_link_libraries(ossimTerraSarXSarSensorModelTest ${OTB_LIBRARIES} ${Boost_LIBRARIES})

#Use ctest
enable_testing()
add_test( ossimSarSensorModelTest ossimSarSensorModelTest )

#S1 tests
file(GLOB s1_files ${CMAKE_SOURCE_DIR}/data/s1/*.xml)

#Loop over all S1 annotation files
foreach(entry ${s1_files})
  get_filename_component(name ${entry} NAME_WE)
  add_test( s1_inverse_${name} ossimSentinel1SarSensorModelTest 1 ${entry})
  add_test( s1_forward_${name} ossimSentinel1SarSensorModelTest 0 ${entry})
endforeach()

#TSX tests
file(GLOB MTDATA ${CMAKE_SOURCE_DIR}/data/tsx/*.xml)
LIST( LENGTH MTDATA count )
#Loop over all S1 annotation files
foreach(index RANGE 0 ${count-1})
  list ( GET MTDATA ${index} entry )
  add_test( ossimTerraSarXSarSensorModelTest_${index} ossimTerraSarXSarSensorModelTest ${entry} )
endforeach()
