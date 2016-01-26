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

# #TSX tests
# file(GLOB tsx_files ${CMAKE_SOURCE_DIR}/data/tsx/*)
# #Loop over all S1 annotation files
# foreach(entry ${tsx_files})
#   file(GLOB tsx_files ${CMAKE_SOURCE_DIR}/data/tsx/*)
#   get_filename_component(name ${entry} NAME_WE)
#   add_test( tsx_inverse_${name} ossimTerraSarXSarSensorModelTest ${entry} )
#    add_test( tsx_forward_${name} ossimSentinel1SarSensorModelTest 0 ${entry})
# endforeach()

add_test( tsx_inverse_PANGKALANBUUN ossimTerraSarXSarSensorModelTest 1 ${CMAKE_SOURCE_DIR}/data/tsx/PANGKALANBUUN/TSX1_SAR__SSC______SM_D_SRA_20080313T221955_20080313T222003.xml ${CMAKE_SOURCE_DIR}/data/tsx/PANGKALANBUUN/GEOREF.xml)

add_test( tsx_forward_PANGKALANBUUN ossimTerraSarXSarSensorModelTest 0 ${CMAKE_SOURCE_DIR}/data/tsx/PANGKALANBUUN/TSX1_SAR__SSC______SM_D_SRA_20080313T221955_20080313T222003.xml ${CMAKE_SOURCE_DIR}/data/tsx/PANGKALANBUUN/GEOREF.xml)

add_test( tsx_inverse_UPSALA ossimTerraSarXSarSensorModelTest 1 ${CMAKE_SOURCE_DIR}/data/tsx/UPSALA_GLACIER/TSX1_SAR__MGD.xml ${CMAKE_SOURCE_DIR}/data/tsx/UPSALA_GLACIER/GEOREF.xml)

add_test( tsx_forward_UPSALA ossimTerraSarXSarSensorModelTest 0 ${CMAKE_SOURCE_DIR}/data/tsx/UPSALA_GLACIER/TSX1_SAR__MGD.xml ${CMAKE_SOURCE_DIR}/data/tsx/UPSALA_GLACIER/GEOREF.xml)
