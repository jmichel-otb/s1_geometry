cmake_minimum_required(VERSION 2.6)
SET(CMAKE_MODULE_PATH "${${PROJECT_NAME}_SOURCE_DIR}/cmake;${CMAKE_MODULE_PATH}")

# Boost
find_package ( Boost
  1.35.0
  REQUIRED
  COMPONENTS date_time
  )

MESSAGE("Boost libraries used and found: ${Boost_LIBRARIES}")

# Ossim
find_package(ossim)
message(STATUS "OSSIM_INSTALL_PREFIX: $ENV{OSSIM_INSTALL_PREFIX}")
if (NOT OSSIM_FOUND)
    message(FATAL_ERROR "Cannot find OSSIM.")
endif()


# s1_inversemodel
include_directories(${CMAKE_CURRENT_SOURCE_DIR} ${OSSIM_INCLUDE_DIR})

add_executable(s1_inversemodel
    s1_inversemodel.cpp ossimSarSensorModel.cpp
    ossimSentinel1SarSensorModel.cpp xmlTools.cpp
    )
target_link_libraries(s1_inversemodel
    ${Boost_LIBRARIES} ${OSSIM_LIBRARIES})


#Test executables
add_executable(ossimSarSensorModelTest ossimSarSensorModelTest.cpp ossimSarSensorModel.cpp xmlTools.cpp)
target_link_libraries(ossimSarSensorModelTest
    ${Boost_LIBRARIES} ${OSSIM_LIBRARIES})

add_executable(ossimSentinel1SarSensorModelTest ossimSentinel1SarSensorModelTest.cpp ossimSarSensorModel.cpp ossimSentinel1SarSensorModel.cpp xmlTools.cpp)
target_link_libraries(ossimSentinel1SarSensorModelTest
    ${Boost_LIBRARIES} ${OSSIM_LIBRARIES})

add_executable(ossimTerraSarXSarSensorModelTest ossimTerraSarXSarSensorModelTest.cpp ossimSarSensorModel.cpp ossimTerraSarXSarSensorModel.cpp xmlTools.cpp)
target_link_libraries(ossimTerraSarXSarSensorModelTest
    ${Boost_LIBRARIES} ${OSSIM_LIBRARIES})

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
