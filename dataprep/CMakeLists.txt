cmake_minimum_required(VERSION 3.6)

project(dataprep CXX)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/")

set( larlite_DIR $ENV{LARLITE_LIBDIR}/cmake/larlite )
set( larcv_DIR $ENV{LARCV_LIBDIR}/cmake/larcv )
find_package( ROOT REQUIRED )
find_package( larlite REQUIRED )
find_package( larcv REQUIRED )

#add_executable
add_executable(collect_testdata collect_testdata.cxx)

# To set the C++ standard to be used, only one of the following two lines should
# be un-commented.
#set_property(TARGET demo PROPERTY CXX_STANDARD 14)
set_property(TARGET collect_testdata PROPERTY CXX_STANDARD 14)

include_directories( ${LARLITE_INCLUDE_DIR} )
include_directories( ${LARCV_INCLUDE_DIR} )
include_directories( ${ROOT_INCLUDE_DIRS} )

set( LINK_LIBS larlite_DataFormat larlite_LArUtil larlite_Base LArCVCoreDataFormat LArCVCoreBase ${ROOT_LIBRARIES} )

target_link_libraries(collect_testdata ${LINK_LIBS})

