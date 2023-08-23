find_package(Cairo)
include_directories(${CAIRO_INCLUDE_DIR})



find_package(CUDA REQUIRED)
add_subdirectory(cuda)
include_directories(${CMAKE_SOURCE_DIR}/src/cuda)

find_package(OpenGL REQUIRED)
find_package(GLUT REQUIRED)
include_directories( ${OPENGL_INCLUDE_DIRS}  ${GLUT_INCLUDE_DIRS} )


find_package(GSL REQUIRED)
link_directories(${GSL_LIBRARY_DIRS})
include_directories(${GSL_INCLUDE_DIRS})

include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-std=c++0x" COMPILER_SUPPORTS_CXX0X)
if(COMPILER_SUPPORTS_CXX0X)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")
else()
        message(STATUS "The compiler ${CMAKE_CXX_COMPILER} has no C++11 support. Please use a different C++ compiler.")
endif()

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS}  -g -O3 -Wall ")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}  -g   -O3 -Wall")
set(CMAKE_EXE_LINKER_FLAGS "${GCC_COVERAGE_LINK_FLAGS} ${CMAKE_EXE_LINKER_FLAGS}")


find_package(Boost REQUIRED COMPONENTS random regex)
link_directories(${BOOST_LIBRARY_DIRS})
include_directories(${BOOST_INCLUDE_DIRS})


find_package(ITK REQUIRED)
include(${ITK_USE_FILE})


find_package(VTK REQUIRED)
include(${VTK_USE_FILE})


cuda_add_executable (proj main.cpp itkScatterplotDrawer.h
    dataset.cpp parsedata.cpp projection2d.cpp skeleton.cpp principalgraph.cpp
    propagate2d.cpp transform2d.cpp plotanalysis.cpp principalcurve.cpp curve2d.cpp draw.cpp pipeline.cpp skelft_core.cxx
    skeleton/flags.cpp skeleton/fmm.cpp skeleton/genrl.cpp skeleton/io.cpp skeleton/mfmm.cpp skeletongenerator.cxx
    )

#cuda_add_executable (proj main.cpp
#    dataset.cpp parsedata.cpp projection2d.cpp skeleton.cpp principalgraph.cpp
#    propagate2d.cpp transform2d.cpp plotanalysis.cpp principalcurve.cpp curve2d.cpp draw.cpp pipeline.cpp skelft_core.cxx
#    skeleton/flags.cpp skeleton/fmm.cpp skeleton/genrl.cpp skeleton/io.cpp skeleton/mfmm.cpp skeletongenerator.cxx
#    )

add_dependencies(proj skelGPU)


target_link_libraries(proj skelGPU ${CAIRO_LIBRARIES} ${ITK_LIBRARIES}   ${VTK_LIBRARIES} ${Boost_LIBRARIES} ${OPENGL_LIBRARIES} ${GLUT_LIBRARY} ${GSL_LIBRARIES} )
set_target_properties(proj PROPERTIES RUNTIME_OUTPUT_DIRECTORY ..)



find_package(PkgConfig)
pkg_search_module(GTKMM REQUIRED gtkmm-3.0)
link_directories(${GTKMM_LIBRARY_DIRS})
include_directories(${GTKMM_INCLUDE_DIRS})



#add_dependencies( skelGPU)


cuda_add_executable (hks guimain.cpp
    dataset.cpp parsedata.cpp projection2d.cpp skeleton.cpp
    propagate2d.cpp transform2d.cpp plotanalysis.cpp principalcurve.cpp principalgraph.cpp curve2d.cpp draw.cpp pipeline.cpp skelft_core.cxx
    randomprojection.cpp mds.cpp tsne/tsne.cpp tsne/sptree.cpp
    gui/projectionmatrix.cpp gui/plist.cpp gui/imagecontroller.cpp
    skeleton/flags.cpp skeleton/fmm.cpp skeleton/genrl.cpp skeleton/io.cpp skeleton/mfmm.cpp
    skeletongenerator.cxx
    )
target_link_libraries(hks skelGPU ${GTKMM_LIBRARIES} ${BOOST_RANDOM_LIBRARY} ${ITK_LIBRARIES} ${GSL_LIBRARIES}  ${OPENGL_LIBRARIES} ${GLUT_LIBRARY} )
set_target_properties(hks PROPERTIES RUNTIME_OUTPUT_DIRECTORY ..)
