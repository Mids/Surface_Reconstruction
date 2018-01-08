cmake_minimum_required(VERSION 3.3 FATAL_ERROR)


PROJECT(Shape_Reconstruction)
set(VTK_DIR "../VTK/cmake-build-debug")
find_package(VTK REQUIRED)
find_package(Threads)
vtk_module_config(VTK
        vtkCommonCore
        vtkCommonDataModel
        vtkFiltersGeometry
        vtkInteractionStyle
        vtkRendering${VTK_RENDERING_BACKEND}
        vtkIOXML
        vtkIOPLY
        )
include(${VTK_USE_FILE})

add_executable(Surface Surface.cxx TetGen.cxx TetGen.h RawDataReader.cxx RawDataReader.h TetGenMesh.cxx TetGenMesh.h TetGenBehavior.cxx TetGenBehavior.h predicates.cxx predicates.h TetGenIO.cxx TetGenIO.h)

target_link_libraries(Surface ${VTK_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})
