cmake_minimum_required(VERSION 3.3 FATAL_ERROR)


PROJECT(Shape_Reconstruction)
set(VTK_DIR "/home/default/Documents/Dev/VTK/build/")
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

add_executable(Surface MACOSX_BUNDLE Surface.cxx RawDataReader.cxx RawDataReader.h)

target_link_libraries(Surface ${VTK_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})
