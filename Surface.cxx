/*=========================================================================

  Program:   Visualization Toolkit
  Module:    Cube.cxx

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
// This example shows how to manually create vtkPolyData.

#include <vtkDelaunay2D.h>
#include "vtkActor.h"
#include "vtkCamera.h"
#include "vtkFloatArray.h"
#include "vtkPointData.h"
#include "vtkPolyDataMapper.h"
#include "vtkRenderWindow.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkRenderer.h"
#include "vtkProperty.h"

#define SCALE 10
#define KINECT_X_RES 640
#define KINECT_Y_RES 480

static const int MESH_X_RES = KINECT_X_RES / SCALE;
static const int MESH_Y_RES = KINECT_Y_RES / SCALE;

float *getDepthData();

vtkFloatArray * getVertices(float *depthData);

int main() {
	int i;

	// Get vertices and triangles from depth data
	float *depthData = getDepthData();
	vtkFloatArray *vertices = getVertices(depthData);
	delete (depthData);


	// Create the building blocks of polydata
	vtkPolyData *surface = vtkPolyData::New();
	vtkPoints *points = vtkPoints::New();
	vtkFloatArray *scalars = vtkFloatArray::New();


	// Load the point, cell, and data attributes.
	points->SetData(vertices);
	for (i = 0; i < MESH_X_RES * MESH_Y_RES; i++) scalars->InsertTuple1(i, i);


	// Assign the pieces to the vtkPolyData.
	surface->SetPoints(points);
	points->Delete();
	surface->GetPointData()->SetScalars(scalars);
	scalars->Delete();


	// Perform a 2D Delaunay triangulation on them.
	vtkDelaunay2D *delaunay = vtkDelaunay2D::New();
	delaunay->SetInputData(surface);
	delaunay->SetTolerance(0.001);


	// Now we'll look at it.
	vtkPolyDataMapper *PolyMapper = vtkPolyDataMapper::New();
	PolyMapper->SetInputConnection(delaunay->GetOutputPort());
	PolyMapper->SetScalarRange(0, MESH_X_RES * MESH_Y_RES - 1);
	vtkActor *surfaceActor = vtkActor::New();
	surfaceActor->SetMapper(PolyMapper);
	surfaceActor->GetProperty()->SetRepresentationToWireframe();

	// The usual rendering stuff.
	vtkCamera *camera = vtkCamera::New();
	camera->SetPosition(1, 1, 1);
	camera->SetFocalPoint(0, 0, 0);

	vtkRenderer *renderer = vtkRenderer::New();
	vtkRenderWindow *renWin = vtkRenderWindow::New();
	renWin->AddRenderer(renderer);

	vtkRenderWindowInteractor *iren = vtkRenderWindowInteractor::New();
	iren->SetRenderWindow(renWin);

	renderer->AddActor(surfaceActor);
	renderer->SetActiveCamera(camera);
	renderer->ResetCamera();
	renderer->SetBackground(1, 1, 1);

	renWin->SetSize(300, 300);

	// interact with data
	renWin->Render();
	iren->Start();

	// Clean up
	vertices->Delete();
	surface->Delete();
	PolyMapper->Delete();
	surfaceActor->Delete();
	camera->Delete();
	renderer->Delete();
	renWin->Delete();
	iren->Delete();

	return 0;
}


// Building Sample data of object
float *getDepthData() {
	const float PLANEDEPTH = 50000;
	float *depthData = new float[KINECT_X_RES * KINECT_Y_RES];
	float depth;

	for (int i = 0; i < KINECT_X_RES; ++i)
		for (int j = 0; j < KINECT_Y_RES; ++j) {
			depth = 0;
			if (i > KINECT_X_RES / 6 && i < 5 * KINECT_X_RES / 6 && j > KINECT_Y_RES / 6 && j < 5 * KINECT_Y_RES / 6) {

				int x2 = (KINECT_X_RES / 2 - i) * (KINECT_X_RES / 2 - i);
				int y2 = (KINECT_Y_RES / 2 - j) * (KINECT_Y_RES / 2 - j);
				depth = PLANEDEPTH - x2 / 3 - y2 / 3;
			}

			depthData[i * KINECT_Y_RES + j] = depth;
		}

	return depthData;
}


// Locate vertices by depth ( x , y , depth )
vtkFloatArray *getVertices(float *depthData) {
	vtkFloatArray *points = vtkFloatArray::New();
	points->SetNumberOfComponents(3);

	for (int i = 0; i < MESH_X_RES; ++i)
		for (int j = 0; j < MESH_Y_RES; ++j) {
			if (depthData[i * SCALE * KINECT_Y_RES + j * SCALE] < 100) continue;
			points->InsertNextValue(i * SCALE * 100);
			points->InsertNextValue(j * SCALE * 100);
			points->InsertNextValue(depthData[i * SCALE * KINECT_Y_RES + j * SCALE]);
		}

	return points;
}
