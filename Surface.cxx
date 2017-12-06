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

#include "vtkActor.h"
#include "vtkCamera.h"
#include "vtkCellArray.h"
#include "vtkFloatArray.h"
#include "vtkPointData.h"
#include "vtkPoints.h"
#include "vtkPolyData.h"
#include "vtkPolyDataMapper.h"
#include "vtkRenderWindow.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkRenderer.h"

#define KINECT_X_RES 640
#define KINECT_Y_RES 480

float *getDepthData();

float *getVertices(float *depthData);

float *getTriangles();

int main() {
	int i;
	static float x[8][3] = {{0, 0, 0},
							{1, 0, 0},
							{1, 1, 0},
							{0, 1, 0},
							{0, 0, 1},
							{1, 0, 1},
							{1, 1, 1},
							{0, 1, 1}};
	static vtkIdType pts[6][4] = {{0, 1, 2, 3},
								  {4, 5, 6, 7},
								  {0, 1, 5, 4},
								  {1, 2, 6, 5},
								  {2, 3, 7, 6},
								  {3, 0, 4, 7}};

	// Get vertices and triangles from depth data
	float *depthData = getDepthData();
	float *vertices = getVertices(depthData); // TODO : Free data
	delete (depthData);
	float *triangles = getTriangles(); // TODO : Free data




	// We'll create the building blocks of polydata including data attributes.
	vtkPolyData *cube = vtkPolyData::New();
	vtkPoints *points = vtkPoints::New();
	vtkCellArray *polys = vtkCellArray::New();
	vtkFloatArray *scalars = vtkFloatArray::New();

	// Load the point, cell, and data attributes.
	for (i = 0; i < 8; i++) points->InsertPoint(i, x[i]);
	for (i = 0; i < 6; i++) polys->InsertNextCell(4, pts[i]);
	for (i = 0; i < 8; i++) scalars->InsertTuple1(i, i);

	// We now assign the pieces to the vtkPolyData.
	cube->SetPoints(points);
	points->Delete();
	cube->SetPolys(polys);
	polys->Delete();
	cube->GetPointData()->SetScalars(scalars);
	scalars->Delete();

	// Now we'll look at it.
	vtkPolyDataMapper *cubeMapper = vtkPolyDataMapper::New();
	cubeMapper->SetInputData(cube);
	cubeMapper->SetScalarRange(0, 7);
	vtkActor *cubeActor = vtkActor::New();
	cubeActor->SetMapper(cubeMapper);

	// The usual rendering stuff.
	vtkCamera *camera = vtkCamera::New();
	camera->SetPosition(1, 1, 1);
	camera->SetFocalPoint(0, 0, 0);

	vtkRenderer *renderer = vtkRenderer::New();
	vtkRenderWindow *renWin = vtkRenderWindow::New();
	renWin->AddRenderer(renderer);

	vtkRenderWindowInteractor *iren = vtkRenderWindowInteractor::New();
	iren->SetRenderWindow(renWin);

	renderer->AddActor(cubeActor);
	renderer->SetActiveCamera(camera);
	renderer->ResetCamera();
	renderer->SetBackground(1, 1, 1);

	renWin->SetSize(300, 300);

	// interact with data
	renWin->Render();
	iren->Start();

	// Clean up

	cube->Delete();
	cubeMapper->Delete();
	cubeActor->Delete();
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
float *getVertices(float *depthData) {
	float *points = new float[KINECT_X_RES * KINECT_Y_RES * 3];

	for (int i = 0; i < KINECT_X_RES; ++i)
		for (int j = 0; j < KINECT_Y_RES; ++j) {
			int coord = i * KINECT_Y_RES + j;
			points[coord] = i;
			points[coord + 1] = j;
			points[coord + 2] = depthData[coord];
		}

	return points;
}


// Build triangles
float *getTriangles() {
	float *triangles = new float[(KINECT_X_RES - 1) * (KINECT_Y_RES - 1) * 2 * 3]; // Build two triangles

	int triCoord = 0;
	for (int i = 1; i < KINECT_X_RES; ++i)
		for (int j = 1; j < KINECT_Y_RES; ++j) {
			// Right Top Triangle
			triangles[triCoord] = (i - 1) * KINECT_Y_RES + (j - 1); // Left Top
			triangles[triCoord + 1] = (i) * KINECT_Y_RES + (j - 1); // Top
			triangles[triCoord + 2] = i * KINECT_Y_RES + j; // Self

			// Left Bottom Triangle
			triangles[triCoord + 3] = (i - 1) * KINECT_Y_RES + (j - 1); // Left Top
			triangles[triCoord + 4] = (i - 1) * KINECT_Y_RES + (j); // Left
			triangles[triCoord + 5] = i * KINECT_Y_RES + j; // Self

			triCoord += 6;
		}

	return triangles;
}