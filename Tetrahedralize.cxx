//
// Created by default on 17. 12. 16.
//

#include "Tetrahedralize.h"

vtkPolyData Tetrahedralize::GetOutput(vtkPolyData *input) {
	float meanX;
	float meanY;
	float meanZ;
	vtkIdType size;
	vtkPoints *points;
	double *point;

	points = input->GetPoints();
	size = points->GetNumberOfPoints();

	// Get sum of variable of each points;
	for (vtkIdType i = 0; i < size; ++i) {
		point = points->GetPoint(i);
		meanX += point[0];
		meanY += point[1];
		meanZ += point[2];
	}

	// Get mean variable of each points
	meanX /= size;
	meanY /= size;
	meanZ /= size;

	printf("(%f, %f, %f)", meanX, meanY, meanZ);


	// Put the point in the point set
	points->InsertNextPoint(meanX, meanY, meanZ);


}