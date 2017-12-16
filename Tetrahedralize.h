//
// Created by default on 17. 12. 16.
//

#ifndef SHAPE_RECONSTRUCTION_TETRAHEDRALIZE_H
#define SHAPE_RECONSTRUCTION_TETRAHEDRALIZE_H

#include <vtkPolyData.h>

class Tetrahedralize {
public:
	vtkPolyData GetOutput(vtkPolyData *input);
};


#endif //SHAPE_RECONSTRUCTION_TETRAHEDRALIZE_H
