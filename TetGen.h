//
// Created by default on 17. 12. 16.
//

#ifndef SHAPE_RECONSTRUCTION_TETRAHEDRALIZE_H
#define SHAPE_RECONSTRUCTION_TETRAHEDRALIZE_H

#include "predicates.h"
#include <vtkPolyData.h>
#include "TetGenMesh.h"


class TetGen {
public:
	TetGenMesh tetGenMesh;

	void tetrahedralize(vtkPoints *iPoints, vtkCellArray *iPolys);


};

#endif //SHAPE_RECONSTRUCTION_TETRAHEDRALIZE_H
