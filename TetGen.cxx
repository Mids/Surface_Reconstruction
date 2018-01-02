//
// Created by default on 17. 12. 16.
//

#include "TetGen.h"

vtkPolyData TetGen::tetrahedralize(vtkPolyData *input) {

	TetGenMesh tetGenMesh;
	clock_t ts[5]; // Timing informations (defined in time.h)

//	tetGenMesh.b = b;
//	tetGenMesh.in = in;
//	tetGenMesh.addin = addin;

	tetGenMesh.initializepools();
	tetGenMesh.transfernodes();

	// 0,0,0,3060,820,688
	exactinit(0, 0, 0, tetGenMesh.xmax - tetGenMesh.xmin, tetGenMesh.ymax - tetGenMesh.ymin, tetGenMesh.zmax - tetGenMesh.zmin);


	tetGenMesh.incrementaldelaunay(ts[0]);


	tetGenMesh.meshsurface();
	tetGenMesh.constraineddelaunay(ts[0]);


//	if (!b->nojettison && ((tetGenMesh.dupverts > 0) || (tetGenMesh.unuverts > 0)
//						   || (b->refine && (in->numberofcorners == 10)))) {
//		tetGenMesh.jettisonnodes();
//	}
//
//	if ((b->order == 2) && !b->convex) {
//		tetGenMesh.highorder();
//	}
}