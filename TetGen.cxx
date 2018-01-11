//
// Created by default on 17. 12. 16.
//

#include "TetGen.h"

void TetGen::tetrahedralize(vtkPoints *iPoints, vtkCellArray *iPolys) {
	clock_t ts[5]; // Timing informations (defined in time.h)

	tetGenMesh.b = new TetGenBehavior;
	tetGenMesh.b->plc = 1;
	tetGenMesh.b->verbose = 1;
	strcpy(tetGenMesh.b->outfilename, "TetGenResult");
	tetGenMesh.vtkToTetGenMesh(iPoints, iPolys);

	tetGenMesh.initializepools();
	tetGenMesh.transfernodes();

	// 0,0,0,3060,820,688
	exactinit(0, 0, 0, tetGenMesh.xmax - tetGenMesh.xmin, tetGenMesh.ymax - tetGenMesh.ymin, tetGenMesh.zmax - tetGenMesh.zmin);


	tetGenMesh.incrementaldelaunay(ts[0]);


	tetGenMesh.meshsurface();
	tetGenMesh.constraineddelaunay(ts[0]);

	tetGenMesh.outnodes(NULL);
	tetGenMesh.outelements(NULL);


//	if (!b->nojettison && ((tetGenMesh.dupverts > 0) || (tetGenMesh.unuverts > 0)
//						   || (b->refine && (in->numberofcorners == 10)))) {
//		tetGenMesh.jettisonnodes();
//	}
//
//	if ((b->order == 2) && !b->convex) {
//		tetGenMesh.highorder();
//	}
}