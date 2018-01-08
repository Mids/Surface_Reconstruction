//
// Created by default on 18. 1. 2.
//

#ifndef SHAPE_RECONSTRUCTION_TETGENIO_H
#define SHAPE_RECONSTRUCTION_TETGENIO_H


//#include "TetGen.h"

#include <cstdio>
#include "predicates.h"

class TetGenIO {


public:

	// A "polygon" describes a simple polygon (no holes). It is not necessarily
	//   convex. Each polygon contains a number of corners (points) and the same
	//   number of sides (edges).  The points of the polygon must be given in
	//   either counterclockwise or clockwise order and they form a ring, so
	//   every two consecutive points forms an edge of the polygon.
	typedef struct {
		int *vertexlist;
		int numberofvertices;
	} polygon;

	// A "facet" describes a polygonal region possibly with holes, edges, and
	//   points floating in it.  Each facet consists of a list of polygons and
	//   a list of hole points (which lie strictly inside holes).
	typedef struct {
		polygon *polygonlist;
		int numberofpolygons;
		REAL *holelist;
		int numberofholes;
	} facet;

	// A "voroedge" is an edge of the Voronoi diagram. It corresponds to a
	//   Delaunay face.  Each voroedge is either a line segment connecting
	//   two Voronoi vertices or a ray starting from a Voronoi vertex to an
	//   "infinite vertex".  'v1' and 'v2' are two indices pointing to the
	//   list of Voronoi vertices. 'v1' must be non-negative, while 'v2' may
	//   be -1 if it is a ray, in this case, the unit normal of this ray is
	//   given in 'vnormal'.
	typedef struct {
		int v1, v2;
		REAL vnormal[3];
	} voroedge;

	// A "vorofacet" is an facet of the Voronoi diagram. It corresponds to a
	//   Delaunay edge.  Each Voronoi facet is a convex polygon formed by a
	//   list of Voronoi edges, it may not be closed.  'c1' and 'c2' are two
	//   indices pointing into the list of Voronoi cells, i.e., the two cells
	//   share this facet.  'elist' is an array of indices pointing into the
	//   list of Voronoi edges, 'elist[0]' saves the number of Voronoi edges
	//   (including rays) of this facet.
	typedef struct {
		int c1, c2;
		int *elist;
	} vorofacet;


	// Additional parameters associated with an input (or mesh) vertex.
	//   These informations are provided by CAD libraries.
	typedef struct {
		REAL uv[2];
		int tag;
		int type; // 0, 1, or 2.
	} pointparam;

	// Callback functions for meshing PSCs.
	typedef REAL (*GetVertexParamOnEdge)(void *, int, int);

	typedef void (*GetSteinerOnEdge)(void *, int, REAL, REAL *);

	typedef void (*GetVertexParamOnFace)(void *, int, int, REAL *);

	typedef void (*GetEdgeSteinerParamOnFace)(void *, int, REAL, int, REAL *);

	typedef void (*GetSteinerOnFace)(void *, int, REAL *, REAL *);

	// A callback function for mesh refinement.
	typedef bool (*TetSizeFunc)(REAL *, REAL *, REAL *, REAL *, REAL *, REAL);

	// Items are numbered starting from 'firstnumber' (0 or 1), default is 0.
	int firstnumber;

	// Dimension of the mesh (2 or 3), default is 3.
	int mesh_dim;

	// Does the lines in .node file contain index or not, default is 1.
	int useindex;

	// 'pointlist':  An array of point coordinates.  The first point's x
	//   coordinate is at index [0] and its y coordinate at index [1], its
	//   z coordinate is at index [2], followed by the coordinates of the
	//   remaining points.  Each point occupies three REALs.
	// 'pointattributelist':  An array of point attributes.  Each point's
	//   attributes occupy 'numberofpointattributes' REALs.
	// 'pointmtrlist': An array of metric tensors at points. Each point's
	//   tensor occupies 'numberofpointmtr' REALs.
	// 'pointmarkerlist':  An array of point markers; one integer per point.
	REAL *pointlist;
	REAL *pointattributelist;
	REAL *pointmtrlist;
	int *pointmarkerlist;
	pointparam *pointparamlist;
	int numberofpoints;
	int numberofpointattributes;
	int numberofpointmtrs;

	// 'tetrahedronlist':  An array of tetrahedron corners.  The first
	//   tetrahedron's first corner is at index [0], followed by its other
	//   corners, followed by six nodes on the edges of the tetrahedron if the
	//   second order option (-o2) is applied. Each tetrahedron occupies
	//   'numberofcorners' ints.  The second order nodes are ouput only.
	// 'tetrahedronattributelist':  An array of tetrahedron attributes.  Each
	//   tetrahedron's attributes occupy 'numberoftetrahedronattributes' REALs.
	// 'tetrahedronvolumelist':  An array of constraints, i.e. tetrahedron's
	//   volume; one REAL per element.  Input only.
	// 'neighborlist':  An array of tetrahedron neighbors; 4 ints per element.
	//   Output only.
	int *tetrahedronlist;
	REAL *tetrahedronattributelist;
	REAL *tetrahedronvolumelist;
	int *neighborlist;
	int numberoftetrahedra;
	int numberofcorners;
	int numberoftetrahedronattributes;

	// 'facetlist':  An array of facets.  Each entry is a structure of facet.
	// 'facetmarkerlist':  An array of facet markers; one int per facet.
	facet *facetlist;
	int *facetmarkerlist;
	int numberoffacets;

	// 'holelist':  An array of holes (in volume).  Each hole is given by a
	//   seed (point) which lies strictly inside it. The first seed's x, y and z
	//   coordinates are at indices [0], [1] and [2], followed by the
	//   remaining seeds.  Three REALs per hole.
	REAL *holelist;
	int numberofholes;

	// 'regionlist': An array of regions (subdomains).  Each region is given by
	//   a seed (point) which lies strictly inside it. The first seed's x, y and
	//   z coordinates are at indices [0], [1] and [2], followed by the regional
	//   attribute at index [3], followed by the maximum volume at index [4].
	//   Five REALs per region.
	// Note that each regional attribute is used only if you select the 'A'
	//   switch, and each volume constraint is used only if you select the
	//   'a' switch (with no number following).
	REAL *regionlist;
	int numberofregions;

	// 'facetconstraintlist':  An array of facet constraints.  Each constraint
	//   specifies a maximum area bound on the subfaces of that facet.  The
	//   first facet constraint is given by a facet marker at index [0] and its
	//   maximum area bound at index [1], followed by the remaining facet con-
	//   straints. Two REALs per facet constraint.  Note: the facet marker is
	//   actually an integer.
	REAL *facetconstraintlist;
	int numberoffacetconstraints;

	// 'segmentconstraintlist': An array of segment constraints. Each constraint
	//   specifies a maximum length bound on the subsegments of that segment.
	//   The first constraint is given by the two endpoints of the segment at
	//   index [0] and [1], and the maximum length bound at index [2], followed
	//   by the remaining segment constraints.  Three REALs per constraint.
	//   Note the segment endpoints are actually integers.
	REAL *segmentconstraintlist;
	int numberofsegmentconstraints;


	// 'trifacelist':  An array of face (triangle) corners.  The first face's
	//   three corners are at indices [0], [1] and [2], followed by the remaining
	//   faces.  Three ints per face.
	// 'trifacemarkerlist':  An array of face markers; one int per face.
	// 'o2facelist':  An array of second order nodes (on the edges) of the face.
	//   It is output only if the second order option (-o2) is applied. The
	//   first face's three second order nodes are at [0], [1], and [2],
	//   followed by the remaining faces.  Three ints per face.
	// 'adjtetlist':  An array of adjacent tetrahedra to the faces. The first
	//   face's two adjacent tetrahedra are at indices [0] and [1], followed by
	//   the remaining faces.  A '-1' indicates outside (no adj. tet). This list
	//   is output when '-nn' switch is used. Output only.
	int *trifacelist;
	int *trifacemarkerlist;
	int *o2facelist;
	int *adjtetlist;
	int numberoftrifaces;

	// 'edgelist':  An array of edge endpoints.  The first edge's endpoints
	//   are at indices [0] and [1], followed by the remaining edges.
	//   Two ints per edge.
	// 'edgemarkerlist':  An array of edge markers; one int per edge.
	// 'o2edgelist':  An array of midpoints of edges. It is output only if the
	//   second order option (-o2) is applied. One int per edge.
	// 'edgeadjtetlist':  An array of adjacent tetrahedra to the edges.  One
	//   tetrahedron (an integer) per edge.
	int *edgelist;
	int *edgemarkerlist;
	int *o2edgelist;
	int *edgeadjtetlist;
	int numberofedges;

	// 'vpointlist':  An array of Voronoi vertex coordinates (like pointlist).
	// 'vedgelist':  An array of Voronoi edges.  Each entry is a 'voroedge'.
	// 'vfacetlist':  An array of Voronoi facets. Each entry is a 'vorofacet'.
	// 'vcelllist':  An array of Voronoi cells.  Each entry is an array of
	//   indices pointing into 'vfacetlist'. The 0th entry is used to store
	//   the length of this array.
	REAL *vpointlist;
	voroedge *vedgelist;
	vorofacet *vfacetlist;
	int **vcelllist;
	int numberofvpoints;
	int numberofvedges;
	int numberofvfacets;
	int numberofvcells;

	// Variable (and callback functions) for meshing PSCs.
	void *geomhandle;
	GetVertexParamOnEdge getvertexparamonedge;
	GetSteinerOnEdge getsteineronedge;
	GetVertexParamOnFace getvertexparamonface;
	GetEdgeSteinerParamOnFace getedgesteinerparamonface;
	GetSteinerOnFace getsteineronface;

	// A callback function.
	TetSizeFunc tetunsuitable;

	// Input & output routines.
	bool load_node_call(FILE *infile, int markers, int uvflag, char *);

	bool load_node(char *);

	bool load_edge(char *);

	bool load_face(char *);

	bool load_tet(char *);

	bool load_vol(char *);

	bool load_var(char *);

	bool load_mtr(char *);

	bool load_pbc(char *);

	bool load_poly(char *);

	bool load_off(char *);

	bool load_ply(char *);

	bool load_stl(char *);

	bool load_vtk(char *);

	bool load_medit(char *, int);

	bool load_plc(char *, int);

	bool load_tetmesh(char *, int);

	void save_nodes(char *);

	void save_elements(char *);

	void save_faces(char *);

	void save_edges(char *);

	void save_neighbors(char *);

	void save_poly(char *);

	void save_faces2smesh(char *);

	// Read line and parse string functions.
	char *readline(char *string, FILE *infile, int *linenumber);

	char *findnextfield(char *string);

	char *readnumberline(char *string, FILE *infile, char *infilename);

	char *findnextnumber(char *string);

	static void init(polygon *p) {
		p->vertexlist = (int *) NULL;
		p->numberofvertices = 0;
	}

	static void init(facet *f) {
		f->polygonlist = (polygon *) NULL;
		f->numberofpolygons = 0;
		f->holelist = (REAL *) NULL;
		f->numberofholes = 0;
	}

	// Initialize routine.
	void initialize() {
		firstnumber = 0;
		mesh_dim = 3;
		useindex = 1;

		pointlist = (REAL *) NULL;
		pointattributelist = (REAL *) NULL;
		pointmtrlist = (REAL *) NULL;
		pointmarkerlist = (int *) NULL;
		pointparamlist = (pointparam *) NULL;
		numberofpoints = 0;
		numberofpointattributes = 0;
		numberofpointmtrs = 0;

		tetrahedronlist = (int *) NULL;
		tetrahedronattributelist = (REAL *) NULL;
		tetrahedronvolumelist = (REAL *) NULL;
		neighborlist = (int *) NULL;
		numberoftetrahedra = 0;
		numberofcorners = 4;
		numberoftetrahedronattributes = 0;

		trifacelist = (int *) NULL;
		trifacemarkerlist = (int *) NULL;
		o2facelist = (int *) NULL;
		adjtetlist = (int *) NULL;
		numberoftrifaces = 0;

		edgelist = (int *) NULL;
		edgemarkerlist = (int *) NULL;
		o2edgelist = (int *) NULL;
		edgeadjtetlist = (int *) NULL;
		numberofedges = 0;

		facetlist = (facet *) NULL;
		facetmarkerlist = (int *) NULL;
		numberoffacets = 0;

		holelist = (REAL *) NULL;
		numberofholes = 0;

		regionlist = (REAL *) NULL;
		numberofregions = 0;

		facetconstraintlist = (REAL *) NULL;
		numberoffacetconstraints = 0;
		segmentconstraintlist = (REAL *) NULL;
		numberofsegmentconstraints = 0;


		vpointlist = (REAL *) NULL;
		vedgelist = (voroedge *) NULL;
		vfacetlist = (vorofacet *) NULL;
		vcelllist = (int **) NULL;
		numberofvpoints = 0;
		numberofvedges = 0;
		numberofvfacets = 0;
		numberofvcells = 0;

		tetunsuitable = NULL;

		geomhandle = NULL;
		getvertexparamonedge = NULL;
		getsteineronedge = NULL;
		getvertexparamonface = NULL;
		getedgesteinerparamonface = NULL;
		getsteineronface = NULL;
	}

	// Free the memory allocated in 'tetgenio'.  Note that it assumes that the
	//   memory was allocated by the "new" operator (C++).
	void deinitialize() {
		int i, j;

		if (pointlist != (REAL *) NULL) {
			delete[] pointlist;
		}
		if (pointattributelist != (REAL *) NULL) {
			delete[] pointattributelist;
		}
		if (pointmtrlist != (REAL *) NULL) {
			delete[] pointmtrlist;
		}
		if (pointmarkerlist != (int *) NULL) {
			delete[] pointmarkerlist;
		}
		if (pointparamlist != (pointparam *) NULL) {
			delete[] pointparamlist;
		}

		if (tetrahedronlist != (int *) NULL) {
			delete[] tetrahedronlist;
		}
		if (tetrahedronattributelist != (REAL *) NULL) {
			delete[] tetrahedronattributelist;
		}
		if (tetrahedronvolumelist != (REAL *) NULL) {
			delete[] tetrahedronvolumelist;
		}
		if (neighborlist != (int *) NULL) {
			delete[] neighborlist;
		}

		if (trifacelist != (int *) NULL) {
			delete[] trifacelist;
		}
		if (trifacemarkerlist != (int *) NULL) {
			delete[] trifacemarkerlist;
		}
		if (o2facelist != (int *) NULL) {
			delete[] o2facelist;
		}
		if (adjtetlist != (int *) NULL) {
			delete[] adjtetlist;
		}

		if (edgelist != (int *) NULL) {
			delete[] edgelist;
		}
		if (edgemarkerlist != (int *) NULL) {
			delete[] edgemarkerlist;
		}
		if (o2edgelist != (int *) NULL) {
			delete[] o2edgelist;
		}
		if (edgeadjtetlist != (int *) NULL) {
			delete[] edgeadjtetlist;
		}

		if (facetlist != (facet *) NULL) {
			facet *f;
			polygon *p;
			for (i = 0; i < numberoffacets; i++) {
				f = &facetlist[i];
				for (j = 0; j < f->numberofpolygons; j++) {
					p = &f->polygonlist[j];
					delete[] p->vertexlist;
				}
				delete[] f->polygonlist;
				if (f->holelist != (REAL *) NULL) {
					delete[] f->holelist;
				}
			}
			delete[] facetlist;
		}
		if (facetmarkerlist != (int *) NULL) {
			delete[] facetmarkerlist;
		}

		if (holelist != (REAL *) NULL) {
			delete[] holelist;
		}
		if (regionlist != (REAL *) NULL) {
			delete[] regionlist;
		}
		if (facetconstraintlist != (REAL *) NULL) {
			delete[] facetconstraintlist;
		}
		if (segmentconstraintlist != (REAL *) NULL) {
			delete[] segmentconstraintlist;
		}
		if (vpointlist != (REAL *) NULL) {
			delete[] vpointlist;
		}
		if (vedgelist != (voroedge *) NULL) {
			delete[] vedgelist;
		}
		if (vfacetlist != (vorofacet *) NULL) {
			for (i = 0; i < numberofvfacets; i++) {
				delete[] vfacetlist[i].elist;
			}
			delete[] vfacetlist;
		}
		if (vcelllist != (int **) NULL) {
			for (i = 0; i < numberofvcells; i++) {
				delete[] vcelllist[i];
			}
			delete[] vcelllist;
		}
	}

	// Constructor & destructor.
	TetGenIO() { initialize(); }

	~TetGenIO() { deinitialize(); }

}; // class tetgenio


#endif //SHAPE_RECONSTRUCTION_TETGENIO_H
