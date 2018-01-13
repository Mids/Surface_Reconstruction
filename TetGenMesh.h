//
// Created by default on 17. 12. 27.
//

#ifndef SHAPE_RECONSTRUCTION_TETGENMESH_H
#define SHAPE_RECONSTRUCTION_TETGENMESH_H

//#include "TetGen.h"
#include <ctime>
#include <vtkPolyData.h>
#include "predicates.h"
#include "TetGenBehavior.h"
#include "TetGenIO.h"

class TetGenMesh {

public:

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// Mesh data structure                                                       //
//                                                                           //
// A tetrahedral mesh T of a 3D piecewise linear complex (PLC) X is a 3D     //
// simplicial complex whose underlying space is equal to the space of X.  T  //
// contains a 2D subcomplex S which is a triangular mesh of the boundary of  //
// X. S contains a 1D subcomplex L which is a linear mesh of the boundary of //
// S. Faces and edges in S and L are respectively called subfaces and segme- //
// nts to distinguish them from others in T.                                 //
//                                                                           //
// TetGen stores the tetrahedra and vertices of T. The basic structure of a  //
// tetrahedron contains pointers to its vertices and adjacent tetrahedra. A  //
// vertex stores its x-, y-, and z-coordinates, and a pointer to a tetrahed- //
// ron containing it. Both tetrahedra and vertices may contain user data.    // 
//                                                                           //
// Each face of T belongs to either two tetrahedra or one tetrahedron. In    //
// the latter case, the face is an exterior boundary face of T.  TetGen adds //
// fictitious tetrahedra (one-to-one) at such faces, and connects them to an //
// "infinite vertex" (which has no geometric coordinates).  One can imagine  //
// such a vertex lies in 4D space and is visible by all exterior boundary    //
// faces.  The extended set of tetrahedra (including the infinite vertex) is //
// a tetrahedralization of a 3-pseudomanifold without boundary.  It has the  //
// property that every face is shared by exactly two tetrahedra.             // 
//                                                                           //
// The current version of TetGen stores explicitly the subfaces and segments //
// (which are in surface mesh S and the linear mesh L), respectively.  Extra //
// pointers are allocated in tetrahedra and subfaces to point each others.   //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

	// The tetrahedron data structure.  It includes the following fields:
	//   - a list of four adjoining tetrahedra;
	//   - a list of four vertices;
	//   - a pointer to a list of four subfaces (optional, for -p switch);
	//   - a pointer to a list of six segments  (optional, for -p switch);
	//   - a list of user-defined floating-point attributes (optional);
	//   - a volume constraint (optional, for -a switch);
	//   - an integer of element marker (and flags);
	// The structure of a tetrahedron is an array of pointers.  Its actual size
	//   (the length of the array) is determined at runtime.

	typedef REAL **tetrahedron;

	// The subface data structure.  It includes the following fields:
	//   - a list of three adjoining subfaces;
	//   - a list of three vertices;
	//   - a list of three adjoining segments;
	//   - two adjoining tetrahedra;
	//   - an area constraint (optional, for -q switch);
	//   - an integer for boundary marker;
	//   - an integer for type, flags, etc.

	typedef REAL **shellface;

	// The point data structure.  It includes the following fields:
	//   - x, y and z coordinates;
	//   - a list of user-defined point attributes (optional);
	//   - u, v coordinates (optional, for -s switch);
	//   - a metric tensor (optional, for -q or -m switch);
	//   - a pointer to an adjacent tetrahedron;
	//   - a pointer to a parent (or a duplicate) point;
	//   - a pointer to an adjacent subface or segment (optional, -p switch);
	//   - a pointer to a tet in background mesh (optional, for -m switch);
	//   - an integer for boundary marker (point index);
	//   - an integer for point type (and flags).
	//   - an integer for geometry tag (optional, for -s switch).
	// The structure of a point is an array of REALs.  Its acutal size is 
	//   determined at the runtime.

	typedef REAL *point;

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// Handles                                                                   //
//                                                                           //
// Navigation and manipulation in a tetrahedralization are accomplished by   //
// operating on structures referred as ``handles". A handle is a pair (t,v), //
// where t is a pointer to a tetrahedron, and v is a 4-bit integer, in the   //
// range from 0 to 11. v is called the ``version'' of a tetrahedron, it rep- //
// resents a directed edge of a specific face of the tetrahedron.            //
//                                                                           //
// There are 12 even permutations of the four vertices, each of them corres- //
// ponds to a directed edge (a version) of the tetrahedron.  The 12 versions //
// can be grouped into 4 distinct ``edge rings'' in 4 ``oriented faces'' of  //
// this tetrahedron.  One can encode each version (a directed edge) into a   //
// 4-bit integer such that the two upper bits encode the index (from 0 to 2) //
// of this edge in the edge ring, and the two lower bits encode the index (  //
// from 0 to 3) of the oriented face which contains this edge.               //  
//                                                                           //
// The four vertices of a tetrahedron are indexed from 0 to 3 (according to  //
// their storage in the data structure).  Give each face the same index as   //
// the node opposite it in the tetrahedron.  Denote the edge connecting face //
// i to face j as i/j. We number the twelve versions as follows:             //
//                                                                           //
//           |   edge 0     edge 1     edge 2                                //
//   --------|--------------------------------                               //
//    face 0 |   0 (0/1)    4 (0/3)    8 (0/2)                               //
//    face 1 |   1 (1/2)    5 (1/3)    9 (1/0)                               //
//    face 2 |   2 (2/3)    6 (2/1)   10 (2/0)                               //
//    face 3 |   3 (3/0)    7 (3/1)   11 (3/2)                               //
//                                                                           //
// Similarly, navigation and manipulation in a (boundary) triangulation are  //
// done by using handles of triangles. Each handle is a pair (s, v), where s //
// is a pointer to a triangle, and v is a version in the range from 0 to 5.  //
// Each version corresponds to a directed edge of this triangle.             //
//                                                                           //
// Number the three vertices of a triangle from 0 to 2 (according to their   //
// storage in the data structure). Give each edge the same index as the node //
// opposite it in the triangle. The six versions of a triangle are:          //
//                                                                           //
//                 | edge 0   edge 1   edge 2                                //
//  ---------------|--------------------------                               //
//   ccw orieation |   0        2        4                                   //
//    cw orieation |   1        3        5                                   //
//                                                                           //
// In the following, a 'triface' is a handle of tetrahedron, and a 'face' is //
// a handle of a triangle.                                                   //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

	class triface {
	public:
		tetrahedron *tet;
		int ver; // Range from 0 to 11.
		triface() : tet(0), ver(0) {}
		triface& operator=(const triface& t) {
			tet = t.tet; ver = t.ver;
			return *this;
		}
	};

	class face {
	public:
		shellface *sh;
		int shver; // Range from 0 to 5.
		face() : sh(0), shver(0) {}
		face& operator=(const face& s) {
			sh = s.sh; shver = s.shver;
			return *this;
		}
	};

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// Arraypool                                                                 //
//                                                                           //
// A dynamic linear array. (It is written by J. Shewchuk)                    //
//                                                                           //
// Each arraypool contains an array of pointers to a number of blocks.  Each //
// block contains the same fixed number of objects.  Each index of the array //
// addresses a particular object in the pool. The most significant bits add- //
// ress the index of the block containing the object. The less significant   //
// bits address this object within the block.                                //
//                                                                           //
// 'objectbytes' is the size of one object in blocks; 'log2objectsperblock'  //
// is the base-2 logarithm of 'objectsperblock'; 'objects' counts the number //
// of allocated objects; 'totalmemory' is the total memory in bytes.         //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

	class arraypool {

	public:

		int objectbytes;
		int objectsperblock;
		int log2objectsperblock;
		int objectsperblockmark;
		int toparraylen;
		char **toparray;
		long objects;
		unsigned long totalmemory;

		void restart();
		void poolinit(int sizeofobject, int log2objperblk);
		char* getblock(int objectindex);
		void* lookup(int objectindex);
		int newindex(void **newptr);

		arraypool(int sizeofobject, int log2objperblk);
		~arraypool();
	};

// fastlookup() -- A fast, unsafe operation. Return the pointer to the object
//   with a given index.  Note: The object's block must have been allocated,
//   i.e., by the function newindex().

#define fastlookup(pool, index) \
  (void *) ((pool)->toparray[(index) >> (pool)->log2objectsperblock] + \
            ((index) & (pool)->objectsperblockmark) * (pool)->objectbytes)

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// Memorypool                                                                //
//                                                                           //
// A structure for memory allocation. (It is written by J. Shewchuk)         //
//                                                                           //
// firstblock is the first block of items. nowblock is the block from which  //
//   items are currently being allocated. nextitem points to the next slab   //
//   of free memory for an item. deaditemstack is the head of a linked list  //
//   (stack) of deallocated items that can be recycled.  unallocateditems is //
//   the number of items that remain to be allocated from nowblock.          //
//                                                                           //
// Traversal is the process of walking through the entire list of items, and //
//   is separate from allocation.  Note that a traversal will visit items on //
//   the "deaditemstack" stack as well as live items.  pathblock points to   //
//   the block currently being traversed.  pathitem points to the next item  //
//   to be traversed.  pathitemsleft is the number of items that remain to   //
//   be traversed in pathblock.                                              //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

	class memorypool {

	public:

		void **firstblock, **nowblock;
		void *nextitem;
		void *deaditemstack;
		void **pathblock;
		void *pathitem;
		int  alignbytes;
		int  itembytes, itemwords;
		int  itemsperblock;
		long items, maxitems;
		int  unallocateditems;
		int  pathitemsleft;

		memorypool();
		memorypool(int, int, int, int);
		~memorypool();

		void poolinit(int, int, int, int);
		void restart();
		void *alloc();
		void dealloc(void*);
		void traversalinit();
		void *traverse();
	};

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// badface                                                                   //
//                                                                           //
// Despite of its name, a 'badface' can be used to represent one of the      //
// following objects:                                                        //
//   - a face of a tetrahedron which is (possibly) non-Delaunay;             //
//   - an encroached subsegment or subface;                                  //
//   - a bad-quality tetrahedron, i.e, has too large radius-edge ratio;      //
//   - a sliver, i.e., has good radius-edge ratio but nearly zero volume;    //
//   - a recently flipped face (saved for undoing the flip later).           //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

	class badface {
	public:
		triface tt;
		face ss;
		REAL key, cent[6];  // circumcenter or cos(dihedral angles) at 6 edges.
		point forg, fdest, fapex, foppo, noppo;
		badface *nextitem;
		badface() : key(0), forg(0), fdest(0), fapex(0), foppo(0), noppo(0),
					nextitem(0) {}
	};

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// insertvertexflags                                                         //
//                                                                           //
// A collection of flags that pass to the routine insertvertex().            //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

	class insertvertexflags {

	public:

		int iloc;  // input/output.
		int bowywat, lawson;
		int splitbdflag, validflag, respectbdflag;
		int rejflag, chkencflag, cdtflag;
		int assignmeshsize;
		int sloc, sbowywat;

		// Used by Delaunay refinement.
		int refineflag; // 0, 1, 2, 3
		triface refinetet;
		face refinesh;
		int smlenflag; // for useinsertradius.
		REAL smlen; // for useinsertradius.
		point parentpt;

		insertvertexflags() {
			iloc = bowywat = lawson = 0;
			splitbdflag = validflag = respectbdflag = 0;
			rejflag = chkencflag = cdtflag = 0;
			assignmeshsize = 0;
			sloc = sbowywat = 0;

			refineflag = 0;
			refinetet.tet = NULL;
			refinesh.sh = NULL;
			smlenflag = 0;
			smlen = 0.0;
		}
	};

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// flipconstraints                                                           //
//                                                                           //
// A structure of a collection of data (options and parameters) which pass   //
// to the edge flip function flipnm().                                       //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

	class flipconstraints {

	public:

		// Elementary flip flags.
		int enqflag; // (= flipflag)
		int chkencflag;

		// Control flags
		int unflip;  // Undo the performed flips.
		int collectnewtets; // Collect the new tets created by flips.
		int collectencsegflag;

		// Optimization flags.
		int remove_ndelaunay_edge; // Remove a non-Delaunay edge.
		REAL bak_tetprism_vol; // The value to be minimized.
		REAL tetprism_vol_sum;
		int remove_large_angle; // Remove a large dihedral angle at edge.
		REAL cosdihed_in; // The input cosine of the dihedral angle (> 0).
		REAL cosdihed_out; // The improved cosine of the dihedral angle.

		// Boundary recovery flags.
		int checkflipeligibility;
		point seg[2];  // A constraining edge to be recovered.
		point fac[3];  // A constraining face to be recovered.
		point remvert; // A vertex to be removed.


		flipconstraints() {
			enqflag = 0;
			chkencflag = 0;

			unflip = 0;
			collectnewtets = 0;
			collectencsegflag = 0;

			remove_ndelaunay_edge = 0;
			bak_tetprism_vol = 0.0;
			tetprism_vol_sum = 0.0;
			remove_large_angle = 0;
			cosdihed_in = 0.0;
			cosdihed_out = 0.0;

			checkflipeligibility = 0;
			seg[0] = NULL;
			fac[0] = NULL;
			remvert = NULL;
		}
	};

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// optparameters                                                             //
//                                                                           //
// Optimization options and parameters.                                      //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

	class optparameters {

	public:

		// The one of goals of optimization.
		int max_min_volume;      // Maximize the minimum volume.
		int max_min_aspectratio; // Maximize the minimum aspect ratio.
		int min_max_dihedangle;  // Minimize the maximum dihedral angle. 

		// The initial and improved value.
		REAL initval, imprval;

		int numofsearchdirs;
		REAL searchstep;
		int maxiter;  // Maximum smoothing iterations (disabled by -1).
		int smthiter; // Performed iterations.


		optparameters() {
			max_min_volume = 0;
			max_min_aspectratio = 0;
			min_max_dihedangle = 0;

			initval = imprval = 0.0;

			numofsearchdirs = 10;
			searchstep = 0.01;
			maxiter = -1;   // Unlimited smoothing iterations.
			smthiter = 0;

		}
	};


///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// Labels (enumeration declarations) used by TetGen.                         //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

	// Labels that signify the type of a vertex. 
	enum verttype {UNUSEDVERTEX, DUPLICATEDVERTEX, RIDGEVERTEX, ACUTEVERTEX,
		FACETVERTEX, VOLVERTEX, FREESEGVERTEX, FREEFACETVERTEX,
		FREEVOLVERTEX, NREGULARVERTEX, DEADVERTEX};

	// Labels that signify the result of triangle-triangle intersection test.
	enum interresult {DISJOINT, INTERSECT, SHAREVERT, SHAREEDGE, SHAREFACE,
		TOUCHEDGE, TOUCHFACE, ACROSSVERT, ACROSSEDGE, ACROSSFACE,
		COLLISIONFACE, ACROSSSEG, ACROSSSUB};

	// Labels that signify the result of point location.
	enum locateresult {UNKNOWN, OUTSIDE, INTETRAHEDRON, ONFACE, ONEDGE, ONVERTEX,
		ENCVERTEX, ENCSEGMENT, ENCSUBFACE, NEARVERTEX, NONREGULAR,
		INSTAR, BADELEMENT};

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// Variables of TetGen                                                       //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

	// Pointer to the input data (a set of nodes, a PLC, or a mesh).
	TetGenIO *in, *addin;

	// Pointer to the switches and parameters.
	TetGenBehavior *b;

	// Pointer to a background mesh (contains size specification map).
	TetGenMesh *bgm;

	// Memorypools to store mesh elements (points, tetrahedra, subfaces, and
	//   segments) and extra pointers between tetrahedra, subfaces, and segments.
	memorypool *tetrahedrons, *subfaces, *subsegs, *points;
	memorypool *tet2subpool, *tet2segpool;

	// Memorypools to store bad-quality (or encroached) elements.
	memorypool *badtetrahedrons, *badsubfacs, *badsubsegs;

	// A memorypool to store faces to be flipped.
	memorypool *flippool;
	arraypool *unflipqueue;
	badface *flipstack;

	// Arrays used for point insertion (the Bowyer-Watson algorithm).
	arraypool *cavetetlist, *cavebdrylist, *caveoldtetlist;
	arraypool *cavetetshlist, *cavetetseglist, *cavetetvertlist;
	arraypool *caveencshlist, *caveencseglist;
	arraypool *caveshlist, *caveshbdlist, *cavesegshlist;

	// Stacks used for CDT construction and boundary recovery.
	arraypool *subsegstack, *subfacstack, *subvertstack;

	// Arrays of encroached segments and subfaces (for mesh refinement).
	arraypool *encseglist, *encshlist;

	// The map between facets to their vertices (for mesh refinement).
	int *idx2facetlist;
	point *facetverticeslist;

	// The map between segments to their endpoints (for mesh refinement).
	point *segmentendpointslist;

	// The infinite vertex.
	point dummypoint;
	// The recently visited tetrahedron, subface.
	triface recenttet;
	face recentsh;

	// PI is the ratio of a circle's circumference to its diameter.
	static REAL PI;

	// Array (size = numberoftetrahedra * 6) for storing high-order nodes of
	//   tetrahedra (only used when -o2 switch is selected).
	point *highordertable;

	// Various variables.
	int numpointattrib;                          // Number of point attributes.
	int numelemattrib;                     // Number of tetrahedron attributes.
	int sizeoftensor;                     // Number of REALs per metric tensor.
	int pointmtrindex;           // Index to find the metric tensor of a point.
	int pointparamindex;       // Index to find the u,v coordinates of a point.
	int point2simindex;         // Index to find a simplex adjacent to a point.
	int pointmarkindex;            // Index to find boundary marker of a point.
	int elemattribindex;          // Index to find attributes of a tetrahedron.
	int volumeboundindex;       // Index to find volume bound of a tetrahedron.
	int elemmarkerindex;              // Index to find marker of a tetrahedron.
	int shmarkindex;             // Index to find boundary marker of a subface.
	int areaboundindex;               // Index to find area bound of a subface.
	int checksubsegflag;   // Are there segments in the tetrahedralization yet?
	int checksubfaceflag;  // Are there subfaces in the tetrahedralization yet?
	int checkconstraints;  // Are there variant (node, seg, facet) constraints?
	int nonconvex;                               // Is current mesh non-convex?
	int autofliplinklevel;        // The increase of link levels, default is 1.
	int useinsertradius;       // Save the insertion radius for Steiner points.
	long samples;               // Number of random samples for point location.
	unsigned long randomseed;                    // Current random number seed.
	REAL cosmaxdihed, cosmindihed;    // The cosine values of max/min dihedral.
	REAL cossmtdihed;     // The cosine value of a bad dihedral to be smoothed.
	REAL cosslidihed;      // The cosine value of the max dihedral of a sliver.
	REAL minfaceang, minfacetdihed;     // The minimum input (dihedral) angles.
	REAL tetprism_vol_sum;   // The total volume of tetrahedral-prisms (in 4D).
	REAL longest;                          // The longest possible edge length.
	REAL xmax, xmin, ymax, ymin, zmax, zmin;         // Bounding box of points.

	// Counters.
	long insegments;                               // Number of input segments.
	long hullsize;                        // Number of exterior boundary faces.
	long meshedges;                                    // Number of mesh edges.
	long meshhulledges;                       // Number of boundary mesh edges.
	long steinerleft;                 // Number of Steiner points not yet used.
	long dupverts;                            // Are there duplicated vertices?
	long unuverts;                                // Are there unused vertices?
	long nonregularcount;                    // Are there non-regular vertices?
	long st_segref_count, st_facref_count, st_volref_count;  // Steiner points.
	long fillregioncount, cavitycount, cavityexpcount;
	long flip14count, flip26count, flipn2ncount;
	long flip23count, flip32count, flip44count, flip41count;
	long flip31count, flip22count;
	unsigned long totalworkmemory;      // Total memory used by working arrays.


///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// Mesh manipulation primitives                                              //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

	// Fast lookup tables for mesh manipulation primitives.
	static int bondtbl[12][12], fsymtbl[12][12];
	static int esymtbl[12], enexttbl[12], eprevtbl[12];
	static int enextesymtbl[12], eprevesymtbl[12];
	static int eorgoppotbl[12], edestoppotbl[12];
	static int facepivot1[12], facepivot2[12][12];
	static int orgpivot[12], destpivot[12], apexpivot[12], oppopivot[12];
	static int tsbondtbl[12][6], stbondtbl[12][6];
	static int tspivottbl[12][6], stpivottbl[12][6];
	static int ver2edge[12], edge2ver[6], epivot[12];
	static int sorgpivot [6], sdestpivot[6], sapexpivot[6];
	static int snextpivot[6];

	void inittables();

	// Primitives for tetrahedra.
	inline tetrahedron encode(triface& t);
	inline tetrahedron encode2(tetrahedron* ptr, int ver);
	inline void decode(tetrahedron ptr, triface& t);
	inline void bond(triface& t1, triface& t2);
	inline void dissolve(triface& t);
	inline void esym(triface& t1, triface& t2);
	inline void esymself(triface& t);
	inline void enext(triface& t1, triface& t2);
	inline void enextself(triface& t);
	inline void eprev(triface& t1, triface& t2);
	inline void eprevself(triface& t);
	inline void enextesym(triface& t1, triface& t2);
	inline void enextesymself(triface& t);
	inline void eprevesym(triface& t1, triface& t2);
	inline void eprevesymself(triface& t);
	inline void eorgoppo(triface& t1, triface& t2);
	inline void eorgoppoself(triface& t);
	inline void edestoppo(triface& t1, triface& t2);
	inline void edestoppoself(triface& t);
	inline void fsym(triface& t1, triface& t2);
	inline void fsymself(triface& t);
	inline void fnext(triface& t1, triface& t2);
	inline void fnextself(triface& t);
	inline point org (triface& t);
	inline point dest(triface& t);
	inline point apex(triface& t);
	inline point oppo(triface& t);
	inline void setorg (triface& t, point p);
	inline void setdest(triface& t, point p);
	inline void setapex(triface& t, point p);
	inline void setoppo(triface& t, point p);
	inline REAL elemattribute(tetrahedron* ptr, int attnum);
	inline void setelemattribute(tetrahedron* ptr, int attnum, REAL value);
	inline REAL volumebound(tetrahedron* ptr);
	inline void setvolumebound(tetrahedron* ptr, REAL value);
	inline int  elemindex(tetrahedron* ptr);
	inline void setelemindex(tetrahedron* ptr, int value);
	inline int  elemmarker(tetrahedron* ptr);
	inline void setelemmarker(tetrahedron* ptr, int value);
	inline void infect(triface& t);
	inline void uninfect(triface& t);
	inline bool infected(triface& t);
	inline void marktest(triface& t);
	inline void unmarktest(triface& t);
	inline bool marktested(triface& t);
	inline void markface(triface& t);
	inline void unmarkface(triface& t);
	inline bool facemarked(triface& t);
	inline void markedge(triface& t);
	inline void unmarkedge(triface& t);
	inline bool edgemarked(triface& t);
	inline void marktest2(triface& t);
	inline void unmarktest2(triface& t);
	inline bool marktest2ed(triface& t);
	inline int  elemcounter(triface& t);
	inline void setelemcounter(triface& t, int value);
	inline void increaseelemcounter(triface& t);
	inline void decreaseelemcounter(triface& t);
	inline bool ishulltet(triface& t);
	inline bool isdeadtet(triface& t);

	// Primitives for subfaces and subsegments.
	inline void sdecode(shellface sptr, face& s);
	inline shellface sencode(face& s);
	inline shellface sencode2(shellface *sh, int shver);
	inline void spivot(face& s1, face& s2);
	inline void spivotself(face& s);
	inline void sbond(face& s1, face& s2);
	inline void sbond1(face& s1, face& s2);
	inline void sdissolve(face& s);
	inline point sorg(face& s);
	inline point sdest(face& s);
	inline point sapex(face& s);
	inline void setsorg(face& s, point pointptr);
	inline void setsdest(face& s, point pointptr);
	inline void setsapex(face& s, point pointptr);
	inline void sesym(face& s1, face& s2);
	inline void sesymself(face& s);
	inline void senext(face& s1, face& s2);
	inline void senextself(face& s);
	inline void senext2(face& s1, face& s2);
	inline void senext2self(face& s);
	inline REAL areabound(face& s);
	inline void setareabound(face& s, REAL value);
	inline int shellmark(face& s);
	inline void setshellmark(face& s, int value);
	inline void sinfect(face& s);
	inline void suninfect(face& s);
	inline bool sinfected(face& s);
	inline void smarktest(face& s);
	inline void sunmarktest(face& s);
	inline bool smarktested(face& s);
	inline void smarktest2(face& s);
	inline void sunmarktest2(face& s);
	inline bool smarktest2ed(face& s);
	inline void smarktest3(face& s);
	inline void sunmarktest3(face& s);
	inline bool smarktest3ed(face& s);
	inline void setfacetindex(face& f, int value);
	inline int  getfacetindex(face& f);

	// Primitives for interacting tetrahedra and subfaces.
	inline void tsbond(triface& t, face& s);
	inline void tsdissolve(triface& t);
	inline void stdissolve(face& s);
	inline void tspivot(triface& t, face& s);
	inline void stpivot(face& s, triface& t);

	// Primitives for interacting tetrahedra and segments.
	inline void tssbond1(triface& t, face& seg);
	inline void sstbond1(face& s, triface& t);
	inline void tssdissolve1(triface& t);
	inline void sstdissolve1(face& s);
	inline void tsspivot1(triface& t, face& s);
	inline void sstpivot1(face& s, triface& t);

	// Primitives for interacting subfaces and segments.
	inline void ssbond(face& s, face& edge);
	inline void ssbond1(face& s, face& edge);
	inline void ssdissolve(face& s);
	inline void sspivot(face& s, face& edge);

	// Primitives for points.
	inline int  pointmark(point pt);
	inline void setpointmark(point pt, int value);
	inline enum verttype pointtype(point pt);
	inline void setpointtype(point pt, enum verttype value);
	inline int  pointgeomtag(point pt);
	inline void setpointgeomtag(point pt, int value);
	inline REAL pointgeomuv(point pt, int i);
	inline void setpointgeomuv(point pt, int i, REAL value);
	inline void pinfect(point pt);
	inline void puninfect(point pt);
	inline bool pinfected(point pt);
	inline void pmarktest(point pt);
	inline void punmarktest(point pt);
	inline bool pmarktested(point pt);
	inline void pmarktest2(point pt);
	inline void punmarktest2(point pt);
	inline bool pmarktest2ed(point pt);
	inline void pmarktest3(point pt);
	inline void punmarktest3(point pt);
	inline bool pmarktest3ed(point pt);
	inline tetrahedron point2tet(point pt);
	inline void setpoint2tet(point pt, tetrahedron value);
	inline shellface point2sh(point pt);
	inline void setpoint2sh(point pt, shellface value);
	inline point point2ppt(point pt);
	inline void setpoint2ppt(point pt, point value);
	inline tetrahedron point2bgmtet(point pt);
	inline void setpoint2bgmtet(point pt, tetrahedron value);
	inline void setpointinsradius(point pt, REAL value);
	inline REAL getpointinsradius(point pt);

	// Advanced primitives.
	inline void point2tetorg(point pt, triface& t);
	inline void point2shorg(point pa, face& s);
	inline point farsorg(face& seg);
	inline point farsdest(face& seg);

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  Memory managment                                                         //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

	void tetrahedrondealloc(tetrahedron*);
	tetrahedron *tetrahedrontraverse();
	tetrahedron *alltetrahedrontraverse();
	void shellfacedealloc(memorypool*, shellface*);
	shellface *shellfacetraverse(memorypool*);
	void pointdealloc(point);
	point pointtraverse();

	void makeindex2pointmap(point*&);
	void makepoint2submap(memorypool*, int*&, face*&);
	void maketetrahedron(triface*);
	void makeshellface(memorypool*, face*);
	void makepoint(point*, enum verttype);

	void initializepools();

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// Advanced geometric predicates and calculations                            //
//                                                                           //
// TetGen uses a simplified symbolic perturbation scheme from Edelsbrunner,  //
// et al [*].  Hence the point-in-sphere test never returns a zero. The idea //
// is to perturb the weights of vertices in the fourth dimension.  TetGen    //
// uses the indices of the vertices decide the amount of perturbation. It is //
// implemented in the routine insphere_s().
//                                                                           //
// The routine tri_edge_test() determines whether or not a triangle and an   //
// edge intersect in 3D. If they intersect, their intersection type is also  //
// reported. This test is a combination of n 3D orientation tests (n is bet- //
// ween 3 and 9). It uses the robust orient3d() test to make the branch dec- //
// isions.  The routine tri_tri_test() determines whether or not two triang- //
// les intersect in 3D. It also uses the robust orient3d() test.             //
//                                                                           //
// There are a number of routines to calculate geometrical quantities, e.g., //
// circumcenters, angles, dihedral angles, face normals, face areas, etc.    //
// They are so far done by the default floating-point arithmetics which are  //
// non-robust. They should be improved in the future.                        //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

	// Symbolic perturbations (robust)
	REAL insphere_s(REAL*, REAL*, REAL*, REAL*, REAL*);
	REAL orient4d_s(REAL*, REAL*, REAL*, REAL*, REAL*,
					REAL, REAL, REAL, REAL, REAL);

	// Triangle-edge intersection test (robust)
	int tri_edge_2d(point, point, point, point, point, point, int, int*, int*);
	int tri_edge_tail(point, point, point, point, point, point, REAL, REAL, int,
					  int*, int*);
	int tri_edge_test(point, point, point, point, point, point, int, int*, int*);

	// Triangle-triangle intersection test (robust)
	int tri_edge_inter_tail(point, point, point, point, point, REAL, REAL);
	int tri_tri_inter(point, point, point, point, point, point);

	// Linear algebra functions
	inline REAL dot(REAL* v1, REAL* v2);
	inline void cross(REAL* v1, REAL* v2, REAL* n);
	bool lu_decmp(REAL lu[4][4], int n, int* ps, REAL* d, int N);
	void lu_solve(REAL lu[4][4], int n, int* ps, REAL* b, int N);

	// An embedded 2-dimensional geometric predicate (non-robust)
	REAL incircle3d(point pa, point pb, point pc, point pd);

	// Geometric calculations (non-robust)
	REAL orient3dfast(REAL *pa, REAL *pb, REAL *pc, REAL *pd);
	inline REAL norm2(REAL x, REAL y, REAL z);
	inline REAL distance(REAL* p1, REAL* p2);
	void facenormal(point pa, point pb, point pc, REAL *n, int pivot, REAL *lav);
	REAL shortdistance(REAL* p, REAL* e1, REAL* e2);
	REAL triarea(REAL* pa, REAL* pb, REAL* pc);
	REAL interiorangle(REAL* o, REAL* p1, REAL* p2, REAL* n);
	void projpt2edge(REAL* p, REAL* e1, REAL* e2, REAL* prj);
	void projpt2face(REAL* p, REAL* f1, REAL* f2, REAL* f3, REAL* prj);
	REAL facedihedral(REAL* pa, REAL* pb, REAL* pc1, REAL* pc2);
	bool tetalldihedral(point, point, point, point, REAL*, REAL*, REAL*);
	void tetallnormal(point, point, point, point, REAL N[4][3], REAL* volume);
	REAL tetaspectratio(point, point, point, point);
	bool circumsphere(REAL*, REAL*, REAL*, REAL*, REAL* cent, REAL* radius);
	bool orthosphere(REAL*,REAL*,REAL*,REAL*,REAL,REAL,REAL,REAL,REAL*,REAL*);
	void planelineint(REAL*, REAL*, REAL*, REAL*, REAL*, REAL*, REAL*);
	int linelineint(REAL*, REAL*, REAL*, REAL*, REAL*, REAL*, REAL*, REAL*);
	REAL tetprismvol(REAL* pa, REAL* pb, REAL* pc, REAL* pd);
	bool calculateabovepoint(arraypool*, point*, point*, point*);
	void calculateabovepoint4(point, point, point, point);

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// Local mesh transformations                                                //
//                                                                           //
// A local transformation replaces a small set of tetrahedra with another    //
// set of tetrahedra which fills the same space and the same boundaries.     //
//   In 3D, the most simplest local transformations are the elementary flips //
// performed within the convex hull of five vertices: 2-to-3, 3-to-2, 1-to-4,//
// and 4-to-1 flips,  where the numbers indicate the number of tetrahedra    //
// before and after each flip.  The 1-to-4 and 4-to-1 flip involve inserting //
// or deleting a vertex, respectively.                                       //
//   There are complex local transformations which can be decomposed as a    //
// combination of elementary flips. For example,a 4-to-4 flip which replaces //
// two coplanar edges can be regarded by a 2-to-3 flip and a 3-to-2 flip.    //
// Note that the first 2-to-3 flip will temporarily create a degenerate tet- //
// rahedron which is removed immediately by the followed 3-to-2 flip.  More  //
// generally, a n-to-m flip, where n > 3, m = (n - 2) * 2, which removes an  //
// edge can be done by first performing a sequence of (n - 3) 2-to-3 flips   //
// followed by a 3-to-2 flip.                                                //
//                                                                           //
// The routines flip23(), flip32(), and flip41() perform the three element-  //
// ray flips. The flip14() is available inside the routine insertpoint().    //
//                                                                           //
// The routines flipnm() and flipnm_post() implement a generalized edge flip //
// algorithm which uses a combination of elementary flips.                   //
//                                                                           //
// The routine insertpoint() implements a variant of Bowyer-Watson's cavity  //
// algorithm to insert a vertex. It works for arbitrary tetrahedralization,  //
// either Delaunay, or constrained Delaunay, or non-Delaunay.                //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

	// The elementary flips.
	void flip23(triface*, int, flipconstraints* fc);
	void flip32(triface*, int, flipconstraints* fc);
	void flip41(triface*, int, flipconstraints* fc);

	// A generalized edge flip.
	int flipnm(triface*, int n, int level, int, flipconstraints* fc);
	int flipnm_post(triface*, int n, int nn, int, flipconstraints* fc);

	// Point insertion.
	int  insertpoint(point, triface*, face*, face*, insertvertexflags*);
	void insertpoint_abort(face*, insertvertexflags*);

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// Delaunay tetrahedralization                                               //
//                                                                           //
// The routine incrementaldelaunay() implemented two incremental algorithms  //
// for constructing Delaunay tetrahedralizations (DTs):  the Bowyer-Watson   //
// (B-W) algorithm and the incremental flip algorithm of Edelsbrunner and    //
// Shah, "Incremental topological flipping works for regular triangulation," //
// Algorithmica, 15:233-241, 1996.                                           //
//                                                                           //
// The routine incrementalflip() implements the flip algorithm of [Edelsbru- //
// nner and Shah, 1996].  It flips a queue of locally non-Delaunay faces (in //
// an arbitrary order).  The success is guaranteed when the Delaunay tetrah- //
// edralization is constructed incrementally by adding one vertex at a time. //
//                                                                           //
// The routine locate() finds a tetrahedron contains a new point in current  //
// DT.  It uses a simple stochastic walk algorithm: starting from an arbitr- //
// ary tetrahedron in DT, it finds the destination by visit one tetrahedron  //
// at a time, randomly chooses a tetrahedron if there are more than one      //
// choices. This algorithm terminates due to Edelsbrunner's acyclic theorem. //
//   Choose a good starting tetrahedron is crucial to the speed of the walk. //
// TetGen originally uses the "jump-and-walk" algorithm of Muecke, E.P.,     //
// Saias, I., and Zhu, B. "Fast Randomized Point Location Without Preproces- //
// sing." In Proceedings of the 12th ACM Symposium on Computational Geometry,//
// 274-283, 1996.  It first randomly samples several tetrahedra in the DT    //
// and then choosing the closet one to start walking.                        //
//   The above algorithm slows download dramatically as the number of points //
// grows -- reported in Amenta, N., Choi, S. and Rote, G., "Incremental      //
// construction con {BRIO}," In Proceedings of 19th ACM Symposium on         //
// Computational Geometry, 211-219, 2003.  On the other hand, Liu and        //
// Snoeyink showed that the point location can be made in constant time if   //
// the points are pre-sorted so that the nearby points in space have nearby  //
// indices, then adding the points in this order. They sorted the points     //
// along the 3D Hilbert curve.                                               //
//                                                                           //
// The routine hilbert_sort3() sorts a set of 3D points along the 3D Hilbert //
// curve. It recursively splits a point set according to the Hilbert indices //
// mapped to the subboxes of the bounding box of the point set.              //
//   The Hilbert indices is calculated by Butz's algorithm in 1971.  A nice  //
// exposition of this algorithm can be found in the paper of Hamilton, C.,   //
// "Compact Hilbert Indices", Technical Report CS-2006-07, Computer Science, //
// Dalhousie University, 2006 (the Section 2). My implementation also refer- //
// enced Steven Witham's implementation of "Hilbert walk" (hopefully, it is  //
// still available at: http://www.tiac.net/~sw/2008/10/Hilbert/).            //
//                                                                           //
// TetGen sorts the points using the method in the paper of Boissonnat,J.-D.,//
// Devillers, O. and Hornus, S. "Incremental Construction of the Delaunay    //
// Triangulation and the Delaunay Graph in Medium Dimension," In Proceedings //
// of the 25th ACM Symposium on Computational Geometry, 2009.                //
//   It first randomly sorts the points into subgroups using the Biased Rand-//
// omized Insertion Ordering (BRIO) of Amenta et al 2003, then sorts the     //
// points in each subgroup along the 3D Hilbert curve.  Inserting points in  //
// this order ensures a randomized "sprinkling" of the points over the       //
// domain, while sorting of each subset ensures locality.                    //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

	void transfernodes();

	// Point sorting.
	int  transgc[8][3][8], tsb1mod3[8];
	void hilbert_init(int n);
	int  hilbert_split(point* vertexarray, int arraysize, int gc0, int gc1,
					   REAL, REAL, REAL, REAL, REAL, REAL);
	void hilbert_sort3(point* vertexarray, int arraysize, int e, int d,
					   REAL, REAL, REAL, REAL, REAL, REAL, int depth);
	void brio_multiscale_sort(point*,int,int threshold,REAL ratio,int* depth);

	// Point location.
	unsigned long randomnation(unsigned int choices);
	void randomsample(point searchpt, triface *searchtet);
	enum locateresult locate(point searchpt, triface *searchtet);

	// Incremental flips.
	void flippush(badface*&, triface*);
	int  incrementalflip(point newpt, int, flipconstraints *fc);

	// Incremental Delaunay construction.
	void initialdelaunay(point pa, point pb, point pc, point pd);
	void incrementaldelaunay(clock_t&);

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// Surface triangulation                                                     //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

	void flipshpush(face*);
	void flip22(face*, int, int);
	void flip31(face*, int);
	long lawsonflip();
	int sinsertvertex(point newpt, face*, face*, int iloc, int bowywat, int);
	int sremovevertex(point delpt, face*, face*, int lawson);

	enum locateresult slocate(point, face*, int, int, int);
	enum interresult sscoutsegment(face*, point);
	void scarveholes(int, REAL*);
	void triangulate(int, arraypool*, arraypool*, int, REAL*);

	void unifysubfaces(face*, face*);
	void unifysegments();
	void mergefacets();
	void identifypscedges(point*);
	void meshsurface();

	void interecursive(shellface** subfacearray, int arraysize, int axis,
					   REAL, REAL, REAL, REAL, REAL, REAL, int* internum);
	void detectinterfaces();

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// Constrained Delaunay tetrahedralization                                   //
//                                                                           //
// A constrained Delaunay tetrahedralization (CDT) is a variation of a Dela- //
// unay tetrahedralization (DT) that is constrained to respect the boundary  //
// of a 3D PLC (domain). In a CDT of a 3D PLC, every vertex or edge of the   //
// PLC is also a vertex or an edge of the CDT, every polygon of the PLC is a //
// union of triangles of the CDT. A crucial difference between a CDT and a   //
// DT is that triangles in the PLC's polygons are not required to be locally //
// Delaunay, which frees the CDT to better respect the PLC's polygons. CDTs  //
// have optimal properties similar to those of DTs.                          //
//                                                                           //
// Steiner Points and Steiner CDTs. It is known that even a simple 3D polyh- //
// edron may not have a tetrahedralization which only uses its own vertices. //
// Some extra points, so-called "Steiner points" are needed in order to form //
// a tetrahedralization of such polyhedron.  It is true for tetrahedralizing //
// a 3D PLC as well. A Steiner CDT of a 3D PLC is a CDT containing Steiner   //
// points. The CDT algorithms of TetGen in general create Steiner CDTs.      //
// Almost all of the Steiner points are added in the edges of the PLC. They  //
// guarantee the existence of a CDT of the modified PLC.                     //
//                                                                           //
// The routine constraineddelaunay() starts from a DT of the vertices of a   //
// PLC and creates a (Steiner) CDT of the PLC (including Steiner points). It //
// is constructed by two steps, (1) segment recovery and (2) facet (polygon) //
// recovery. Each step is accomplished by its own algorithm.                 //
//                                                                           //
// The routine delaunizesegments() implements the segment recovery algorithm //
// of Si, H. and Gaertner, K. "Meshing Piecewise Linear Complexes by Constr- //
// ained Delaunay Tetrahedralizations," In Proceedings of the 14th Internat- //
// ional Meshing Roundtable, 147--163, 2005.  It adds Steiner points into    //
// non-Delaunay segments until all subsegments appear together in a DT. The  //
// running time of this algorithm is proportional to the number of added     //
// Steiner points.                                                           //
//                                                                           //
// There are two incremental facet recovery algorithms: the cavity re-trian- //
// gulation algorithm of Si, H. and Gaertner, K. "3D Boundary Recovery by    //
// Constrained Delaunay Tetrahedralization," International Journal for Numer-//
// ical Methods in Engineering, 85:1341-1364, 2011, and the flip algorithm   //
// of Shewchuk, J. "Updating and Constructing Constrained Delaunay and       //
// Constrained Regular Triangulations by Flips." In Proceedings of the 19th  //
// ACM Symposium on Computational Geometry, 86-95, 2003.                     //
//                                                                           //
// It is guaranteed in theory, no Steiner point is needed in both algorithms //
// However, a facet with non-coplanar vertices might cause the  additions of //
// Steiner points. It is discussed in the paper of Si, H., and  Shewchuk, J.,//
// "Incrementally Constructing and Updating Constrained Delaunay             //
// Tetrahedralizations with Finite Precision Coordinates." In Proceedings of //
// the 21th International Meshing Roundtable, 2012.                          //
//                                                                           //
// Our implementation of the facet recovery algorithms recover a "missing    //
// region" at a time. Each missing region is a subset of connected interiors //
// of a polygon. The routine formcavity() creates the cavity of crossing     //
// tetrahedra of the missing region.                                         //
//                                                                           //
// The cavity re-triangulation algorithm is implemented by three subroutines,//
// delaunizecavity(), fillcavity(), and carvecavity(). Since it may fail due //
// to non-coplanar vertices, the subroutine restorecavity() is used to rest- //
// ore the original cavity.                                                  //
//                                                                           //
// The routine flipinsertfacet() implements the flip algorithm. The subrout- //
// ine flipcertify() is used to maintain the priority queue of flips.        // 
//                                                                           //
// The routine refineregion() is called when the facet recovery algorithm    //
// fail to recover a missing region. It inserts Steiner points to refine the //
// missing region. In order to avoid inserting Steiner points very close to  //
// existing segments.  The classical encroachment rules of the Delaunay      //
// refinement algorithm are used to choose the Steiner points.               //
//                                                                           //
// The routine constrainedfacets() does the facet recovery by using either   //
// the cavity re-triangulation algorithm (default) or the flip algorithm. It //
// results a CDT of the (modified) PLC (including Steiner points).           //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

	void makesegmentendpointsmap();

	enum interresult finddirection(triface* searchtet, point endpt);
	enum interresult scoutsegment(point, point, triface*, point*, arraypool*);
	int  getsteinerptonsegment(face* seg, point refpt, point steinpt);
	void delaunizesegments();

	enum interresult scoutsubface(face* searchsh, triface* searchtet);
	void formregion(face*, arraypool*, arraypool*, arraypool*);
	int  scoutcrossedge(triface& crosstet, arraypool*, arraypool*);
	bool formcavity(triface*, arraypool*, arraypool*, arraypool*, arraypool*,
					arraypool*, arraypool*);

	// Facet recovery by cavity re-triangulation [Si and Gaertner 2011].
	void delaunizecavity(arraypool*, arraypool*, arraypool*, arraypool*,
						 arraypool*, arraypool*);
	bool fillcavity(arraypool*, arraypool*, arraypool*, arraypool*,
					arraypool*, arraypool*, triface* crossedge);
	void carvecavity(arraypool*, arraypool*, arraypool*);
	void restorecavity(arraypool*, arraypool*, arraypool*, arraypool*);

	// Facet recovery by flips [Shewchuk 2003].
	void flipcertify(triface *chkface, badface **pqueue, point, point, point);
	void flipinsertfacet(arraypool*, arraypool*, arraypool*, arraypool*);

	bool fillregion(arraypool* missingshs, arraypool*, arraypool* newshs);

	int  insertpoint_cdt(point, triface*, face*, face*, insertvertexflags*,
						 arraypool*, arraypool*, arraypool*, arraypool*,
						 arraypool*, arraypool*);
	void refineregion(face&, arraypool*, arraypool*, arraypool*, arraypool*,
					  arraypool*, arraypool*);

	void constrainedfacets();

	void constraineddelaunay(clock_t&);

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// Constrained tetrahedralizations.                                          //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

	int checkflipeligibility(int fliptype, point, point, point, point, point,
							 int level, int edgepivot, flipconstraints* fc);

	int removeedgebyflips(triface*, flipconstraints*);
	int removefacebyflips(triface*, flipconstraints*);

	int recoveredgebyflips(point, point, triface*, int fullsearch);
	int add_steinerpt_in_schoenhardtpoly(triface*, int, int chkencflag);
	int add_steinerpt_in_segment(face*, int searchlevel);
	int addsteiner4recoversegment(face*, int);
	int recoversegments(arraypool*, int fullsearch, int steinerflag);

	int recoverfacebyflips(point, point, point, face*, triface*);
	int recoversubfaces(arraypool*, int steinerflag);

	int getvertexstar(int, point searchpt, arraypool*, arraypool*, arraypool*);
	int getedge(point, point, triface*);
	int reduceedgesatvertex(point startpt, arraypool* endptlist);
	int removevertexbyflips(point steinerpt);

	int suppressbdrysteinerpoint(point steinerpt);
	int suppresssteinerpoints();

	void recoverboundary(clock_t&);

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// Mesh reconstruction                                                       //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

	void carveholes();

	void reconstructmesh();

	int  scoutpoint(point, triface*, int randflag);
	REAL getpointmeshsize(point, triface*, int iloc);
	void interpolatemeshsize();

	void insertconstrainedpoints(point *insertarray, int arylen, int rejflag);
	void insertconstrainedpoints(TetGenIO *addio);

	void collectremovepoints(arraypool *remptlist);
	void meshcoarsening();

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// Mesh refinement                                                           //
//                                                                           //
// The purpose of mesh refinement is to obtain a tetrahedral mesh with well- //
// -shaped tetrahedra and appropriate mesh size.  It is necessary to insert  //
// new Steiner points to achieve this property. The questions are (1) how to //
// choose the Steiner points? and (2) how to insert them?                    //
//                                                                           //
// Delaunay refinement is a technique first developed by Chew [1989] and     //
// Ruppert [1993, 1995] to generate quality triangular meshes in the plane.  //
// It provides guarantee on the smallest angle of the triangles.  Rupper's   //
// algorithm guarantees that the mesh is size-optimal (to within a constant  //
// factor) among all meshes with the same quality.                           //
//   Shewchuk generalized Ruppert's algorithm into 3D in his PhD thesis      //
// [Shewchuk 1997]. A short version of his algorithm appears in "Tetrahedral //
// Mesh Generation by Delaunay Refinement," In Proceedings of the 14th ACM   //
// Symposium on Computational Geometry, 86-95, 1998.  It guarantees that all //
// tetrahedra of the output mesh have a "radius-edge ratio" (equivalent to   //
// the minimal face angle) bounded. However, it does not remove slivers, a   //
// type of very flat tetrahedra which can have no small face angles but have //
// very small (and large) dihedral angles. Moreover, it may not terminate if //
// the input PLC contains "sharp features", e.g., two edges (or two facets)  //
// meet at an acute angle (or dihedral angle).                               //
//                                                                           //
// TetGen uses the basic Delaunay refinement scheme to insert Steiner points.//
// While it always maintains a constrained Delaunay mesh.  The algorithm is  //
// described in Si, H., "Adaptive Constrained Delaunay Mesh Generation,"     //
// International Journal for Numerical Methods in Engineering, 75:856-880.   //
// This algorithm always terminates and sharp features are easily preserved. //
// The mesh has good quality (same as Shewchuk's Delaunay refinement algori- //
// thm) in the bulk of the mesh domain. Moreover, it supports the generation //
// of adaptive mesh according to a (isotropic) mesh sizing function.         //   
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

	void makefacetverticesmap();
	int segsegadjacent(face *, face *);
	int segfacetadjacent(face *checkseg, face *checksh);
	int facetfacetadjacent(face *, face *);

	int checkseg4encroach(point pa, point pb, point checkpt);
	int checkseg4split(face *chkseg, point&, int&);
	int splitsegment(face *splitseg, point encpt, REAL, point, point, int, int);
	void repairencsegs(int chkencflag);

	void enqueuesubface(memorypool*, face*);
	int checkfac4encroach(point, point, point, point checkpt, REAL*, REAL*);
	int checkfac4split(face *chkfac, point& encpt, int& qflag, REAL *ccent);
	int splitsubface(face *splitfac, point, point, int qflag, REAL *ccent, int);
	void repairencfacs(int chkencflag);

	void enqueuetetrahedron(triface*);
	int checktet4split(triface *chktet, int& qflag, REAL *ccent);
	int splittetrahedron(triface* splittet,int qflag,REAL *ccent, int);
	void repairbadtets(int chkencflag);

	void delaunayrefinement();

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// Mesh optimization                                                         //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

	long lawsonflip3d(flipconstraints *fc);
	void recoverdelaunay();

	int  gettetrahedron(point, point, point, point, triface *);
	long improvequalitybyflips();

	int  smoothpoint(point smtpt, arraypool*, int ccw, optparameters *opm);
	long improvequalitybysmoothing(optparameters *opm);

	int  splitsliver(triface *, REAL, int);
	long removeslivers(int);

	void optimizemesh();

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// Mesh check and statistics                                                 //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

	// Mesh validations.
	int checkmesh(int topoflag);
	int checkshells();
	int checksegments();
	int checkdelaunay();
	int checkregular(int);
	int checkconforming(int);

	//  Mesh statistics.
	void printfcomma(unsigned long n);
	void qualitystatistics();
	void memorystatistics();
	void statistics();

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// Mesh output                                                               //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

	void jettisonnodes();
	void highorder();
	void numberedges();
	void outnodes(TetGenIO*);
	void outmetrics(TetGenIO*);
	void outelements(TetGenIO*);
	void outfaces(TetGenIO*);
	void outhullfaces(TetGenIO*);
	void outsubfaces(TetGenIO*);
	void outedges(TetGenIO*);
	void outsubsegments(TetGenIO*);
	void outneighbors(TetGenIO*);
	void outvoronoi(TetGenIO*);
	void outsmesh(char*);
	void outmesh2medit(char*);
	void outmesh2vtk(char*);



///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// Constructor & destructor                                                  //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

	TetGenMesh()
	{
		in  = addin = NULL;
		b   = NULL;
		bgm = NULL;

		tetrahedrons = subfaces = subsegs = points = NULL;
		badtetrahedrons = badsubfacs = badsubsegs = NULL;
		tet2segpool = tet2subpool = NULL;
		flippool = NULL;

		dummypoint = NULL;
		flipstack = NULL;
		unflipqueue = NULL;

		cavetetlist = cavebdrylist = caveoldtetlist = NULL;
		cavetetshlist = cavetetseglist = cavetetvertlist = NULL;
		caveencshlist = caveencseglist = NULL;
		caveshlist = caveshbdlist = cavesegshlist = NULL;

		subsegstack = subfacstack = subvertstack = NULL;
		encseglist = encshlist = NULL;
		idx2facetlist = NULL;
		facetverticeslist = NULL;
		segmentendpointslist = NULL;

		highordertable = NULL;

		numpointattrib = numelemattrib = 0;
		sizeoftensor = 0;
		pointmtrindex = 0;
		pointparamindex = 0;
		pointmarkindex = 0;
		point2simindex = 0;
		elemattribindex = 0;
		volumeboundindex = 0;
		shmarkindex = 0;
		areaboundindex = 0;
		checksubsegflag = 0;
		checksubfaceflag = 0;
		checkconstraints = 0;
		nonconvex = 0;
		autofliplinklevel = 1;
		useinsertradius = 0;
		samples = 0l;
		randomseed = 1l;
		minfaceang = minfacetdihed = PI;
		tetprism_vol_sum = 0.0;
		longest = 0.0;
		xmax = xmin = ymax = ymin = zmax = zmin = 0.0;

		insegments = 0l;
		hullsize = 0l;
		meshedges = meshhulledges = 0l;
		steinerleft = -1;
		dupverts = 0l;
		unuverts = 0l;
		nonregularcount = 0l;
		st_segref_count = st_facref_count = st_volref_count = 0l;
		fillregioncount = cavitycount = cavityexpcount = 0l;
		flip14count = flip26count = flipn2ncount = 0l;
		flip23count = flip32count = flip44count = flip41count = 0l;
		flip22count = flip31count = 0l;
		totalworkmemory = 0l;


	} // TetGenMesh()

	void freememory()
	{
		if (bgm != NULL) {
			delete bgm;
		}

		if (points != (memorypool *) NULL) {
			delete points;
			delete [] dummypoint;
		}

		if (tetrahedrons != (memorypool *) NULL) {
			delete tetrahedrons;
		}

		if (subfaces != (memorypool *) NULL) {
			delete subfaces;
			delete subsegs;
		}

		if (tet2segpool != NULL) {
			delete tet2segpool;
			delete tet2subpool;
		}

		if (flippool != NULL) {
			delete flippool;
			delete unflipqueue;
		}

		if (cavetetlist != NULL) {
			delete cavetetlist;
			delete cavebdrylist;
			delete caveoldtetlist;
			delete cavetetvertlist;
		}

		if (caveshlist != NULL) {
			delete caveshlist;
			delete caveshbdlist;
			delete cavesegshlist;
			delete cavetetshlist;
			delete cavetetseglist;
			delete caveencshlist;
			delete caveencseglist;
		}

		if (subsegstack != NULL) {
			delete subsegstack;
			delete subfacstack;
			delete subvertstack;
		}

		if (idx2facetlist != NULL) {
			delete [] idx2facetlist;
			delete [] facetverticeslist;
		}

		if (segmentendpointslist != NULL) {
			delete [] segmentendpointslist;
		}

		if (highordertable != NULL) {
			delete [] highordertable;
		}
	}

	~TetGenMesh()
	{
		freememory();
	} // ~TetGenMesh()

	void vtkToTetGenMesh(vtkPoints *iPoints, vtkCellArray *iPolys);
};                                               // End of class TetGenMesh.

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// tetrahedralize()    Interface for using TetGen's library to generate      //
//                     Delaunay tetrahedralizations, constrained Delaunay    //
//                     tetrahedralizations, quality tetrahedral meshes.      //
//                                                                           //
// 'in' is an object of 'TetGenIO' which contains a PLC you want to tetrahed-//
// ralize or a previously generated tetrahedral mesh you want to refine.  It //
// must not be a NULL. 'out' is another object of 'TetGenIO' for storing the //
// generated tetrahedral mesh. It can be a NULL. If so, the output will be   //
// saved to file(s). If 'bgmin' != NULL, it contains a background mesh which //
// defines a mesh size function.                                             //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

void tetrahedralize(TetGenBehavior *b, TetGenIO *in, TetGenIO *out,
					TetGenIO *addin = NULL, TetGenIO *bgmin = NULL);

#ifdef TETLIBRARY
void tetrahedralize(char *switches, TetGenIO *in, TetGenIO *out,
                    TetGenIO *addin = NULL, TetGenIO *bgmin = NULL);
#endif // #ifdef TETLIBRARY

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// terminatetetgen()    Terminate TetGen with a given exit code.             //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

inline void terminatetetgen(TetGenMesh *m, int x)
{
	// Release the allocated memory.
	if (m) {
		m->freememory();
	}
#ifdef TETLIBRARY
	throw x;
#else
	switch (x) {
		case 1: // Out of memory.
			printf("Error:  Out of memory.\n");
			break;
		case 2: // Encounter an internal error.
			printf("Please report this bug to Hang.Si@wias-berlin.de. Include\n");
			printf("  the message above, your input data set, and the exact\n");
			printf("  command line you used to run this program, thank you.\n");
			break;
		case 3:
			printf("A self-intersection was detected. Program stopped.\n");
			printf("Hint: use -d option to detect all self-intersections.\n");
			break;
		case 4:
			printf("A very small input feature size was detected. Program stopped.\n");
			printf("Hint: use -T option to set a smaller tolerance.\n");
			break;
		case 5:
			printf("Two very close input facets were detected. Program stopped.\n");
			printf("Hint: use -Y option to avoid adding Steiner points in boundary.\n");
			break;
		case 10:
			printf("An input error was detected. Program stopped.\n");
			break;
	} // switch (x)
	exit(x);
#endif // #ifdef TETLIBRARY
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// Primitives for tetrahedra                                                 //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

// encode()  compress a handle into a single pointer.  It relies on the 
//   assumption that all addresses of tetrahedra are aligned to sixteen-
//   byte boundaries, so that the last four significant bits are zero.

inline TetGenMesh::tetrahedron TetGenMesh::encode(triface& t) {
	return (tetrahedron) ((uintptr_t) (t).tet | (uintptr_t) (t).ver);
}

inline TetGenMesh::tetrahedron TetGenMesh::encode2(tetrahedron* ptr, int ver) {
	return (tetrahedron) ((uintptr_t) (ptr) | (uintptr_t) (ver));
}

// decode()  converts a pointer to a handle. The version is extracted from
//   the four least significant bits of the pointer.

inline void TetGenMesh::decode(tetrahedron ptr, triface& t) {
	(t).ver = (int) ((uintptr_t) (ptr) & (uintptr_t) 15);
	(t).tet = (tetrahedron *) ((uintptr_t) (ptr) ^ (uintptr_t) (t).ver);
}

// bond()  connects two tetrahedra together. (t1,v1) and (t2,v2) must 
//   refer to the same face and the same edge. 

inline void TetGenMesh::bond(triface& t1, triface& t2) {
	t1.tet[t1.ver & 3] = encode2(t2.tet, bondtbl[t1.ver][t2.ver]);
	t2.tet[t2.ver & 3] = encode2(t1.tet, bondtbl[t2.ver][t1.ver]);
}


// dissolve()  a bond (from one side).

inline void TetGenMesh::dissolve(triface& t) {
	t.tet[t.ver & 3] = NULL;
}

// enext()  finds the next edge (counterclockwise) in the same face.

inline void TetGenMesh::enext(triface& t1, triface& t2) {
	t2.tet = t1.tet;
	t2.ver = enexttbl[t1.ver];
}

inline void TetGenMesh::enextself(triface& t) {
	t.ver = enexttbl[t.ver];
}

// eprev()   finds the next edge (clockwise) in the same face.

inline void TetGenMesh::eprev(triface& t1, triface& t2) {
	t2.tet = t1.tet;
	t2.ver = eprevtbl[t1.ver];
}

inline void TetGenMesh::eprevself(triface& t) {
	t.ver = eprevtbl[t.ver];
}

// esym()  finds the reversed edge.  It is in the other face of the
//   same tetrahedron.

inline void TetGenMesh::esym(triface& t1, triface& t2) {
	(t2).tet = (t1).tet;
	(t2).ver = esymtbl[(t1).ver];
}

inline void TetGenMesh::esymself(triface& t) {
	(t).ver = esymtbl[(t).ver];
}

// enextesym()  finds the reversed edge of the next edge. It is in the other
//   face of the same tetrahedron. It is the combination esym() * enext(). 

inline void TetGenMesh::enextesym(triface& t1, triface& t2) {
	t2.tet = t1.tet;
	t2.ver = enextesymtbl[t1.ver];
}

inline void TetGenMesh::enextesymself(triface& t) {
	t.ver = enextesymtbl[t.ver];
}

// eprevesym()  finds the reversed edge of the previous edge.

inline void TetGenMesh::eprevesym(triface& t1, triface& t2) {
	t2.tet = t1.tet;
	t2.ver = eprevesymtbl[t1.ver];
}

inline void TetGenMesh::eprevesymself(triface& t) {
	t.ver = eprevesymtbl[t.ver];
}

// eorgoppo()    Finds the opposite face of the origin of the current edge.
//               Return the opposite edge of the current edge.

inline void TetGenMesh::eorgoppo(triface& t1, triface& t2) {
	t2.tet = t1.tet;
	t2.ver = eorgoppotbl[t1.ver];
}

inline void TetGenMesh::eorgoppoself(triface& t) {
	t.ver = eorgoppotbl[t.ver];
}

// edestoppo()    Finds the opposite face of the destination of the current 
//                edge. Return the opposite edge of the current edge.

inline void TetGenMesh::edestoppo(triface& t1, triface& t2) {
	t2.tet = t1.tet;
	t2.ver = edestoppotbl[t1.ver];
}

inline void TetGenMesh::edestoppoself(triface& t) {
	t.ver = edestoppotbl[t.ver];
}

// fsym()  finds the adjacent tetrahedron at the same face and the same edge.

inline void TetGenMesh::fsym(triface& t1, triface& t2) {
	decode((t1).tet[(t1).ver & 3], t2);
	t2.ver = fsymtbl[t1.ver][t2.ver];
}


#define fsymself(t) \
  t1ver = (t).ver; \
  decode((t).tet[(t).ver & 3], (t));\
  (t).ver = fsymtbl[t1ver][(t).ver]

// fnext()  finds the next face while rotating about an edge according to
//   a right-hand rule. The face is in the adjacent tetrahedron.  It is
//   the combination: fsym() * esym().

inline void TetGenMesh::fnext(triface& t1, triface& t2) {
	decode(t1.tet[facepivot1[t1.ver]], t2);
	t2.ver = facepivot2[t1.ver][t2.ver];
}


#define fnextself(t) \
  t1ver = (t).ver; \
  decode((t).tet[facepivot1[(t).ver]], (t)); \
  (t).ver = facepivot2[t1ver][(t).ver]


// The following primtives get or set the origin, destination, face apex,
//   or face opposite of an ordered tetrahedron.

inline TetGenMesh::point TetGenMesh::org(triface& t) {
	return (point) (t).tet[orgpivot[(t).ver]];
}

inline TetGenMesh::point TetGenMesh:: dest(triface& t) {
	return (point) (t).tet[destpivot[(t).ver]];
}

inline TetGenMesh::point TetGenMesh:: apex(triface& t) {
	return (point) (t).tet[apexpivot[(t).ver]];
}

inline TetGenMesh::point TetGenMesh:: oppo(triface& t) {
	return (point) (t).tet[oppopivot[(t).ver]];
}

inline void TetGenMesh:: setorg(triface& t, point p) {
	(t).tet[orgpivot[(t).ver]] = (tetrahedron) (p);
}

inline void TetGenMesh:: setdest(triface& t, point p) {
	(t).tet[destpivot[(t).ver]] = (tetrahedron) (p);
}

inline void TetGenMesh:: setapex(triface& t, point p) {
	(t).tet[apexpivot[(t).ver]] = (tetrahedron) (p);
}

inline void TetGenMesh:: setoppo(triface& t, point p) {
	(t).tet[oppopivot[(t).ver]] = (tetrahedron) (p);
}

#define setvertices(t, torg, tdest, tapex, toppo) \
  (t).tet[orgpivot[(t).ver]] = (tetrahedron) (torg);\
  (t).tet[destpivot[(t).ver]] = (tetrahedron) (tdest); \
  (t).tet[apexpivot[(t).ver]] = (tetrahedron) (tapex); \
  (t).tet[oppopivot[(t).ver]] = (tetrahedron) (toppo)

// Check or set a tetrahedron's attributes.

inline REAL TetGenMesh::elemattribute(tetrahedron* ptr, int attnum) {
	return ((REAL *) (ptr))[elemattribindex + attnum];
}

inline void TetGenMesh::setelemattribute(tetrahedron* ptr, int attnum,
										 REAL value) {
	((REAL *) (ptr))[elemattribindex + attnum] = value;
}

// Check or set a tetrahedron's maximum volume bound.

inline REAL TetGenMesh::volumebound(tetrahedron* ptr) {
	return ((REAL *) (ptr))[volumeboundindex];
}

inline void TetGenMesh::setvolumebound(tetrahedron* ptr, REAL value) {
	((REAL *) (ptr))[volumeboundindex] = value;
}

// Get or set a tetrahedron's index (only used for output).
//    These two routines use the reserved slot ptr[10].

inline int TetGenMesh::elemindex(tetrahedron* ptr) {
	int *iptr = (int *) &(ptr[10]);
	return iptr[0];
}

inline void TetGenMesh::setelemindex(tetrahedron* ptr, int value) {
	int *iptr = (int *) &(ptr[10]);
	iptr[0] = value;
}

// Get or set a tetrahedron's marker. 
//   Set 'value = 0' cleans all the face/edge flags.

inline int TetGenMesh::elemmarker(tetrahedron* ptr) {
	return ((int *) (ptr))[elemmarkerindex];
}

inline void TetGenMesh::setelemmarker(tetrahedron* ptr, int value) {
	((int *) (ptr))[elemmarkerindex] = value;
}

// infect(), infected(), uninfect() -- primitives to flag or unflag a
//   tetrahedron. The last bit of the element marker is flagged (1)
//   or unflagged (0).

inline void TetGenMesh::infect(triface& t) {
	((int *) (t.tet))[elemmarkerindex] |= 1;
}

inline void TetGenMesh::uninfect(triface& t) {
	((int *) (t.tet))[elemmarkerindex] &= ~1;
}

inline bool TetGenMesh::infected(triface& t) {
	return (((int *) (t.tet))[elemmarkerindex] & 1) != 0;
}

// marktest(), marktested(), unmarktest() -- primitives to flag or unflag a
//   tetrahedron.  Use the second lowerest bit of the element marker.

inline void TetGenMesh::marktest(triface& t) {
	((int *) (t.tet))[elemmarkerindex] |= 2;
}

inline void TetGenMesh::unmarktest(triface& t) {
	((int *) (t.tet))[elemmarkerindex] &= ~2;
}

inline bool TetGenMesh::marktested(triface& t) {
	return (((int *) (t.tet))[elemmarkerindex] & 2) != 0;
}

// markface(), unmarkface(), facemarked() -- primitives to flag or unflag a
//   face of a tetrahedron.  From the last 3rd to 6th bits are used for
//   face markers, e.g., the last third bit corresponds to loc = 0. 

inline void TetGenMesh::markface(triface& t) {
	((int *) (t.tet))[elemmarkerindex] |= (4 << (t.ver & 3));
}

inline void TetGenMesh::unmarkface(triface& t) {
	((int *) (t.tet))[elemmarkerindex] &= ~(4 << (t.ver & 3));
}

inline bool TetGenMesh::facemarked(triface& t) {
	return (((int *) (t.tet))[elemmarkerindex] & (4 << (t.ver & 3))) != 0;
}

// markedge(), unmarkedge(), edgemarked() -- primitives to flag or unflag an
//   edge of a tetrahedron.  From the last 7th to 12th bits are used for
//   edge markers, e.g., the last 7th bit corresponds to the 0th edge, etc. 
//   Remark: The last 7th bit is marked by 2^6 = 64.

inline void TetGenMesh::markedge(triface& t) {
	((int *) (t.tet))[elemmarkerindex] |= (int) (64 << ver2edge[(t).ver]);
}

inline void TetGenMesh::unmarkedge(triface& t) {
	((int *) (t.tet))[elemmarkerindex] &= ~(int) (64 << ver2edge[(t).ver]);
}

inline bool TetGenMesh::edgemarked(triface& t) {
	return (((int *) (t.tet))[elemmarkerindex] &
			(int) (64 << ver2edge[(t).ver])) != 0;
}

// marktest2(), unmarktest2(), marktest2ed() -- primitives to flag and unflag
//   a tetrahedron. The 13th bit (2^12 = 4096) is used for this flag.

inline void TetGenMesh::marktest2(triface& t) {
	((int *) (t.tet))[elemmarkerindex] |= (int) (4096);
}

inline void TetGenMesh::unmarktest2(triface& t) {
	((int *) (t.tet))[elemmarkerindex] &= ~(int) (4096);
}

inline bool TetGenMesh::marktest2ed(triface& t) {
	return (((int *) (t.tet))[elemmarkerindex] & (int) (4096)) != 0;
}

// elemcounter(), setelemcounter() -- primitives to read or ser a (small)
//   integer counter in this tet. It is saved from the 16th bit. On 32 bit
//   system, the range of the counter is [0, 2^15 = 32768]. 

inline int TetGenMesh::elemcounter(triface& t) {
	return (((int *) (t.tet))[elemmarkerindex]) >> 16;
}

inline void TetGenMesh::setelemcounter(triface& t, int value) {
	int c = ((int *) (t.tet))[elemmarkerindex];
	// Clear the old counter while keep the other flags.
	c &= 65535; // sum_{i=0^15} 2^i
	c |= (value << 16);
	((int *) (t.tet))[elemmarkerindex] = c;
}

inline void TetGenMesh::increaseelemcounter(triface& t) {
	int c = elemcounter(t);
	setelemcounter(t, c + 1);
}

inline void TetGenMesh::decreaseelemcounter(triface& t) {
	int c = elemcounter(t);
	setelemcounter(t, c - 1);
}

// ishulltet()  tests if t is a hull tetrahedron.

inline bool TetGenMesh::ishulltet(triface& t) {
	return (point) (t).tet[7] == dummypoint;
}

// isdeadtet()  tests if t is a tetrahedron is dead.

inline bool TetGenMesh::isdeadtet(triface& t) {
	return ((t.tet == NULL) || (t.tet[4] == NULL));
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// Primitives for subfaces and subsegments                                   //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

// Each subface contains three pointers to its neighboring subfaces, with
//   edge versions.  To save memory, both information are kept in a single
//   pointer. To make this possible, all subfaces are aligned to eight-byte
//   boundaries, so that the last three bits of each pointer are zeros. An
//   edge version (in the range 0 to 5) is compressed into the last three
//   bits of each pointer by 'sencode()'.  'sdecode()' decodes a pointer,
//   extracting an edge version and a pointer to the beginning of a subface.

inline void TetGenMesh::sdecode(shellface sptr, face& s) {
	s.shver = (int) ((uintptr_t) (sptr) & (uintptr_t) 7);
	s.sh = (shellface *) ((uintptr_t) (sptr) ^ (uintptr_t) (s.shver));
}

inline TetGenMesh::shellface TetGenMesh::sencode(face& s) {
	return (shellface) ((uintptr_t) s.sh | (uintptr_t) s.shver);
}

inline TetGenMesh::shellface TetGenMesh::sencode2(shellface *sh, int shver) {
	return (shellface) ((uintptr_t) sh | (uintptr_t) shver);
}

// sbond() bonds two subfaces (s1) and (s2) together. s1 and s2 must refer
//   to the same edge. No requirement is needed on their orientations.

inline void TetGenMesh::sbond(face& s1, face& s2)
{
	s1.sh[s1.shver >> 1] = sencode(s2);
	s2.sh[s2.shver >> 1] = sencode(s1);
}

// sbond1() bonds s1 <== s2, i.e., after bonding, s1 is pointing to s2,
//   but s2 is not pointing to s1.  s1 and s2 must refer to the same edge.
//   No requirement is needed on their orientations.

inline void TetGenMesh::sbond1(face& s1, face& s2)
{
	s1.sh[s1.shver >> 1] = sencode(s2);
}

// Dissolve a subface bond (from one side).  Note that the other subface
//   will still think it's connected to this subface.

inline void TetGenMesh::sdissolve(face& s)
{
	s.sh[s.shver >> 1] = NULL;
}

// spivot() finds the adjacent subface (s2) for a given subface (s1).
//   s1 and s2 share at the same edge.

inline void TetGenMesh::spivot(face& s1, face& s2)
{
	shellface sptr = s1.sh[s1.shver >> 1];
	sdecode(sptr, s2);
}

inline void TetGenMesh::spivotself(face& s)
{
	shellface sptr = s.sh[s.shver >> 1];
	sdecode(sptr, s);
}

// These primitives determine or set the origin, destination, or apex
//   of a subface with respect to the edge version.

inline TetGenMesh::point TetGenMesh::sorg(face& s)
{
	return (point) s.sh[sorgpivot[s.shver]];
}

inline TetGenMesh::point TetGenMesh::sdest(face& s)
{
	return (point) s.sh[sdestpivot[s.shver]];
}

inline TetGenMesh::point TetGenMesh::sapex(face& s)
{
	return (point) s.sh[sapexpivot[s.shver]];
}

inline void TetGenMesh::setsorg(face& s, point pointptr)
{
	s.sh[sorgpivot[s.shver]] = (shellface) pointptr;
}

inline void TetGenMesh::setsdest(face& s, point pointptr)
{
	s.sh[sdestpivot[s.shver]] = (shellface) pointptr;
}

inline void TetGenMesh::setsapex(face& s, point pointptr)
{
	s.sh[sapexpivot[s.shver]] = (shellface) pointptr;
}

#define setshvertices(s, pa, pb, pc)\
  setsorg(s, pa);\
  setsdest(s, pb);\
  setsapex(s, pc)

// sesym()  reserves the direction of the lead edge.

inline void TetGenMesh::sesym(face& s1, face& s2)
{
	s2.sh = s1.sh;
	s2.shver = (s1.shver ^ 1);  // Inverse the last bit.
}

inline void TetGenMesh::sesymself(face& s)
{
	s.shver ^= 1;
}

// senext()  finds the next edge (counterclockwise) in the same orientation
//   of this face.

inline void TetGenMesh::senext(face& s1, face& s2)
{
	s2.sh = s1.sh;
	s2.shver = snextpivot[s1.shver];
}

inline void TetGenMesh::senextself(face& s)
{
	s.shver = snextpivot[s.shver];
}

inline void TetGenMesh::senext2(face& s1, face& s2)
{
	s2.sh = s1.sh;
	s2.shver = snextpivot[snextpivot[s1.shver]];
}

inline void TetGenMesh::senext2self(face& s)
{
	s.shver = snextpivot[snextpivot[s.shver]];
}


// Check or set a subface's maximum area bound.

inline REAL TetGenMesh::areabound(face& s)
{
	return ((REAL *) (s.sh))[areaboundindex];
}

inline void TetGenMesh::setareabound(face& s, REAL value)
{
	((REAL *) (s.sh))[areaboundindex] = value;
}

// These two primitives read or set a shell marker.  Shell markers are used
//   to hold user boundary information.

inline int TetGenMesh::shellmark(face& s)
{
	return ((int *) (s.sh))[shmarkindex];
}

inline void TetGenMesh::setshellmark(face& s, int value)
{
	((int *) (s.sh))[shmarkindex] = value;
}



// sinfect(), sinfected(), suninfect() -- primitives to flag or unflag a
//   subface. The last bit of ((int *) ((s).sh))[shmarkindex+1] is flagged.

inline void TetGenMesh::sinfect(face& s)
{
	((int *) ((s).sh))[shmarkindex+1] =
			(((int *) ((s).sh))[shmarkindex+1] | (int) 1);
}

inline void TetGenMesh::suninfect(face& s)
{
	((int *) ((s).sh))[shmarkindex+1] =
			(((int *) ((s).sh))[shmarkindex+1] & ~(int) 1);
}

// Test a subface for viral infection.

inline bool TetGenMesh::sinfected(face& s)
{
	return (((int *) ((s).sh))[shmarkindex+1] & (int) 1) != 0;
}

// smarktest(), smarktested(), sunmarktest() -- primitives to flag or unflag
//   a subface. The last 2nd bit of the integer is flagged.

inline void TetGenMesh::smarktest(face& s)
{
	((int *) ((s).sh))[shmarkindex+1] =
			(((int *)((s).sh))[shmarkindex+1] | (int) 2);
}

inline void TetGenMesh::sunmarktest(face& s)
{
	((int *) ((s).sh))[shmarkindex+1] =
			(((int *)((s).sh))[shmarkindex+1] & ~(int)2);
}

inline bool TetGenMesh::smarktested(face& s)
{
	return ((((int *) ((s).sh))[shmarkindex+1] & (int) 2) != 0);
}

// smarktest2(), smarktest2ed(), sunmarktest2() -- primitives to flag or 
//   unflag a subface. The last 3rd bit of the integer is flagged.

inline void TetGenMesh::smarktest2(face& s)
{
	((int *) ((s).sh))[shmarkindex+1] =
			(((int *)((s).sh))[shmarkindex+1] | (int) 4);
}

inline void TetGenMesh::sunmarktest2(face& s)
{
	((int *) ((s).sh))[shmarkindex+1] =
			(((int *)((s).sh))[shmarkindex+1] & ~(int)4);
}

inline bool TetGenMesh::smarktest2ed(face& s)
{
	return ((((int *) ((s).sh))[shmarkindex+1] & (int) 4) != 0);
}

// The last 4th bit of ((int *) ((s).sh))[shmarkindex+1] is flagged.

inline void TetGenMesh::smarktest3(face& s)
{
	((int *) ((s).sh))[shmarkindex+1] =
			(((int *)((s).sh))[shmarkindex+1] | (int) 8);
}

inline void TetGenMesh::sunmarktest3(face& s)
{
	((int *) ((s).sh))[shmarkindex+1] =
			(((int *)((s).sh))[shmarkindex+1] & ~(int)8);
}

inline bool TetGenMesh::smarktest3ed(face& s)
{
	return ((((int *) ((s).sh))[shmarkindex+1] & (int) 8) != 0);
}


// Each facet has a unique index (automatically indexed). Starting from '0'.
// We save this index in the same field of the shell type. 

inline void TetGenMesh::setfacetindex(face& s, int value)
{
	((int *) (s.sh))[shmarkindex + 2] = value;
}

inline int TetGenMesh::getfacetindex(face& s)
{
	return ((int *) (s.sh))[shmarkindex + 2];
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// Primitives for interacting between tetrahedra and subfaces                //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

// tsbond() bond a tetrahedron (t) and a subface (s) together.
// Note that t and s must be the same face and the same edge. Moreover,
//   t and s have the same orientation. 
// Since the edge number in t and in s can be any number in {0,1,2}. We bond
//   the edge in s which corresponds to t's 0th edge, and vice versa.

inline void TetGenMesh::tsbond(triface& t, face& s)
{
	if ((t).tet[9] == NULL) {
		// Allocate space for this tet.
		(t).tet[9] = (tetrahedron) tet2subpool->alloc();
		// Initialize.
		for (int i = 0; i < 4; i++) {
			((shellface *) (t).tet[9])[i] = NULL;
		}
	}
	// Bond t <== s.
	((shellface *) (t).tet[9])[(t).ver & 3] =
			sencode2((s).sh, tsbondtbl[t.ver][s.shver]);
	// Bond s <== t.
	s.sh[9 + ((s).shver & 1)] =
			(shellface) encode2((t).tet, stbondtbl[t.ver][s.shver]);
}

// tspivot() finds a subface (s) abutting on the given tetrahdera (t).
//   Return s.sh = NULL if there is no subface at t. Otherwise, return
//   the subface s, and s and t must be at the same edge wth the same
//   orientation.

inline void TetGenMesh::tspivot(triface& t, face& s)
{
	if ((t).tet[9] == NULL) {
		(s).sh = NULL;
		return;
	}
	// Get the attached subface s.
	sdecode(((shellface *) (t).tet[9])[(t).ver & 3], (s));
	(s).shver = tspivottbl[t.ver][s.shver];
}

// Quickly check if the handle (t, v) is a subface.
#define issubface(t) \
  ((t).tet[9] && ((t).tet[9])[(t).ver & 3])

// stpivot() finds a tetrahedron (t) abutting a given subface (s).
//   Return the t (if it exists) with the same edge and the same
//   orientation of s.

inline void TetGenMesh::stpivot(face& s, triface& t)
{
	decode((tetrahedron) s.sh[9 + (s.shver & 1)], t);
	if ((t).tet == NULL) {
		return;
	}
	(t).ver = stpivottbl[t.ver][s.shver];
}

// Quickly check if this subface is attached to a tetrahedron.

#define isshtet(s) \
  ((s).sh[9 + ((s).shver & 1)])

// tsdissolve() dissolve a bond (from the tetrahedron side).

inline void TetGenMesh::tsdissolve(triface& t)
{
	if ((t).tet[9] != NULL) {
		((shellface *) (t).tet[9])[(t).ver & 3] = NULL;
	}
}

// stdissolve() dissolve a bond (from the subface side).

inline void TetGenMesh::stdissolve(face& s)
{
	(s).sh[9] = NULL;
	(s).sh[10] = NULL;
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// Primitives for interacting between subfaces and segments                  //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

// ssbond() bond a subface to a subsegment.

inline void TetGenMesh::ssbond(face& s, face& edge)
{
	s.sh[6 + (s.shver >> 1)] = sencode(edge);
	edge.sh[0] = sencode(s);
}

inline void TetGenMesh::ssbond1(face& s, face& edge)
{
	s.sh[6 + (s.shver >> 1)] = sencode(edge);
	//edge.sh[0] = sencode(s);
}

// ssdisolve() dissolve a bond (from the subface side)

inline void TetGenMesh::ssdissolve(face& s)
{
	s.sh[6 + (s.shver >> 1)] = NULL;
}

// sspivot() finds a subsegment abutting a subface.

inline void TetGenMesh::sspivot(face& s, face& edge)
{
	sdecode((shellface) s.sh[6 + (s.shver >> 1)], edge);
}

// Quickly check if the edge is a subsegment.

#define isshsubseg(s) \
  ((s).sh[6 + ((s).shver >> 1)])

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// Primitives for interacting between tetrahedra and segments                //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

inline void TetGenMesh::tssbond1(triface& t, face& s)
{
	if ((t).tet[8] == NULL) {
		// Allocate space for this tet.
		(t).tet[8] = (tetrahedron) tet2segpool->alloc();
		// Initialization.
		for (int i = 0; i < 6; i++) {
			((shellface *) (t).tet[8])[i] = NULL;
		}
	}
	((shellface *) (t).tet[8])[ver2edge[(t).ver]] = sencode((s));
}

inline void TetGenMesh::sstbond1(face& s, triface& t)
{
	((tetrahedron *) (s).sh)[9] = encode(t);
}

inline void TetGenMesh::tssdissolve1(triface& t)
{
	if ((t).tet[8] != NULL) {
		((shellface *) (t).tet[8])[ver2edge[(t).ver]] = NULL;
	}
}

inline void TetGenMesh::sstdissolve1(face& s)
{
	((tetrahedron *) (s).sh)[9] = NULL;
}

inline void TetGenMesh::tsspivot1(triface& t, face& s)
{
	if ((t).tet[8] != NULL) {
		sdecode(((shellface *) (t).tet[8])[ver2edge[(t).ver]], s);
	} else {
		(s).sh = NULL;
	}
}

// Quickly check whether 't' is a segment or not.

#define issubseg(t) \
  ((t).tet[8] && ((t).tet[8])[ver2edge[(t).ver]])

inline void TetGenMesh::sstpivot1(face& s, triface& t)
{
	decode((tetrahedron) s.sh[9], t);
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// Primitives for points                                                     //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

inline int TetGenMesh::pointmark(point pt) {
	return ((int *) (pt))[pointmarkindex];
}

inline void TetGenMesh::setpointmark(point pt, int value) {
	((int *) (pt))[pointmarkindex] = value;
}


// These two primitives set and read the type of the point.

inline enum TetGenMesh::verttype TetGenMesh::pointtype(point pt) {
	return (enum verttype) (((int *) (pt))[pointmarkindex + 1] >> (int) 8);
}

inline void TetGenMesh::setpointtype(point pt, enum verttype value) {
	((int *) (pt))[pointmarkindex + 1] =
			((int) value << 8) + (((int *) (pt))[pointmarkindex + 1] & (int) 255);
}

// Read and set the geometry tag of the point (used by -s option).

inline int TetGenMesh::pointgeomtag(point pt) {
	return ((int *) (pt))[pointmarkindex + 2];
}

inline void TetGenMesh::setpointgeomtag(point pt, int value) {
	((int *) (pt))[pointmarkindex + 2] = value;
}

// Read and set the u,v coordinates of the point (used by -s option).

inline REAL TetGenMesh::pointgeomuv(point pt, int i) {
	return pt[pointparamindex + i];
}

inline void TetGenMesh::setpointgeomuv(point pt, int i, REAL value) {
	pt[pointparamindex + i] = value;
}

// pinfect(), puninfect(), pinfected() -- primitives to flag or unflag
//   a point. The last bit of the integer '[pointindex+1]' is flagged.

inline void TetGenMesh::pinfect(point pt) {
	((int *) (pt))[pointmarkindex + 1] |= (int) 1;
}

inline void TetGenMesh::puninfect(point pt) {
	((int *) (pt))[pointmarkindex + 1] &= ~(int) 1;
}

inline bool TetGenMesh::pinfected(point pt) {
	return (((int *) (pt))[pointmarkindex + 1] & (int) 1) != 0;
}

// pmarktest(), punmarktest(), pmarktested() -- more primitives to 
//   flag or unflag a point. 

inline void TetGenMesh::pmarktest(point pt) {
	((int *) (pt))[pointmarkindex + 1] |= (int) 2;
}

inline void TetGenMesh::punmarktest(point pt) {
	((int *) (pt))[pointmarkindex + 1] &= ~(int) 2;
}

inline bool TetGenMesh::pmarktested(point pt) {
	return (((int *) (pt))[pointmarkindex + 1] & (int) 2) != 0;
}

inline void TetGenMesh::pmarktest2(point pt) {
	((int *) (pt))[pointmarkindex + 1] |= (int) 4;
}

inline void TetGenMesh::punmarktest2(point pt) {
	((int *) (pt))[pointmarkindex + 1] &= ~(int) 4;
}

inline bool TetGenMesh::pmarktest2ed(point pt) {
	return (((int *) (pt))[pointmarkindex + 1] & (int) 4) != 0;
}

inline void TetGenMesh::pmarktest3(point pt) {
	((int *) (pt))[pointmarkindex + 1] |= (int) 8;
}

inline void TetGenMesh::punmarktest3(point pt) {
	((int *) (pt))[pointmarkindex + 1] &= ~(int) 8;
}

inline bool TetGenMesh::pmarktest3ed(point pt) {
	return (((int *) (pt))[pointmarkindex + 1] & (int) 8) != 0;
}

// These following primitives set and read a pointer to a tetrahedron
//   a subface/subsegment, a point, or a tet of background mesh.

inline TetGenMesh::tetrahedron TetGenMesh::point2tet(point pt) {
	return ((tetrahedron *) (pt))[point2simindex];
}

inline void TetGenMesh::setpoint2tet(point pt, tetrahedron value) {
	((tetrahedron *) (pt))[point2simindex] = value;
}

inline TetGenMesh::point TetGenMesh::point2ppt(point pt) {
	return (point) ((tetrahedron *) (pt))[point2simindex + 1];
}

inline void TetGenMesh::setpoint2ppt(point pt, point value) {
	((tetrahedron *) (pt))[point2simindex + 1] = (tetrahedron) value;
}

inline TetGenMesh::shellface TetGenMesh::point2sh(point pt) {
	return (shellface) ((tetrahedron *) (pt))[point2simindex + 2];
}

inline void TetGenMesh::setpoint2sh(point pt, shellface value) {
	((tetrahedron *) (pt))[point2simindex + 2] = (tetrahedron) value;
}


inline TetGenMesh::tetrahedron TetGenMesh::point2bgmtet(point pt) {
	return ((tetrahedron *) (pt))[point2simindex + 3];
}

inline void TetGenMesh::setpoint2bgmtet(point pt, tetrahedron value) {
	((tetrahedron *) (pt))[point2simindex + 3] = value;
}


// The primitives for saving and getting the insertion radius.
inline void TetGenMesh::setpointinsradius(point pt, REAL value)
{
	pt[pointmtrindex + sizeoftensor - 1] = value;
}

inline REAL TetGenMesh::getpointinsradius(point pt)
{
	return pt[pointmtrindex + sizeoftensor - 1];
}

// point2tetorg()    Get the tetrahedron whose origin is the point.

inline void TetGenMesh::point2tetorg(point pa, triface& searchtet)
{
	decode(point2tet(pa), searchtet);
	if ((point) searchtet.tet[4] == pa) {
		searchtet.ver = 11;
	} else if ((point) searchtet.tet[5] == pa) {
		searchtet.ver = 3;
	} else if ((point) searchtet.tet[6] == pa) {
		searchtet.ver = 7;
	} else {
		assert((point) searchtet.tet[7] == pa); // SELF_CHECK
		searchtet.ver = 0;
	}
}

// point2shorg()    Get the subface/segment whose origin is the point.

inline void TetGenMesh::point2shorg(point pa, face& searchsh)
{
	sdecode(point2sh(pa), searchsh);
	if ((point) searchsh.sh[3] == pa) {
		searchsh.shver = 0;
	} else if ((point) searchsh.sh[4] == pa) {
		searchsh.shver = (searchsh.sh[5] != NULL ? 2 : 1);
	} else {
		assert((point) searchsh.sh[5] == pa); // SELF_CHECK
		searchsh.shver = 4;
	}
}

// farsorg()    Return the origin of the subsegment.
// farsdest()   Return the destination of the subsegment.

inline TetGenMesh::point TetGenMesh::farsorg(face& s)
{
	face travesh, neighsh;

	travesh = s;
	while (1) {
		senext2(travesh, neighsh);
		spivotself(neighsh);
		if (neighsh.sh == NULL) break;
		if (sorg(neighsh) != sorg(travesh)) sesymself(neighsh);
		senext2(neighsh, travesh);
	}
	return sorg(travesh);
}

inline TetGenMesh::point TetGenMesh::farsdest(face& s)
{
	face travesh, neighsh;

	travesh = s;
	while (1) {
		senext(travesh, neighsh);
		spivotself(neighsh);
		if (neighsh.sh == NULL) break;
		if (sdest(neighsh) != sdest(travesh)) sesymself(neighsh);
		senext(neighsh, travesh);
	}
	return sdest(travesh);
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// Linear algebra operators.                                                 //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

// dot() returns the dot product: v1 dot v2.
inline REAL TetGenMesh::dot(REAL* v1, REAL* v2)
{
	return v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2];
}

// cross() computes the cross product: n = v1 cross v2.
inline void TetGenMesh::cross(REAL* v1, REAL* v2, REAL* n)
{
	n[0] =   v1[1] * v2[2] - v2[1] * v1[2];
	n[1] = -(v1[0] * v2[2] - v2[0] * v1[2]);
	n[2] =   v1[0] * v2[1] - v2[0] * v1[1];
}

// distance() computes the Euclidean distance between two points.
inline REAL TetGenMesh::distance(REAL* p1, REAL* p2)
{
	return sqrt((p2[0] - p1[0]) * (p2[0] - p1[0]) +
				(p2[1] - p1[1]) * (p2[1] - p1[1]) +
				(p2[2] - p1[2]) * (p2[2] - p1[2]));
}

inline REAL TetGenMesh::norm2(REAL x, REAL y, REAL z)
{
	return (x) * (x) + (y) * (y) + (z) * (z);
}



#endif //SHAPE_RECONSTRUCTION_TETGENMESH_H
