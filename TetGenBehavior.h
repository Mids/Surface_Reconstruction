//
// Created by default on 17. 12. 27.
//

#ifndef SHAPE_RECONSTRUCTION_TETGENBEHAVIOR_H
#define SHAPE_RECONSTRUCTION_TETGENBEHAVIOR_H



///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// tetgenbehavior                                                            //
//                                                                           //
// A structure for maintaining the switches and parameters used by TetGen's  //
// mesh data structure and algorithms.                                       //
//                                                                           //
// All switches and parameters are initialized with default values. They can //
// be set by the command line arguments (a list of strings) of TetGen.       //
//                                                                           //
// NOTE: Some of the switches are incompatible. While some may depend on     //
// other switches.  The routine parse_commandline() sets the switches from   //
// the command line (a list of strings) and checks the consistency of the    //
// applied switches.                                                         //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "predicates.h"

class TetGenBehavior {

public:

	// Switches of TetGen.
	int plc;                                                         // '-p', 0.
	int psc;                                                         // '-s', 0.
	int refine;                                                      // '-r', 0.
	int quality;                                                     // '-q', 0.
	int nobisect;                                                    // '-Y', 0.
	int coarsen;                                                     // '-R', 0.
	int weighted;                                                    // '-w', 0.
	int brio_hilbert;                                                // '-b', 1.
	int incrflip;                                                    // '-l', 0.
	int flipinsert;                                                  // '-L', 0.
	int metric;                                                      // '-m', 0.
	int varvolume;                                                   // '-a', 0.
	int fixedvolume;                                                 // '-a', 0.
	int regionattrib;                                                // '-A', 0.
	int conforming;                                                  // '-D', 0.
	int insertaddpoints;                                             // '-i', 0.
	int diagnose;                                                    // '-d', 0.
	int convex;                                                      // '-c', 0.
	int nomergefacet;                                                // '-M', 0.
	int nomergevertex;                                               // '-M', 0.
	int noexact;                                                     // '-X', 0.
	int nostaticfilter;                                              // '-X', 0.
	int zeroindex;                                                   // '-z', 0.
	int facesout;                                                    // '-f', 0.
	int edgesout;                                                    // '-e', 0.
	int neighout;                                                    // '-n', 0.
	int voroout;                                                     // '-v', 0.
	int meditview;                                                   // '-g', 0.
	int vtkview;                                                     // '-k', 0.
	int nobound;                                                     // '-B', 0.
	int nonodewritten;                                               // '-N', 0.
	int noelewritten;                                                // '-E', 0.
	int nofacewritten;                                               // '-F', 0.
	int noiterationnum;                                              // '-I', 0.
	int nojettison;                                                  // '-J', 0.
	int reversetetori;                                               // '-R', 0.
	int docheck;                                                     // '-C', 0.
	int quiet;                                                       // '-Q', 0.
	int verbose;                                                     // '-V', 0.

	// Parameters of TetGen.
	int vertexperblock;                                           // '-x', 4092.
	int tetrahedraperblock;                                       // '-x', 8188.
	int shellfaceperblock;                                        // '-x', 2044.
	int nobisect_param;                                              // '-Y', 2.
	int addsteiner_algo;                                            // '-Y/', 1.
	int coarsen_param;                                               // '-R', 0.
	int weighted_param;                                              // '-w', 0.
	int fliplinklevel;                                                    // -1.
	int flipstarsize;                                                     // -1.
	int fliplinklevelinc;                                                 //  1.
	int reflevel;                                                    // '-D', 3.
	int optlevel;                                                    // '-O', 2.
	int optscheme;                                                   // '-O', 7.
	int delmaxfliplevel;                                                   // 1.
	int order;                                                       // '-o', 1.
	int steinerleft;                                                 // '-S', 0.
	int no_sort;                                                           // 0.
	int hilbert_order;                                           // '-b///', 52.
	int hilbert_limit;                                             // '-b//'  8.
	int brio_threshold;                                              // '-b' 64.
	REAL brio_ratio;                                             // '-b/' 0.125.
	REAL facet_ang_tol;                                          // '-p', 179.9.
	REAL maxvolume;                                               // '-a', -1.0.
	REAL minratio;                                                 // '-q', 0.0.
	REAL mindihedral;                                              // '-q', 5.0.
	REAL optmaxdihedral;                                               // 165.0.
	REAL optminsmtdihed;                                               // 179.0.
	REAL optminslidihed;                                               // 179.0.
	REAL epsilon;                                               // '-T', 1.0e-8.
	REAL minedgelength;                                                  // 0.0.
	REAL coarsen_percent;                                         // -R1/#, 1.0.

	// Strings of command line arguments and input/output file names.
	char commandline[1024];
	char infilename[1024];
	char outfilename[1024];
	char addinfilename[1024];
	char bgmeshfilename[1024];

	// The input object of TetGen. They are recognized by either the input
	//   file extensions or by the specified options.
	// Currently the following objects are supported:
	//   - NODES, a list of nodes (.node);
	//   - POLY, a piecewise linear complex (.poly or .smesh);
	//   - OFF, a polyhedron (.off, Geomview's file format);
	//   - PLY, a polyhedron (.ply, file format from gatech, only ASCII);
	//   - STL, a surface mesh (.stl, stereolithography format);
	//   - MEDIT, a surface mesh (.mesh, Medit's file format);
	//   - MESH, a tetrahedral mesh (.ele).
	// If no extension is available, the imposed command line switch
	//   (-p or -r) implies the object.
	enum objecttype {NODES, POLY, OFF, PLY, STL, MEDIT, VTK, MESH} object;


	void syntax();
	void usage();

	// Command line parse routine.
	bool parse_commandline(int argc, char **argv);
	bool parse_commandline(char *switches) {
		return parse_commandline(0, &switches);
	}

	// Initialize all variables.
	TetGenBehavior()
	{
		plc = 0;
		psc = 0;
		refine = 0;
		quality = 0;
		nobisect = 0;
		coarsen = 0;
		metric = 0;
		weighted = 0;
		brio_hilbert = 1;
		incrflip = 0;
		flipinsert = 0;
		varvolume = 0;
		fixedvolume = 0;
		noexact = 0;
		nostaticfilter = 0;
		insertaddpoints = 0;
		regionattrib = 0;
		conforming = 0;
		diagnose = 0;
		convex = 0;
		zeroindex = 0;
		facesout = 0;
		edgesout = 0;
		neighout = 0;
		voroout = 0;
		meditview = 0;
		vtkview = 0;
		nobound = 0;
		nonodewritten = 0;
		noelewritten = 0;
		nofacewritten = 0;
		noiterationnum = 0;
		nomergefacet = 0;
		nomergevertex = 0;
		nojettison = 0;
		reversetetori = 0;
		docheck = 0;
		quiet = 0;
		verbose = 0;

		vertexperblock = 4092;
		tetrahedraperblock = 8188;
		shellfaceperblock = 4092;
		nobisect_param = 2;
		addsteiner_algo = 1;
		coarsen_param = 0;
		weighted_param = 0;
		fliplinklevel = -1; // No limit on linklevel.
		flipstarsize = -1;  // No limit on flip star size.
		fliplinklevelinc = 1;
		reflevel = 3;
		optscheme = 7;  // 1 & 2 & 4, // min_max_dihedral.
		optlevel = 2;
		delmaxfliplevel = 1;
		order = 1;
		steinerleft = -1;
		no_sort = 0;
		hilbert_order = 52; //-1;
		hilbert_limit = 8;
		brio_threshold = 64;
		brio_ratio = 0.125;
		facet_ang_tol = 179.9;
		maxvolume = -1.0;
		minratio = 2.0;
		mindihedral = 0.0; // 5.0;
		optmaxdihedral = 165.00; // without -q, default is 179.0
		optminsmtdihed = 179.00; // without -q, default is 179.999
		optminslidihed = 179.00; // without -q, default is 179.999
		epsilon = 1.0e-8;
		minedgelength = 0.0;
		coarsen_percent = 1.0;
		object = NODES;

		commandline[0] = '\0';
		infilename[0] = '\0';
		outfilename[0] = '\0';
		addinfilename[0] = '\0';
		bgmeshfilename[0] = '\0';

	}

}; // class tetgenbehavior



#endif //SHAPE_RECONSTRUCTION_TETGENBEHAVIOR_H
