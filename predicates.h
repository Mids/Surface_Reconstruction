//
// Created by default on 18. 1. 2.
//

#ifndef SHAPE_RECONSTRUCTION_CONSTANTS_H
#define SHAPE_RECONSTRUCTION_CONSTANTS_H


// TetGen default uses the double precision (64 bit) for a real number.
//   Alternatively, one can use the single precision (32 bit) 'float' if the
//   memory is limited.

#define REAL double  // #define REAL float

// Maximum number of characters in a file name (including the null).

#define FILENAMESIZE 1024

// Maximum number of chars in a line read from a file (including the null).

#define INPUTLINESIZE 2048


void exactinit(int, int, int, REAL, REAL, REAL);
REAL orient3d(REAL *pa, REAL *pb, REAL *pc, REAL *pd);
REAL insphere(REAL *pa, REAL *pb, REAL *pc, REAL *pd, REAL *pe);
REAL orient4d(REAL *pa, REAL *pb, REAL *pc, REAL *pd, REAL *pe,
			  REAL ah, REAL bh, REAL ch, REAL dh, REAL eh);


#endif //SHAPE_RECONSTRUCTION_CONSTANTS_H
