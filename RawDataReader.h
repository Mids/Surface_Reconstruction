//
// Created by default on 17. 12. 16.
//

#ifndef SHAPE_RECONSTRUCTION_RAWDATAREADER_H
#define SHAPE_RECONSTRUCTION_RAWDATAREADER_H


#include <fstream>

class RawDataReader {
public:
	/**
	 *
	 * @param _arg Filename
	 * @return 640*480 float array
	 */
	float * ReadData(const char *_arg);

	/**
	 *
	 * @param input	Input filestream
	 * @return	640*480 float array of current frame
	 */
	float *ReadFrame(std::ifstream &input);
};


#endif //SHAPE_RECONSTRUCTION_RAWDATAREADER_H
