//
// Created by default on 17. 12. 16.
//

#ifndef SHAPE_RECONSTRUCTION_RAWDATAREADER_H
#define SHAPE_RECONSTRUCTION_RAWDATAREADER_H


#include <fstream>

class RawDataReader {
public:
	~RawDataReader();

	/**
	 *
	 * @param _arg Filename
	 * @return 640*480 float array
	 */
	void SetFileName(const char *_arg);

	void SetFileName(const char *_arg, int startFrame);

	/**
	 *
	 * @return 640*480 float array of next frame
	 */
	float *ReadNextFrame();

private:
	int index;
	std::ifstream file;
	float *bufferFrame;
};


#endif //SHAPE_RECONSTRUCTION_RAWDATAREADER_H
