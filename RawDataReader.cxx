//
// Created by default on 17. 12. 16.
//

#include "RawDataReader.h"
#include <iostream>

#define KINECT_X_RES 640
#define KINECT_Y_RES 480
#define DUMPSIZE 150

// Read Data from the file
void RawDataReader::SetFileName(const char *_arg) {
	SetFileName(_arg, 0);
}

void RawDataReader::SetFileName(const char *_arg, int startFrame) {
	file.open(_arg, std::ios::binary);
	index = DUMPSIZE;
	if (file.is_open()) {
		// The first address of the frame
		index += KINECT_X_RES * KINECT_Y_RES * startFrame * 2;

		// Move pointer to the exact location
		file.seekg(index, std::ios_base::beg);
	}
}

// Read 640*480 data
float *RawDataReader::ReadNextFrame() {
	if(bufferFrame == nullptr)
		bufferFrame = new float[KINECT_X_RES * KINECT_Y_RES];
	unsigned short buffer;

	// Check whether opened or not
	if (!file.is_open()) {
		std::cout << "File is not opened\n";
		return nullptr;
	}

	// Read data
	for (int i = 0; i < KINECT_Y_RES; ++i) {
		for (int j = 0; j < KINECT_X_RES; ++j) {
			// If it's EOF, return null
			if (!file.read(reinterpret_cast<char *>(&buffer), sizeof(unsigned short))) return nullptr;

			bufferFrame[i * KINECT_X_RES + j] = buffer;
		}
	}

	return bufferFrame;
}

RawDataReader::~RawDataReader() {
	file.close();
}
