//
// Created by default on 17. 12. 16.
//

#include "RawDataReader.h"
#include <iostream>

#define KINECT_X_RES 640
#define KINECT_Y_RES 480
#define FRAME_NUM 100
#define DUMPSIZE 150

// Read Data from the file
float *RawDataReader::ReadData(const char *_arg) {
	std::ifstream inputFile(_arg, std::ios::binary);
	float *result = nullptr;
	if (inputFile.is_open()) {
		// The first address of the frame
		int startIndex = KINECT_X_RES * KINECT_Y_RES * FRAME_NUM * 2;

		// Move pointer to the exact location
		inputFile.seekg(startIndex + DUMPSIZE, std::ios_base::beg);

		// Read frame
		result = ReadFrame(inputFile);

		inputFile.close();
	}
	return result;
}

// Read 640*480 data
float *RawDataReader::ReadFrame(std::ifstream &input) {
	float *array = new float[KINECT_X_RES * KINECT_Y_RES];
	unsigned short buffer;

	// Check whether opened or not
	if (!input.is_open()) return nullptr;

	// Read data
	for (int i = 0; i < KINECT_Y_RES; ++i) {
		for (int j = 0; j < KINECT_X_RES; ++j) {
			// If it's EOF, return null
			if (!input.read(reinterpret_cast<char *>(&buffer), sizeof(unsigned short))) return nullptr;

			array[i * KINECT_X_RES + j] = buffer;
		}
	}

	return array;
}
