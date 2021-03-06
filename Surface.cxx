#include <vtkDelaunay2D.h>
#include <vtkDelaunay3D.h>
#include <vtkXMLDataSetWriter.h>
#include <thread>
#include "vtkXMLPolyDataWriter.h"
#include "vtkActor.h"
#include "vtkCamera.h"
#include "vtkFloatArray.h"
#include "vtkPointData.h"
#include "vtkPolyDataMapper.h"
#include "vtkRenderWindow.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkRenderer.h"
#include "vtkProperty.h"
#include "vtkPLYWriter.h"
#include "vtkGeometryFilter.h"
#include "RawDataReader.h"
#include "TetGen.h"

#define SCALE 16
#define KINECT_X_RES 640
#define KINECT_Y_RES 480

static const int MESH_X_RES = KINECT_X_RES / SCALE;
static const int MESH_Y_RES = KINECT_Y_RES / SCALE;

float *getDepthData();

void CallNextFrames();

vtkSmartPointer<vtkFloatArray> getVertices(float *depthData);

RawDataReader *rawDataReader;
vtkPolyData *surface;
vtkRenderWindow *renWin;


int main() {
	int i;

	// Get vertices and triangles from depth data
	rawDataReader = new RawDataReader();
	rawDataReader->SetFileName("kinectsdk_depth.data", 63);
	float *depthData = getDepthData();
//	float *depthData = rawDataReader->ReadNextFrame();
	vtkSmartPointer<vtkFloatArray> vertices = getVertices(depthData);


	// Create the building blocks of polydata
	surface = vtkPolyData::New();
	vtkPoints *points = vtkPoints::New();
	vtkFloatArray *scalars = vtkFloatArray::New();


	// Load the point, cell, and data attributes.
	points->SetData(vertices);
	for (i = 0; i < MESH_X_RES * MESH_Y_RES; i++) scalars->InsertTuple1(i, i);


	// Assign the pieces to the vtkPolyData.
	surface->SetPoints(points);
	surface->GetPointData()->SetScalars(scalars);
	scalars->Delete();


	// Perform a 2D Delaunay triangulation on them.
	vtkDelaunay2D *delaunay = vtkDelaunay2D::New();
	delaunay->SetInputData(surface);
	delaunay->SetTolerance(0.001);


	// Now we'll look at it.
	vtkPolyDataMapper *PolyMapper = vtkPolyDataMapper::New();
	PolyMapper->SetInputConnection(delaunay->GetOutputPort());
	PolyMapper->SetScalarRange(0, MESH_X_RES * MESH_Y_RES - 1);

	vtkActor *surfaceActor = vtkActor::New();
	surfaceActor->SetMapper(PolyMapper);
	surfaceActor->GetProperty()->SetRepresentationToWireframe();

	// Write the file
	vtkSmartPointer<vtkXMLPolyDataWriter> writer =
			vtkSmartPointer<vtkXMLPolyDataWriter>::New();
	writer->SetFileName("test.vtp");
#if VTK_MAJOR_VERSION <= 5
	writer->SetInput(polydata);
#else
	writer->SetInputConnection(delaunay->GetOutputPort());
#endif

	// Optional - set the mode. The default is binary.
	//writer->SetDataModeToBinary();
	//writer->SetDataModeToAscii();
	writer->Write();


	// Generate a tetrahedral mesh from the input points. By
	// default, the generated volume is the convex hull of the points.
	vtkSmartPointer<vtkDelaunay3D> delaunay3D =
			vtkSmartPointer<vtkDelaunay3D>::New();
	delaunay3D->SetTolerance(0.01);
	delaunay3D->BoundingTriangulationOff();
	delaunay3D->SetInputConnection(delaunay->GetOutputPort());

//	// Write the mesh as an unstructured grid
//	vtkSmartPointer<vtkXMLDataSetWriter> writer3D =
//			vtkSmartPointer<vtkXMLDataSetWriter>::New();
//	writer3D->SetFileName("test.vtu");
//	writer3D->SetInputConnection(delaunay3D->GetOutputPort());
//	writer3D->Write();


	// Convert UnstructuredGrid to Polydata to create .ply file.
	vtkSmartPointer<vtkGeometryFilter> geometryFilter =
			vtkSmartPointer<vtkGeometryFilter>::New();
#if VTK_MAJOR_VERSION <= 5
	geometryFilter->SetInput(delaunay3D->GetOutputPort());
#else
	geometryFilter->SetInputConnection(delaunay3D->GetOutputPort());
#endif
	geometryFilter->Update();

	// Write PLY
	vtkPolyData* filteredData = geometryFilter->GetOutput();
	vtkSmartPointer<vtkPLYWriter> plyWriter = vtkSmartPointer<vtkPLYWriter>::New();
	plyWriter->SetFileName("test.ply");
	plyWriter->SetFileTypeToASCII(); // tetgen only can read ASCII
	plyWriter->SetInputData(filteredData);
	plyWriter->Write();

	TetGen tetGen;
	tetGen.tetrahedralize(filteredData->GetPoints(), filteredData->GetPolys());


	// The usual rendering stuff.
	vtkCamera *camera = vtkCamera::New();
	camera->SetPosition(1, 1, 1);
	camera->SetFocalPoint(0, 0, 0);

	vtkRenderer *renderer = vtkRenderer::New();
	renWin = vtkRenderWindow::New();
	renWin->AddRenderer(renderer);

	vtkRenderWindowInteractor *iren = vtkRenderWindowInteractor::New();
	iren->SetRenderWindow(renWin);

	renderer->AddActor(surfaceActor);
	renderer->SetActiveCamera(camera);
	renderer->ResetCamera();
	renderer->SetBackground(1, 1, 1);

	renWin->SetSize(300, 300);

	// interact with data
	renWin->Render();
	iren->Initialize();

	// Play the frame
//	CallNextFrames();

	// Stopped
	iren->Start();

	// Clean up
	delete (depthData);
	points->Delete();
	delete rawDataReader;
	surface->Delete();
	PolyMapper->Delete();
	surfaceActor->Delete();
	camera->Delete();
	renderer->Delete();
	renWin->Delete();
	iren->Delete();

	return 0;
}


// Building Sample data of object
float *getDepthData() {
	const float PLANEDEPTH = 50000;
	float *depthData = new float[KINECT_X_RES * KINECT_Y_RES];
	float depth;

	for (int i = 0; i < KINECT_X_RES; ++i)
		for (int j = 0; j < KINECT_Y_RES; ++j) {
			depth = 0;
			if (i > KINECT_X_RES / 6 && i < 5 * KINECT_X_RES / 6 && j > KINECT_Y_RES / 6 && j < 5 * KINECT_Y_RES / 6) {

				int x2 = (KINECT_X_RES / 2 - i) * (KINECT_X_RES / 2 - i);
				int y2 = (KINECT_Y_RES / 2 - j) * (KINECT_Y_RES / 2 - j);
				depth = PLANEDEPTH - x2 / 3 - y2 / 3;
			}

			depthData[i * KINECT_Y_RES + j] = depth;
		}

	return depthData;
}


// Locate vertices by depth ( x , y , depth )
vtkSmartPointer<vtkFloatArray> getVertices(float *depthData) {
	vtkSmartPointer<vtkFloatArray> points = vtkSmartPointer<vtkFloatArray>::New();
	points->SetNumberOfComponents(3);

	for (int i = 0; i < MESH_X_RES; ++i)
		for (int j = 0; j < MESH_Y_RES; ++j) {
			if (depthData[i * SCALE * KINECT_Y_RES + j * SCALE] < 200) continue;
			if (depthData[i * SCALE * KINECT_Y_RES + j * SCALE] > 50000) continue;
			points->InsertNextValue(i * SCALE * 100);
			points->InsertNextValue(j * SCALE * 100);
			points->InsertNextValue(depthData[i * SCALE * KINECT_Y_RES + j * SCALE]);
		}

	return points;
}

// Render every frames until the end
void CallNextFrames() {
	while (true) {
		std::this_thread::sleep_for(std::chrono::milliseconds(10));
		float *nextFrameDepth = rawDataReader->ReadNextFrame();

		// If it's finished
		if (nextFrameDepth == nullptr) return;

		// Reset points
		surface->GetPoints()->SetData(getVertices(nextFrameDepth));

		renWin->Render();
	}
}
