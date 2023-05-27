
#include "curved-triangulation.h"
#include "utils.h"

void CurvedTriangulation::compute_approx(double shrink)
{
	int fnb = (int)_mesh.n_faces();

#pragma omp parallel for
	for (int f = 0; f < fnb; ++f)
	{
		FaceHandle fh = _mesh.face_handle(f);
		compute_approx(fh, shrink);
	}
}

void CurvedTriangulation::compute_approx(const std::vector<int> &faces, double shrink)
{
	int fnb = (int)faces.size();

#pragma omp parallel for
	for (int k = 0; k < fnb; ++k)
	{
		int f = faces[k];
		FaceHandle fh = _mesh.face_handle(f);
		compute_approx(fh, shrink);
	}
}

void CurvedTriangulation::compute_approx(const FaceHandle &fh, double shrink)
{
	std::vector<double> polygon;
	meshface_polygon(fh, shrink, polygon);

	PixelSet pixels;
	MyRasterizer ras(_width, _height);
	ras.compute(&polygon[0], (int)polygon.size() / 2, pixels);
	
	MyFaceData *fp = &(_mesh.property(FPropData, fh));
	fp->approx.set_degree(_degree);
	fp->approx.compute_factors(_image, _width, _height, _channel, &pixels);
	fp->approxEnergy = fp->approx.compute_energy(_image, _width, _height, _channel, &pixels);
}
