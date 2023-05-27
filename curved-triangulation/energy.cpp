
#include "curved-triangulation.h"

double CurvedTriangulation::compute_regenergy(const VertexHandle &vh) const
{
	double vreg = 0.0;
	int count = 0;

	const Point& vp = _mesh.point(vh);

	for (auto vvit = _mesh.cvv_begin(vh); vvit != _mesh.cvv_end(vh); ++vvit)
	{
		Normal temp = _mesh.point(*vvit) - vp;
		vreg += OpenMesh::dot(temp, temp);
		++count;
	}

	if (count > 0)
		vreg /= count;

	return 0.5 * vreg;
}

double CurvedTriangulation::compute_regenergy(const FaceHandle &fh) const
{
	double freg = 0.0;

	HalfedgeHandle fhe = _mesh.halfedge_handle(fh);
	EdgeHandle eh = _mesh.edge_handle(fhe);
	double sqrlen = _mesh.calc_edge_sqr_length(eh);
	
	if (_mesh.is_boundary(eh))
		freg += sqrlen;
	else
		freg += 0.5 * sqrlen;

	fhe = _mesh.next_halfedge_handle(fhe);
	eh = _mesh.edge_handle(fhe);
	sqrlen = _mesh.calc_edge_sqr_length(eh);

	if (_mesh.is_boundary(eh))
		freg += sqrlen;
	else
		freg += 0.5 * sqrlen;

	fhe = _mesh.next_halfedge_handle(fhe);
	eh = _mesh.edge_handle(fhe);
	sqrlen = _mesh.calc_edge_sqr_length(eh);

	if (_mesh.is_boundary(eh))
		freg += sqrlen;
	else
		freg += 0.5 * sqrlen;

	return freg;
}

double CurvedTriangulation::compute_regenergy() const
{
	double sum = 0.0;

	int vnb = (int)_mesh.n_vertices();

#pragma omp parallel for reduction(+:sum)
	for (int v = 0; v < vnb; ++v)
	{
		VertexHandle vh = _mesh.vertex_handle(v);
		sum += compute_regenergy(vh);
	}

	return sum;
}

double CurvedTriangulation::compute_approx_energy() const
{
	double sum = 0.0;

	int fnb = (int)_mesh.n_faces();

#pragma omp parallel for reduction(+:sum)
	for (int f = 0; f < fnb; ++f)
	{
		FaceHandle fh = _mesh.face_handle(f);
		sum += _mesh.property(FPropData, fh).approxEnergy;
	}

	return sum;
}

double CurvedTriangulation::compute_energy() const
{
	double approxEnergy = compute_approx_energy();
	double regEnergy = compute_regenergy();

	return approxEnergy + _regweight * regEnergy;
}
