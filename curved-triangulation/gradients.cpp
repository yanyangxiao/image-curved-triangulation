
#include "curved-triangulation.h"
#include "bezier2.h"

void CurvedTriangulation::compute_vertex_gradient(const VertexHandle &vh, double *vgrad) const
{
	vgrad[0] = 0.0;
	vgrad[1] = 0.0;

	const Point& vp = _mesh.point(vh);

	// approx
	for (auto vohit = _mesh.cvoh_begin(vh); vohit != _mesh.cvoh_end(vh); ++vohit)
	{
		VertexHandle nvh = _mesh.to_vertex_handle(*vohit);
		const Point& nvp = _mesh.point(nvh);

		EdgeHandle eh = _mesh.edge_handle(*vohit);
		const Point& econ = _mesh.property(EPropControl, eh);

		FaceHandle leftFace = _mesh.face_handle(*vohit);
		const MyPolynomial *leftApprox = NULL;
		if (leftFace.idx() != -1)
			leftApprox = &(_mesh.property(FPropData, leftFace).approx);

		FaceHandle rightFace = _mesh.opposite_face_handle(*vohit);
		const MyPolynomial *rightApprox = NULL;
		if (rightFace.idx() != -1)
			rightApprox = &(_mesh.property(FPropData, rightFace).approx);

		double temp[2] = { 0.0, 0.0 };
		vertex_integral(&vp[0], &econ[0], &nvp[0], leftApprox, rightApprox, 0.0, 1.0, temp);

		vgrad[0] += temp[0];
		vgrad[1] += temp[1];
	}

	// reg
	if (_regweight > 1e-10)
	{
		double regGrad[2] = { 0.0, 0.0 };
		int count = 0;
		for (auto vvit = _mesh.cvv_begin(vh); vvit != _mesh.cvv_end(vh); ++vvit)
		{
			const Point& nvp = _mesh.point(*vvit);
			regGrad[0] += vp[0] - nvp[0];
			regGrad[1] += vp[1] - nvp[1];
			++count;
		}

		if (count > 0)
		{
			regGrad[0] /= count;
			regGrad[1] /= count;
		}

		vgrad[0] += _regweight * regGrad[0];
		vgrad[1] += _regweight * regGrad[1];
	}

	if (_mesh.is_boundary(vh))
	{
		for (auto vohit = _mesh.cvoh_begin(vh); vohit != _mesh.cvoh_end(vh); ++vohit)
		{
			if (_mesh.is_boundary(_mesh.edge_handle(*vohit)))
			{
				Normal dir = _mesh.calc_edge_vector(*vohit);
				dir.normalize();

				double dot = (dir[0] * vgrad[0] + dir[1] * vgrad[1]);
				if (dot >= 0.0)
				{
					dir *= dot;
					vgrad[0] = dir[0];
					vgrad[1] = dir[1];
					break;
				}
			}
		}
	}
}

void CurvedTriangulation::vertex_integral(
	const double *pi,
	const double *control,
	const double *pj,
	const MyPolynomial *leftApprox,
	const MyPolynomial *rightApprox,
	double t0,
	double t1,
	double *result) const
{
	double t = (t0 + t1) * 0.5;
	double dt = fabs(t1 - t0);

	double p[2], q[2];
	bezier2_point(pi, control, pj, t0, p);
	bezier2_point(pi, control, pj, t1, q);

	double dx = q[0] - p[0];
	double dy = q[1] - p[1];
	double ds = std::sqrt(dx * dx + dy * dy);
	if (ds < _pixwidth)
	{
		double sample[2], fn[2];
		bezier2_point(pi, control, pj, t, sample);

		dx = -2 * (1 - t) * pi[0] + (2 - 4 * t) * control[0] + 2 * t * pj[0];
		dy = -2 * (1 - t) * pi[1] + (2 - 4 * t) * control[1] + 2 * t * pj[1];

		fn[0] = dy;
		fn[1] = -dx;

		int i = int(sample[0] * _width);
		if (i < 0) i = 0;
		if (i >= _width) i = _width - 1;
		int j = int(sample[1] * _width);
		if (j < 0) j = 0;
		if (j >= _height) j = _height - 1;

		int pid = j * _width + i;
		const float *pcolor = &_image[_channel * pid];

		double leftError = 0.0;
		double rightError = 0.0;
		for (int k = 0; k < _channel; ++k)
		{
			double diff = leftApprox ? ((double)pcolor[k] - leftApprox->evaluate(k, sample)) : 0.0;
			leftError += diff * diff;

			diff = rightApprox ? ((double)pcolor[k] - rightApprox->evaluate(k, sample)) : 0.0;
			rightError += diff * diff;
		}

		double temp = leftError - rightError;

		temp *= ((1 - t) * (1 - t) * dt);
		result[0] += temp * fn[0];
		result[1] += temp * fn[1];
	}
	else
	{
		vertex_integral(pi, control, pj, leftApprox, rightApprox, t0, t, result);
		vertex_integral(pi, control, pj, leftApprox, rightApprox, t, t1, result);
	}
}

void CurvedTriangulation::compute_control_gradient(const EdgeHandle &eh, double *cgrad) const
{
	cgrad[0] = 0.0;
	cgrad[1] = 0.0;

	const Point& econ = _mesh.property(EPropControl, eh);

	HalfedgeHandle heha = _mesh.halfedge_handle(eh, 0);
	HalfedgeHandle hehb = _mesh.halfedge_handle(eh, 1);

	FaceHandle fha = _mesh.face_handle(heha);
	FaceHandle fhb = _mesh.face_handle(hehb);

	MyFaceData fpa = _mesh.property(FPropData, fha);
	MyFaceData fpb = _mesh.property(FPropData, fhb);

	const Point& src = _mesh.point(_mesh.from_vertex_handle(heha));
	const Point& tgt = _mesh.point(_mesh.to_vertex_handle(heha));

	control_integral(&src[0], &econ[0], &tgt[0], &fpa.approx, &fpb.approx, 0, 1, cgrad);
}

void CurvedTriangulation::control_integral(
	const double *pi,
	const double *control,
	const double *pj,
	const MyPolynomial *leftApprox,
	const MyPolynomial *rightApprox,
	double t0,
	double t1,
	double *result) const
{
	double t = (t0 + t1) * 0.5;
	double dt = fabs(t1 - t0);

	double p[2], q[2];
	bezier2_point(pi, control, pj, t0, p);
	bezier2_point(pi, control, pj, t1, q);

	double dx = q[0] - p[0];
	double dy = q[1] - p[1];
	double ds = std::sqrt(dx * dx + dy * dy);
	if (ds < _pixwidth)
	{
		double sample[2], fn[2];
		bezier2_point(pi, control, pj, t, sample);

		dx = -2 * (1 - t) * pi[0] + (2 - 4 * t) * control[0] + 2 * t * pj[0];
		dy = -2 * (1 - t) * pi[1] + (2 - 4 * t) * control[1] + 2 * t * pj[1];

		fn[0] = dy;
		fn[1] = -dx;

		int i = int(sample[0] * _width);
		if (i < 0) i = 0;
		if (i >= _width) i = _width - 1;
		int j = int(sample[1] * _width);
		if (j < 0) j = 0;
		if (j >= _height) j = _height - 1;

		int pid = j * _width + i;
		const float *pcolor = &_image[_channel * pid];

		double leftError = 0.0;
		double rightError = 0.0;
		for (int k = 0; k < _channel; ++k)
		{
			double diff = leftApprox ? ((double)pcolor[k] - leftApprox->evaluate(k, sample)) : 0.0;
			leftError += diff * diff;

			diff = rightApprox ? ((double)pcolor[k] - rightApprox->evaluate(k, sample)) : 0.0;
			rightError += diff * diff;
		}

		double temp = leftError - rightError;

		temp *= (2 * t * (1 - t) * dt);
		result[0] += temp * fn[0];
		result[1] += temp * fn[1];
	}
	else
	{
		control_integral(pi, control, pj, leftApprox, rightApprox, t0, t, result);
		control_integral(pi, control, pj, leftApprox, rightApprox, t, t1, result);
	}
}

