
#include <unordered_set>
#include "curved-triangulation.h"
#include "Sutherland-Hodgeman2.h"
#include "utils.h"

void CurvedTriangulation::optimize(int maxiter, bool optallcontrols)
{
	if (!_image || _mesh.faces_empty())
		return;
	
	double err = compute_energy();
	xlog("init status: %d vertices, %d faces, energy = %f", 
		(int)_mesh.n_vertices(), (int)_mesh.n_faces(), err);
	
	for (int it = 1; it <= 5; ++it)
	{
		optimize_edges();
		optimize_controls(optallcontrols);
	}

	for (int it = 1; it <= maxiter; ++it)
	{
		optimize_vertices();
		process_boundary(true);
#ifdef _DEBUG
		svg_mesh("iter-mesh-vopt.svg");
#endif
		optimize_edges();
#ifdef _DEBUG
		svg_mesh("iter-mesh-eopt.svg");
#endif
		optimize_controls(optallcontrols);
#ifdef _DEBUG
		svg_mesh("iter-mesh-copt.svg");
#endif

		err = compute_energy();
		xlog("iter = %d, %d vertices, %d faces, energy = %f",
			it, (int)_mesh.n_vertices(), (int)_mesh.n_faces(), err);
	}
}

void CurvedTriangulation::optimize_vertices()
{
	xlog_debug("begin");
	int vnb = (int)_mesh.n_vertices();
		
	for (int v = 4; v < vnb; ++v)
	{
		VertexHandle vh = _mesh.vertex_handle(v);
		
		double vgrad[2];
		compute_vertex_gradient(vh, vgrad);		
		optimize_vertex(vh, vgrad);
#ifdef _DEBUG
		svg_mesh("iter-mesh-vopt.svg");
#endif
	}
	xlog_debug("end");
}

#define VERTEX_OPTIMIZATION_MAXANGLE (M_PI * 90.0 / 180.0)
void CurvedTriangulation::optimize_vertex(const VertexHandle &vh, const double *vgrad)
{
	double norm = vgrad[0] * vgrad[0] + vgrad[1] * vgrad[1];
	if (norm < 1e-16)
		return;

	Point vp = _mesh.point(vh);
	Normal dir(-vgrad[0], -vgrad[1], 0);

	// too thin-long
	for (auto vohit = _mesh.voh_begin(vh); vohit != _mesh.voh_end(vh); ++vohit)
	{
		if (_mesh.is_boundary(*vohit))
			continue;

		HalfedgeHandle heh = _mesh.next_halfedge_handle(*vohit);
		EdgeHandle eh = _mesh.edge_handle(heh);
		if (!_mesh.property(EPropFeature, eh))
			continue;

		VertexHandle vhb = _mesh.from_vertex_handle(heh);
		VertexHandle vhc = _mesh.to_vertex_handle(heh);

		const Point& vpb = _mesh.point(vhb);
		const Point& vpc = _mesh.point(vhc);

		Normal dira = vpb - vp;
		Normal dirb = vpc - vp;

		double across = dira[0] * dir[1] - dira[1] * dir[0];
		double bcross = dir[0] * dirb[1] - dir[1] * dirb[0];
		if (across * bcross < 0)
			continue;

		heh = _mesh.next_halfedge_handle(heh);
		double angle = _mesh.calc_sector_angle(heh);
		if (angle > VERTEX_OPTIMIZATION_MAXANGLE)
			return;

		break;
	}

	std::vector<double> region;
	vertex_constraint(vh, region);
	if (region.size() < 6)
		return;

#ifdef _DEBUG
	svg_polygon(&vp[0], &region[0], (int)region.size() / 2, "vertex-constraint.svg");
#endif

	bool in = inside_convex(&region[0], (int)region.size() / 2, &vp[0]);
	if (!in)
		return;
	
	norm = std::sqrt(norm);
	dir[0] /= norm;
	dir[1] /= norm;

	double safestep = compute_safestep(&region[0], (int)region.size() / 2, &vp[0], &dir[0]);
	//xlog_debug("v = %d, safestep = %.10f", vh.idx(), safestep);
	if (safestep > 1)
	{		
		return;
	}

	double oldError = 0.0;
	for (auto vfit = _mesh.cvf_begin(vh); vfit != _mesh.cvf_end(vh); ++vfit)
	{
		oldError += _mesh.property(FPropData, *vfit).approxEnergy;
	}
	if (_regweight > 1e-10)
	{
		double vreg = compute_regenergy(vh);
		oldError += _regweight * vreg;
	}

	int count = 6;
	bool success = false;

	double step = safestep * 0.5;
	while (--count)
	{
		Point newp = vp + step * dir;

		double newError = move_vertex(vh, newp);
		if (_regweight > 1e-10)
		{
			double vreg = compute_regenergy(vh);
			newError += _regweight * vreg;
		}

		if (newError < oldError)
		{
			success = true;
			break;
		}

		step *= 0.2;
	}
	
	if (!success)
	{
		double scale = double(rand()) / RAND_MAX;
		step = _pixwidth * 0.2 * scale;
		if (step > safestep)
			step = safestep * 0.2;

		Point newp = vp + step * dir;
		move_vertex(vh, newp);
		//move_vertex(vh, vp);
	}

#ifdef _DEBUG
	if (check_self_intersection() > 0)
		xlog("v = %d", vh.idx());
#endif
}

// maybe concave
void CurvedTriangulation::vertex_constraint(const VertexHandle &vh, std::vector<double> &region) const
{
	region.clear();

	const Point& vp = _mesh.point(vh);
	const double radius = 1.0;

	region.push_back(vp[0] - radius);
	region.push_back(vp[1] - radius);
	region.push_back(vp[0] + radius);
	region.push_back(vp[1] - radius);
	region.push_back(vp[0] + radius);
	region.push_back(vp[1] + radius);
	region.push_back(vp[0] - radius);
	region.push_back(vp[1] + radius);

	std::vector<double> tempRegion;
	std::vector<double> *ping = &region, *pong = &tempRegion;

	std::vector<double> planes;

	for (auto vohit = _mesh.cvoh_begin(vh); vohit != _mesh.cvoh_end(vh); ++vohit)
	{
		EdgeHandle eh = _mesh.edge_handle(*vohit);
		const Point& econ = _mesh.property(EPropControl, eh);

		if (!_mesh.is_boundary(*vohit))
		{
			HalfedgeHandle nextheh = _mesh.next_halfedge_handle(*vohit);
			VertexHandle vhb = _mesh.from_vertex_handle(nextheh);
			VertexHandle vhc = _mesh.to_vertex_handle(nextheh);

			const Point& vpb = _mesh.point(vhb);
			const Point& vpc = _mesh.point(vhc);

			Normal dir = vpc - econ;

			planes.push_back(econ[0]);
			planes.push_back(econ[1]);
			planes.push_back(dir[1]);
			planes.push_back(-dir[0]);

			EdgeHandle nexteh = _mesh.edge_handle(nextheh);
			const Point& nextecon = _mesh.property(EPropControl, nexteh);

			Normal dira = vpc - vpb;
			Normal dirb = nextecon - vpb;
			double cross = dira[0] * dirb[1] - dira[1] * dirb[0];
			if (cross > 0)
			{
				planes.push_back(nextecon[0]);
				planes.push_back(nextecon[1]);
				planes.push_back(dirb[1]);
				planes.push_back(-dirb[0]);

				dirb = vpc - nextecon;
				planes.push_back(nextecon[0]);
				planes.push_back(nextecon[1]);
				planes.push_back(dirb[1]);
				planes.push_back(-dirb[0]);
			}

			HalfedgeHandle prevheh = _mesh.prev_halfedge_handle(*vohit);
			EdgeHandle preveh = _mesh.edge_handle(prevheh);
			const Point& prevecon = _mesh.property(EPropControl, preveh);

			dir = prevecon - econ;

			planes.push_back(econ[0]);
			planes.push_back(econ[1]);
			planes.push_back(dir[1]);
			planes.push_back(-dir[0]);
		}
		
		HalfedgeHandle oppheh = _mesh.opposite_halfedge_handle(*vohit);
		if (!_mesh.is_boundary(oppheh))
		{
			HalfedgeHandle tempheh = _mesh.next_halfedge_handle(oppheh);
			VertexHandle tempvh = _mesh.to_vertex_handle(tempheh);
			const Point& tempvp = _mesh.point(tempvh);

			Normal dir = econ - tempvp;

			planes.push_back(econ[0]);
			planes.push_back(econ[1]);
			planes.push_back(dir[1]);
			planes.push_back(-dir[0]);
		}
	}

	int nb = (int)planes.size() / 4;
	for (int k = 0; k < nb; ++k)
	{
		Sutherland_Hodgeman2(&planes[4 * k], &planes[4 * k + 2], *ping, *pong);
		if (pong->empty())
			return;

		std::vector<double> *ptr = ping;
		ping = pong;
		pong = ptr;
	}

	if (ping != &region)
		region = *ping;
}

#define FLIPPING_MAX_ANGLE (M_PI * 150.0 / 180.0)
#define FLIPPING_MIN_ANGLE (M_PI * 30.0 / 180.0)
void CurvedTriangulation::optimize_edges()
{
	xlog_debug("begin");
	int flips = 0;
	for (auto eit = _mesh.edges_begin(); eit != _mesh.edges_end(); ++eit)
	{
		if (_mesh.property(EPropFeature, *eit))
			continue;

		if (!_mesh.is_flippable(*eit))
			continue;

		// angle test
		HalfedgeHandle heha = _mesh.halfedge_handle(*eit, 0);
		HalfedgeHandle hehb = _mesh.halfedge_handle(*eit, 1);

		HalfedgeHandle neighs[4];
		neighs[0] = _mesh.next_halfedge_handle(heha);
		neighs[1] = _mesh.next_halfedge_handle(neighs[0]);
		neighs[2] = _mesh.next_halfedge_handle(hehb);
		neighs[3] = _mesh.next_halfedge_handle(neighs[2]);

		double angleT = _mesh.calc_sector_angle(neighs[0]);
		double angleB = _mesh.calc_sector_angle(neighs[2]);
		if (angleT < FLIPPING_MIN_ANGLE && angleB < FLIPPING_MIN_ANGLE)
			continue;

		double maxangle = angleT;
		HalfedgeHandle maxheh = neighs[0];
		if (maxangle < angleB)
		{
			maxangle = angleB;
			maxheh = neighs[2];
		}

		VertexHandle maxvh = _mesh.to_vertex_handle(maxheh);
		if (_mesh.property(VPropFeature, maxvh) &&
			_mesh.property(EPropFeature, *eit))
			continue;

		double angleR = _mesh.calc_sector_angle(heha) + _mesh.calc_sector_angle(neighs[3]);
		if (angleR > FLIPPING_MAX_ANGLE && angleR > maxangle)
			continue;

		double angleL = _mesh.calc_sector_angle(neighs[1]) + _mesh.calc_sector_angle(hehb);
		if (angleL > FLIPPING_MAX_ANGLE && angleL > maxangle)
			continue;

		if (!legal_flipping(*eit)) // will self-intersect
			continue;

		++flips;

		FaceHandle fha = _mesh.face_handle(heha);
		FaceHandle fhb = _mesh.face_handle(hehb);

		if (maxangle > FLIPPING_MAX_ANGLE)
		{
			bool isfeat = _mesh.property(EPropFeature, *eit);
			VertexHandle vha = _mesh.from_vertex_handle(heha);
			VertexHandle vhb = _mesh.from_vertex_handle(hehb);
			VertexHandle maxvh = _mesh.to_vertex_handle(maxheh);
			
			flip_edge(*eit);
			compute_approx(fha);
			compute_approx(fhb);

			if (isfeat)
			{
				_mesh.property(EPropFeature, *eit) = false;

				for (auto vohit = _mesh.voh_begin(maxvh); vohit != _mesh.voh_end(maxvh); ++vohit)
				{
					if (_mesh.to_vertex_handle(*vohit) == vha ||
						_mesh.to_vertex_handle(*vohit) == vhb)
					{
						_mesh.property(EPropFeature, _mesh.edge_handle(*vohit)) = true;
					}
				}
			}

			continue;
		}

		//if (_mesh.property(EPropFeature, *eit))
		//	continue;

		Point econ = _mesh.property(EPropControl, *eit);

		MyFaceData fpa = _mesh.property(FPropData, fha);
		MyFaceData fpb = _mesh.property(FPropData, fhb);

		double oldError = fpa.approxEnergy + fpb.approxEnergy;
		if (_regweight > 1e-10)
		{
			double frega = compute_regenergy(fha);
			double fregb = compute_regenergy(fhb);
			oldError += _regweight * (frega + fregb);
		}

		flip_edge(*eit);

		compute_approx(fha);
		compute_approx(fhb);

		const MyFaceData& newfpa = _mesh.property(FPropData, fha);
		const MyFaceData& newfpb = _mesh.property(FPropData, fhb);

		double newError = newfpa.approxEnergy + newfpb.approxEnergy;
		if (_regweight > 1e-10)
		{
			double frega = compute_regenergy(fha);
			double fregb = compute_regenergy(fhb);
			newError += _regweight * (frega + fregb);
		}

		if (newError > oldError)
		{
			// flip back
			_mesh.flip(*eit);
			_mesh.property(EPropControl, *eit) = econ;
			_mesh.property(FPropData, fha) = fpb;
			_mesh.property(FPropData, fhb) = fpa;
			--flips;
		}
	}
	xlog_debug("end, flips = %d", flips);
}

bool CurvedTriangulation::legal_flipping(const EdgeHandle &eh) const
{
	HalfedgeHandle heha = _mesh.halfedge_handle(eh, 0);
	HalfedgeHandle hehb = _mesh.halfedge_handle(eh, 1);

	HalfedgeHandle neighs[4];
	neighs[0] = _mesh.next_halfedge_handle(heha);
	neighs[1] = _mesh.next_halfedge_handle(neighs[0]);
	neighs[2] = _mesh.next_halfedge_handle(hehb);
	neighs[3] = _mesh.next_halfedge_handle(neighs[2]);

	// avoid overlapping
	Point quad[4];
	for (int i = 0; i < 4; ++i)
	{
		quad[i] = _mesh.point(_mesh.from_vertex_handle(neighs[i]));
	}

	for (int i = 0; i < 4; ++i)
	{
		int next = (i + 1) % 4;

		EdgeHandle eeh = _mesh.edge_handle(neighs[i]);
		Point econ = _mesh.property(EPropControl, eeh);

		Normal dira = quad[next] - quad[i];
		Normal dirb = econ - quad[i];

		double cross = dira[0] * dirb[1] - dira[1] * dirb[0];
		if (cross <= 0.0)
			continue;

		Point emid = (quad[next] + quad[i]) * 0.5;

		Normal ab = quad[1] - emid;
		Normal ad = quad[3] - emid;
		double midcross = ab[0] * ad[1] - ab[1] * ad[0];

		ab = quad[1] - econ;
		ad = quad[3] - econ;
		double concross = ab[0] * ad[1] - ab[1] * ad[0];

		if (midcross * concross <= 0.0)
		{
			return false;
		}
	}

	return true;
}

void CurvedTriangulation::optimize_controls(bool optall)
{
	xlog_debug("begin");
	
	for (auto eit = _mesh.edges_begin(); eit != _mesh.edges_end(); ++eit)
	{
		if (_mesh.is_boundary(*eit))
			continue;

		if (!_mesh.property(EPropFeature, *eit))
			continue;

		double cgrad[2];
		compute_control_gradient(*eit, cgrad);
		optimize_control(*eit, cgrad);
	}

	if (!optall)
		return;

	for (auto eit = _mesh.edges_begin(); eit != _mesh.edges_end(); ++eit)
	{
		if (_mesh.is_boundary(*eit))
			continue;

		if (_mesh.property(EPropFeature, *eit))
			continue;

		double cgrad[2];
		compute_control_gradient(*eit, cgrad);
		optimize_control(*eit, cgrad);
	}

	xlog_debug("end");
}

void CurvedTriangulation::optimize_control(const EdgeHandle &eh, const double *cgrad)
{
	double norm = cgrad[0] * cgrad[0] + cgrad[1] * cgrad[1];
	if (norm < 1e-16)
		return;

	std::vector<double> constraint;
	control_constraint(eh, constraint);
	if (constraint.size() < 6)
		return;

	Point econ = _mesh.property(EPropControl, eh);

#ifdef _DEBUG
	svg_polygon(&econ[0], &constraint[0], (int)constraint.size() / 2, "control-constraint.svg");
#endif

	Normal dir(-cgrad[0], -cgrad[1], 0);
	norm = std::sqrt(norm);
	dir[0] /= norm;
	dir[1] /= norm;

	double safestep = compute_safestep(&constraint[0], (int)constraint.size() / 2, &econ[0], &dir[0]);
	//xlog_debug("e = %d, safestep = %.10f", eh.idx(), safestep);
	if (safestep < 0 || safestep > 1)
	{
		return;
	}

	if (!_mesh.property(EPropFeature, eh))
		safestep /= 10;

	HalfedgeHandle heha = _mesh.halfedge_handle(eh, 0);
	HalfedgeHandle hehb = _mesh.halfedge_handle(eh, 1);

	FaceHandle fha = _mesh.face_handle(heha);
	FaceHandle fhb = _mesh.face_handle(hehb);

	MyFaceData fpa = _mesh.property(FPropData, fha);
	MyFaceData fpb = _mesh.property(FPropData, fhb);

	double oldError = fpa.approxEnergy + fpb.approxEnergy;

	int count = 6;
	bool success = false;

	double step = 0.8 * safestep;
	while (--count)
	{
		Normal move = step * dir;
		_mesh.property(EPropControl, eh) = econ + move;

		compute_approx(fha);
		compute_approx(fhb);

		const MyFaceData& newfpa = _mesh.property(FPropData, fha);
		const MyFaceData& newfpb = _mesh.property(FPropData, fhb);

		double newError = newfpa.approxEnergy + newfpb.approxEnergy;

		if (newError < oldError)
		{
			success = true;
			break;
		}

		step *= 0.2;
	}

	if (!success)
	{
		double scale = double(rand()) / RAND_MAX;
		step = _pixwidth * 0.2 * scale;
		if (step > safestep)
			step = safestep * 0.2;

		Normal move = step * dir;
		_mesh.property(EPropControl, eh) = econ + move;

		//_mesh.property(EPropControl, eh) = econ;

		compute_approx(fha);
		compute_approx(fhb);
	}
}

void CurvedTriangulation::control_constraint(const EdgeHandle &eh, std::vector<double> &region) const
{
	region.clear();

	const Point& econ = _mesh.property(EPropControl, eh);
	const double radius = 1.0;

	region.push_back(econ[0] - radius);
	region.push_back(econ[1] - radius);
	region.push_back(econ[0] + radius);
	region.push_back(econ[1] - radius);
	region.push_back(econ[0] + radius);
	region.push_back(econ[1] + radius);
	region.push_back(econ[0] - radius);
	region.push_back(econ[1] + radius);

	std::vector<double> tempRegion;
	std::vector<double> *ping = &region, *pong = &tempRegion;

	HalfedgeHandle heha = _mesh.halfedge_handle(eh, 0);
	HalfedgeHandle hehb = _mesh.halfedge_handle(eh, 1);

	std::vector<double> planes;

	HalfedgeHandle tempheh = _mesh.next_halfedge_handle(heha);
	VertexHandle vha = _mesh.from_vertex_handle(tempheh);
	VertexHandle vhb = _mesh.to_vertex_handle(tempheh);
	Point vpa = _mesh.point(vha);
	Point vpb = _mesh.point(vhb);

	EdgeHandle tempeh = _mesh.edge_handle(tempheh);
	Point tempecon = _mesh.property(EPropControl, tempeh);;

	Normal dira = vpb - vpa;
	Normal dirb = tempecon - vpa;
	double cross = dira[0] * dirb[1] - dira[1] * dirb[0];

	if (cross > 0)
	{
		planes.push_back(vpa[0]);
		planes.push_back(vpa[1]);
		planes.push_back(dirb[1]);
		planes.push_back(-dirb[0]);
	}
	else
	{
		planes.push_back(vpa[0]);
		planes.push_back(vpa[1]);
		planes.push_back(dira[1]);
		planes.push_back(-dira[0]);
	}

	tempheh = _mesh.prev_halfedge_handle(heha);
	vha = _mesh.from_vertex_handle(tempheh);
	vhb = _mesh.to_vertex_handle(tempheh);
	vpa = _mesh.point(vha);
	vpb = _mesh.point(vhb);

	tempeh = _mesh.edge_handle(tempheh);
	tempecon = _mesh.property(EPropControl, tempeh);

	dira = vpa - vpb;
	dirb = tempecon - vpb;
	cross = dira[0] * dirb[1] - dira[1] * dirb[0];

	if (cross < 0)
	{
		planes.push_back(vpb[0]);
		planes.push_back(vpb[1]);
		planes.push_back(-dirb[1]);
		planes.push_back(dirb[0]);
	}
	else
	{
		planes.push_back(vpb[0]);
		planes.push_back(vpb[1]);
		planes.push_back(-dira[1]);
		planes.push_back(dira[0]);
	}

	tempheh = _mesh.next_halfedge_handle(hehb);
	vha = _mesh.from_vertex_handle(tempheh);
	vhb = _mesh.to_vertex_handle(tempheh);
	vpa = _mesh.point(vha);
	vpb = _mesh.point(vhb);

	tempeh = _mesh.edge_handle(tempheh);
	tempecon = _mesh.property(EPropControl, tempeh);

	dira = vpb - vpa;
	dirb = tempecon - vpa;
	cross = dira[0] * dirb[1] - dira[1] * dirb[0];

	if (cross > 0)
	{
		planes.push_back(vpa[0]);
		planes.push_back(vpa[1]);
		planes.push_back(dirb[1]);
		planes.push_back(-dirb[0]);
	}
	else
	{
		planes.push_back(vpa[0]);
		planes.push_back(vpa[1]);
		planes.push_back(dira[1]);
		planes.push_back(-dira[0]);
	}

	tempheh = _mesh.prev_halfedge_handle(hehb);
	vha = _mesh.from_vertex_handle(tempheh);
	vhb = _mesh.to_vertex_handle(tempheh);
	vpa = _mesh.point(vha);
	vpb = _mesh.point(vhb);

	tempeh = _mesh.edge_handle(tempheh);
	tempecon = _mesh.property(EPropControl, tempeh);

	dira = vpa - vpb;
	dirb = tempecon - vpb;
	cross = dira[0] * dirb[1] - dira[1] * dirb[0];

	if (cross < 0)
	{
		planes.push_back(vpb[0]);
		planes.push_back(vpb[1]);
		planes.push_back(-dirb[1]);
		planes.push_back(dirb[0]);
	}
	else
	{
		planes.push_back(vpb[0]);
		planes.push_back(vpb[1]);
		planes.push_back(-dira[1]);
		planes.push_back(dira[0]);
	}

	int nb = (int)planes.size() / 4;
	for (int k = 0; k < nb; ++k)
	{
		Sutherland_Hodgeman2(&planes[4 * k], &planes[4 * k + 2], *ping, *pong);
		if (pong->empty())
			return;

		std::vector<double> *ptr = ping;
		ping = pong;
		pong = ptr;
	}

	if (ping != &region)
		region = *ping;
}

