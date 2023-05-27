
#include "curved-triangulation.h"
#include "bezier2.h"

#define CURVESEGMENT_NUMBER 10
void CurvedTriangulation::meshface_polygon(
	const FaceHandle &fh, 
	double shrink, 
	std::vector<double> &polygon) const
{
	polygon.clear();

	HalfedgeHandle fhe = _mesh.halfedge_handle(fh);

	Point pts[3];
	pts[0] = _mesh.point(_mesh.from_vertex_handle(fhe));
	pts[1] = _mesh.point(_mesh.to_vertex_handle(fhe));
	pts[2] = _mesh.point(_mesh.to_vertex_handle(_mesh.next_halfedge_handle(fhe)));

	Point fcent = (pts[0] + pts[1] + pts[2]) / 3;
	for (int k = 0; k < 3; ++k)
	{
		pts[k] = fcent + (pts[k] - fcent) * shrink;
	}

	polygon.resize(CURVESEGMENT_NUMBER * 3 * 2);

	for (int k = 0; k < 3; ++k)
	{
		int next = (k + 1) % 3;

		EdgeHandle eh = _mesh.edge_handle(fhe);
		const Point& econ = _mesh.property(EPropControl, eh);

		for (int s = 0; s < CURVESEGMENT_NUMBER; ++s)
		{
			double t = double(s) / CURVESEGMENT_NUMBER;
			bezier2_point(&pts[k][0], &econ[0], &pts[next][0], t, &polygon[2 * (k * CURVESEGMENT_NUMBER + s)]);
		}

		fhe = _mesh.next_halfedge_handle(fhe);
	}
}

#define BOUNDARY_ANGLE_THRESHOLD (M_PI * 160.0 / 180.0)
int CurvedTriangulation::process_boundary(bool updateApprox)
{
	int vnb = (int)_mesh.n_vertices();
	if (vnb < 3)
		return 0;

	int count = 0;
	std::vector<bool> visited(vnb, false);

	int bv = 0;
	while (1)
	{
		if (visited[bv])
			break;

		visited[bv] = true;

		VertexHandle bvh = _mesh.vertex_handle(bv);
		HalfedgeHandle *bhe = NULL;
		for (auto vohit = _mesh.voh_begin(bvh); vohit != _mesh.voh_end(bvh); ++vohit)
		{
			if (_mesh.is_boundary(*vohit))
			{
				bhe = &(*vohit);
				break;
			}
		}

		if (!bhe)
			break;

		HalfedgeHandle heh = _mesh.next_halfedge_handle(_mesh.opposite_halfedge_handle(*bhe));
		VertexHandle vh = _mesh.to_vertex_handle(heh);

		double threshold = BOUNDARY_ANGLE_THRESHOLD;
		if (_mesh.property(VPropFeature, vh))
			threshold = M_PI * 175.0 / 180.0;

		double angle = _mesh.calc_sector_angle(heh);
		if (angle > threshold)
		{
			bool editable = true;
			
			Normal dir = _mesh.calc_edge_vector(*bhe);
			dir.normalize();

			for (auto vohit = _mesh.voh_begin(bvh); vohit != _mesh.voh_end(bvh); ++vohit)
			{
				HalfedgeHandle
					temp = _mesh.opposite_halfedge_handle(_mesh.next_halfedge_handle(*vohit));
				if (temp != *bhe && _mesh.is_boundary(temp))
				{
					Normal dir2 = _mesh.calc_edge_vector(temp);
					dir2.normalize();
					if (OpenMesh::dot(dir, dir2) > 0.5)
					{
						editable = false;
						break;
					}
				}
			}

			if (editable)
			{
				double len = _mesh.calc_edge_length(*bhe);
				double area = _mesh.calc_sector_area(heh);
				double height = area * 2.0 / len;

				Normal perp(-dir[1], dir[0], 0.0);
				perp.normalize();

				FaceHandle fh = _mesh.opposite_face_handle(*bhe);

				_mesh.delete_face(fh);
				_mesh.garbage_collection();
				_mesh.set_point(vh, _mesh.point(vh) + height * perp);

				for (auto vohit = _mesh.voh_begin(vh); vohit != _mesh.voh_end(vh); ++vohit)
				{
					EdgeHandle eh = _mesh.edge_handle(*vohit);
					_mesh.property(EPropControl, eh) = _mesh.calc_edge_midpoint(eh);

					if (_mesh.is_boundary(*vohit))
						continue;

					if (updateApprox)
					{
						FaceHandle fh = _mesh.face_handle(*vohit);
						compute_approx(fh);
					}
				}

				visited[bv] = false;

				++count;
			}
		}
		else
		{
			bv = _mesh.to_vertex_handle(*bhe).idx();
		}
	}

	return count;
}

double CurvedTriangulation::move_vertex(const VertexHandle &vh, const Point &p)
{
	_mesh.set_point(vh, p);

	double newError = 0.0;
	for (auto vohit = _mesh.voh_begin(vh); vohit != _mesh.voh_end(vh); ++vohit)
	{
		EdgeHandle teh = _mesh.edge_handle(*vohit);
		if (_mesh.is_boundary(teh))
			_mesh.property(EPropControl, teh) = _mesh.calc_edge_midpoint(teh);

		if (_mesh.is_boundary(*vohit))
			continue;

		FaceHandle tfh = _mesh.face_handle(*vohit);
		compute_approx(tfh);

		newError += _mesh.property(FPropData, tfh).approxEnergy;
	}

	return newError;
}

void CurvedTriangulation::flip_edge(const EdgeHandle &eh)
{
	_mesh.flip(eh);
	_mesh.property(EPropControl, eh) = _mesh.calc_edge_midpoint(eh);
}

void CurvedTriangulation::delaunay(const VertexHandle &vh, std::unordered_set<int> &fset)
{
	fset.clear();

	std::list<EdgeHandle> edgeList;
	std::vector<bool> inList(_mesh.n_edges(), false);

	for (auto voheit = _mesh.voh_begin(vh); voheit != _mesh.voh_end(vh); ++voheit)
	{
		EdgeHandle eh = _mesh.edge_handle(*voheit);

		if (!_mesh.property(EPropFeature, eh) && !_mesh.is_delaunay(eh))
		{
			edgeList.push_back(eh);
			inList[eh.idx()] = true;
		}

		eh = _mesh.edge_handle(_mesh.next_halfedge_handle(*voheit));
		if (!_mesh.property(EPropFeature, eh) && !_mesh.is_delaunay(eh))
		{
			edgeList.push_back(eh);
			inList[eh.idx()] = true;
		}
	}

	while (!edgeList.empty())
	{
		EdgeHandle eh = edgeList.front();
		edgeList.pop_front();
		inList[eh.idx()] = false;

		// self-defined flip
		flip_edge(eh);

		HalfedgeHandle heh = _mesh.halfedge_handle(eh, 0);

		// faces required to update
		fset.insert(_mesh.face_handle(heh).idx());
		fset.insert(_mesh.opposite_face_handle(heh).idx());

		EdgeHandle temp = _mesh.edge_handle(_mesh.next_halfedge_handle(heh));

		if (!inList[temp.idx()] && !_mesh.property(EPropFeature, temp) && !_mesh.is_delaunay(temp))
		{
			edgeList.push_front(temp);
			inList[temp.idx()] = true;
		}

		temp = _mesh.edge_handle(_mesh.next_halfedge_handle(_mesh.next_halfedge_handle(heh)));
		if (!inList[temp.idx()] && !_mesh.property(EPropFeature, temp) && !_mesh.is_delaunay(temp))
		{
			edgeList.push_front(temp);
			inList[temp.idx()] = true;
		}

		temp = _mesh.edge_handle(_mesh.next_halfedge_handle(_mesh.opposite_halfedge_handle(heh)));
		if (!inList[temp.idx()] && !_mesh.property(EPropFeature, temp) && !_mesh.is_delaunay(temp))
		{
			edgeList.push_front(temp);
			inList[temp.idx()] = true;
		}

		temp = _mesh.edge_handle(_mesh.next_halfedge_handle(_mesh.next_halfedge_handle(_mesh.opposite_halfedge_handle(heh))));
		if (!inList[temp.idx()] && !_mesh.property(EPropFeature, temp) && !_mesh.is_delaunay(temp))
		{
			edgeList.push_front(temp);
			inList[temp.idx()] = true;
		}
	}
}

bool CurvedTriangulation::self_intersect(const EdgeHandle &eh) const
{
	if (_mesh.is_boundary(eh))
		return false;

	const Point& econ = _mesh.property(EPropControl, eh);

	HalfedgeHandle heh = _mesh.halfedge_handle(eh, 0);
	VertexHandle vha = _mesh.from_vertex_handle(heh);
	VertexHandle vhb = _mesh.to_vertex_handle(heh);

	const Point& vpa = _mesh.point(vha);
	const Point& vpb = _mesh.point(vhb);

	Normal base = vpb - vpa;
	Normal dira = econ - vpa;
	Normal dirb = econ - vpb;

	double across = base[0] * dira[1] - base[1] * dira[0];
	if (abs(across) < 1e-16)
		return false;
	else if (across > 0)
	{
		// left-left
		HalfedgeHandle tempheh = _mesh.prev_halfedge_handle(heh);
		VertexHandle vhc = _mesh.from_vertex_handle(tempheh);
		const Point& vpc = _mesh.point(vhc);

		EdgeHandle tempeh = _mesh.edge_handle(tempheh);
		Point tempecon = _mesh.property(EPropControl, tempeh);

		Normal veca = vpc - vpa;
		Normal vecb = tempecon - vpa;

		double bcross = veca[0] * vecb[1] - veca[1] * vecb[0];
		if (bcross >= 0)
		{
			double ccross = dira[0] * veca[1] - dira[1] * veca[0];
			if (ccross < -1e-16)
				return true;
		}
		else
		{
			double ccross = dira[0] * vecb[1] - dira[1] * vecb[0];
			if (ccross < -1e-16)
				return true;
		}

		// right-left
		tempheh = _mesh.next_halfedge_handle(heh);
		tempeh = _mesh.edge_handle(tempheh);
		tempecon = _mesh.property(EPropControl, tempeh);

		veca = vpc - vpb;
		vecb = tempecon - vpb;

		bcross = veca[0] * vecb[1] - veca[1] * vecb[0];
		if (bcross <= 0)
		{
			double ccross = dirb[0] * veca[1] - dirb[1] * veca[0];
			if (ccross > 1e-16)
				return true;
		}
		else
		{
			double ccross = dirb[0] * vecb[1] - dirb[1] * vecb[0];
			if (ccross > 1e-16)
				return true;
		}
	}
	else
	{
		// left-right
		HalfedgeHandle tempheh = _mesh.next_halfedge_handle(_mesh.opposite_halfedge_handle(heh));
		VertexHandle vhd = _mesh.to_vertex_handle(tempheh);
		const Point& vpd = _mesh.point(vhd);

		EdgeHandle tempeh = _mesh.edge_handle(tempheh);
		Point tempecon = _mesh.property(EPropControl, tempeh);

		Normal veca = vpd - vpa;
		Normal vecb = tempecon - vpa;

		double bcross = veca[0] * vecb[1] - veca[1] * vecb[0];
		if (bcross <= 0)
		{
			double ccross = dira[0] * veca[1] - dira[1] * veca[0];
			if (ccross > 1e-16)
				return true;
		}
		else
		{
			double ccross = dira[0] * vecb[1] - dira[1] * vecb[0];
			if (ccross > 1e-16)
				return true;
		}

		// right-right
		tempheh = _mesh.prev_halfedge_handle(_mesh.opposite_halfedge_handle(heh));
		tempeh = _mesh.edge_handle(tempheh);
		tempecon = _mesh.property(EPropControl, tempeh);

		veca = vpd - vpb;
		vecb = tempecon - vpb;

		bcross = veca[0] * vecb[1] - veca[1] * vecb[0];
		if (bcross >= 0)
		{
			double ccross = dirb[0] * veca[1] - dirb[1] * veca[0];
			if (ccross < -1e-16)
				return true;
		}
		else
		{
			double ccross = dirb[0] * vecb[1] - dirb[1] * vecb[0];
			if (ccross < -1e-16)
				return true;
		}
	}

	return false;
}

int CurvedTriangulation::check_self_intersection() const
{
	int count = 0;
	for (auto eit = _mesh.edges_begin(); eit != _mesh.edges_end(); ++eit)
	{
		if (self_intersect(*eit))
			++count;
	}

	return count;
}

bool CurvedTriangulation::self_intersect(const VertexHandle &vh) const
{
	for (auto vohit = _mesh.cvoh_begin(vh); vohit != _mesh.cvoh_end(vh); ++vohit)
	{
		EdgeHandle eh = _mesh.edge_handle(*vohit);
		if (self_intersect(eh))
			return true;
	}
	return false;
}
