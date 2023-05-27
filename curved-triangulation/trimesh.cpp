
#include <list>
#include "trimesh.h"

void Trimesh::face_triangle(const FaceHandle &fh, double *tri) const
{
	HalfedgeHandle heh = halfedge_handle(fh);
	const Point& pa = point(from_vertex_handle(heh));
	const Point& pb = point(to_vertex_handle(heh));

	heh = next_halfedge_handle(heh);
	const Point& pc = point(to_vertex_handle(heh));

	tri[0] = pa[0];
	tri[1] = pa[1];
	tri[2] = pb[0];
	tri[3] = pb[1];
	tri[4] = pc[0];
	tri[5] = pc[1];
}

bool Trimesh::is_flippable(const EdgeHandle &eh) const
{
	if (is_boundary(eh) || !is_flip_ok(eh))
		return false;

	HalfedgeHandle heh = halfedge_handle(eh, 0);
	Normal v1v2 = -calc_edge_vector(prev_halfedge_handle(heh));
	Normal v1v4 = calc_edge_vector(next_halfedge_handle(opposite_halfedge_handle(heh)));
	Normal v3v4 = -calc_edge_vector(prev_halfedge_handle(opposite_halfedge_handle(heh)));
	Normal v3v2 = calc_edge_vector(next_halfedge_handle(heh));

	double cross1 = v1v2[0] * v1v4[1] - v1v4[0] * v1v2[1];
	double cross2 = v3v4[0] * v3v2[1] - v3v2[0] * v3v4[1];

	return (cross1 * cross2 > 0.0);
}

bool Trimesh::is_delaunay(const EdgeHandle &eh) const
{
	if (!is_flippable(eh))
		return true;

	HalfedgeHandle heh = halfedge_handle(eh, 0);
	Point p1 = point(to_vertex_handle(heh));
	Point p2 = point(to_vertex_handle(next_halfedge_handle(heh)));
	Point p3 = point(from_vertex_handle(heh));
	Point p4 = point(to_vertex_handle(next_halfedge_handle(opposite_halfedge_handle(heh))));

	Normal p2p3 = p3 - p2;
	Normal p2p1 = p1 - p2;
	Normal p4p1 = p1 - p4;
	Normal p4p3 = p3 - p4;

	double len1 = p2p3.length();
	if (len1 < DBL_EPSILON)
		return true;

	double len2 = p2p1.length();
	if (len2 < DBL_EPSILON)
		return true;

	double len3 = p4p1.length();
	if (len3 < DBL_EPSILON)
		return true;

	double len4 = p4p3.length();
	if (len4 < DBL_EPSILON)
		return true;

	double dot1 = OpenMesh::dot(p2p3, p2p1);
	double dot2 = OpenMesh::dot(p4p1, p4p3);

	double angle1 = std::acos(dot1 / (len1 * len2));
	double angle2 = std::acos(dot2 / (len3 * len4));

	if (angle1 + angle2 > M_PI + DBL_EPSILON)
		return false;

	return true;
}

void Trimesh::delaunay(const VertexHandle &vh)
{
	std::list<EdgeHandle> edgeList;
	for (VertexOHalfedgeIter voheit = voh_begin(vh);
		voheit != voh_end(vh); 
		++voheit)
	{
		EdgeHandle eh = edge_handle(*voheit);
		
		if (!is_delaunay(eh))
			edgeList.push_back(eh);

		eh = edge_handle(next_halfedge_handle(*voheit));
		if (!is_delaunay(eh))
			edgeList.push_back(eh);
	}

	while (!edgeList.empty())
	{
		EdgeHandle eh = edgeList.front();
		edgeList.pop_front();

		flip(eh);

		HalfedgeHandle heh = halfedge_handle(eh, 0);
		EdgeHandle temp = edge_handle(next_halfedge_handle(heh));
		
		if (!is_delaunay(temp))
			edgeList.push_front(temp);

		temp = edge_handle(next_halfedge_handle(next_halfedge_handle(heh)));
		if (!is_delaunay(temp))
			edgeList.push_front(temp);

		temp = edge_handle(next_halfedge_handle(opposite_halfedge_handle(heh)));
		if (!is_delaunay(temp))
			edgeList.push_front(temp);

		temp = edge_handle(next_halfedge_handle(next_halfedge_handle(opposite_halfedge_handle(heh))));
		if (!is_delaunay(temp))
			edgeList.push_front(temp);
	}
}

void Trimesh::delaunay(const VertexHandle &vh, std::unordered_set<int> &faces)
{
	faces.clear();

	std::list<EdgeHandle> edgeList;
	for (VertexOHalfedgeIter voheit = voh_begin(vh);
		voheit != voh_end(vh);
		++voheit)
	{
		EdgeHandle eh = edge_handle(*voheit);

		if (!is_delaunay(eh))
			edgeList.push_back(eh);

		eh = edge_handle(next_halfedge_handle(*voheit));
		if (!is_delaunay(eh))
			edgeList.push_back(eh);
	}

	while (!edgeList.empty())
	{
		EdgeHandle eh = edgeList.front();
		edgeList.pop_front();

		flip(eh);

		HalfedgeHandle heh = halfedge_handle(eh, 0);

		faces.insert(face_handle(heh).idx());
		faces.insert(opposite_face_handle(heh).idx());

		EdgeHandle temp = edge_handle(next_halfedge_handle(heh));

		if (!is_delaunay(temp))
			edgeList.push_front(temp);

		temp = edge_handle(next_halfedge_handle(next_halfedge_handle(heh)));
		if (!is_delaunay(temp))
			edgeList.push_front(temp);

		temp = edge_handle(next_halfedge_handle(opposite_halfedge_handle(heh)));
		if (!is_delaunay(temp))
			edgeList.push_front(temp);

		temp = edge_handle(next_halfedge_handle(next_halfedge_handle(opposite_halfedge_handle(heh))));
		if (!is_delaunay(temp))
			edgeList.push_front(temp);
	}
}

void Trimesh::delaunay()
{
	int enb = (int)n_edges();

	std::list<int> edgeList;
	for (int e = 0; e < enb; ++e)
	{
		edgeList.push_back(e);
	}

	std::vector<bool> inList(enb, true);

	while (!edgeList.empty())
	{
		int e = edgeList.front();
		edgeList.pop_front();
		inList[e] = false;

		EdgeHandle eh = edge_handle(e);
		if (!is_flippable(eh))
			continue;

		HalfedgeHandle heh1 = halfedge_handle(eh, 0);
		HalfedgeHandle heh2 = next_halfedge_handle(heh1);
		HalfedgeHandle heh3 = opposite_halfedge_handle(heh1);
		HalfedgeHandle heh4 = next_halfedge_handle(heh3);

		double angle1 = calc_sector_angle(heh2);
		double angle2 = calc_sector_angle(heh4);

		if (angle1 + angle2 > M_PI)
		{
			flip(eh);

			int nes[4];
			nes[0] = edge_handle(heh2).idx();
			nes[1] = edge_handle(prev_halfedge_handle(heh1)).idx();
			nes[2] = edge_handle(heh4).idx();
			nes[3] = edge_handle(prev_halfedge_handle(heh3)).idx();

			for (int k = 0; k < 4; ++k)
			{
				if (!inList[nes[k]])
				{
					edgeList.push_back(nes[k]);
					inList[nes[k]] = true;
				}
			}
		}
	}
}

int Trimesh::make_vertices_boundary(double maxangle)
{
	int vnb = (int)n_vertices();
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

		VertexHandle bvh = vertex_handle(bv);
		HalfedgeHandle *bhe = NULL;

		for (auto vohit = voh_begin(bvh); vohit != voh_end(bvh); ++vohit)
		{
			if (is_boundary(*vohit))
			{
				bhe = &(*vohit);
				break;
			}
		}

		if (!bhe)
			break;

		HalfedgeHandle heh = next_halfedge_handle(opposite_halfedge_handle(*bhe));
		double angle = calc_sector_angle(heh);
		if (angle > maxangle)
		{
			bool editable = true;
			VertexHandle vh = to_vertex_handle(heh);
			Normal dir = calc_edge_vector(*bhe);
			dir.normalize();

			for (auto vohit = voh_begin(bvh); vohit != voh_end(bvh); ++vohit)
			{
				HalfedgeHandle
					temp = opposite_halfedge_handle(next_halfedge_handle(*vohit));
				if (temp != *bhe && is_boundary(temp))
				{
					Normal dir2 = calc_edge_vector(temp);
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
				double len = calc_edge_length(*bhe);
				double area = calc_sector_area(heh);
				double height = area * 2.0 / len;

				Normal perp(-dir[1], dir[0], 0);
				perp.normalize();

				FaceHandle fh(opposite_face_handle(*bhe));

				delete_face(fh);
				garbage_collection();
				set_point(vh, point(vh) + height * perp);

				visited[bv] = false;

				++count;
			}
		}
		else
		{
			bv = to_vertex_handle(*bhe).idx();
		}
	}

	return count;
}

int Trimesh::remove_small_triangles(double minarea)
{
	if (n_faces() < 2)
		return 0;

	int count = 0;
	int fnb = (int)n_faces();

	for (int f = fnb - 1; f > -1; --f)
	{
		if (f >= (int)n_faces())
			f = (int)n_faces() - 1;

		FaceHandle fh = face_handle(f);
		double farea = calc_face_area(fh);
		if (farea > minarea)
			continue;

		double minlen = 1e10;
		HalfedgeHandle *heh = NULL;
		for (auto fhit = fh_begin(fh); fhit != fh_end(fh); ++fhit)
		{
			double len = calc_edge_length(*fhit);
			if (len < minlen)
			{
				minlen = len;
				heh = &(*fhit);
			}
		}

		if (!heh)
			continue;

		if (to_vertex_handle(*heh).idx() > from_vertex_handle(*heh).idx())
			heh = &(opposite_halfedge_handle(*heh));

		if (is_boundary(from_vertex_handle(*heh)) && !is_boundary(to_vertex_handle(*heh)))
			heh = &(opposite_halfedge_handle(*heh));

		if (!is_collapse_ok(*heh))
			continue;

		collapse(*heh);
		garbage_collection();

		++count;
	}

	return count;
}

