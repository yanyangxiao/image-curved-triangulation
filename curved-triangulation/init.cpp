
#include <map>
#include <unordered_set>
#include "curved-triangulation.h"
#include "Douglas-Peucker.h"
#include "utils.h"
#include "../extern/triangle/triangle.h"

#define INIT_MAX_ANGLE (M_PI * 100.0 / 180.0)
#define DOUGLAS_PEUCKER_THREASHOLD 3.0

void CurvedTriangulation::random_init(int vnb)
{
	if (!_image)
		return;

	double ratio = double(_height) / _width;

	_mesh.clean();
	VertexHandle vh0 = _mesh.add_vertex(Point(0, 0, 0));
	VertexHandle vh1 = _mesh.add_vertex(Point(1, 0, 0));
	VertexHandle vh2 = _mesh.add_vertex(Point(1, ratio, 0));
	VertexHandle vh3 = _mesh.add_vertex(Point(0, ratio, 0));

	_mesh.property(VPropFeature, vh0) = false;
	_mesh.property(VPropFeature, vh1) = false;
	_mesh.property(VPropFeature, vh2) = false;
	_mesh.property(VPropFeature, vh3) = false;

	_mesh.add_face(vh0, vh1, vh2);
	_mesh.add_face(vh0, vh2, vh3);

	const double shrink = 0.9;

	srand((unsigned int)time(NULL));
	for (int i = 4; i < vnb; ++i)
	{
		double rdm = double(rand()) / RAND_MAX;

		int f = int(rdm * _mesh.n_faces());
		if (f == _mesh.n_faces())
			--f;
		
		//int f = rand() % _mesh.n_faces();
		FaceHandle fh = _mesh.face_handle(f);

		HalfedgeHandle fhe = _mesh.halfedge_handle(fh);
		double angle = _mesh.calc_sector_angle(fhe);
		double maxangle = angle;
		HalfedgeHandle maxhe = _mesh.prev_halfedge_handle(fhe);

		fhe = _mesh.next_halfedge_handle(fhe);
		angle = _mesh.calc_sector_angle(fhe);
		if (angle > maxangle)
		{
			maxangle = angle;
			maxhe = _mesh.prev_halfedge_handle(fhe);
		}

		fhe = _mesh.next_halfedge_handle(fhe);
		angle = _mesh.calc_sector_angle(fhe);
		if (angle > maxangle)
		{
			maxangle = angle;
			maxhe = _mesh.prev_halfedge_handle(fhe);
		}

		VertexHandle vh;

		if (maxangle > INIT_MAX_ANGLE)
		{
			EdgeHandle eh = _mesh.edge_handle(maxhe);
			vh = _mesh.split(eh, _mesh.calc_edge_midpoint(eh));
		}
		else
		{
			Point p0 = _mesh.point(_mesh.from_vertex_handle(fhe));
			Point p1 = _mesh.point(_mesh.to_vertex_handle(fhe));
			fhe = _mesh.next_halfedge_handle(fhe);
			Point p2 = _mesh.point(_mesh.to_vertex_handle(fhe));

			Point pc = (p0 + p1 + p2) / 3;
			p0 = shrink * p0 + (1.0 - shrink) * pc;
			p1 = shrink * p1 + (1.0 - shrink) * pc;
			p2 = shrink * p2 + (1.0 - shrink) * pc;

			double alpha = double(rand()) / RAND_MAX;
			double beta = double(rand()) / RAND_MAX * (1 - alpha);
			double gama = 1 - alpha - beta;
			Point p = p0 * alpha + p1 * beta + p2 * gama;

			vh = _mesh.split(fh, p);
		}

		_mesh.property(VPropFeature, vh) = false;

		_mesh.delaunay(vh);
	}

	process_boundary(false);

	for (auto eit = _mesh.edges_begin(); eit != _mesh.edges_end(); ++eit)
	{
		_mesh.property(EPropFeature, *eit) = false;
		_mesh.property(EPropControl, *eit) = _mesh.calc_edge_midpoint(*eit);
	}

	xlog("vnb = %d, fnb = %d", (int)_mesh.n_vertices(), (int)_mesh.n_faces());
}

void CurvedTriangulation::greedy_init(int vnb)
{
	if (!_image)
		return;

	random_init(4);

	compute_approx();
	insert_vertices(vnb - 4, false);

	process_boundary(true);

	xlog("vnb = %d, fnb = %d", (int)_mesh.n_vertices(), (int)_mesh.n_faces());
}

void CurvedTriangulation::insert_vertices(int vnb, bool considerArea)
{
	if (vnb < 1)
		return;

	std::multimap<double, int, std::greater<double>> errorMap;
	for (int f = 0; f < (int)_mesh.n_faces(); ++f)
	{
		FaceHandle fh = _mesh.face_handle(f);
		MyFaceData *fp = &(_mesh.property(FPropData, fh));
		
		fp->regEnergy = 0.0;
		if (_regweight > 1e-10)
			fp->regEnergy = compute_regenergy(fh);

		double key = fp->approxEnergy + _regweight * fp->regEnergy;

		if (considerArea)
		{
			double area = _mesh.calc_face_area(fh);
			key = std::pow(key, area);
			//key *= area;
		}

		fp->totalEnergy = key;
		errorMap.insert(std::make_pair(key, f));
	}

	std::unordered_set<int> faceList;

	while (vnb-- && !errorMap.empty())
	{
		int fnb = (int)_mesh.n_faces();

		auto topIt = errorMap.begin();
		int f = topIt->second;
		errorMap.erase(topIt);

		FaceHandle fh = _mesh.face_handle(f);
		HalfedgeHandle heh = _mesh.halfedge_handle(fh);

		double angle = _mesh.calc_sector_angle(heh);
		double maxangle = angle;
		HalfedgeHandle maxheh = _mesh.prev_halfedge_handle(heh);

		heh = _mesh.next_halfedge_handle(heh);
		angle = _mesh.calc_sector_angle(heh);
		if (angle > maxangle)
		{
			maxangle = angle;
			maxheh = _mesh.prev_halfedge_handle(heh);
		}

		heh = _mesh.next_halfedge_handle(heh);
		angle = _mesh.calc_sector_angle(heh);
		if (angle > maxangle)
		{
			maxangle = angle;
			maxheh = _mesh.prev_halfedge_handle(heh);
		}

		VertexHandle newvh, fromvh(0), tovh(0);

		EdgeHandle eh = _mesh.edge_handle(maxheh);
		bool isFeature = _mesh.property(EPropFeature, eh);

		bool edgeSplit = false;
		if (!isFeature && maxangle > INIT_MAX_ANGLE)
		{
			edgeSplit = true;

			fromvh = _mesh.from_vertex_handle(maxheh);
			tovh = _mesh.to_vertex_handle(maxheh);

			newvh = _mesh.split(eh, _mesh.calc_edge_midpoint(eh));
		}

		if (!edgeSplit)
		{
			newvh = _mesh.split(fh, _mesh.calc_face_centroid(fh));
		}

		_mesh.property(VPropFeature, newvh) = false;

		for (auto vohit = _mesh.voh_begin(newvh); vohit != _mesh.voh_end(newvh); ++vohit)
		{
			EdgeHandle eh = _mesh.edge_handle(*vohit);
			_mesh.property(EPropControl, eh) = _mesh.calc_edge_midpoint(eh);

			if (edgeSplit &&
				(_mesh.to_vertex_handle(*vohit) == fromvh ||
					_mesh.to_vertex_handle(*vohit) == tovh))
			{
				_mesh.property(EPropFeature, eh) = isFeature;
			}
			else
			{
				_mesh.property(EPropFeature, eh) = false;
			}
		}

		faceList.clear();
		for (auto vfit = _mesh.vf_begin(newvh); vfit != _mesh.vf_end(newvh); ++vfit)
		{
			faceList.insert(vfit->idx());
		}

		// delaunay
		std::unordered_set<int> fset;
		delaunay(newvh, fset);

		faceList.insert(fset.begin(), fset.end());

		for (auto sit = faceList.begin(); sit != faceList.end(); ++sit)
		{
			if (*sit < fnb)
			{
				const MyFaceData& fp = _mesh.property(FPropData, _mesh.face_handle(*sit));

				auto fbegin(errorMap.lower_bound(fp.totalEnergy));
				auto fend(errorMap.upper_bound(fp.totalEnergy));

				for (auto mit = fbegin; mit != fend; ++mit)
				{
					if (mit->second == *sit)
					{
						errorMap.erase(mit);
						break;
					}
				}
			}
		}

		for (auto sit = faceList.begin(); sit != faceList.end(); ++sit)
		{
			FaceHandle tempfh = _mesh.face_handle(*sit);
			
			compute_approx(tempfh);

			MyFaceData *fp = &(_mesh.property(FPropData, tempfh));

			fp->regEnergy = 0.0;
			if (_regweight > 1e-10)
				fp->regEnergy = compute_regenergy(tempfh);

			double key = fp->approxEnergy + _regweight * fp->regEnergy;

			if (considerArea)
			{
				double area = _mesh.calc_face_area(tempfh);
				key = std::pow(key, area);
				//key *= area;
			}

			fp->totalEnergy = key;
			errorMap.insert(std::make_pair(key, *sit));
		}
	}
}

void CurvedTriangulation::features_init(const std::vector<std::vector<int>> &paths, int vnb)
{
	if (!_image)
		return;

	if (paths.empty())
		return greedy_init(vnb);

	int featPixNumber = 0;
	std::multimap<int, int, std::greater<int>> orderedFeats;

	int nb = (int)paths.size();
	for (int k = 0; k < nb; ++k)
	{
		int len = (int)paths[k].size() / 2;
		featPixNumber += len;
		orderedFeats.insert(std::make_pair(len, k));
	}

	std::vector<double> featPoints;
	std::vector<int> featSegments;
	std::vector<std::vector<int>> featMap;

	init_featuremap(paths, featMap, featPoints, featSegments);

	std::set<std::pair<int, int>> segset;

	int step = int(featPixNumber / (vnb * 0.9));
	if (step < 5)
		step = 5;
	
	for (auto mit = orderedFeats.begin(); mit != orderedFeats.end(); ++mit)
	{
		int len = mit->first;
		int k = mit->second;
		
		bool *picked = (bool *)calloc(len, sizeof(bool));
		Douglas_Peucker(&paths[k][0], 0, len - 1, DOUGLAS_PEUCKER_THREASHOLD, step, picked);

		int prev = -1;
		for (int i = 0; i < len; ++i)
		{
			if (!picked[i])
				continue;

			const int *p = &paths[k][2 * i];
			int id = visit_featuremap(featMap, featPoints, p[0], p[1]);

			if (prev > -1 && prev != id)
			{
				// avoid duplication
				std::pair<int, int> temp(prev, id);
				if (prev > id)
				{
					temp.first = id;
					temp.second = prev;
				}

				if (segset.find(temp) == segset.end())
				{
					segset.insert(temp);

					featSegments.push_back(prev);
					featSegments.push_back(id);
				}
			}

			prev = id;
		}

		free(picked);
	}

	struct triangulateio in, out;
	init_triangulateio(&in);

	in.pointlist = &featPoints[0];
	in.numberofpoints = (int)featPoints.size() / 2;
	in.segmentlist = &featSegments[0];
	in.numberofsegments = (int)featSegments.size() / 2;

	init_triangulateio(&out);
	triangulate("pz", &in, &out, (struct triangulateio *)NULL);

	_mesh.clean();
	for (int i = 0; i < out.numberofpoints; ++i)
	{
		VertexHandle vh = _mesh.add_vertex(Point(out.pointlist[2 * i], out.pointlist[2 * i + 1], 0));
		_mesh.property(VPropFeature, vh) = true;
	}

	for (int k = 0; k < out.numberoftriangles; ++k)
	{
		VertexHandle vha = _mesh.vertex_handle(out.trianglelist[3 * k]);
		VertexHandle vhb = _mesh.vertex_handle(out.trianglelist[3 * k + 1]);
		VertexHandle vhc = _mesh.vertex_handle(out.trianglelist[3 * k + 2]);
		_mesh.add_face(vha, vhb, vhc);
	}

	for (auto eit = _mesh.edges_begin(); eit != _mesh.edges_end(); ++eit)
	{
		_mesh.property(EPropFeature, *eit) = false;
	}

	for (int i = 0; i < (int)featSegments.size() / 2; ++i)
	{
		int a = featSegments[2 * i];
		int b = featSegments[2 * i + 1];

		VertexHandle vha = _mesh.vertex_handle(a);
		VertexHandle vhb = _mesh.vertex_handle(b);
		for (auto vohit = _mesh.voh_begin(vha); vohit != _mesh.voh_end(vha); ++vohit)
		{
			if (_mesh.to_vertex_handle(*vohit) == vhb)
			{
				_mesh.property(EPropFeature, _mesh.edge_handle(*vohit)) = true;
				break;
			}
		}
	}

	free_triangulateio(&out);

	for (auto eit = _mesh.edges_begin(); eit != _mesh.edges_end(); ++eit)
	{
		_mesh.property(EPropControl, *eit) = _mesh.calc_edge_midpoint(*eit);
	}

	process_boundary(false);

	compute_approx();
	insert_vertices(vnb - (int)_mesh.n_vertices(), true);
	process_boundary(true);

	xlog("vnb = %d, fnb = %d", (int)_mesh.n_vertices(), (int)_mesh.n_faces());

	svg_selected_features(featPoints, featSegments, "3-selected-features.svg");
}

void CurvedTriangulation::init_featuremap(
	const std::vector<std::vector<int>> &paths,
	std::vector<std::vector<int>> &featMap,
	std::vector<double> &featPoints,
	std::vector<int> &featSegments) const
{
	if (paths.empty())
		return;

	featPoints.clear();
	featSegments.clear();

	featMap.clear();
	featMap.resize(_width + 1, std::vector<int>(_height + 1, -1));
	
	int corners[4];
	corners[0] = visit_featuremap(featMap, featPoints, 0, 0);
	corners[1] = visit_featuremap(featMap, featPoints, _width, 0);
	corners[2] = visit_featuremap(featMap, featPoints, _width, _height);
	corners[3] = visit_featuremap(featMap, featPoints, 0, _height);

	std::vector<std::vector<int>> pathmap(_width, std::vector<int>(_height, -1));
	for (int k = 0; k < (int)paths.size(); ++k)
	{
		for (int i = 0; i < (int)paths[k].size() / 2; ++i)
		{
			const int *p = &paths[k][2 * i];
			pathmap[p[0]][p[1]] = k;
		}
	}

	std::vector<int> boundaryIndices;
	for (int i = 2; i < _width - 2; ++i)
	{
		if (pathmap[i][0] > -1)
			boundaryIndices.push_back(visit_featuremap(featMap, featPoints, i, 0));
	}
	if (boundaryIndices.empty())
		boundaryIndices.push_back(visit_featuremap(featMap, featPoints, _width / 2, 0));

	featSegments.push_back(corners[0]);
	featSegments.push_back(boundaryIndices[0]);

	int bsize = (int)boundaryIndices.size();
	for (int k = 0; k < bsize - 1; ++k)
	{
		featSegments.push_back(boundaryIndices[k]);
		featSegments.push_back(boundaryIndices[k + 1]);
	}

	featSegments.push_back(boundaryIndices[bsize - 1]);
	featSegments.push_back(corners[1]);

	boundaryIndices.clear();
	for (int j = 2; j < _height - 2; ++j)
	{
		if (pathmap[_width - 1][j] > -1)
			boundaryIndices.push_back(visit_featuremap(featMap, featPoints, _width, j));
	}
	if (boundaryIndices.empty())
		boundaryIndices.push_back(visit_featuremap(featMap, featPoints, _width, _height / 2));

	featSegments.push_back(corners[1]);
	featSegments.push_back(boundaryIndices[0]);

	bsize = (int)boundaryIndices.size();
	for (int k = 0; k < bsize - 1; ++k)
	{
		featSegments.push_back(boundaryIndices[k]);
		featSegments.push_back(boundaryIndices[k + 1]);
	}

	featSegments.push_back(boundaryIndices[bsize - 1]);
	featSegments.push_back(corners[2]);

	boundaryIndices.clear();
	for (int i = _width - 3; i > 1; --i)
	{
		if (pathmap[i][_height - 1] > -1)
			boundaryIndices.push_back(visit_featuremap(featMap, featPoints, i, _height));
	}
	if (boundaryIndices.empty())
		boundaryIndices.push_back(visit_featuremap(featMap, featPoints, _width / 2, _height));

	featSegments.push_back(corners[2]);
	featSegments.push_back(boundaryIndices[0]);

	bsize = (int)boundaryIndices.size();
	for (int k = 0; k < bsize - 1; ++k)
	{
		featSegments.push_back(boundaryIndices[k]);
		featSegments.push_back(boundaryIndices[k + 1]);
	}

	featSegments.push_back(boundaryIndices[bsize - 1]);
	featSegments.push_back(corners[3]);

	boundaryIndices.clear();
	for (int j = _height - 3; j > 1; --j)
	{
		if (pathmap[0][j] > -1)
			boundaryIndices.push_back(visit_featuremap(featMap, featPoints, 0, j));
	}
	if (boundaryIndices.empty())
		boundaryIndices.push_back(visit_featuremap(featMap, featPoints, 0, _height / 2));

	featSegments.push_back(corners[3]);
	featSegments.push_back(boundaryIndices[0]);

	bsize = (int)boundaryIndices.size();
	for (int k = 0; k < bsize - 1; ++k)
	{
		featSegments.push_back(boundaryIndices[k]);
		featSegments.push_back(boundaryIndices[k + 1]);
	}

	featSegments.push_back(boundaryIndices[bsize - 1]);
	featSegments.push_back(corners[0]);
}

int CurvedTriangulation::visit_featuremap(
	std::vector<std::vector<int>> &featMap,
	std::vector<double> &featPoints,
	int i, int j, int r) const
{
	if (-1 == featMap[i][j])
	{
		int id = (int)featPoints.size() / 2;
		for (int a = -r; a <= r; ++a)
		{
			int y = j + a;
			if (y < 0 || y > _height)
				continue;

			for (int b = -r; b <= r; ++b)
			{
				int x = i + b;
				if (x < 0 || x > _width)
					continue;

				featMap[x][y] = id;
			}
		}

		if (0 == i || _width == i)
			featPoints.push_back(i / _width);
		else
			featPoints.push_back((double(i) + 0.5) / _width);

		if (0 == j || _height == j)
			featPoints.push_back(double(j) / _width);
		else
			featPoints.push_back((double(j) + 0.5) / _width);
	}

	return featMap[i][j];
}

