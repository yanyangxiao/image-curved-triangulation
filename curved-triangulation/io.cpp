
#include "curved-triangulation.h"
#include "svg-utils.h"

void CurvedTriangulation::save_mesh(const char *meshname) const
{
	if (_mesh.vertices_empty())
		return;

	if (!OpenMesh::IO::write_mesh(_mesh, meshname))
	{
		xlog("failed");
		return;
	}

	xlog("success, %s", meshname);
}

bool CurvedTriangulation::load_mesh(const char *meshname)
{
	_mesh.clean();

	if (OpenMesh::IO::read_mesh(_mesh, meshname))
	{
		xlog("success, vnb = %d, fnb = %d", (int)_mesh.n_vertices(), (int)_mesh.n_faces());
		return true;
	}

	xlog("failed, %s", meshname);
	return false;
}

void CurvedTriangulation::svg_mesh(const char *svgname) const
{
	if (_mesh.faces_empty())
		return;

	FILE *svgFile = fopen(svgname, "wb");
	if (!svgFile)
		return;

	int margin = 20;
	int svgw = _width + margin;
	int svgh = _height + margin;
	fprintf(svgFile, "<svg xmlns=\"http://www.w3.org/2000/svg\" xmlns:xlink=\"http://www.w3.org/1999/xlink\"\n  width=\"%d\" height=\"%d\" viewBox=\"0 0 %d %d\">\n",
		svgw, svgh, svgw, svgh);
	fprintf(svgFile, "<metadata>\n");
	fprintf(svgFile, "  Created by Xiao, Yanyang (yanyangxiaoxyy@gmail.com)\n");
	fprintf(svgFile, "  All copyrights reserved\n");
	fprintf(svgFile, "</metadata>\n");

	// curves
	double linew = 2.0;
	fprintf(svgFile, "<g fill=\"none\" stroke=\"rgb(%f, %f, %f)\" stroke-width=\"%f\" stroke-linecap=\"round\">\n",
		50.0, 50.0, 250.0, linew);
	for (auto eit = _mesh.edges_begin(); eit != _mesh.edges_end(); ++eit)
	{
		HalfedgeHandle heh = _mesh.halfedge_handle(*eit, 0);

		const Point& p1 = _mesh.point(_mesh.from_vertex_handle(heh));
		const Point& p2 = _mesh.property(EPropControl, *eit);
		const Point& p3 = _mesh.point(_mesh.to_vertex_handle(heh));
		
		double x1 = p1[0] * _width + margin / 2;
		double y1 = svgh - (p1[1] * _width + margin / 2);
		double x2 = p2[0] * _width + margin / 2;
		double y2 = svgh - (p2[1] * _width + margin / 2);
		double x3 = p3[0] * _width + margin / 2;
		double y3 = svgh - (p3[1] * _width + margin / 2);

		fprintf(svgFile, "<!--eid=%d-->\n", eit->idx());
		fprintf(svgFile, "<path d=\"M %f %f ", x1, y1);
		fprintf(svgFile, "Q %f %f %f %f\"/>\n", x2, y2, x3, y3);
	}

	fprintf(svgFile, "</g>\n");

	// control segments
	linew = 1.0;
	fprintf(svgFile, "<g fill=\"none\" stroke=\"rgb(%f, %f, %f)\" stroke-width=\"%f\" stroke-linecap=\"round\">\n",
		50.0, 250.0, 50.0, linew);
	for (auto eit = _mesh.edges_begin(); eit != _mesh.edges_end(); ++eit)
	{
		if (self_intersect(*eit))
			continue;

		HalfedgeHandle heh = _mesh.halfedge_handle(*eit, 0);

		const Point& p1 = _mesh.point(_mesh.from_vertex_handle(heh));
		const Point& p2 = _mesh.property(EPropControl, *eit);
		const Point& p3 = _mesh.point(_mesh.to_vertex_handle(heh));
		
		double x1 = p1[0] * _width + margin / 2;
		double y1 = svgh - (p1[1] * _width + margin / 2);
		double x2 = p2[0] * _width + margin / 2;
		double y2 = svgh - (p2[1] * _width + margin / 2);
		double x3 = p3[0] * _width + margin / 2;
		double y3 = svgh - (p3[1] * _width + margin / 2);

		fprintf(svgFile, "<!--eid=%d-->\n", eit->idx());
		fprintf(svgFile, "<path d=\"M %f %f ", x1, y1);
		fprintf(svgFile, "L %f %f L %f %f\"/>\n", x2, y2, x3, y3);
	}

	fprintf(svgFile, "</g>\n");

	// control segments
	linew = 1.0;
	fprintf(svgFile, "<g fill=\"none\" stroke=\"rgb(%f, %f, %f)\" stroke-width=\"%f\" stroke-linecap=\"round\">\n",
		250.0, 50.0, 250.0, linew);
	for (auto eit = _mesh.edges_begin(); eit != _mesh.edges_end(); ++eit)
	{
		if (!self_intersect(*eit))
			continue;

#ifdef _DEBUG
		xlog("self-intersect: e = %d", eit->idx());
#endif

		HalfedgeHandle heh = _mesh.halfedge_handle(*eit, 0);

		const Point& p1 = _mesh.point(_mesh.from_vertex_handle(heh));
		const Point& p2 = _mesh.property(EPropControl, *eit);
		const Point& p3 = _mesh.point(_mesh.to_vertex_handle(heh));

		double x1 = p1[0] * _width + margin / 2;
		double y1 = svgh - (p1[1] * _width + margin / 2);
		double x2 = p2[0] * _width + margin / 2;
		double y2 = svgh - (p2[1] * _width + margin / 2);
		double x3 = p3[0] * _width + margin / 2;
		double y3 = svgh - (p3[1] * _width + margin / 2);

		fprintf(svgFile, "<!--eid=%d-->\n", eit->idx());
		fprintf(svgFile, "<path d=\"M %f %f ", x1, y1);
		fprintf(svgFile, "L %f %f L %f %f\"/>\n", x2, y2, x3, y3);
	}

	fprintf(svgFile, "</g>\n");

	// vertices
	double radius = 4.0;
	fprintf(svgFile, "<g fill=\"rgb(%f, %f, %f)\">\n", 250.0, 50.0, 50.0);
	for (auto vit = _mesh.vertices_begin(); vit != _mesh.vertices_end(); ++vit)
	{
		const Point& p = _mesh.point(*vit);
		double x = p[0] * _width + margin / 2;
		double y = svgh - (p[1] * _width + margin / 2);

		fprintf(svgFile, "<!--vid=%d-->\n", vit->idx());
		fprintf(svgFile, "<circle cx=\"%f\" cy=\"%f\" r=\"%f\"/>\n", x, y, radius);
	}
	fprintf(svgFile, "</g>\n");

	fprintf(svgFile, "</svg>\n");
	fclose(svgFile);
}

void CurvedTriangulation::svg_approx(const char *svgname) const
{
	if (_mesh.faces_empty())
		return;

	FILE *svgFile = fopen(svgname, "wb");
	if (!svgFile)
		return;

	fprintf(svgFile, "<svg xmlns=\"http://www.w3.org/2000/svg\" xmlns:xlink=\"http://www.w3.org/1999/xlink\"\n  width=\"%d\" height=\"%d\" viewBox=\"0 0 %d %d\">\n",
		_width, _height, _width, _height);
	fprintf(svgFile, "<metadata>\n");
	fprintf(svgFile, "  Created by Xiao, Yanyang (yanyangxiaoxyy@gmail.com)\n");
	fprintf(svgFile, "  All copyrights reserved\n");
	fprintf(svgFile, "</metadata>\n");

	for (auto fit = _mesh.faces_begin(); fit != _mesh.faces_end(); ++fit)
	{
		HalfedgeHandle fhe = _mesh.halfedge_handle(*fit);
		const MyFaceData *fp = &(_mesh.property(FPropData, *fit));
		const std::vector<double>* coeff = fp->approx.coefficients();

		if (0 == _degree)
		{
			if (1 == _channel)
				fprintf(svgFile, "<g fill=\"rgb(%f, %f, %f)\">\n",
				(*coeff)[0] * 255, (*coeff)[0] * 255, (*coeff)[0] * 255);
			else if (3 == _channel)
				fprintf(svgFile, "<g fill=\"rgb(%f, %f, %f)\">\n",
				(*coeff)[0] * 255, (*coeff)[1] * 255, (*coeff)[2] * 255);
		}
		else if (1 == _degree)
		{
			std::vector<double> polygon;
			meshface_polygon(*fit, 1.0, polygon);

			double dir[2] = { 0.0 };
			for (int k = 0; k < _channel; ++k)
			{
				dir[0] += (*coeff)[3 * k];
				dir[1] += (*coeff)[3 * k + 1];
			}
			dir[0] /= _channel;
			dir[1] /= _channel;

			double stopa[5], stopb[5];
			svg_linestops(&polygon[0], (int)polygon.size() / 2, dir, stopa, stopb);

			for (int k = 0; k < _channel; ++k)
			{
				stopa[2 + k] = (*coeff)[3 * k] * stopa[0] + (*coeff)[3 * k + 1] * stopa[1] + (*coeff)[3 * k + 2];
				stopa[2 + k] *= 255;
				if (stopa[2 + k] < 0.0)
					stopa[2 + k] = 0.0;
				if (stopa[2 + k] > 255.0)
					stopa[2 + k] = 255.0;

				stopb[2 + k] = (*coeff)[3 * k] * stopb[0] + (*coeff)[3 * k + 1] * stopb[1] + (*coeff)[3 * k + 2];
				stopb[2 + k] *= 255;
				if (stopb[2 + k] < 0.0)
					stopb[2 + k] = 0.0;
				if (stopb[2 + k] > 255.0)
					stopb[2 + k] = 255.0;
			}

			stopa[0] *= _width;
			stopa[1] = _height - stopa[1] * _width;
			stopb[0] *= _width;
			stopb[1] = _height - stopb[1] * _width;

			fprintf(svgFile, "<linearGradient id=\"grad%d\" gradientUnits=\"userSpaceOnUse\" x1=\"%f\" y1=\"%f\" x2=\"%f\" y2=\"%f\">\n",
				fit->idx(), stopa[0], stopa[1], stopb[0], stopb[1]);
			fprintf(svgFile, "  <stop offset=\"0%%\" stop-color=\"rgb(%f, %f, %f)\"/>\n",
				stopa[2], stopa[3], stopa[4]);
			fprintf(svgFile, "  <stop offset=\"100%%\" stop-color=\"rgb(%f, %f, %f)\"/>\n",
				stopb[2], stopb[3], stopb[4]);
			fprintf(svgFile, "</linearGradient>\n");
			fprintf(svgFile, "<g fill=\"url(#grad%d)\">\n", fit->idx());
		}

		EdgeHandle eh = _mesh.edge_handle(fhe);

		const Point& p1 = _mesh.point(_mesh.from_vertex_handle(fhe));
		const Point& p2 = _mesh.property(EPropControl, eh);
		const Point& p3 = _mesh.point(_mesh.to_vertex_handle(fhe));

		fhe = _mesh.next_halfedge_handle(fhe);
		eh = _mesh.edge_handle(fhe);
		const Point& p4 = _mesh.property(EPropControl, eh);
		const Point& p5 = _mesh.point(_mesh.to_vertex_handle(fhe));

		fhe = _mesh.next_halfedge_handle(fhe);
		eh = _mesh.edge_handle(fhe);
		const Point& p6 = _mesh.property(EPropControl, eh);

		fprintf(svgFile, "<!--fid=%d-->\n", fit->idx());
		fprintf(svgFile, "<path d=\"M %f %f ", p1[0] * _width, _height - p1[1] * _width);
		fprintf(svgFile, "Q %f %f %f %f ", p2[0] * _width, _height - p2[1] * _width, p3[0] * _width, _height - p3[1] * _width);
		fprintf(svgFile, "Q %f %f %f %f ", p4[0] * _width, _height - p4[1] * _width, p5[0] * _width, _height - p5[1] * _width);
		fprintf(svgFile, "Q %f %f %f %f Z\" />\n", p6[0] * _width, _height - p6[1] * _width, p1[0] * _width, _height - p1[1] * _width);
		fprintf(svgFile, "</g>\n");
	}

	fprintf(svgFile, "</svg>\n");
	fclose(svgFile);
}

void CurvedTriangulation::svg_selected_features(
	const std::vector<double> &points, 
	const std::vector<int> &segments, 
	const char *svgname) const
{
	FILE *svgFile = fopen(svgname, "wb");
	if (!svgFile)
		return;

	int margin = 20;
	int svgw = _width + margin;
	int svgh = _height + margin;
	fprintf(svgFile, "<svg xmlns=\"http://www.w3.org/2000/svg\" xmlns:xlink=\"http://www.w3.org/1999/xlink\"\n  width=\"%d\" height=\"%d\" viewBox=\"0 0 %d %d\">\n",
		svgw, svgh, svgw, svgh);
	fprintf(svgFile, "<metadata>\n");
	fprintf(svgFile, "  Created by Xiao, Yanyang (yanyangxiaoxyy@gmail.com)\n");
	fprintf(svgFile, "  All copyrights reserved\n");
	fprintf(svgFile, "</metadata>\n");

	// segments
	double linew = 3.0;
	fprintf(svgFile, "<g fill=\"none\" stroke=\"rgb(%f, %f, %f)\" stroke-width=\"%f\" stroke-linecap=\"round\">\n",
		50.0, 250.0, 50.0, linew);
	for (int k = 0; k < (int)segments.size() / 2; ++k)
	{
		int a = segments[2 * k];
		int b = segments[2 * k + 1];

		const double *pa = &points[2 * a];
		const double *pb = &points[2 * b];

		double x1 = pa[0] * _width + margin / 2;
		double y1 = svgh - (pa[1] * _width + margin / 2);
		double x2 = pb[0] * _width + margin / 2;
		double y2 = svgh - (pb[1] * _width + margin / 2);

		fprintf(svgFile, "<line x1=\"%f\" y1=\"%f\" x2=\"%f\" y2=\"%f\"/>\n", x1, y1, x2, y2);
	}

	fprintf(svgFile, "</g>\n");

	// vertices
	double radius = 5.0;
	fprintf(svgFile, "<g fill=\"rgb(%f, %f, %f)\">\n", 250.0, 50.0, 50.0);
	for (int k = 0; k < (int)points.size() / 2; ++k)
	{
		const double *p = &points[2 * k];

		double x = p[0] * _width + margin / 2;
		double y = svgh - (p[1] * _width + margin / 2);

		fprintf(svgFile, "<circle cx=\"%f\" cy=\"%f\" r=\"%f\"/>\n", x, y, radius);
	}
	fprintf(svgFile, "</g>\n");

	fprintf(svgFile, "</svg>\n");
	fclose(svgFile);
}

void CurvedTriangulation::svg_polygon(
	const double *cent, 
	const double *polygon, 
	int nb, 
	const char *svgname) const
{
	FILE *svgFile = fopen(svgname, "wb");
	if (!svgFile)
		return;

	int margin = 20;
	int svgw = _width + margin;
	int svgh = _height + margin;

	fprintf(svgFile, "<svg xmlns=\"http://www.w3.org/2000/svg\" xmlns:xlink=\"http://www.w3.org/1999/xlink\"\n  width=\"%d\" height=\"%d\" viewBox=\"0 0 %d %d\">\n",
		svgw, svgh, svgw, svgh);
	fprintf(svgFile, "<metadata>\n");
	fprintf(svgFile, "  Created by Xiao, Yanyang (yanyangxiaoxyy@gmail.com)\n");
	fprintf(svgFile, "  All copyrights reserved\n");
	fprintf(svgFile, "</metadata>\n");

	fprintf(svgFile, "<g fill=\"rgb(%f, %f, %f)\">\n", 100.0, 200.0, 100.0);
	fprintf(svgFile, "<path d=\"");
	for (int k = 0; k < nb; ++k)
	{
		double x = _width * polygon[2 * k] + margin / 2;
		double y = svgh - (_width * polygon[2 * k + 1] + margin / 2);

		if (0 == k)
			fprintf(svgFile, "M %f %f ", x, y);
		else
			fprintf(svgFile, "L %f %f ", x, y);
	}
	fprintf(svgFile, "Z\"/>\n");
	fprintf(svgFile, "</g>\n");

	// vertices
	double radius = 5.0;
	fprintf(svgFile, "<g fill=\"rgb(%f, %f, %f)\">\n", 250.0, 50.0, 50.0);
	double x = cent[0] * _width + margin / 2;
	double y = svgh - (cent[1] * _width + margin / 2);
	fprintf(svgFile, "<circle cx=\"%f\" cy=\"%f\" r=\"%f\"/>\n", x, y, radius);
	fprintf(svgFile, "</g>\n");

	fprintf(svgFile, "</svg>\n");
	fclose(svgFile);
}
