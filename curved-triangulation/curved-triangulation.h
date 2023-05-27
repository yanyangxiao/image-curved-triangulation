
#ifndef CURVED_TRIANGULATION_H
#define CURVED_TRIANGULATION_H

#include "trimesh.h"
#include "rasterizer.h"
#include "polynomial.h"
#include "../xlog.h"

#define DEFAULT_SHRINK 0.99
#define OUTPUT_SHRINK 0.9

class CurvedTriangulation
{
	typedef Trimesh::VertexHandle   VertexHandle;
	typedef Trimesh::HalfedgeHandle HalfedgeHandle;
	typedef Trimesh::EdgeHandle     EdgeHandle;
	typedef Trimesh::FaceHandle     FaceHandle;
	typedef Trimesh::Point          Point;
	typedef Trimesh::Normal         Normal;

	typedef Rasterizer<double> MyRasterizer;
	typedef Polynomial<double> MyPolynomial;

	struct MyFaceData
	{
		MyPolynomial approx;
		double       approxEnergy;
		double       regEnergy;
		double       totalEnergy;

		MyFaceData()
			: approx(MyPolynomial(0)), approxEnergy(0.0), regEnergy(0.0), totalEnergy(0.0)
		{ }

		MyFaceData& operator= (const MyFaceData &rhs)
		{
			approx = rhs.approx;
			approxEnergy = rhs.approxEnergy;
			regEnergy = rhs.regEnergy;
			totalEnergy = rhs.totalEnergy;

			return *this;
		}
	};

protected:
	// image
	float *_image; // [0, 1]
	int _width;
	int _height;
	int _channel;
	double _pixwidth;

	// mesh
	Trimesh _mesh;
	OpenMesh::FPropHandleT<MyFaceData> FPropData;
	OpenMesh::EPropHandleT<Point>      EPropControl;
	OpenMesh::EPropHandleT<bool>       EPropFeature;
	OpenMesh::VPropHandleT<bool>       VPropFeature;

	// parameters
	int    _degree;
	double _regweight;

public:
	CurvedTriangulation();
	~CurvedTriangulation();

	// setter
	void set_degree(int d) { _degree = d; }
	void set_regweight(double regw) { _regweight = regw; }
	void set_image(const unsigned char *image, int width, int height, int channel);
	
	// mesh
	void meshface_polygon(const FaceHandle &fh, double shrink, std::vector<double> &polygon) const;
	int process_boundary(bool updateApprox);
	void flip_edge(const EdgeHandle &eh);
	void delaunay(const VertexHandle &vh, std::unordered_set<int> &fset);
	bool self_intersect(const EdgeHandle &eh) const;
	int check_self_intersection() const;
	bool self_intersect(const VertexHandle &vh) const;
	double move_vertex(const VertexHandle &vh, const Point &p);

	// approx	
	void compute_approx(double shrink = DEFAULT_SHRINK);
	void compute_approx(const std::vector<int> &faces, double shrink = DEFAULT_SHRINK);
	void compute_approx(const FaceHandle &fh, double shrink = DEFAULT_SHRINK);

	// energy
	double compute_regenergy(const VertexHandle &vh) const;
	double compute_regenergy(const FaceHandle &fh) const;
	double compute_regenergy() const;
	double compute_approx_energy() const;
	double compute_energy() const;
	
	// init
	void random_init(int vnb);
	void greedy_init(int vnb);
	void features_init(const std::vector<std::vector<int>> &features, int vnb);
	void insert_vertices(int vnb, bool considerArea);
	void init_featuremap(
		const std::vector<std::vector<int>> &paths,
		std::vector<std::vector<int>> &featMap,
		std::vector<double> &featPoints,
		std::vector<int> &featSegments) const;
	int visit_featuremap(
		std::vector<std::vector<int>> &featMap,
		std::vector<double> &featPoints,
		int i, int j, int r = 2) const;

	// gradient
	void compute_vertex_gradient(const VertexHandle &vh, double *vgrad) const;
	void vertex_integral(
		const double *pi, const double *control, const double *pj,
		const MyPolynomial *leftApprox, const MyPolynomial *rightApprox,
		double t0, double t1, double *result) const;
	void compute_control_gradient(const EdgeHandle &eh, double *cgrad) const;
	void control_integral(
		const double *pi, const double *control, const double *pj,
		const MyPolynomial *leftApprox, const MyPolynomial *rightApprox,
		double t0, double t1, double *result) const;

	// optimize
	void optimize(int iter, bool optallcontrols);
	void optimize_edges();
	bool legal_flipping(const EdgeHandle &eh) const;
	void optimize_vertices();
	void optimize_vertex(const VertexHandle &vh, const double *vgrad);
	void vertex_constraint(const VertexHandle &vh, std::vector<double> &region) const;
	void optimize_controls(bool optall);
	void optimize_control(const EdgeHandle &eh, const double *cgrad);
	void control_constraint(const EdgeHandle &eh, std::vector<double> &region) const;

	// io
	void save_mesh(const char *meshname) const;
	bool load_mesh(const char *meshname);
	void svg_mesh(const char *svgname) const;
	void svg_approx(const char *svgname) const;
	void svg_selected_features(
		const std::vector<double> &points, 
		const std::vector<int> &segments, 
		const char *svgname) const;
	void svg_polygon(
		const double *cent, 
		const double *polygon, 
		int nb, 
		const char *svgname) const;
};

#endif