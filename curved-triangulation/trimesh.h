
#ifndef OPENMESH_TRIMESH_H
#define OPENMESH_TRIMESH_H

#include <unordered_set>

#include "OpenMesh/Core/Mesh/Traits.hh"
#include "OpenMesh/Core/Geometry/VectorT.hh"
#include "OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh"
#include "OpenMesh/Core/Mesh/TriConnectivity.hh"
#include "OpenMesh/Core/IO/MeshIO.hh"
#include "OpenMesh/Core/IO/Options.hh"

// define traits
struct MyTraits : public OpenMesh::DefaultTraits
{
	typedef OpenMesh::Vec3d Point;
	typedef OpenMesh::Vec3d Normal;

	VertexAttributes(OpenMesh::Attributes::Status);
	FaceAttributes(OpenMesh::Attributes::Status);
	EdgeAttributes(OpenMesh::Attributes::Status);
};
// Select mesh type (TriMesh) and kernel (ArrayKernel)
// and define my personal mesh type (Trimesh)
typedef OpenMesh::TriMesh_ArrayKernelT<MyTraits>  OpenMeshTrimesh;

class Trimesh : public OpenMeshTrimesh
{
public:
	Trimesh() : OpenMeshTrimesh() {}

	void face_triangle(const FaceHandle &fh, double *tri) const;
	bool is_flippable(const EdgeHandle &eh) const;
	bool is_delaunay(const EdgeHandle &eh) const;

	void delaunay(const VertexHandle &vh);
	void delaunay(const VertexHandle &vh, std::unordered_set<int> &faces);
	void delaunay();

	int make_vertices_boundary(double maxangle);
	int remove_small_triangles(double minarea);
};

#endif
