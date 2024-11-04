/*

Data structures for learnply

Eugene Zhang 2005

*/

#ifndef __LEARNPLY_H__
#define __LEARNPLY_H__


#include <vector>
#include "ply.h"
#include "icVector.H"
#define M_PI 3.14159265358979323846
/* forward declarations */
class Triangle;
class Vertex {
public:
  int index;
  icVector3 pos;
  icVector3 normal;
  icVector3 color;

  std::vector<Triangle*> tris;

  void *other_props;

public:
	Vertex(icVector3& p) : index(-1), pos(p), normal(icVector3(0.0, 1.0, 0.0)) {  }
	Vertex(double xx, double yy, double zz): index(-1), pos(icVector3(xx,yy,zz)), normal(icVector3(0.0, 1.0, 0.0)) {  }

	int ntris() { return (int)tris.size(); }
};

class Edge 
{
public:
  int index;

  Vertex *verts[2];

  std::vector<Triangle*> tris;

  double length;

public:
	Edge() : index(-1), length(0.0)
	{ 
		verts[0] = NULL; verts[1] = NULL;
	}
	int ntris() { return (int)tris.size(); }
};

class Triangle {
public:
  int index;

  icVector3 normal;
  double  area;

  int nverts;
  Vertex *verts[3];
  Edge   *edges[3];

  void *other_props; 
  std::vector<icVector3> inter_pts;//interpolated points

public:
	Triangle() : index(-1), nverts(3), normal(icVector3(0.0, 1.0, 0.0)), area(0.0)
	{  
		verts[0] = NULL;
		verts[1] = NULL;
		verts[2] = NULL;
		edges[0] = NULL;
		edges[1] = NULL;
		edges[2] = NULL;
	}
};

class Weights {
public:
	double value;
	Vertex* vertex;
public:
	Weights() : value(0.0)
	{
		vertex = NULL;
	}
};

class Polyhedron 
{
public:
	PlyFile* in_ply;
	std::vector<Triangle*> tlist;/* list of triangles */
	std::vector<Vertex*>   vlist;/* list of vertices */
	std::vector<Edge*>     elist;/* list of edges */

	unsigned char orientation;  // 0=ccw, 1=cw
	icVector3 center;
	double radius;
	double area;

	int seed;//selection
public:

	Polyhedron();
	Polyhedron(std::vector<Vertex*>& verts, std::vector<Triangle*>& tris, bool re_index = true);
	Polyhedron(FILE*);
	//
	int ntris() { return (int)tlist.size(); }
	int nverts() { return (int)vlist.size(); }
	int nedges() { return (int)elist.size(); }
	void generateColorControlMesh();
	void InterpolateColor(const Polyhedron& Control_Mesh, Polyhedron& Model_Mesh, const std::vector<std::vector<Weights>>& Weights_Map);
	void InterpolatePos(const Polyhedron& Control_Mesh, Polyhedron& Model_Mesh, const std::vector<std::vector<Weights>>& Weights_Map);

	// initialization and finalization
	void initialize();
	void finalize();
	std::vector<std::vector<Weights>> getMeanValueCoordinates(Polyhedron& controlMesh) const;
	std::vector<Weights> getVertexMeanValues(const Vertex* modelVertex, Polyhedron& controlMesh) const;
	double determinant(const icVector3& u0, const icVector3& u1, const icVector3& u2) const;

private:
	  PlyOtherProp *vert_other, *face_other;
	  void write_file(FILE*);

	  void calc_vert_normals();

	  void create_edge(Vertex *, Vertex *);
	  void create_edges();

	  void mean_value_colorinter();
	  int face_to_vertex_ref(Triangle *, Vertex *);
	  void order_vertex_to_tri_ptrs(Vertex *);
	  void vertex_to_tri_ptrs();

	  Triangle *find_common_edge(Triangle *, Vertex *, Vertex *);
	  Triangle *other_triangle(Edge *, Triangle *);

	  void calc_bounding_sphere();
	  void calc_face_normals_and_area();
	  void calc_edge_length();


	  void create_pointers();

};

#endif /* __LEARNPLY_H__ */

