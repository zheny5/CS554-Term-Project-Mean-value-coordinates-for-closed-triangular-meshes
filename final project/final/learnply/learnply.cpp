/*

Functions for learnply

Eugene Zhang, 2005
*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <fstream>
#include "ply.h"
#include "icVector.H"
#include "icMatrix.H"
#include "learnply.h"
#include "learnply_io.h"
#include <iostream>
#include <cmath>

Polyhedron::Polyhedron() : orientation(0), seed(-1), center(0.0), radius(0.0), area(0.0), in_ply(NULL)
{

}

Polyhedron::Polyhedron(std::vector<Vertex*>& verts, std::vector<Triangle*>& tris, bool re_index)
{
    vlist = verts;
    tlist = tris;
    if (re_index)
    {
        /* index the vertices and triangles */
        for (int i = 0; i < nverts(); i++) { vlist[i]->index = i; }
        for (int i = 0; i < ntris(); i++) { tlist[i]->index = i; }
    }
}

/******************************************************************************
Read in a polyhedron from a file.
******************************************************************************/
Polyhedron::Polyhedron(FILE *file) : orientation(0), seed(-1), center(0.0), radius(0.0), area(0.0), in_ply(NULL)
{
  int i,j;
  int elem_count;
  char *elem_name;

  /*** Read in the original PLY object ***/
  in_ply = read_ply (file);

  for (i = 0; i < in_ply->num_elem_types; i++) 
  {
    /* prepare to read the i'th list of elements */
    elem_name = setup_element_read_ply (in_ply, i, &elem_count);

    if (equal_strings ("vertex", elem_name)) 
	{
      /* create a vertex list to hold all the vertices */
	  vlist.resize(elem_count, NULL);
      /* set up for getting vertex elements */
      setup_property_ply (in_ply, &vert_props[0]);
      setup_property_ply (in_ply, &vert_props[1]);
      setup_property_ply (in_ply, &vert_props[2]);
      vert_other = get_other_properties_ply (in_ply, offsetof(Vertex_io,other_props));
      /* grab all the vertex elements */
      for (j = 0; j < nverts(); j++) 
	  {
        Vertex_io vert;
        get_element_ply (in_ply, (void *) &vert);
        /* copy info from the "vert" structure */
        vlist[j] = new Vertex (vert.x, vert.y, vert.z);
        vlist[j]->other_props = vert.other_props;
      }
    }
    else if (equal_strings ("face", elem_name)) 
	{
      /* create a list to hold all the face elements */
	  tlist.resize(elem_count, NULL);
      /* set up for getting face elements */
      setup_property_ply (in_ply, &face_props[0]);
      face_other = get_other_properties_ply (in_ply, offsetof(Face_io,other_props));

      /* grab all the face elements */
      for (j = 0; j < elem_count; j++) 
	  {
        Face_io face;
        get_element_ply (in_ply, (void *) &face);

        if (face.nverts != 3) 
		{
          fprintf (stderr, "Face has %d vertices (should be three).\n", face.nverts);
          exit (-1);
        }

        /* copy info from the "face" structure */
        tlist[j] = new Triangle;
        tlist[j]->verts[0] = (Vertex *) face.verts[0];
        tlist[j]->verts[1] = (Vertex *) face.verts[1];
        tlist[j]->verts[2] = (Vertex *) face.verts[2];
        tlist[j]->other_props = face.other_props;
      }
    }
    else
      get_other_element_ply (in_ply);
  }

  /* close the file */
  close_ply (in_ply);

  /* fix up vertex pointers in triangles */
  for (i = 0; i < ntris(); i++) 
  {
    tlist[i]->verts[0] = vlist[(int) tlist[i]->verts[0]];
    tlist[i]->verts[1] = vlist[(int) tlist[i]->verts[1]];
    tlist[i]->verts[2] = vlist[(int) tlist[i]->verts[2]];
  }

  /* get rid of triangles that use the same vertex more than once */
  for (std::vector<Triangle*>::iterator it = tlist.begin(); it != tlist.end();)
  {
	Triangle *tri = *it;
	Vertex *v0 = tri->verts[0];
	Vertex *v1 = tri->verts[1];
	Vertex *v2 = tri->verts[2];

	if (v0 == v1 || v1 == v2 || v2 == v0) 
	{
		delete tri;
		it = tlist.erase(it);
	}
	else { it++; }
  }

  /* index the vertices and triangles */
  for (int i = 0; i < nverts(); i++) { vlist[i]->index = i; }
  for (int i = 0; i < ntris(); i++)  { tlist[i]->index = i; }
}


/******************************************************************************
Write out a polyhedron to a file.
******************************************************************************/
void Polyhedron::write_file(FILE *file)
{
  int i;
  PlyFile *ply;
  char **elist;
  int num_elem_types;

  /*** Write out the transformed PLY object ***/
  elist = get_element_list_ply (in_ply, &num_elem_types);
  ply = write_ply (file, num_elem_types, elist, in_ply->file_type);

  /* describe what properties go into the vertex elements */
  describe_element_ply (ply, "vertex", nverts());
  describe_property_ply (ply, &vert_props[0]);
  describe_property_ply (ply, &vert_props[1]);
  describe_property_ply (ply, &vert_props[2]);

  describe_element_ply (ply, "face", ntris());
  describe_property_ply (ply, &face_props[0]);

  copy_comments_ply (ply, in_ply);
  char mm[1024];
  sprintf(mm, "modified by learnply");
  append_comment_ply (ply, mm);
  copy_obj_info_ply (ply, in_ply);

  header_complete_ply (ply);

  /* set up and write the vertex elements */
  put_element_setup_ply (ply, "vertex");
  for (i = 0; i < nverts(); i++) 
  {
    Vertex_io vert;
    /* copy info to the "vert" structure */
    vert.x = (float)vlist[i]->pos.x;
    vert.y = (float)vlist[i]->pos.y;
    vert.z = (float)vlist[i]->pos.z;
    vert.other_props = vlist[i]->other_props;

    put_element_ply (ply, (void *) &vert);
  }

  /* index all the vertices */
  for (i = 0; i < nverts(); i++)
    vlist[i]->index = i;

  /* set up and write the face elements */
  put_element_setup_ply (ply, "face");

  Face_io face;
  face.verts = new int[3];
  
  for (i = 0; i < ntris(); i++) 
  {
    /* copy info to the "face" structure */
    face.nverts = 3;
    face.verts[0] = tlist[i]->verts[0]->index;
    face.verts[1] = tlist[i]->verts[1]->index;
    face.verts[2] = tlist[i]->verts[2]->index;
    face.other_props = tlist[i]->other_props;
    put_element_ply (ply, (void *) &face);
  }
  put_other_elements_ply (ply);

  close_ply (ply);
  free_ply (ply);
}

void Polyhedron::initialize()
{
	create_pointers();
	calc_edge_length();
    calc_bounding_sphere();
    calc_face_normals_and_area();
    calc_vert_normals();
	seed = -1;
}

void Polyhedron::finalize(){
	int i;

	for (i=0; i<ntris(); i++){
		free(tlist[i]->other_props);
		free(tlist[i]);
	}
	for (i=0; i<nedges(); i++) {
		elist[i]->tris.clear();
		free(elist[i]);
	}
	for (i=0; i<nverts(); i++) {
        vlist[i]->tris.clear();
		free(vlist[i]->other_props);
		free(vlist[i]);
	}
	tlist.clear();
	vlist.clear();
	elist.clear();

	if (!vert_other)
		free(vert_other);
	if (!face_other)
		free(face_other);
}

/******************************************************************************
Find out if there is another face that shares an edge with a given face.

Entry:
  f1    - face that we're looking to share with
  v1,v2 - two vertices of f1 that define edge

Exit:
  return the matching face, or NULL if there is no such face
******************************************************************************/
Triangle *Polyhedron::find_common_edge(Triangle *f1, Vertex *v1, Vertex *v2)
{
  int i,j;
  Triangle *f2;
  Triangle *adjacent = NULL;

  /* look through all faces of the first vertex */

  for (i = 0; i < v1->ntris(); i++) {
    f2 = v1->tris[i];
    if (f2 == f1)
      continue;
    /* examine the vertices of the face for a match with the second vertex */
    for (j = 0; j < f2->nverts; j++) {

      /* look for a match */
      if (f2->verts[j] == v2) 
	  {
		/* if we've got a match, return this face */
        return (f2);
      }
    }
  }
  return (adjacent);
}

/******************************************************************************
Create an edge.

Entry:
  v1,v2 - two vertices of f1 that define edge
******************************************************************************/
void Polyhedron::create_edge(Vertex *v1, Vertex *v2)
{
  int i,j;
  Triangle *f;

  /* create the edge */
  Edge *e  = new Edge;
  e->index = nedges();
  e->verts[0] = v1;
  e->verts[1] = v2;
  elist.push_back(e);

  /* count all triangles that will share the edge, and do this */
  /* by looking through all faces of the first vertex */
  int ntris = 0;
  for (i = 0; i < v1->ntris(); i++)
  {
    f = v1->tris[i];
    /* examine the vertices of the face for a match with the second vertex */
    for (j = 0; j < 3; j++) 
    { /* look for a match */
        if (f->verts[j] == v2)  { ntris++; break; }
    }
  }

  /* make room for the face pointers (at least two) */
  e->tris.resize(ntris);

  /* create pointers from edges to faces and vice-versa */
  int tidx = 0;
  for (i = 0; i < v1->ntris(); i++) 
  {
    f = v1->tris[i];
    /* examine the vertices of the face for a match with the second vertex */
    for (j = 0; j < 3; j++)
    {
        if (f->verts[j] == v2)
        {
            e->tris[tidx] = f; tidx++;

            if      (f->verts[(j + 1) % 3] == v1) { f->edges[j] = e; }
            else if (f->verts[(j + 2) % 3] == v1) { f->edges[(j + 2) % 3] = e; }
            else 
            {
                fprintf(stderr, "Non-recoverable inconsistancy in create_edge()\n");
                exit(-1);
            }
            break;  /* we'll only find one instance of v2 */
        }
    }
  }
}


/******************************************************************************
Create edges.
******************************************************************************/
void Polyhedron::create_edges()
{
    /* create all the edges by examining all the triangles */
    for (int i = 0; i < ntris(); i++) 
    {
        Triangle* f = tlist[i];
        for (int j = 0; j < 3; j++) 
        {
            if (f->edges[j]) { continue; }/* skip over edges that we've already created */
            Vertex* v1 = f->verts[j];
            Vertex* v2 = f->verts[(j+1) % f->nverts];
            create_edge (v1, v2);
        }
    }
}

/******************************************************************************
Mean value coordinates color interpolation .
******************************************************************************/
void Polyhedron::mean_value_colorinter()
{




}
/******************************************************************************
Create pointers from vertices to faces.
******************************************************************************/
void Polyhedron::vertex_to_tri_ptrs()
{
    for (int i = 0; i < ntris(); i++)
    {
        Triangle* f = tlist[i];
        for (int j = 0; j < f->nverts; j++) 
        {
            Vertex* v = f->verts[j];
            v->tris.push_back(f);
        }
    }
}

/******************************************************************************
Find the other triangle that is incident on an edge, or NULL if there is
no other.
******************************************************************************/
Triangle *Polyhedron::other_triangle(Edge *edge, Triangle *tri)
{
  /* search for any other triangle */
  for (int i = 0; i < edge->ntris(); i++)
    if (edge->tris[i] != tri)
      return (edge->tris[i]);

  /* there is no such other triangle if we get here */
  return (NULL);
}

/******************************************************************************
Order the pointers to faces that are around a given vertex.

Entry:
  v - vertex whose face list is to be ordered
******************************************************************************/
void Polyhedron::order_vertex_to_tri_ptrs(Vertex *v)
{
  int i,j;
  Triangle *f;
  Triangle *fnext;
  int nf;
  int vindex;
  int boundary;
  int count;

  nf = v->ntris();
  f = v->tris[0];

  /* go backwards (clockwise) around faces that surround a vertex */
  /* to find out if we reach a boundary */

  boundary = 0;

  for (i = 1; i <= nf; i++) {

    /* find reference to v in f */
    vindex = -1;
    for (j = 0; j < f->nverts; j++)
      if (f->verts[j] == v) {
	vindex = j;
	break;
      }

    /* error check */
    if (vindex == -1) {
      fprintf (stderr, "can't find vertex #1\n");
      exit (-1);
    }

    /* corresponding face is the previous one around v */
    fnext = other_triangle (f->edges[vindex], f);

    /* see if we've reached a boundary, and if so then place the */
    /* current face in the first position of the vertice's face list */

    if (fnext == NULL) {
      /* find reference to f in v */
      for (j = 0; j < v->ntris(); j++)
        if (v->tris[j] == f) 
        {
	        v->tris[j] = v->tris[0];
	        v->tris[0] = f;
	        break;
        }
        boundary = 1;
        break;
    }

    f = fnext;
  }

  /* now walk around the faces in the forward direction and place */
  /* them in order */

  f = v->tris[0];
  count = 0;

  for (i = 1; i < nf; i++) {

    /* find reference to vertex in f */
    vindex = -1;
    for (j = 0; j < f->nverts; j++)
      if (f->verts[(j+1) % f->nverts] == v) {
	vindex = j;
	break;
      }

    /* error check */
    if (vindex == -1) {
      fprintf (stderr, "can't find vertex #2\n");
      exit (-1);
    }

    /* corresponding face is next one around v */
    fnext = other_triangle (f->edges[vindex], f);

    /* break out of loop if we've reached a boundary */
    count = i;
    if (fnext == NULL) { break; }

    /* swap the next face into its proper place in the face list */
    for (j = 0; j < v->ntris(); j++)
      if (v->tris[j] == fnext) {
	v->tris[j] = v->tris[i];
	v->tris[i] = fnext;
	break;
      }

    f = fnext;
  }
}


/******************************************************************************
Find the index to a given vertex in the list of vertices of a given face.

Entry:
  f - face whose vertex list is to be searched
  v - vertex to return reference to

Exit:
  returns index in face's list, or -1 if vertex not found
******************************************************************************/
int Polyhedron::face_to_vertex_ref(Triangle *f, Vertex *v)
{
  int j;
  int vindex = -1;

  for (j = 0; j < f->nverts; j++)
    if (f->verts[j] == v) 
    {
      vindex = j;
      break;
    }

  return (vindex);
}

/******************************************************************************
Create various face and vertex pointers.
******************************************************************************/
void Polyhedron::create_pointers()
{
	/* create pointers from vertices to triangles */
	vertex_to_tri_ptrs();
	/* make edges */
	create_edges();
	/* order the pointers from vertices to faces */
	for (int i = 0; i < nverts(); i++)
	{
		order_vertex_to_tri_ptrs(vlist[i]);
	}
}

void Polyhedron::calc_bounding_sphere()
{
  icVector3 minV, maxV;
  for (int i=0; i<nverts(); i++) 
  {
    if (i==0)  
	{
		minV.set(vlist[i]->pos);
		maxV.set(vlist[i]->pos);
    }
    else 
	{
		if (vlist[i]->pos.x < minV.entry[0]) { minV.entry[0] = vlist[i]->pos.x; }
		if (vlist[i]->pos.x > maxV.entry[0]) { maxV.entry[0] = vlist[i]->pos.x; }
		if (vlist[i]->pos.y < minV.entry[1]) { minV.entry[1] = vlist[i]->pos.y; }
		if (vlist[i]->pos.y > maxV.entry[1]) { maxV.entry[1] = vlist[i]->pos.y; }
		if (vlist[i]->pos.z < minV.entry[2]) { minV.entry[2] = vlist[i]->pos.z; }
		if (vlist[i]->pos.z > maxV.entry[2]) { maxV.entry[2] = vlist[i]->pos.z; }
	}
  }
  center = (minV + maxV) * 0.5;
  radius = length(center - minV);
}

void Polyhedron::generateColorControlMesh() {
    // Calculate the bounding box
    icVector3 minV, maxV;

    for (int i = 0; i < nverts(); i++)
    {
        if (i == 0)
        {
            minV.set(vlist[i]->pos);
            maxV.set(vlist[i]->pos);
        }
        else
        {
            if (vlist[i]->pos.x < minV.entry[0]) { minV.entry[0] = vlist[i]->pos.x; }
            if (vlist[i]->pos.x > maxV.entry[0]) { maxV.entry[0] = vlist[i]->pos.x; }
            if (vlist[i]->pos.y < minV.entry[1]) { minV.entry[1] = vlist[i]->pos.y; }
            if (vlist[i]->pos.y > maxV.entry[1]) { maxV.entry[1] = vlist[i]->pos.y; }
            if (vlist[i]->pos.z < minV.entry[2]) { minV.entry[2] = vlist[i]->pos.z; }
            if (vlist[i]->pos.z > maxV.entry[2]) { maxV.entry[2] = vlist[i]->pos.z; }
        }
    }
    double epsilon = std::numeric_limits<double>::epsilon();
    double dx = maxV.x - minV.x + epsilon;
    double dy = maxV.y - minV.y + epsilon;
    double dz = maxV.z - minV.z + epsilon;


    // Generate colors based on vertex positions
    for (int i = 0; i < nverts(); i++) {
        Vertex* vertex = vlist[i];

        double r = (vertex->pos.x - minV.x) / dx;
        double g = (vertex->pos.y - minV.y) / dy;
        double b = (vertex->pos.z - minV.z) / dz;

        vertex->color = icVector3(r, g, b);
        //std::cout << vertex->color.x<<","<< vertex->color.y << ","<<vertex->color.z << "," << std::endl;
    }
}

double Polyhedron::determinant(const icVector3& u0, const icVector3& u1, const icVector3& u2) const {
    return u0.entry[0] * (u1.entry[1] * u2.entry[2] - u2.entry[1] * u1.entry[2]) -
        u0.entry[1] * (u1.entry[0] * u2.entry[2] - u2.entry[0] * u1.entry[2]) +
        u0.entry[2] * (u1.entry[0] * u2.entry[1] - u2.entry[0] * u1.entry[1]);
}

std::vector<std::vector<Weights>> Polyhedron::getMeanValueCoordinates(Polyhedron& controlMesh) const {
    std::vector<std::vector<Weights>> weights;
    for (const auto& modelVertex : vlist) {
        weights.push_back(getVertexMeanValues(modelVertex, controlMesh));
    }
    return weights;
}

std::vector<Weights> Polyhedron::getVertexMeanValues(const Vertex* modelVertex, Polyhedron& controlMesh) const {
    double tolerance = 0.01;
    std::vector<Weights> weights(controlMesh.nverts());

    // Initialize the vertex for each weight
    for (size_t i = 0; i < controlMesh.nverts(); ++i) {
        weights[i].vertex = controlMesh.vlist[i];
    }

    for (size_t i = 0; i < controlMesh.nverts(); ++i) {
        icVector3 p = controlMesh.vlist[i]->pos;
        weights[i].vertex = controlMesh.vlist[i];

        double d = length(p - modelVertex->pos);
        if (d <= tolerance) {
            weights[i].value = 1.0;
            return weights;
        }
    }

    for (auto& controlTriangle: controlMesh.tlist) {
        double w0, w1, w2;
        icVector3 p0 = controlTriangle->verts[0]->pos;
        icVector3 p1 = controlTriangle->verts[1]->pos;
        icVector3 p2 = controlTriangle->verts[2]->pos;
        
        size_t index0 = controlTriangle->verts[0]->index;
        size_t index1 = controlTriangle->verts[1]->index;
        size_t index2 = controlTriangle->verts[2]->index;

        double d0 = length(p0 - modelVertex->pos);
        double d1 = length(p1 - modelVertex->pos);
        double d2 = length(p2 - modelVertex->pos);

        icVector3 u0 = (p0 - modelVertex->pos); u0 /= d0;
        icVector3 u1 = (p1 - modelVertex->pos); u1 /= d1;
        icVector3 u2 = (p2 - modelVertex->pos); u2 /= d2;

        double l0 = length(u1 - u2), l1 = length(u0 - u2), l2 = length(u0 - u1);
        double theta_0 = 2 * asin(l0 / 2), theta_1 = 2 * asin(l1 / 2), theta_2 = 2 * asin(l2 / 2);
        double h = (theta_0 + theta_1 + theta_2) / 2;

        if ((M_PI - h) <= tolerance) {
            for (size_t i = 0; i < controlMesh.nverts(); ++i) {
                weights[i].value = 0.0;
            }
            double w0 = sin(theta_0) * l1 * l2;
            double w1 = sin(theta_1) * l0 * l2;
            double w2 = sin(theta_2) * l0 * l1;


            weights[index0].value = w0; weights[index0].vertex = controlTriangle->verts[0];
            weights[index1].value = w1; weights[index1].vertex = controlTriangle->verts[1];
            weights[index2].value = w2; weights[index2].vertex = controlTriangle->verts[2];

            return weights;
        }

        double c0 = (2 * sin(h) * sin(h - theta_0)) / (sin(theta_1) * sin(theta_2)) - 1;
        double c1 = (2 * sin(h) * sin(h - theta_1)) / (sin(theta_0) * sin(theta_2)) - 1;
        double c2 = (2 * sin(h) * sin(h - theta_2)) / (sin(theta_0) * sin(theta_1)) - 1;

        double det = determinant(u0, u1, u2);
        double sign;
        if (det < 0)
        {
            sign = -1;
        }
        else
        {
            sign = 1;
        }
        double s0 = sign * det * (sqrt(1 - c0 * c0)), s1 = sign * det * (sqrt(1 - c1 * c1)), s2 = sign * det * (sqrt(1 - c2 * c2));

        //If the model vertex lies outside T in the same plane, then ignore
        if (s0 <= tolerance || s1 <= tolerance || s2 <= tolerance)
        {
            continue;
        }
        
        w0 = (theta_0 - c1 * theta_2 - c2 * theta_1) / (2 * d0 * sin(theta_2) * sqrt(1 - c1 * c1));
        w1 = (theta_1 - c0 * theta_2 - c2 * theta_0) / (2 * d1 * sin(theta_0) * sqrt(1 - c2 * c2));
        w2 = (theta_2 - c1 * theta_0 - c0 * theta_1) / (2 * d2 * sin(theta_1) * sqrt(1 - c0 * c0));

        weights[index0].value += w0; weights[index0].vertex = controlTriangle->verts[0];
        weights[index1].value += w1; weights[index1].vertex = controlTriangle->verts[1];
        weights[index2].value += w2; weights[index2].vertex = controlTriangle->verts[2];
    }
    return weights;
}

void Polyhedron::InterpolateColor(const Polyhedron& Control_Mesh, Polyhedron& Model_Mesh, const std::vector<std::vector<Weights>>& Weights_Map)
{
    // Use mean value coordinates to interpolate color
    for (size_t i = 0; i < Model_Mesh.vlist.size(); ++i)
    {
        Vertex* MM_Vertex = Model_Mesh.vlist[i];
        const auto& Weight = Weights_Map[i];

        icVector3 TotalF(0.0, 0.0, 0.0);
        double TotalW = 0;
        for (const Weights& w : Weight)
        {
            const Vertex* CM_Vertex = w.vertex;
            icVector3 f(CM_Vertex->color);
            TotalW += w.value;
            TotalF += w.value * f;
        }
        icVector3 color(TotalF); color /= TotalW; // interpolate the color for the model mesh vertex
        MM_Vertex->color = color;
        //model_vertexColors[i] = color; // Insert the color for the model mesh vertex in the vector
    }
}

void Polyhedron::InterpolatePos(const Polyhedron& Control_Mesh, Polyhedron& Model_Mesh, const std::vector<std::vector<Weights>>& Weights_Map)
{
    //std::vector<icVector3> model_vertexColors(Model_Mesh.vlist.size()); // color vector for the vertices of the model mesh

    // Use mean value coordinates to interpolate position for deformation
    for (size_t i = 0; i < Model_Mesh.vlist.size(); ++i)
    {
        Vertex* MM_Vertex = Model_Mesh.vlist[i];
        const auto& Weight = Weights_Map[i];

        icVector3 TotalF(0.0, 0.0, 0.0);
        double TotalW = 0;
        for (const Weights& w : Weight)
        {
            const Vertex* CM_Vertex = w.vertex;
            icVector3 f(CM_Vertex->pos);
            TotalW += w.value;
            TotalF += w.value * f;
        }
        icVector3 pos(TotalF); pos /= TotalW; // interpolate the color for the model mesh vertex
        MM_Vertex->pos = pos;
        //model_vertexColors[i] = color; // Insert the color for the model mesh vertex in the vector
    }
}

void Polyhedron::calc_edge_length()
{
	icVector3 v1, v2;
	for (int i=0; i<nedges(); i++) 
	{
		v1.set(elist[i]->verts[0]->pos.x, elist[i]->verts[0]->pos.y, elist[i]->verts[0]->pos.z);
		v2.set(elist[i]->verts[1]->pos.x, elist[i]->verts[1]->pos.y, elist[i]->verts[1]->pos.z);
		elist[i]->length = length(v1-v2);
	}
}

void Polyhedron::calc_face_normals_and_area()
{
	icVector3 v0, v1, v2;
	Triangle *temp_t;
	double length[3];

	area = 0.0;
	for (int i=0; i<ntris(); i++)
	{
		for (int j = 0; j < 3; j++) { length[j] = tlist[i]->edges[j]->length; }
		double temp_s = (length[0] + length[1] + length[2]) * 0.5;
		tlist[i]->area = sqrt(temp_s*(temp_s-length[0])*(temp_s-length[1])*(temp_s-length[2]));

		area += tlist[i]->area;
		temp_t = tlist[i];
		v1.set(vlist[tlist[i]->verts[0]->index]->pos.x, vlist[tlist[i]->verts[0]->index]->pos.y, vlist[tlist[i]->verts[0]->index]->pos.z);
		v2.set(vlist[tlist[i]->verts[1]->index]->pos.x, vlist[tlist[i]->verts[1]->index]->pos.y, vlist[tlist[i]->verts[1]->index]->pos.z);
		v0.set(vlist[tlist[i]->verts[2]->index]->pos.x, vlist[tlist[i]->verts[2]->index]->pos.y, vlist[tlist[i]->verts[2]->index]->pos.z);
		tlist[i]->normal = cross(v0-v1, v2-v1);
		normalize(tlist[i]->normal);
	}

	double signedvolume = 0.0;
	icVector3 test = center;
	for (int i=0; i<ntris(); i++)
	{
		const icVector3& cent = vlist[tlist[i]->verts[0]->index]->pos;
		signedvolume += dot(test-cent, tlist[i]->normal)*tlist[i]->area;
	}
	signedvolume /= area;
	if (signedvolume < 0)
	{
		orientation = 0;
	}
	else 
	{
		orientation = 1;
		for (int i = 0; i < ntris(); i++) { tlist[i]->normal *= -1.0; }
	}
}

void Polyhedron::calc_vert_normals()
{
	for (int i = 0; i < nverts(); i++)
	{
		vlist[i]->normal = icVector3(0.0);
		for (int j = 0; j < vlist[i]->ntris(); j++)
		{
			vlist[i]->normal += vlist[i]->tris[j]->normal;
		}
		normalize(vlist[i]->normal);
	}
}