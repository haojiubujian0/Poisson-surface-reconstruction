/*
Copyright (c) 2006, Michael Kazhdan and Matthew Bolitho
All rights reserved.

Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:

Redistributions of source code must retain the above copyright notice, this list of
conditions and the following disclaimer. Redistributions in binary form must reproduce
the above copyright notice, this list of conditions and the following disclaimer
in the documentation and/or other materials provided with the distribution.

Neither the name of the Johns Hopkins University nor the names of its contributors
may be used to endorse or promote products derived from this software without specific
prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY
EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO THE IMPLIED WARRANTIES
OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT
SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
TO, PROCUREMENT OF SUBSTITUTE  GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR
BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
DAMAGE.
*/

#include "ply.h"
#include "mesh.h"
namespace poisson = poisson_reconstruction;

//
// PLY data structures
//
char* elem_names[] = { (char*)"vertex", (char*)"face" };

typedef struct PlyVertex {
	float x, y, z;
} PlyVertex;

typedef struct PlyOrientedVertex {
	float x, y, z, nx, ny, nz;
} PlyOrientedVertex;

typedef struct PlyFace {
	unsigned char nr_vertices;
	int* vertices;
	int segment;
} PlyFace;

static PlyProperty vert_props[] = {
	{(char*)"x", PLY_FLOAT, PLY_FLOAT, offsetof(PlyVertex,x), 0, 0, 0, 0},
	{(char*)"y", PLY_FLOAT, PLY_FLOAT, offsetof(PlyVertex,y), 0, 0, 0, 0},
	{(char*)"z", PLY_FLOAT, PLY_FLOAT, offsetof(PlyVertex,z), 0, 0, 0, 0}
};
static PlyProperty oriented_vert_props[] = {
	{(char*)"x",  PLY_FLOAT, PLY_FLOAT, offsetof(PlyOrientedVertex,x), 0, 0, 0, 0},
	{(char*)"y",  PLY_FLOAT, PLY_FLOAT, offsetof(PlyOrientedVertex,y), 0, 0, 0, 0},
	{(char*)"z",  PLY_FLOAT, PLY_FLOAT, offsetof(PlyOrientedVertex,z), 0, 0, 0, 0},
	{(char*)"nx", PLY_FLOAT, PLY_FLOAT, offsetof(PlyOrientedVertex,nx), 0, 0, 0, 0},
	{(char*)"ny", PLY_FLOAT, PLY_FLOAT, offsetof(PlyOrientedVertex,ny), 0, 0, 0, 0},
	{(char*)"nz", PLY_FLOAT, PLY_FLOAT, offsetof(PlyOrientedVertex,nz), 0, 0, 0, 0}
};

// List of property information for a vertex
static PlyProperty face_props[] = {
	{(char*)"vertex_indices", PLY_INT, PLY_INT, offsetof(PlyFace,vertices),
		1, PLY_UCHAR, PLY_UCHAR, offsetof(PlyFace,nr_vertices)},
};

int PlyDefaultFileType(void) { return PLY_ASCII; }

int PlyWriteTriangles(char* fileName, const poisson::Mesh_data& mesh, int file_type, double center[3], double max_stretch, double scale_ratio, char** comments, const int& commentNum) {
	int i;
	// 顶点的数目和三角形的数目
	int nr_vertices = mesh.intersection_points.size();
	int nr_faces = mesh.triangles.size();
	float version;
	PlyFile* ply = ply_open_for_writing(fileName, 2, elem_names, file_type, &version);
	if (!ply) { return 0; }

	//
	// describe vertex and face properties
	//
	ply_element_count(ply, (char*)"vertex", nr_vertices);
	ply_describe_property(ply, (char*)"vertex", &vert_props[0]);
	ply_describe_property(ply, (char*)"vertex", &vert_props[1]);
	ply_describe_property(ply, (char*)"vertex", &vert_props[2]);

	ply_element_count(ply, (char*)"face", nr_faces);
	ply_describe_property(ply, (char*)"face", &face_props[0]);

	// Write in the comments
	for (i = 0; i < commentNum; i++) { ply_put_comment(ply, comments[i]); }

	ply_header_complete(ply);

	// write vertices
	ply_put_element_setup(ply, (char*)"vertex");
	poisson::Position<float> p;
	for (i = 0; i < nr_vertices; i++) {
		PlyVertex ply_vertex;
		p = mesh.intersection_points[i];
		// 逆变换，将经过变换后的坐标进行逆变换
		ply_vertex.x = (p.x - 0.5) * (max_stretch * scale_ratio) + center[0];
		ply_vertex.y = (p.y - 0.5) * (max_stretch * scale_ratio) + center[1];
		ply_vertex.z = (p.z - 0.5) * (max_stretch * scale_ratio) + center[2];
		ply_put_element(ply, (void*)&ply_vertex);
	}

	// write faces
	ply_put_element_setup(ply, (char*)"face");
	for (i = 0; i < nr_faces; i++) {
		//
		// create and fill a struct that the ply code can handle
		//
		PlyFace ply_face;
		ply_face.nr_vertices = 3;
		ply_face.vertices = new int[3];

		for (int j = 0; j < 3; j++) { ply_face.vertices[j] = mesh.triangles[i][j]; }
		ply_put_element(ply, (void*)&ply_face);
		delete[] ply_face.vertices;
	}  // for, write faces

	ply_close(ply);
	return 1;
}
