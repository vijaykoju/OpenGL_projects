/*
---------------------------------------------------------------------------------
FILENAME   : createMesh.cpp
PROGRAMMER : Vijay Koju
CLASS      : CSCI 7300 (Scientific Visualization and Databases)
DUE DATE   : 11/16/2013
INSTRUCTOR : Dr. Li
---------------------------------------------------------------------------------
This program gereates mesh files for a half sphere and a tapered cylinder. It is a 
supplement program for "Project3B.cpp" which uses these mesh files
"halfSphereMesh.txt" and "taperedCylinderMesh.txt" to draw them. The format of these
mesh files as given below.

To run this program:
$ make -f makeCreateMeshMac or makeCreateMeshLinux
$ ./createMesh

---------------------------------------------------------------------------------
data fields                     example
numVertics numFaces numNorms | 901 900 900
                             |
list of x,y,z components of  | 0 1 2   0 3 2   2 4 2   0.3 3 9 ...
vertices                     | ...
                             |
list or x,y,z components of  | 0.22 0.4 0.5   0 0 1   0.1 0.4 0.33 ...
normals to each face         | ...
                             |
list of r,g,b color for each | 0.2 0.5 0.5   0.4 0.55 0.1   0 0.43 0.2 ...
face                         | ...
                             |
num of vertices on face 1    | 4
VertexIndex and normIndex of | 0 0
each vertex on face 1        | 1 0
                             | 9 0
                             | 8 0
                             | 
num of vertices on face 2    | ...
---------------------------------------------------------------------------------
*/
#include <iostream>
#include <fstream>
#include "mesh.h"
#include <cmath>
using namespace std;

const float PI = 2*asin(1); // constant pi

// function prototypes
// vertex list and color-value for each vertex
void HalfSphereVertexList(float radius, int nStacks, int nSlices, Point3 *pt, Point3 *color);
void TaperedCylinderVertexList(float s, int nStacks, int nSlices, Point3 *pt, Point3 *color);
// face list with number of vertices and normals associted with each face and color
// color-value for each face 
void HalfSphereFaceList(int nStacks, int nSlices, Point3 *pt, Face *face, Vector3 *norm, Point3 *v_color, Point3 *f_color);
void TaperedCylinderFaceList(int nStacks, int nSlices, Point3 *pt, Face *face, Vector3 *norm, Point3 *v_color, Point3 *f_color);

// main function
int main()
{
	// variables and data fields for half sphere
	float    r         = 1.2;                 // radius
	int      nStacks   = 30;                  // # of stacks
	int      nSlices   = 30;                  // # of slices
	int      numVertex = nStacks*nSlices+1;   // # of vertices
	int      numFaces  = nStacks*nSlices;     // # of faces. 
	Point3   pt[numVertex];                   // array of vertex points
	Face     face[numFaces];                  // array of faces
	Vector3  norm[numFaces];                  // array of norms
	Point3   v_color[numVertex];              // array of vertex color
	Point3   f_color[numFaces];               // array of face color
	ofstream meshData;                        // output to a file
	// variables and data fields for tapered cylinder
	float    s         = 0.4;                 // radius of tapered end
	int      cStacks   = 60;                  // # of stacks
	int      cSlices   = 30;                  // # of slices
	int      cVertex   = cSlices*(cStacks+1); // # of vertices
	int      cFaces    = cSlices*cStacks;     // # of of faces
	Point3   cPt[cVertex];                    // array of vertex points
	Face     cFace[cFaces];                   // array of faces
	Vector3  cNorm[cFaces];                   // array of normals
	Point3   cv_color[cVertex];               // array of vetex color
	Point3   cf_color[cFaces];                // array of face color
	ofstream cMeshData;                       // output to a file

	// create Vertex list for half sphere
	HalfSphereVertexList(r, nStacks, nSlices, pt, v_color);
	// create Face, normal, and color list for half sphere
	HalfSphereFaceList(nStacks, nSlices, pt, face, norm, v_color, f_color);
	// create Vertex list for tapered cylinder 
	TaperedCylinderVertexList(s,cStacks,cSlices,cPt, cv_color);
	// create Face, normal, and color list for tapered cylinder 
	TaperedCylinderFaceList(cStacks,cSlices,cPt,cFace, cNorm, cv_color, cf_color);
	
	//-------------------------------------------------------------------//
	meshData.open("halfSphereMesh.txt"); // open output file for half sphere
	meshData << endl;
	// number of vertices   number of faces   number of normals
	meshData << numVertex << " " << numFaces << " " << numFaces << endl;
	//-----------------------Start Vertex List---------------------------//
	meshData << pt[0].x << " " << pt[0].y << " " << pt[0].z <<  "   ";
	for (int i=0; i<nStacks; ++i)
	{
		for (int j=0; j<nSlices; ++j)
			meshData << pt[(j+1)+i*nSlices].x << " " << pt[(j+1)+i*nSlices].y << " " << pt[(j+1)+i*nSlices].z << "   ";
		meshData << endl;
	}
	meshData << endl;
	//-------------------------End Vertex List---------------------------//
	//------------------------Start Normal List--------------------------//
	for (int i=0; i<nStacks; ++i)
	{
		for (int j=0; j<nSlices; ++j)
			meshData << norm[(j)+i*nSlices].x << " " << norm[(j)+i*nSlices].y << " " << norm[(j)+i*nSlices].z << "   ";
		meshData << endl;
	}
	meshData << endl;
	//--------------------------End Normal List--------------------------//
	//-------------------------Start Color List--------------------------//
	for (int i=0; i<nStacks; ++i)
	{
		for (int j=0; j<nSlices; ++j)
			meshData << f_color[(j)+i*nSlices].x << " " << f_color[(j)+i*nSlices].y << " " << f_color[(j)+i*nSlices].z << "   ";
		meshData << endl;
	}
	meshData << endl;
	//---------------------------End Color List--------------------------//
	//--------------------------Start Face List--------------------------//
	for (int f=0; f<numFaces; ++f)
	{
		meshData << face[f].nVerts << endl;
		for (int i=0; i<face[f].nVerts; ++i)
			meshData << face[f].vert[i].vertIndex << " " << face[f].vert[i].normIndex << endl;
		meshData << endl;
	}
	//----------------------------End Face List--------------------------//
	meshData.close(); // close output file for half sphere

	//-------------------------------------------------------------------//
	cMeshData.open("taperedCylinderMesh.txt"); // create output file for
	cMeshData << endl;                         // tapered cylinder
	// number of vertices   number of faces   number of normals
 	cMeshData << cVertex << " " << cFaces << " " << cFaces << endl;
	//-----------------------Start Vertex List---------------------------//
  for (int i=0; i<cStacks+1; ++i)
  { 
    for (int j=0; j<cSlices; ++j)
      cMeshData << cPt[j+i*cSlices].x << " " << cPt[j+i*cSlices].y << " " << cPt[j+i*cSlices].z << "   ";
    cMeshData << endl;
  }
  cMeshData << endl;
	//-------------------------End Vertex List---------------------------//
	//------------------------Start Normal List--------------------------//
	for (int i=0; i<cStacks; ++i)
  {
    for (int j=0; j<cSlices; ++j)
      cMeshData << cNorm[j+i*cSlices].x << " " << cNorm[j+i*cSlices].y << " " << cNorm[j+i*cSlices].z << "   ";
    cMeshData << endl;
  }
  cMeshData << endl; 
	//--------------------------End Normal List--------------------------//
	//-------------------------Start Color List--------------------------//
	for (int i=0; i<cStacks; ++i)
  {
    for (int j=0; j<cSlices; ++j)
      cMeshData << cf_color[j+i*cSlices].x << " " << cf_color[j+i*cSlices].y << " " << cf_color[j+i*cSlices].z << "   ";
    cMeshData << endl;
  }
  cMeshData << endl; 
	//---------------------------End Color List--------------------------//
	//--------------------------Start Face List--------------------------//
	for (int f=0; f<cFaces; ++f)
  {
    cMeshData << cFace[f].nVerts << endl;
    for (int i=0; i<cFace[f].nVerts; ++i)
      cMeshData << cFace[f].vert[i].vertIndex << " " << cFace[f].vert[i].normIndex << endl;
    cMeshData << endl;
  }
	//----------------------------End Face List--------------------------//
 	cMeshData.close(); // close output file for tapered cylinder
 	return 0;
 }
 
void HalfSphereVertexList(float radius, int nStacks, int nSlices, Point3 *pt, Point3 *color)
{
	int   i, j, indx;       // counters
	float u;                // angle phi   -> angle in x-y plane
	float v;                // angle theta -> angle in x-z plane
	float x, y, z;          // x, y, and z components of vertices
	float r, g, b;          // r, g, and b color faces of vertices
	float c, d;             // color function value and sort of distance

	// set (x,y,z) coordinates of each vertex for half sphere
	pt[0].set(0,0,1+radius);      // pole of the half sphere
	for (i=1; i<nStacks+1; ++i)
	{
		v = i*((PI/2)/nStacks);     // step-size for v
		for (j=0; j<nSlices; ++j)
		{
			u = (j)*(2*PI/(nSlices)); // step-size for u
			x = radius*cos(u)*sin(v); // x-component
			y = radius*sin(u)*sin(v); // y-component
			z = 1+radius*cos(v);      // z-component
			indx = 1+j+(i-1)*nStacks; // index of each vertex
			pt[indx].set(x,y,z);      // set vertex
		}
	}

	// set (r,g,b) color value for each vertex	
	d = sqrt(pt[0].x*pt[0].x+pt[0].y*pt[0].y+pt[0].z*pt[0].z);
	color[0].set(0,0,(1+radius)/d); // color value for the pole
	for (i=1; i<nStacks+1; ++i)
	{	
		v = i*((PI/2)/nStacks);       // step-size for u
		for (j=0; j<nSlices; ++j)
		{
			u = j*(2*PI/nSlices);       // step-size for v
			c = u*sin(1+v)+1;           // color function
			indx = 1+j+(i-1)*nStacks;   // index of a vertex
			d = sqrt(c*c*pt[indx].x*pt[indx].x+c*c*pt[indx].y*pt[indx].y+c*c*pt[indx].z*pt[indx].z);
			r = abs(c*pt[indx].x/d);    // r-color value
			g = abs(c*pt[indx].y/d);    // g-color value
			b = abs(c*pt[indx].z/d);    // b-color value
			color[indx].set(r,g,b);     // set color
		}
	}	
}

void HalfSphereFaceList(int nStacks, int nSlices, Point3 *pt, Face *face, Vector3 *norm, Point3 *v_color, Point3 *f_color)
{
	int f, i, j, nextJ;     // counters
	int N;                  // # of vertices on each face
	int index, nextIndex;   // index of current and next vertex

	// set vertexIndex and normIndex associated with each face
	// each face in the first stack of half sphere only has three
	// vertices, so the first handled separately than other stacks
	for (f=0; f<nSlices; ++f)
	{
		face[f].nVerts = 3; // # of vertices on each face
		face[f].vert = new VertexID[face[f].nVerts]; // allocate memory
		// set vertexIndex
		face[f].vert[0].vertIndex = 0; 
		face[f].vert[1].vertIndex = f+1;
		face[f].vert[2].vertIndex = 1+(f+1)%(nSlices);
		// set normIndex
		for (i=0; i<face[f].nVerts; ++i)
			face[f].vert[i].normIndex = f;
	}

	// set vertexIndex and normIndex for the rest of faces
	for (i=1; i<nStacks; ++i)
	{
		for (j=0; j<nSlices; ++j)
		{
			f = i*nSlices+j; // current face index
			face[f].nVerts = 4; // # of vetices on each face
			face[f].vert = new VertexID[face[f].nVerts]; // allocate memory
			// set vertexIndex
			face[f].vert[0].vertIndex = (i-1)*nSlices+(j+1);
			face[f].vert[1].vertIndex = ((i-1)*nSlices+(j+1))+nSlices;
			nextJ = (((i)*nSlices+(j+1)))%nSlices+1;
			face[f].vert[2].vertIndex = nextJ + i*nSlices;
			face[f].vert[3].vertIndex = nextJ + (i-1)*nSlices;
			// set normIndex
			for (int k=0; k<face[f].nVerts; ++k)
				face[f].vert[k].normIndex = f;
		}
	}

	// compute normal and color-value for each face
	for (f=0; f<nStacks*nSlices; ++f)
	{
		N = face[f].nVerts; // # of vertices on each face
		// compute normal using Newell method
		for ( i=0; i<N; ++i)
		{
			index = face[f].vert[i].vertIndex;
			nextIndex = face[f].vert[(i+1)%N].vertIndex;
			norm[f].x += (pt[index].y-pt[nextIndex].y)*(pt[index].z+pt[nextIndex].z);
			norm[f].y += (pt[index].z-pt[nextIndex].z)*(pt[index].x+pt[nextIndex].x);
			norm[f].z += (pt[index].x-pt[nextIndex].x)*(pt[index].y+pt[nextIndex].y);
		}	
		norm[f].normalize(); // normalize the normal vector
		// compute color-value for each face as an average of color-values
		// of vertices associated to that face
		for (i=0; i<N; ++i)
		{
			index = face[f].vert[i].vertIndex;
			f_color[f].x += v_color[index].x/N; // red
			f_color[f].y += v_color[index].y/N; // green
			f_color[f].z += v_color[index].z/N;	// blue
		}
	}
}

void TaperedCylinderVertexList(float s, int nStacks, int nSlices, Point3 *pt, Point3 *color)
{
	int   i, j, indx;       // counters
	float u;                // angle phi   -> angle in x-y plane
	float v;                // angle theta -> angle in x-z plane
	float x, y, z;          // x, y, and z components of vertices
	float r, g, b;          // r, g, and b color faces of vertices
	float c, d;             // color function value and sort of distance

	// set (x,y,z) coordinates of each vertex for tapered cylinder
	for (i=0; i<nStacks+1; ++i)
	{
		v = i*1.0/(nStacks);       // step-size for v
		for (j=0; j<nSlices; ++j) 
		{
			u = j*(2*PI/nSlices);    // step-size for u
			x = (1+(s-1)*v)*cos(u);  // x-component
			y = (1+(s-1)*v)*sin(u);  // y-component
			z = v;                   // z-component
			pt[j+i*nSlices].set(x,y,z); // set vertex
		}
	}

	// set (r,g,b) color value for each vertex
	for (i=0; i<nStacks+1; ++i)
	{	
		v = i*1.0/nStacks;         // step-size for v
		for (j=0; j<nSlices; ++j)
		{
			u = j*(2*PI/nSlices);    // step-size for u
			c = u*sin(1+v)+1;        // color function
			indx = j+i*nStacks;      // index of current vertex
			d = sqrt(c*c*pt[indx].x*pt[indx].x+c*c*pt[indx].y*pt[indx].y+c*c*pt[indx].z*pt[indx].z);
			r = abs(c*pt[indx].x/d); // red
			g = abs(c*pt[indx].y/d); // green
			b = abs(c*pt[indx].z/d); // blue
			color[indx].set(r,g,b);  // set color
		}
	}	
}

void TaperedCylinderFaceList(int nStacks, int nSlices, Point3 *pt, Face *face, Vector3 *norm, Point3 *v_color, Point3 *f_color)
{
	int f, i, j, nextJ;     // counters
	int N;                  // # of vertices on each face
	int index, nextIndex;   // index of current and next vertex

	// set vertexIndex and normIndex associated with each face
	for (i=0; i<nStacks; ++i)
	{
		for (j=0; j<nSlices; ++j)
		{
			f = i*nSlices+j;    // index of the current face
			face[f].nVerts = 4; // # of vertices on each face
			face[f].vert = new VertexID[face[f].nVerts]; // allocate memory
			// set vetexIndex
			face[f].vert[0].vertIndex = f;
			face[f].vert[1].vertIndex = (f+1)%nSlices+i*nSlices;
			nextJ = (f+1)%nSlices;
			face[f].vert[2].vertIndex = nextJ+(i+1)*nSlices;
			face[f].vert[3].vertIndex = f+nSlices;
			// set normIndex
			for (int k=0; k<face[f].nVerts; ++k)
				face[f].vert[k].normIndex = f;
		}
	}

	// compute normal and color-value for each face
	for (f=0; f<nStacks*nSlices; ++f)
  {
    N = face[f].nVerts; // # of vertices on each face
		// compute normal using Newell's method
    for ( i=0; i<N; ++i)
    {
      index = face[f].vert[i].vertIndex; // current vertex index
      nextIndex = face[f].vert[(i+1)%N].vertIndex; // next vertex index
      norm[f].x += (pt[index].y-pt[nextIndex].y)*(pt[index].z+pt[nextIndex].z);
      norm[f].y += (pt[index].z-pt[nextIndex].z)*(pt[index].x+pt[nextIndex].x);
      norm[f].z += (pt[index].x-pt[nextIndex].x)*(pt[index].y+pt[nextIndex].y);
    }
    norm[f].normalize(); // normalize the normal vector
		// compute color-value for each face as an average of colo-values
		// of all the vertices associated to that face
		for (i=0; i<N; ++i)
		{
			index = face[f].vert[i].vertIndex;  // current vertex index
			f_color[f].x += v_color[index].x/N; // red
			f_color[f].y += v_color[index].y/N; // green
			f_color[f].z += v_color[index].z/N;	// blue
		}
  }
}
