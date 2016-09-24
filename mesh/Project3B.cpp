#include <iostream>
#include <fstream>
#include "mesh.h"
#include <cmath>
using namespace std;

const float PI = 2*asin(1);
void HalfSphereVertexList(float radius, int nStacks, int nSlices, Point3 *pt);

int main()
{
	float r = 1;
	int nStacks = 10, nSlices = 5;
	Point3 pt[nStacks*nSlices+1];
	Face face[nStacks*nSlices];

	HalfSphereVertexList(r, nStacks, nSlices, pt);
	cout << pt[10].x << " " << pt[10].y << " " << pt[10].z << endl;
	return 0;
}

void HalfSphereVertexList(float radius, int nStacks, int nSlices, Point3 *pt)
{
	float u; // angle phi   -> angle in x-y plane
	float v; // angle theta -> angle in x-z plane
	int i, j; // loop counters
	float x, y, z;
	int nVert, nNorm, nFace;

	pt[0].set(0,0,0);
	for (i=0; i<nStacks; ++i)
	{
		v = i*(PI/nStacks);
		for (j=0; j<nSlices; ++j)
		{
			u = j*(2*PI/nSlices);
			x = radius*cos(u)*sin(v);
			y = radius*sin(u)*sin(v);
			z = radius*cos(v);
			pt[1+j+i*nSlices].set(x,y,z);
		}
	}
}

void HalfSphereFaceList(int nStacks, int nSlices, Point3 *pt, Face *face)
{
	int f, i, j;

	for (f=0; f<nSlices; ++f)
	{
		face[f].nVerts = 3;
		face[f].vert = new VertexID[face[f].nVerts];
		for (i=0; i<face[f].nVerts; ++i)
		{
			face[f].vert[i].vertIndex = i%(face[f].nVerts)+i*f;
			face[f].vert[i].normIndex = f;
		}
	}

	for (f=nSlices; f<nStacks*nSlices; ++f)
		face[f].nVerts = 4;
		face[f].vert = new VertexID[face[f].nVerts];
		face[f].vert[0].vertIndex = f-(f-1);
		face[f].vert[1].vertIndex = f+1;
		face[f].vert[2].vertIndex = f+2;
		face[f].vert[3].vertIndex = f-(f-2);
		for (i=0; i<face[f].nVerts; ++i)
			face[f].vert[i].normIndex = f;
}
