
//This is the header file for the Point3, Vector3, VertexID, Face, 
//and Mesh classes.
#define MAC

#ifdef WINDOWS
#include <Windows.h>
#include <gl/GL.h>
#include <gl/GLU.h>
#include <gl/glut.h>
#endif

#ifdef LINUX
#include <GL/glut.h>
#endif

#ifdef MAC
#include <OpenGL/gl.h>      
#include <OpenGL/glu.h>     
#include <GLUT/glut.h>      
#endif

#include <cmath>
#include <fstream>
#include <iostream>
using namespace std;

#ifndef MeshClasses
#define MeshClasses
//This is the definition of the point class
class Point3
{
public:
	//the x, y, and z coordinates of the point
	float x, y, z;
	
	//constructors for the Point3 class
	Point3 (float xx, float yy, float zz){x = xx; y = yy; z = zz;}
	Point3() {x = y = z = 0;}
	Point3(const Point3& other) {x=other.x; y=other.y; z=other.z;}
	
	//set the point either by passing x, y, & z or by passing a point
	void set (float dx, float dy, float dz){x = dx; y = dy; z = dz;}
	void set (Point3& p){x = p.x; y = p.y; z = p.z;}
	
	//build a homogeneous point from the point
	void build4tuple(float v[])
	{
		//load 4-tuple with this color: v[3]=1 for homogeneous
		v[0] = x; 
		v[1] = y;
		v[2] = z;
		v[3] = 1.0f;
	}
};

//This is the definition for Vector3 -- a three dimensional vector
//class
class Vector3
{
public:
	//the x, y, and z direction of the vector
	float x, y, z;
	
	//constructors for the vector class
	Vector3() {x = 0; y = 0; z = 0;}
	Vector3(float xx, float yy, float zz) { x = xx; y = yy; z = zz;}
	Vector3(const Vector3& v) { x = v.x; y = v.y; z = v.z;}
	
	//set a vector's values
	void set(float dx, float dy, float dz) {x = dx; y = dy; z = dz;}
	void set (Vector3 & v) { x = v.x; y = v.y; z = v.z;}
	
	//reverse the vector
	void flip() { x = -x; y = -y; z = -z;} 
	
	//determine a vector between two points
	void setDiff (Point3& a, Point3& b)    
	{
		x = a.x - b.x;
		y = a.y - b.y;
		z = a.z - b.z;
	}
	
	//normalize a vector to unit length
	void normalize();	
	
	//determine the cross product of a vector
	Vector3 cross (Vector3 b); 
	
	//determine the dot product of a vector
	float dot (Vector3 b);
};

//################# VertexID ###################
class VertexID{
public:
	int vertIndex;		// index of this vert in the vertex list
	int normIndex;		// index of this vertex's normal
};

//#################### Face ##################
class Face{
public:
	int nVerts;			// number of vertices in this face
	VertexID * vert;	// the list of vertex and normal indices
	Face(){nVerts = 0; vert = NULL;}	// constructor
	~Face(){delete[] vert; nVerts = 0;} // destructor
};

//###################### Mesh #######################
class Mesh{
protected:
	int numVerts;		// number of vertices in the mesh
	Point3* pt;			// array of 3D vertices
	int numNorms;		// number of normal vectors for the mesh
	Vector3 *norm;		// array of normals 
	int numFaces; 		// number of faces in the mesh
	Face* face;			// array of face data
	
public:
	Mesh(); 			// constructor
	~Mesh();			// destructor
	
	// to read in a filed mesh
	int readmesh(char * fileName);
	
	void draw();		// use OpenGL to draw this mesh
	
	//For the ith face, determine a normal
	//using Newell's method.
	Vector3 newellMethod(int i);
	
};	

class ExtrudedMesh: public Mesh
{
public:
	ExtrudedMesh();
	ExtrudedMesh(Point3[], int numPts, int height);
	void draw();
};
#endif


