//This is the header file for the VertexID, Face, 
//and Mesh classes.

#include <fstream>
#include <string>
//#include "Camera.h"
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
using namespace std;
//This is the definition of the point class
class Point3
{
public:
        //the x, y, and z coordinates of the point
        float x, y, z;

        //constructors for the Point3 class
        Point3 (float xx, float yy, float zz){x = xx; y = yy; z = zz;}
        Point3() {x = y = z = 0;}

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

//////////////////////////////////////
//The Vector3 class
//////////////////////////////////////

class Vector3                
{
public:
    float x, y, z; 

    Vector3() {}   // constructor
    Vector3(float new_x, float new_y, float new_z) // initialize constructor   
    	{x = new_x; y = new_y; z = new_z;}

    void set(float dx, float dy, float dz) 
        {x = dx; y = dy; z = dz;}

//    void set(Vector3 & v)  
 //       {x = v.x; y = v.y; z = v.z;}

    void set(Vector3 v)  
        {x = v.x; y = v.y; z = v.z;}

    // overload + operator so that we easier can add vectors
    Vector3 operator+(Vector3 vVector) 
        {return Vector3(vVector.x+x, vVector.y+y, vVector.z+z);}

    // overload - operator that we easier can subtract vectors
    Vector3 operator-(Vector3 vVector) 
        {return Vector3(x-vVector.x, y-vVector.y, z-vVector.z);}

    // overload * operator that we easier can multiply by scalars
    Vector3 operator*(float number)     
        {return Vector3(x*number, y*number, z*number);}

    // overload / operator that we easier can divide by a scalar
    Vector3 operator/(float number)     
        {return Vector3(x/number, y/number, z/number);}

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
private:
	int numVerts;		// number of vertices in the mesh
	Point3* pt;			// array of 3D vertices
	int numNorms;		// number of normal vectors for the mesh
	Vector3 *norm;		// array of normals 
	int numFaces; 		// number of faces in the mesh
	Face* face;			// array of face data
	Point3* f_color;
	// ... others to be added later
public:
	Mesh(); 			// constructor
	~Mesh();			// destructor
	
	// to read in a filed mesh
	int readmesh(string fileName);
	void draw();		// use OpenGL to draw this mesh
	//.. others ..
};	
