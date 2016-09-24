
//This is the implementation file for the Vector3 and Mesh classes.

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

#include "extrudedMesh.h"

//Determine and return the cross product of two vectors.
Vector3 Vector3::cross (Vector3 b) // return this cross b
{
	Vector3 c (y*b.z - z*b.y, z*b.x - x*b.z, x*b.y - y*b.x);
	return c;
}

//Determine and return the dot product of two vectors.
float Vector3::dot (Vector3 b)     //return this dotted with b
{
	return x * b.x + y * b.y + z * b.z;
}

//Normalize the vector so that the vector's length is one.
void Vector3::normalize()          
{
	double sizeSq = x * x + y * y + z * z;  //the current magnitude
	
	//only normalize if the vector isn't the zero vector
	if (sizeSq < 0.0000001)
	{
		cerr << "\nnormalize() sees vector (0,0,0)!";
		return;
	}
	
	//normalize each component of the vector
	float scaleFactor = 1.0 / (float) sqrt (sizeSq);
	x *= scaleFactor;
	y *= scaleFactor;
	z *= scaleFactor;
}

//The default constructor for the mesh class
//sets the number of vertices, normals, and
//faces to 0.
Mesh::Mesh() 			
{
	numVerts=0;
	numNorms=0;
	numFaces=0;
	norm=NULL;
	pt=NULL;
	face=NULL;
}

//The Mesh destructor releases all the space
//allocated to the mesh and sets the number
//of vertices, normals, and faces back to 0.
Mesh::~Mesh()			
{
	delete[] pt;
	numVerts=0;
	delete[] norm;
	numNorms=0;
	delete[] face;
	numFaces=0;
}

//Determine the normal for the ith
//face using Newell's Method.
Vector3 Mesh::newellMethod(int j)
{
	int N = face[j].nVerts;
	int index, nextIndex;
	Vector3 m;
	
	for (int i = 0; i < N; i++)
	{
		//get the next two indices of vertices
		//in the face
		index=face[j].vert[i].vertIndex;
		nextIndex = face[j].vert[(i+1)% N].vertIndex;
		m.x += (pt[index].y - pt[nextIndex].y)*(pt[index].z + pt[nextIndex].z);
		m.y += (pt[index].z - pt[nextIndex].z)*(pt[index].x + pt[nextIndex].x);
		m.z += (pt[index].x - pt[nextIndex].x)*(pt[index].y + pt[nextIndex].y);
	}
	m.normalize();
	return m;
}

//Draw the mesh.  Each face of the object is drawn
//using a different material property.
void Mesh:: draw() 
{
	//set up the beginning material properties
	GLfloat mat_diffuse[] = {0.6, 0.6, 0.6, 1.0};
        GLfloat mat_specular[] = {1.0, 1.0, 1.0, 1.0};
        GLfloat mat_shininess[] = {50.0};
        glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, mat_specular);
        glMaterialfv(GL_FRONT_AND_BACK, GL_SHININESS, mat_shininess);   

	//loop through the faces of the object
	for(int f = 0; f < numFaces; f++) 
	{ 
		//adjust the diffuse material property slightly for each face
		mat_diffuse[2]=0.0 + f * .1; mat_diffuse[1] = 0.0 + .02 * f;
                glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, mat_diffuse);
		
		//draw the face
		glBegin(GL_POLYGON);
		for(int v = 0; v < face[f].nVerts; v++) // for each one..
		{
			//find the next normal and vertex
			int in = face[f].vert[v].normIndex ; // index of this normal
			int iv =  face[f].vert[v].vertIndex ; // index of this vertex
			
			//inform OpenGL of the normal and vertex
			glNormal3f(norm[in].x, norm[in].y, norm[in].z);
			glVertex3f(pt[iv].x, pt[iv].y, pt[iv].z);
		}
		glEnd();
	}
}

//This function reads face information from a data file.
//The name of the file is passed to the function through 
//the argument list
int Mesh:: readmesh(char * fileName)
{
	//open the file and check for file failure
	fstream infile;
	infile.open(fileName, ios::in);
	if(infile.fail()) return -1; // error - can't open file
	if(infile.eof())  return -1; // error - empty file
	
	//the file is OK so read the number of vertices,
	//normals, and faces.
	infile >> numVerts >> numNorms >> numFaces;
	
	//create arrays to hold the vertices, nomrmals,
	//and faces.
	pt = new Point3[numVerts];
	norm = new Vector3[numNorms];
	face = new Face[numFaces];
	
	//check that enough memory was found:
	if( !pt || !norm || !face)return -1; 
	
	//read the vertices
	for(int p = 0; p < numVerts; p++) 
	{
		infile >> pt[p].x >> pt[p].y >> pt[p].z;
		
	}
	
	//read the normals
	for(int n = 0; n < numNorms; n++) 
	{
		infile >> norm[n].x >> norm[n].y >> norm[n].z;	
	}
	
	//read the faces
	for(int f = 0; f < numFaces; f++)
	{
		infile >> face[f].nVerts;
		face[f].vert = new VertexID[face[f].nVerts];
		for(int i = 0; i < face[f].nVerts; i++)
			infile >> face[f].vert[i].vertIndex 
			>> face[f].vert[i].normIndex;
	} 
	return 0; // success
}

//This is the default constructor of the Extruded Mesh.
ExtrudedMesh::ExtrudedMesh():Mesh()
{}

// This assumes the vertices in the base shape are specified in 
// CCW order
ExtrudedMesh::ExtrudedMesh(Point3 array[], int N, int height)
{
	//set the total number of vertices to N
	numVerts = 2*N;
	
	//allocate enough space for the vertex array
	pt = new Point3[2*N];
	
	//create the vertex array
	//the first N vertices come from the base
	for (int i = 0; i < N; i++)
		pt[i] = array[i];
	
	//the next N vertices come from the cap
	for (int i = 0; i < N ; i++)
	{
		pt[i+N] = array[i];
		pt[i+N].z += height;
	}
	
	//set the number of normals and faces and
	//allocate space for the normal vectors
	numNorms = N + 2;
	numFaces = N + 2;
	norm = new Vector3[N + 2];
	face = new Face[N + 2];
	
	//Create the face list
	//create side faces first.  N quadrilaterals
	//make up the side faces.
	int nextJ;
	
	for (int j = 0; j < N; j++)
	{
		face[j].nVerts = 4;
		
		face[j].vert = new VertexID[4];
		face[j].vert[0].vertIndex = j;
		face[j].vert[1].vertIndex = j+N;
		nextJ = (j+1)%N;
		face[j].vert[2].vertIndex = nextJ + N;
		face[j].vert[3].vertIndex = nextJ;
		
		//use Newell's method to determine the
		//normal vector to each face
		norm[j] = newellMethod(j);
		for (int i = 0; i <= 3; i++)
			face[j].vert[i].normIndex = j;
	}

	//add the base face to the face list as the Nth face
	face[N].nVerts = N;
	face[N].vert = new VertexID[N];
	for (int i = 0; i < N; i++)
		face[N].vert[i].vertIndex = i;
	norm[N] = newellMethod(N);
	for (int i = 0; i < N; i++)
		face[N].vert[i].normIndex = N;
	
	//add the cap face to the face list as the N+1st face
	face[N + 1].nVerts = N;
	face[N + 1].vert = new VertexID[N];
	for (int i = 0; i < N; i++)
		face[N + 1].vert[i].vertIndex = N + i;
	norm[N + 1] = newellMethod(N+1);
	for (int i = 0; i < N; i++)
		face[N + 1].vert[i].normIndex = N + 1;
}

void ExtrudedMesh::draw()
{
	Mesh::draw();
}
