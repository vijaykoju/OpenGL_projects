//This is the implementation file for the Vector3 and Mesh classes.

#include "mesh.h"
#include <cassert>
#include <iostream>
using namespace std;

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
Mesh::Mesh() 			
{
	numVerts=0;
	numNorms=0;
	numFaces=0;
	norm=NULL;
	pt=NULL;
	face=NULL;
	f_color=NULL;
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

//Draw the mesh.  Each face of the object is drawn
//using a different material property.
void Mesh:: draw() 
{
    //set up the beginning material properties
    GLfloat mat_diffuse[] = {0, 0, 0, 1.0};
		//GLfloat b = 0.1;
		GLfloat blue[] = {1,0,0};
    GLfloat mat_specular[] = {1.0, 1.0, 1.0, 1.0};
    GLfloat mat_shininess[] = {50.0};
    glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, mat_specular);
    glMaterialfv(GL_FRONT_AND_BACK, GL_SHININESS, mat_shininess);   

    //loop through the faces of the object
    for(int f = 0; f < numFaces; f++) 
    { 
	//adjust the diffuse material property slightly for each face
		//mat_diffuse[0] = f_color[f].x;
	 	//mat_diffuse[1] = f_color[f].y;
		//mat_diffuse[2] = f_color[f].z;
	 //glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, mat_diffuse);
	//blue[2] += 0.2; 	
	//glMaterialfv(GL_FRONT,GL_AMBIENT,blue);
	//draw the face
	glColor3f(f_color[f].x,f_color[f].y,f_color[f].z);
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
int Mesh:: readmesh(string fileName)
{
	//open the file and check for file failure
	ifstream infile;
	infile.open(fileName.c_str());
	assert(infile);
	if(infile.eof())  return -1; // error - empty file
	
	//the file is OK so read the number of vertices,
	//normals, and faces.
	infile >> numVerts >> numNorms >> numFaces;
	
	//create arrays to hold the vertices, nomrmals,
	//and faces.
	pt = new Point3[numVerts];
	norm = new Vector3[numNorms];
	face = new Face[numFaces];
	f_color = new Point3[numFaces];
	
	//check that enough memory was found:
	if( !pt || !norm || !face || !f_color)return -1; 
	
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
	
	for (int c =0; c < numFaces; c++)
	{
		infile >> f_color[c].x >> f_color[c].y >> f_color[c].z;
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


