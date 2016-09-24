
//This program demonstrates the use of Hill's mesh class.
//The mesh class is used to read a mesh from a file and 
//draw it.
//
//If the user clicks the left mouse button, the barn will
//be rotated by 10 degrees each time.

#include "mesh.h"

int spin = 0;           //allow the user to rotate the barn

//This is the display function.  Each time it draws a barn,
//it is rotated by an angle, spin.
void display()
{
	//set the background to black and clear the buffers
	glClearColor(0.5, 0.5, 0.5, 1.0);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	
	//create a mesh for the barn
	Mesh mymesh;
	
	//read the mesh from the file
	//mymesh.readmesh("BarnMeshFile.txt");
	mymesh.readmesh("BarnMeshFile.txt");
	
	//rotate the barn by spin degrees
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	glRotatef(spin, 1.0, 0.0, 0.0);
	
	//draw the mesh
	mymesh.draw();
	
	//make it appear on the screen
	glFlush();
}

//This function initializes the lighting and material properties.
//There is one positional light in the scene.
void myinit(void)
{
	//set up material and lighting arrays
	GLfloat ambient[] = {1.0, 1.0, 1.0, 1.0};
	GLfloat position[] = {2.0, 2.0, 2.0, 1.0};
	GLfloat mat_diffuse[] = {0.6, 0.6, 0.6, 1.0};
	GLfloat mat_specular[] = {1.0, 1.0, 1.0, 1.0};
	GLfloat mat_shininess[] = {128.0};
	
	//enable lighting
	glEnable(GL_LIGHTING);
	glEnable(GL_LIGHT0);
	
	//set the light properties
	glLightfv(GL_LIGHT0, GL_AMBIENT, ambient);
	glLightfv (GL_LIGHT0, GL_DIFFUSE, mat_diffuse);
	glLightfv(GL_LIGHT0, GL_POSITION, position);
	glLightf(GL_LIGHT0,GL_SPOT_EXPONENT,128);
	
	//set the material properties
	mat_diffuse[2]=0.0; mat_diffuse[1] = 0.0;
	glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, mat_diffuse);
	glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, mat_specular);
	glMaterialfv(GL_FRONT_AND_BACK, GL_SHININESS, mat_shininess);
	
	//enable hidden surface removals, etc.
	glEnable(GL_DEPTH_TEST);
	glEnable(GL_NORMALIZE);
	glShadeModel (GL_SMOOTH);
}

//This is the mouse callback function.  Clicking the left
//mouse button causes the angle contained in the variable
//spin to be adjusted by 10 degrees and the window to
//be redisplayed.
void mouse(int button, int state, int x, int y)
{
    switch (button) {
		case GLUT_LEFT_BUTTON:
			if (state == GLUT_DOWN) {
				spin = (spin + 10) % 360;
				glutPostRedisplay();
			}
			break;
		default:
			break;
	}
}
//void myKeyboard(unsigned char theKey, int x, int y)
//{
//	switch(theKey)
//	{
//		case 'q':   // end display
//			exit(0);
//		default:
//			if (theKey == 27)   // ASCII for escape character
//				exit(0);
//	}
//}

//This is the reshape callback.
void myReshape(int w, int h)
{
	//set the viewport to the entire window
	glViewport(0, 0, w, h);
	
	//use an orthographic view volume of a cube
	//centered about the origin
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	if (w <= h)
		glOrtho(-13.0, 13.0, -13.0*(GLfloat)h/(GLfloat)w, 
				13.0*(GLfloat)h/(GLfloat)w, -13.0, 13.0);
	else
		glOrtho(-13.0*(GLfloat)w/(GLfloat)h, 
				13.0*(GLfloat)w/(GLfloat)h, -13.0, 13.0, -13.0, 13.0);
	
}

int main(int argc, char** argv)
{
	//create the window
	glutInit(&argc, argv);
	glutInitDisplayMode (GLUT_SINGLE | GLUT_RGB | GLUT_DEPTH);
	glutInitWindowSize (500, 500); 
	glutInitWindowPosition (100, 100);
	glutCreateWindow (argv[0]);
	
	//initialize lighting & material properties
	myinit();
	
	//set up callbacks
	glutMouseFunc(mouse);
	glutReshapeFunc (myReshape);
	glutDisplayFunc(display);
  //glutKeyboardFunc(myKeyboard);
	
	//start the main loop
	glutMainLoop();
	return 0;             
}
