/*
-----------------------------------------------------------------------------
FILENAME   : halfSpheres.cpp
PROGRAMMER : Vijay Koju
-----------------------------------------------------------------------------
This program draws two objects (a half sphere and a tapered cylinder) using
mesh approximation. It reads in the mesh files for a half sphere and tapered
cylinder separately, and draws them together.
The color value for each face of the objects are set using a color function
f(u,v) = u*sin(1+v)+1.
The mesh files "halfSphereMesh.txt" and "taperedCylinderMesh.txt" are created
using another program "createMesh.cpp".

To create mesh files:
$ make makeCreateMeshMac or makeCreateMeshLinux
$ ./createMesh

To run this program:
$ make makeProject3BMac or makeProject3BLinux
$ ./Project3B
-----------------------------------------------------------------------------
Keyboard interactions:
m or M  --> set move mode
        left  arrow -> move left
				right arrow -> move right
				up arrow    -> move up
				down arrow  -> move down
r or R  --> set rotation mode
				left arrow  -> rotate leftward
				right arrow -> rotate rightward
				up arrow    -> rotate upward
				down arrow  -> rotate downward
a or A  --> prompt user to define point of rotation
						and start animation
s or S  --> stop/pause animation
f or F  --> move forward  / zoom in
b or B  --> move backward / zoom out
-----------------------------------------------------------------------------
*/

#include "Camera.h"
#include <cstdlib>
#include <iostream>
using namespace std;

#define WIDTH 640          // window width
#define HEIGHT 480         // window height

#define startupCameraX 0   // initial x-position of camera
#define startupCameraY 1.1 // initial y-position of camera
#define startupCameraZ 10  // initial z-position of camera

// rotation amount
GLdouble xRot;             // roll angle for animation
GLdouble  x_r, y_r, z_r;   // variable to store user define point of
                           // rotation for animation

// starting camera position
GLdouble  cameraX=startupCameraX;
GLdouble  cameraY=startupCameraY;
GLdouble  cameraZ=startupCameraZ;

// aspect ratio
GLdouble aspectRatio=WIDTH/HEIGHT;
GLdouble ww_x = 14.0*aspectRatio;
GLdouble ww_y = 14.0;

// boolean variables to control different functions
bool animation=false; // control animation
bool moveLRTB=false;  // control left right top bottom motion
bool rotate=false;    // control left right top bottom rotation

Camera myCamera;  // camera object

// function prototypes
void MyInit();                           // set up camera and depth
void Draw();                           	 // draw half sphere and tapered cylinder
void SpecialKeys(int key, int x, int y); // special key functions
void TimerFunction(int value);           // timer function
void myKeyboard(unsigned char theKey, int x, int y); // keyboard functions

// main function
int main(int argc, char **argv)
{
	glutInit(&argc, argv);               // initialize glut varialbles
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH); // set displaymode
	glutInitWindowSize(WIDTH, HEIGHT);   // set window size
	glutInitWindowPosition(100, 100);    // set window position
	glutCreateWindow("Project 3B");      // set window title
	glutDisplayFunc(Draw);               // callback to Draw() function
	glutSpecialFunc(SpecialKeys);        // callback to SpecialKey()
	glutTimerFunc(40, TimerFunction, 1); // callback to TimerFunction()
	glutKeyboardFunc(myKeyboard);        // callback to myKeyboard()
	glViewport(0, 0, WIDTH, HEIGHT);     // set the whole window as viewport
	MyInit();                            // set camera and depth test
	glutMainLoop();                      // enter an infinite loop and wait for events
	return 0;
}

// set initial camera setting and enable depth test
void MyInit()
{
  myCamera.setShape(45, 1, 1, 100);
  myCamera.set(cameraX, cameraY, cameraZ,   0, 0, 0,   0, 1, 0);
	glClearColor(0.4f, 0.4f, 0.4f, 0.0f);
	glEnable(GL_DEPTH_TEST);
}	

//<<<<<<<<<<<<<<<<<<<<<<<<<<<<< draw scene >>>>>>>>>>>>>>>>>>>>>>
void Draw()
{
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT); // clear the screen
	
	Mesh mymesh;  // Mesh object for tapered cylinder
	Mesh mymesh1; // Mesh object for half sphere
	mymesh.readmesh("taperedCylinderMesh.txt"); // read mesh data
	mymesh1.readmesh("halfSphereMesh.txt");     // read mesh data
	
	glPushMatrix();       // push matrix
	glRotatef(-90,1,0,0); // rotate 90 degrees to make the object stand upright
	mymesh.draw();        // draw tapered cylinder
	mymesh1.draw();       // draw half sphere
	glPopMatrix();        // push matrix
	glutSwapBuffers();    // flush the display to the screen
}


// Respond to arrow keys
void SpecialKeys(int key, int x, int y)
{
	if (moveLRTB)                  // if moveLRTB is true
	{
	if(key == GLUT_KEY_UP)         // up arrow    -> move up
  	myCamera.pitch(1.0); 
 	else if(key == GLUT_KEY_DOWN)  // down arrow  -> move down
  	myCamera.pitch(-1.0);
	else if(key == GLUT_KEY_LEFT)  // left arrow  -> move left
  	myCamera.yaw(-1.0);
	else if(key == GLUT_KEY_RIGHT) // right arrow -> move right
  	myCamera.yaw(1.0);
	}

	if (rotate)                    // if rotate is true
	{
	if(key == GLUT_KEY_UP)         // up arrow    -> rotate upward
  	myCamera.slide(0, -0.1, 0);
 	else if(key == GLUT_KEY_DOWN)  // down arrow  -> rotate downward
  	myCamera.slide(0, 0.1, 0);
	else if(key == GLUT_KEY_LEFT)  // left arrow  -> rotate leftward
  	myCamera.slide(0.1, 0, 0); 
	else if(key == GLUT_KEY_RIGHT) // right arrow -> rotate rightward
  	myCamera.slide(-0.1, 0, 0);
	}
	glutPostRedisplay();           // flush the display to the screen
}

// keyboard functions
void myKeyboard(unsigned char theKey, int x, int y)
{
	switch(theKey)
	{
		case 'f': case 'F':		// move camera forward, zoom in
   		myCamera.slide(0, 0, -0.2);
			break;
		case 'b': case 'B':		// move camera backward, zoom out
    	myCamera.slide(0, 0, 0.1);
			break;
		case 'r': case 'R':   // set rotate mode
			moveLRTB = false; 
			rotate = true;
			break;
		case 'm': case 'M':   // set move mode
			rotate = false;
			moveLRTB = true;
			break;
		case 'i':	case 'I':   // set everything to the initial setup
    	myCamera.setShape(45, 1, 1, 100);
			myCamera.set(startupCameraX,startupCameraY,startupCameraZ,0, 0, 0,0, 1, 0);
			break;
		case 'a':	case 'A':   // animation mode
			cout << "Enter point of rotation : " ; // user input for point of rotation
			cin >> x_r >> y_r >> z_r;              // set point of rataion
			animation = true;                      // start animation
			glutTimerFunc(50, TimerFunction, 1);   // timer function
			break;
		case 's':	case 'S':   // stop animation
			animation = false;
			break;
		case 'q':             // quit program
			exit (0);
		default:
			if (theKey == 27)   // ASCII for escape character
				exit(0);
	}

	glutPostRedisplay();    // flush the display to the screen
}

// timer function for animation
void TimerFunction(int value)
{
	if (animation)                             // if animation is true
	{
		xRot -= 15.0;                            // rotation angle increment
		xRot = (GLfloat)((const int)xRot % 360); // reset rot angle angle after 360 deg
		myCamera.set(x_r,y_r,z_r,0,0,0,0,1,0);   // set camera to user defined point
		myCamera.roll(xRot);                     // rotate
		glutPostRedisplay();                     // redisplay the screen
		glutTimerFunc(70, TimerFunction, 1);     // recursive call to timer function
	}
}

