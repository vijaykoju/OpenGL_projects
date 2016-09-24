/*
FILENAME   : fern.cpp
PROGRAMMER : Vijay Koju

This programs displays the pattern of two different ferns using a simple
mathematical structure, along with probability.
When the program is run, it displays a pleasing fern.
Some mouse interactivities are as follow.
left-click  -> Display the second fern
right-click -> Change the color of fern
Keyboard interactions:
q   -> quit the program
esc -> quit the program
*/
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

#include <cstdlib>
#include <time.h>

// set the width and height of the screen window
#define Width 500
#define Height 500

// set the width and height of the world window
#define worldWidth 10
#define worldHeight 11

float R = worldWidth/worldHeight;
bool changeColor=false; // flag to change the color of fern
bool fern2=false;       // flag to change the fern

#define NUM_OF_POINTS 100000

// define a point data type
typedef GLfloat point2[2];
point2 points[NUM_OF_POINTS]; // this array holds all 100000 points to draw
point2 points2[NUM_OF_POINTS];// the ferns


// function prototypes
void ComputeFern();
void MyInit();
void setWindow(double left, double right, double bottom, double top);
void setViewport(double left, double right, double bottom, double top);
void DrawFern();
void myKeyboard(unsigned char theKey, int x, int y);
void myMouse(int button, int state, int x, int y);


///////////////////////////////////////////////////////////////////////////
// Function name  : main()
// Preconditions  : There are none.
// Postconditions : This function initiates glut (and thus OpenGL).
// It requires appropriate callbacks to diplay a fern in the window
// created by glut.
///////////////////////////////////////////////////////////////////////////
int main(int argc, char** argv)
{
	// initialize the OpenGL Utility Toolkit
	glutInit(&argc, argv);

	// sets the display mode
	glutInitDisplayMode(GLUT_DOUBLE|GLUT_RGB);

	// request a screen window Width pixels wide by Height pixels high
	glutInitWindowSize(Width,Height);

	// specify the window position
	glutInitWindowPosition(100,100);

	// open and display the window putting "Fern" on the title bar
	glutCreateWindow("Fern");

	// set up the initial state of some of OpenGL's variables
	MyInit();
	
	// register the DrawFern() function as the function to activate when
	// a redraw event occurs
	glutDisplayFunc(DrawFern);

	// register myKeyboard() function as the function to activate keyboard
	// interactions
	glutKeyboardFunc(myKeyboard);

	// register myMouse() function as the function to activate mouse
	// interactions
	glutMouseFunc(myMouse);

	// enter an unending loop waiting from events to occur
	glutMainLoop();
	return 0;
}

//----------------------------- setWindow ---------------------------------
void setWindow(double left, double right, double bottom, double top)
{
	//glMatrixMode(GL_PROJECTION);
	//glLoadIdentity();
	gluOrtho2D(left,right,bottom,top);
}

//--------------------------- setViewport ---------------------------------
void setViewport(double left, double right, double bottom, double top)
{
	glViewport(left,bottom,right-left,top-bottom);
}

///////////////////////////////////////////////////////////////////////////
// Function name  : MyInit()
// Preconditions  : OpenGL and GLUT must be initialized/
// Postconditions : This function intializes the backgroud color and sets 
// the world to window transformations.
///////////////////////////////////////////////////////////////////////////
void MyInit()
{
	// initialize backgound color
	glClearColor(0.0,0.0,0.0,1.0);
	// set point size
	glPointSize(1.0);
	
	// computes all the points needed to draw a fern
	ComputeFern();
	
	// coordinate transformations
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	setWindow(-5,5,0,11);
	setViewport(0,Width,0,Height);
}

///////////////////////////////////////////////////////////////////////////
// Function name  : ComputeFern()
// This function does all the computations required to draw a fern. 
///////////////////////////////////////////////////////////////////////////
void ComputeFern()
{
	//srand(time(0));

	// define variables
	GLint k;
	GLfloat random_num1,random_num2;
	GLfloat a, b, c, d, e, f, a2, b2, c2, d2, e2, f2;

	// define initial points for the first fern
	points[0][0] = ((double) rand()/RAND_MAX);
	points[0][1] = ((double) rand()/RAND_MAX);

	// define initial points for the second fern
	points2[0][0] = ((double) rand()/RAND_MAX);
	points2[0][1] = ((double) rand()/RAND_MAX);

	// compute the points using the following form
	// x_(n+1) = a*x_(n) + b*y_(n) + e
	// y_(n+1) = c*x_(n) + d*y_(n) + f
	for (k=0;k<NUM_OF_POINTS-1;k++)
	{
		// random numbers between 0 and 1.
		random_num1 = ((double) rand()/RAND_MAX); // for the first fern
		random_num2 = ((double) rand()/RAND_MAX); // for the second fern

		// a, b, c, d, e, and f for the first fern according to the 
    // difined probability
		if (random_num1>=0 && random_num1<=0.1)
		{
			a=0.0; b=0.0; c=0.0; d=0.16; e=0.0; f=0.0;
		}
		else if (random_num1>0.1 && random_num1<=0.18)
		{
			a=0.2; b=-0.26; c=0.23; d=0.22; e=0.0; f=1.6;
		}
		else if (random_num1>0.18 && random_num1<=0.26)
		{
			a=-0.15; b=0.28; c=0.26; d=0.24; e=0.0; f=0.44;
		}
		else if (random_num1>0.26 && random_num1<=1.0)
		{
			a=0.75; b=0.04; c=-0.04; d=0.85; e=0.0; f=1.6;
		}

		// a, b, c, d, e, and f for the second fern according to the 
    // defined probability
		if (random_num2>=0 && random_num2<=0.01)
		{
			a2=0.0; b2=0.0; c2=0.0; d2=0.16; e2=0.0; f2=0.0;
		}
		else if (random_num2>0.01 && random_num2<=0.08)
		{
			a2=0.2; b2=-0.26; c2=0.23; d2=0.22; e2=0.0; f2=1.6;
		}
		else if (random_num2>0.08 && random_num2<=0.15)
		{
			a2=-0.15; b2=0.28; c2=0.26; d2=0.24; e2=0.0; f2=0.44;
		}
		else if (random_num2>0.15 && random_num2<=1.0)
		{
			a2=0.85; b2=0.04; c2=-0.04; d2=0.85; e2=0.0; f2=1.6;
		}
	
		// compute points for the first fern
		points[k+1][0] = a*points[k][0] + b*points[k][1] + e;
		points[k+1][1] = c*points[k][0] + d*points[k][1] + f;

		// compute points for the second fern
		points2[k+1][0] = a2*points2[k][0] + b2*points2[k][1] + e2;
		points2[k+1][1] = c2*points2[k][0] + d2*points2[k][1] + f2;
	}
}

///////////////////////////////////////////////////////////////////////////
// Function name  : DrawFern()
// Draws the fern. This function should be registered as a callback for a
// redraw evnent.
// All the points needed to draw the fern should already be computed.
///////////////////////////////////////////////////////////////////////////
void DrawFern()
{
	glClear(GL_COLOR_BUFFER_BIT);
	// if changeColor is true, set the color to a different shade of green
	// else set it to pure green
	if (changeColor)
		glColor3f(0.2,1.0,0.4);
	else
		glColor3f(0.0,1.0,0.0);

	// if fern2 is true, display the second fern
	// else display the first fern
	if (fern2)
	{
	// begin plotting the points
	glBegin(GL_POINTS);
		for (GLint i=0;i<NUM_OF_POINTS;i++)
			glVertex2fv(points2[i]);
	glEnd();
	}
	else
	{
	glBegin(GL_POINTS);
		for (GLint i=0;i<NUM_OF_POINTS;i++)
			glVertex2fv(points[i]);
	glEnd();
	}

	glutSwapBuffers();
}

///////////////////////////////////////////////////////////////////////////
// Function name  : myKeyboard()
// Preconditions  : OpenGL and GLUT must be initialized and this function
// must be registered as the keyboard callback.
// Postconditions : This function hands the keyboard callback. On pressing
// 'q' or 'esc' key on the keyboard, it quits the program.
///////////////////////////////////////////////////////////////////////////
void myKeyboard(unsigned char theKey, int x, int y)
{
	// define the keyboard functionality
	// quit the program when 'q' or 'esc' key is pressed
	if (theKey=='q' or theKey==27)
		exit(0);
}

///////////////////////////////////////////////////////////////////////////
// Function name  : myMouse()
// Preconditions  : OpenGL and GLUT must be initialized and this function
// must be registered as the mouse callback.
// Postconditions : This function handles the mouse callback.
// The left-button switches between the two ferns, and the right-button
// changes the fern color to a different shade of green.
///////////////////////////////////////////////////////////////////////////
void myMouse(int button, int state, int x, int y)
{
	// define the right-button functionalities
	if (button==GLUT_RIGHT_BUTTON && state==GLUT_DOWN && changeColor==false)
		changeColor=true;
	else if (button==GLUT_RIGHT_BUTTON && state==GLUT_DOWN && changeColor==true)
		changeColor=false;

	// define the left-button functioanlities
	else if (button==GLUT_LEFT_BUTTON && state==GLUT_DOWN && fern2==false)
		fern2=true;
	else if (button==GLUT_LEFT_BUTTON && state==GLUT_DOWN && fern2==true)
		fern2=false;
	glutPostRedisplay();
}
