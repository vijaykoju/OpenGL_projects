/*
---------------------------------------------------------------------------
FILENAME   : heat.cpp
PROGRAMMER : Vijay Koju

---------------------------------------------------------------------------
This programs a heat map of 50 sets of data. For each set of data, the
color value is determined from the first data point. Along with the heat
map, it also plots the data for a specific data set when clicked on the
respective region on the heat map. A detailed statistic of the specified
data set is also shown.
---------------------------------------------------------------------------
Mouse interaction:
left-click on any region of the heat map --> Displays the plot and statis-
tics of the respective data set.
---------------------------------------------------------------------------
Keyboard interaction:
q   --> quit the program
exc --> quit the program
w   --> change the color tone of the heat map to warm
c   --> change the color tone of the heat map to cool
g   --> change the color tone of the heat map to grey (default)
f   --> flip the color tome of the heat map
---------------------------------------------------------------------------
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

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <iomanip>
#include <cstring>
#include <sstream>
#include "Plot.h"
using namespace std;


// function prototypes
void setWindow(double left, double right, double bottom, double top);
void setViewport(double left, double right, double bottom, double top);
void compute();
void MyInit();
void Draw();
void myKeyboard(unsigned char theKey, int x, int y);
double findMin(double list[], int list_size);
double findMax(double list[], int list_size);
void drawBitmapText(const char* str,float x,float y,float z,int n);
void myMouse(int button, int state, int x, int y);
//char* doubleToString(double num);

//-----------------------------------------------------------------------//
// global variables declarations
#define Width 600
#define Height 600
int num_Rows, num_Columns, data_size;
Plot plot_object[50];
double first_val[50], max_val[50], min_val[50], range_val[50], mean_val[50], stdDev_val[50];
double data_set[20];
double max_first_val, min_first_val;
bool warm=false;
bool cool=false;
bool grey=false;
bool flip=false;
bool flip1=false;
char color;
int o_num=0;
//typedef GLfloat points[2];
//points point[800*400];

//-----------------------------------------------------------------------//
// Function name  : main()
// Preconditions  : There are none.
// Postconditions : This function initiates glut (and thus OpenGL).
// It requires appropriate callbacks to diplay the heat map, plot of the
// data, and its statistics in a window created by glut.
int main(int argc, char** argv)
{
	// initialize OpenGL Utility Toolkit
	glutInit(&argc, argv);
	// sets the display mode
	glutInitDisplayMode(GLUT_DOUBLE|GLUT_RGB);
	// requests a screen window Width pixels wide by Height pixels high
	glutInitWindowSize(Width,Height);
	// specify the window position
	glutInitWindowPosition(100,100);
	// open and display the window putting "Heat Map" on the title bar
	glutCreateWindow("Heat Map");
	// set up the initial state of some of OpenGL's variables
	MyInit();
	// register the Draw() function as the function to activate when a redraw
	// event occurs
	glutDisplayFunc(Draw);
	// register myKeyboard() function as the function to activate keyboard
	// interactions
	glutKeyboardFunc(myKeyboard);
	// register myMouse() function as the function to activate mouse
	// interactions
	glutMouseFunc(myMouse);
	// enter an unending loop waiting for events to occur
	glutMainLoop();
	return 0;
}

//-------------------------------setWindow-------------------------------//
void setWindow(double left, double right, double bottom, double top)
{
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluOrtho2D(left,right,bottom,top);
}

//------------------------------setViewport------------------------------//
void setViewport(double left, double right, double bottom, double top)
{
	glViewport(left,bottom,right-left,top-bottom);
}

//-----------------------------------------------------------------------//
// Function name   : MyInit()
// Preconditions   : OpenGL and GLUT must be initialized
// Postconditions  : This function initialized the backgroud color and sets
// the world window
void MyInit()
{
	// initialize background color to white
	glClearColor(1,1,1,1);
	// all the major computations such as calculating mean, std dev, range,
  //  min, max are done in compute()
	compute();
	// cordinate transformations
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	setWindow(0,Width,0,Height);
}

//-----------------------------------------------------------------------//
// Function name   : Draw()
// This is a callback function to display the heat map, plot and statistics
// on the display window
void Draw()
{
	double mini, maxi;
	double c_val;
	char val_min[10],val_max[10],val_mean[10],val_range[10],val_stdDev[10];

	// viewport 1 ------------- left side ------------------
	glViewport(10,10,300,600);
	glClear(GL_COLOR_BUFFER_BIT);

	// Generate heat map with 50 boxes, each box representing one object.
	// The color value for the heat map corresponds to the first value of
	// each object. 
	for (int i=0; i<10; ++i)
	{
		for (int j=0; j<5; ++j)
		{
			// normalized color value
			c_val = first_val[49-(j+5*i)]/range_val[49-(j+5*i)];
			//if (flip==true)
			//	c_val = 1-c_val;
			if (warm==true) 
				glColor3f(c_val,0,0); // warm tone
			else if (cool==true)
				glColor3f(0,0,c_val); // cool tone
			else
				glColor3f(c_val,c_val,c_val); // grey tone
			glRecti(116*j,58*i,116*(j+1),58*(i+1)); // generate separate rectan-
																							// gle for each object
		}
	}

	// viewport 2 ----------------top- right side ----------------
	glViewport(310,58*8,290,200);
	drawBitmapText("Heatmap Program",0,150,0,1);	
	drawBitmapText("low",0,0,0,1);	
	drawBitmapText("high",500,0,0,1);	

	// normalized mimimum and maximum value for the color scale
	mini = min_first_val/(max_first_val-min_first_val);
	maxi = max_first_val/(max_first_val-min_first_val);
	//if (flip1==true)
	//{
	//	mini=1-mini;
	//	maxi=1-maxi;
	//}

	// section for generating the color bar with smooth blend effect
	glShadeModel(GL_SMOOTH);
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);
	// make a rectangular box for the color bar
	glBegin(GL_POLYGON);
	if (warm==true)
		glColor4f(mini,0,0,1); // warm minimum
	else if (cool==true)
		glColor4f(0,0,mini,1); // cool minimum
	else
		glColor4f(mini,mini,mini,1); // grey
	glVertex3f(0,50,0); // coordinate of lower-left corner
	glVertex3f(0,50+80,0); // coordinate of upper-left corner

	if (warm==true)
		glColor4f(maxi,0,0,1); // warm maximum
	else if (cool==true) 
		glColor4f(0,0,maxi,1); // cool maximum
	else
		glColor4f(maxi,maxi,maxi,1); // grey
	glVertex3f(580,50+80,0); // coordinate of upper-right corner
	glVertex3f(580,50,0); // coordinate of lower-right corner
	glEnd();

	glDisable(GL_BLEND); // end blending

	warm=false;
	cool=false;
	flip=false;
	flip1=false;
	
	// viewport 3 ----------------lower-right side----------------------
	glViewport(310,0,290,400);

	// get the data set of the requested object and offset it appropriately
	// so that the minimum value lies at the bottom of the plot
	for (int i=0; i<20; ++i)
		data_set[i] = plot_object[49-o_num].GetData(20)[i]-min_val[49-o_num];

	// convert numerical values to characters to display on the window
	sprintf(val_min,"%2.5f",min_val[49-o_num]); // min
	sprintf(val_max,"%2.5f",max_val[49-o_num]); // max
	sprintf(val_mean,"%2.5f",mean_val[49-o_num]); // mean
	sprintf(val_range,"%2.5f",range_val[49-o_num]); // range
	sprintf(val_stdDev,"%2.5f",stdDev_val[49-o_num]); // std-dev

	glColor3f(0,0,0);
	glBegin(GL_LINE_STRIP);
	glVertex3f(0,580,0); // top point of y-axix
	glVertex3f(0,350,0); // origin of the axis
	glVertex3f(580,350,0); // right point of x-axis
	glEnd();

	// vertex of 20 data points associated with the requested object
	glBegin(GL_LINE_STRIP);
	for (int i=0; i<20; ++i)
		glVertex3f(i*(580.0/19),350+data_set[i]*((580-350)/(max_val[49-o_num]-min_val[49-o_num])),0);
	glEnd();

	// text display section
	drawBitmapText(val_max,0,583,0,2);
	drawBitmapText(val_min,0,335,0,2);
	drawBitmapText("mean     = ",0,250,0,2);	
	drawBitmapText(val_mean,150,250,0,2);
	drawBitmapText("min      = ",0,250-24,0,2);	
	drawBitmapText(val_min,150,250-24,0,2);
	drawBitmapText("max      = ",0,250-2*24,0,2);	
	drawBitmapText(val_max,150,250-2*24,0,2);
	drawBitmapText("range    = ",0,250-3*24,0,2);	
	drawBitmapText(val_range,150,250-3*24,0,2);
	drawBitmapText("std dev  = ",0,250-4*24,0,2);	
	drawBitmapText(val_stdDev,150,250-4*24,0,2);
	drawBitmapText("KEY OPTIONS",0,250-6*24,0,2);	
	drawBitmapText("w = warm",0,250-7*24,0,2);	
	drawBitmapText("g = grey",0,250-8*24,0,2);	
	drawBitmapText("c = cool",0,250-9*24,0,2);	
	drawBitmapText("f = flip color scale",0,250-10*24,0,2);	

	glutSwapBuffers();	
}

//-----------------------------------------------------------------------//
// Function name  : compute()
// This function does all the major computations. The data file is opened
// only once to save the data into an array. 50 plot_objects are created
// and each of them is assigned a set of 20 data points.
void compute()
{
	// variables declarations
	ifstream dataFile;
	int x_data, index1,n;
	double y_data;
	double* array_of_y_data;
	dataFile.open("heat.txt"); // open txt file for reading
	// read data from the first line of dataFile
	dataFile >> num_Rows >> num_Columns;
	//cout << num_Rows << " " << num_Columns << endl;
	n = num_Rows*num_Columns;
	// dynamic allocation of array to store the first data of each dataset
	// read data from the second line of dataFile
	dataFile >> data_size;
	//cout << data_size << endl;
	// dynamic allocation of array to store all data
	array_of_y_data = new double[n*data_size];
	
	index1 = -1; // index counters for arrays
	// read and store data the remaing data from dataset
	while (dataFile)
	{
		dataFile >> x_data >> y_data; // read
		array_of_y_data[++index1] = y_data; // store it in array
		if (!dataFile) break; // prevents reading the last line twice
	}
	dataFile.close(); // close dataFile
	
	// set data for each Plot object and compute statistics
	for (int i=0; i<n; ++i)
	{
		plot_object[i].SetObjectNumber(i); // object number
		plot_object[i].SetData(array_of_y_data); // set data
		first_val[i] = plot_object[i].GetFirstVal(); // first value
		max_val[i] = plot_object[i].GetMax(); // maximum value
		min_val[i] = plot_object[i].GetMin(); // minimum value
		range_val[i] = plot_object[i].GetRange(); // range
		mean_val[i] = plot_object[i].GetMean(); // mean
		stdDev_val[i] = plot_object[i].GetStdDev(); // standard deviation
	}
	
	max_first_val = findMax(first_val,50);
	min_first_val = findMin(first_val,50);

	delete[] array_of_y_data;
}

//-----------------------------------------------------------------------//
// Function name  : myKeyboard()
// Preconditions  : OpenGL and GLUT must be initialized and this function
// must be registered as the keyboard callback.
// Postconditions : This function hands the keyboard callback.
// Keyboard function:
// q      --> quit
// esc    --> quit
// w      --> warm color tone
// c      --> cool color tone
// g      --> grey color tone
// f      --> flip color value
void myKeyboard(unsigned char theKey, int x, int y)
{
	switch (theKey)
	{
		case 'q':
		case 27:
			exit(0);
		case 'w':
			color = 'w';
		case 'c':
			color = 'c';
		case 'g':
			color = 'g';
		case 'f':
			color = 'f';
		default:
			color = 'g';
	}
	glutPostRedisplay();
}

// ------------------------------------------------------------------------
// Function name   : findMin()
// Input parameter : double list[] --> array
//                 : int list_size --> arry size
// Output          : minimum of lis[]
double findMin(double list[],int list_size)
{
	double min = list[0];
	for (int i=1; i<list_size; ++i)
	{
		if (list[i]<min)
			min = list[i];
	}
	return min;
}

// ------------------------------------------------------------------------
// Function name   : findMax()
// Input parameter : double list[] --> array
//                 : int list_size --> arry size
// Output          : maximum of lis[]
double findMax(double list[],int list_size)
{
	double max = list[0];
	for (int i=1; i<list_size; ++i)
	{
		if (list[i]>max)
			max = list[i];
	}
	return max;
}

// ------------------------------------------------------------------------
// Function name : drawBitmapText()
// This function is a wrapper function to glutBitmapCharacter() function.
// Input parameter : const char* str --> pointer to a character string
//                 : float x, y, z   --> display postion
//                 : int n           --> font type (1 or 2)
void drawBitmapText(const char* str,float x,float y,float z,int n)
{  
	int len = strlen(str); // string length
	glRasterPos3f(x, y,z); // display position
	if (n==1)
	{
		for (int i=0; i<len; ++i)
			glutBitmapCharacter(GLUT_BITMAP_HELVETICA_18, str[i]); // font 1
	}
	else if (n==2)
	{
		for (int i=0; i<len; ++i)
			glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12, str[i]); // font 2
	}
}

// ------------------------------------------------------------------------
// Function name  : myMouse()
// Preconditions  : OpenGL and GLUT must be initialized and this function
// must be registered as the mouse callback.
// Postconditions : This function handles the mouse callback.
// Left-click on any box on the heat map generates the associated plot of
// the data linked to that box, and also display its statistics.
void myMouse(int button, int state, int x, int y)
{
	if (button==GLUT_LEFT_BUTTON && state==GLUT_DOWN)
	{
		for (int i=0; i<10; ++i)
		{
			for (int j=0; j<5; ++j)
			{
				// check on which box the click position lies 
				if (x>10+(58*j)&&x<=10+(58*(j+1))&&y>10+(58*i)&&y<=10+(58*(i+1)))
					o_num=5*i+j; // update the o_num to respective object number
			}
		}
	}
	glutPostRedisplay();
}
