/*
-------------------------------------------------------------------------------
FILENAME   : FinalProject.cpp
PROGRAMMER : Vijay Koju
CLASS      : CSCI 7300 (Scientific Visualization and Databases)
DUE DATE   : 12/06/2013
INSTRUCTOR : Dr. Li
-------------------------------------------------------------------------------
Required files to run this program :
1) FinalPro_GalSim.cpp
2) intialConditionData.dat
3) CoolWarmFloat257.dat
4) Jet257.dat
5) makefile

To run this program (if all the files listed above are present):
$ make
$ ./galSim

For a brief description of the problem and to view the keyboard and mouse options
available for user interaction during simulation, please run the program using
commands listed above in the teminal and press 'h'.
*/
#define LINUX 
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
#include <iostream>
#include <fstream>
#include <cmath>
#include <cassert>
#include <vector>
#include <ctime>
#include <cstring>
#include <algorithm>
using namespace std;

#define WIDTH  1000 // window width
#define HEIGHT 700  // window height
#define GLX    0    // galaxy mode
#define KE     1    // energy mode
GLint display=GLX;  // set default to galaxy mode

// initial camera position
GLfloat cX  = 0;   // x-position
GLfloat cY  = 0;   // y-position
GLfloat cZ  = 77;  // z-position

// variables for zooming and moving in l,r,u,d direction
GLfloat s_f = 1;   // scale factor for zooming in/out
GLfloat xT  = 0;   // translate factor in x-direction
GLfloat yT  = 0;   // translate factor in y-direction

// variables for rotation
int spin_x = 0, spin_y = 0;
int old_x  = 0, old_y  = 0;

// boolean variables
bool animate     = false; // start/stop animation     - default stop
bool o           = true;  // boundary (box) on/off    - default on
bool changeColor = true;  // change color cool/warm   - default warm
bool fullScreen  = true;  // fullscreen mode on/off   - default off
bool t           = false; // trace galaxy path on/off - default off

// numbers for this simulation
int   N      = 32768; // total number of stars in the simulation
int   M      = 2;     // number of galaxies
float tt     = 0;     // initial time
float ep     = 45;    // thresshold distance of interaction (units arbitrary)
float p_fact = 7;     // scale factor for distance
float v_fact = 13;    // scale factor for velocity
float dt     = 0.1;   // time step-size to march forward
float th1    = 20.0;   // tilt of galaxy 1
float th2    = 30.0;   // tilt of galaxy 2

// arrays for storing data
float S_m[32768];     // star-mass array
float S_r[32768][3];  // star-postion array        - <x,y,z> components
float S_v[32768][3];  // star-velocity array       - <x,y,z> components
float S_a[32768][3];  // star-acceleration array   - <x,y,z> components
float S_jk[32768][3]; // star-jerk array           - <x,y,z> components
float G_m[2];         // galaxy-mass array
float G_r[2][3];      // galaxy-position array     - <x,y,z> components
float G_v[2][3];      // galaxy-velocity array     - <x,y,z> components
float G_a[2][3];      // galaxy-acceleration array - <x,y,z> components
float G_jk[2][3];     // galaxy-jerk array         - <x,y,z> components
char  val_min[10];    // minimum velocity value to display below color bar
char  val_max[10];    // maximum velocity value to display above color bar
char  str_v_fact[10]; // velocity scale factor to display. changes when 'r' is pressed

// the size of these vector lists increase at every time step as they get
// appropriate updates 
// C.O.M. = Center of Mass
vector <float> tm1;   	// time-array
vector <float> KE_G1; 	// kinetic energy of galaxy 1 core
vector <float> KE_G2; 	// kenetic energy of galaxy 2 core
vector <float> PE_G1G2; // potential energy between the two galaxy cores 
vector <float> G1_x;  	// x-position of galaxy 1 core
vector <float> G1_y;  	// y-position of galaxy 1 core
vector <float> G1_z;  	// z-position of galaxy 1 core
vector <float> G2_x;  	// x-position of galaxy 2 core
vector <float> G2_y;  	// y-position of galaxy 2 core
vector <float> G2_z;  	// z-position of galaxy 2 core
vector <float> COM[3];	// coordinates of C.O.M. of the two galaxy cores

// struct type of color-info
typedef struct
{
	float fval;         // color-function value
	float rgb[3];       // color-values - <r,g,b>
}color;
color cVal[257];      // struct for color

// user difined function
void MyInit();          // initialize required opengl variables
void drawBoundary();    // draw 3-d cube
void drawAxis();        // draw 3-d axis at the bottom-left corner
void initCondition();   // initialize galaxy simulation
void evolveGalaxies();  // evolve galaxies - C.O.M. of the galaxies
void evolveStars();     // evolve all the stars
void drawGalaxy();      // draw galaxies
void drawEnergy();      // draw kinetic energy plot of the galaxies
void draw();            // wrapper to drawGalaxies() and drawEnergy()
void drawColorBar();    // draw colorbar
void processMenuEvents(int option);                  // process menu
void myKeyboard(unsigned char theKey, int x, int y); // keyboard functions
void freeMemory();                                   // deallocate memory
void specialKeys(int key, int x, int y);             // spcial keys
void mouse(int button, int state, int x, int y);     // mouse function
void motion(int x,int y);                            // handle mouse motion
void TimerFunc(int value);                           // timer function (animate)
void loadColor();                                    // laod color info
void drawBitmapText(const char* str,float x,float y,float z,int n,char c); //text wrapper
void displayInfo();      // display info and all the options

// main program
int main(int argc, char **argv)
{
	srand(time(0));
	glutInit(&argc, argv); // initialize glut variables
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH); // set displaymode
	glutInitWindowSize(WIDTH,HEIGHT); // window size
	glutInitWindowPosition(0,0);      // window position
	glutCreateWindow("Galaxy Collision Simulation"); // window title
	glutDisplayFunc(draw); 				// draw galaxies
	initCondition();       				// initialize galaxies and stars
	glutCreateMenu(processMenuEvents);               // process menu
	glutAddMenuEntry("Galaxy Collision",GLX);        // add item to menu
	glutAddMenuEntry("Energy",KE);           // add item to menu
	glutAttachMenu(GLUT_RIGHT_BUTTON);               // attach menu
	glutSpecialFunc(specialKeys);  // callback to special keys functions
	glutKeyboardFunc(myKeyboard);  // callback to keyboard functions
	glutMouseFunc(mouse);          // callback to mouse functions
	glutMotionFunc(motion);        // handle mouse motion
	glutTimerFunc(50,TimerFunc,1); // callback to timer function
	glViewport(0,0,WIDTH,HEIGHT);  // set the whole widow as viewport
	MyInit();                      // set required matrix modes
	glutMainLoop();                // enter and infinte loop and wait for events
	return 0;
}

// set the required projection and modelview matrix modes
void MyInit()
{
	glMatrixMode(GL_PROJECTION);    // load projection matrix
	glLoadIdentity();               // load identity matrix
	gluPerspective(60,(float)WIDTH/(float)HEIGHT,0.1,10000); // set perspective view
	
	glMatrixMode(GL_MODELVIEW);     // load modelview matrix
	glLoadIdentity();               // load identity matrix
	gluLookAt(cX,cY,p_fact*cZ,0,0,0,0.0,1.0,0.0); // camera and lookat positions

	glClearColor(0.0f,0.0f,0.0f,1.0f);            // set background color to black

	glEnable(GL_DEPTH_TEST);                      // enable depth
}

// draw 3-d cube as a reference boundary
void drawBoundary()
{
		glPushMatrix();                          // save current tranformation matrix
		glScaled(p_fact*50,p_fact*50,p_fact*50); // scaling
			glutWireCube(1.0);                     // draw a wireframe cube
		glPopMatrix();	                         
}

// set initial conditions (position,velocity,acceleration, and jerk) of stars.
// load the initial condition data from "initialConditionData.dat" file.
// this data is acquired from "http://bima.astro.umd.edu/nemo/archive/#dubinski".
// for the initial condition of galaxies take the average values of position,
// and velocities of their respective stars.
// initial acceleration and jerk are just set to zero
void initCondition()
{
	ifstream myIn;                         // file stream
	myIn.open("initialConditionData.dat"); // open data file
	assert(myIn);                          // check for error while opening
	float avg_Gr[2][3], avg_Gv[2][3];// average position and velocities for galaxies

	for (int i=0; i<N; ++i)
	{
		myIn>>S_m[i];            // read mass 
		for (int k=0; k<3; ++k)
			myIn >> S_r[i][k];     // read position <x,y,z> components
		for (int k=0; k<3; ++k)
			myIn >> S_v[i][k];     // read velocity <x,y,z> components
	}

	for (int i=0; i<N; ++i)
	{
		for (int k=0; k<3; ++k)
		{
			//if (i<N/2)
			//	S_r[i][k] = p_fact*S_r[i][k]*cos(th1*3.14159/180);
			//else 
			//	S_r[i][k] = p_fact*S_r[i][k]*cos(th2*3.14159/180);
			S_r[i][k] = p_fact*S_r[i][k]; // scale position data
			S_v[i][k] = v_fact*S_v[i][k]; // scale velocity data
			S_a[i][k] = S_jk[i][k] = 0;   // set acceleration and jerk to zero
		}
	}

	for (int i=0; i<M; ++i)
	{
		for (int k=0; k<3; ++k)
			avg_Gr[i][k]=avg_Gv[i][k]=0;  // initialize average position and velocity 
	}                                 // of galaxy cores to zero

	for (int i=0; i<N/2; ++i)
	{
		for (int k=0; k<3; ++k)
		{
			// add-up postion data for averaging
			avg_Gr[0][k]  += S_r[i][k];  avg_Gr[1][k]  += S_r[i+N/2][k];
			// add-up velocity data for averaging
			avg_Gv[0][k]  += S_v[i][k];  avg_Gv[1][k]  += S_v[i+N/2][k];
		}
	}
	for (int i=0; i<M; ++i)
	{
		for (int k=0; k<3; ++k)
		{
			avg_Gr[i][k] /= N/2; // positions of galaxy cores (averaged from respective stars)
			avg_Gv[i][k] /= N/2; // velociteis of galxy cores (averaged)
			G_r[i][k] = avg_Gr[i][k]; // set initial positions for galaxies
			G_v[i][k] = avg_Gv[i][k]; // set initial velocities for galaxies
			G_a[i][k] = G_jk[i][k] = 0; // set initial accelration and jerk to zero
		}
	}
	G1_x.push_back(G_r[0][0]); // store initial x-position of galaxy 1 core
	G1_y.push_back(G_r[0][1]); // store initial y-position of galaxy 1 core
	G1_z.push_back(G_r[0][2]); // store initial z-position of galaxy 1 core
	G2_x.push_back(G_r[1][0]); // store initial x-position of galaxy 2 core
	G2_y.push_back(G_r[1][1]); // store initial y-position of galaxy 2 core
	G2_z.push_back(G_r[1][2]); // store initial z-position of galaxy 2 core

	G_m[0]=4500;                  // set mass of galaxy 1 core
	G_m[1]=6400;                  // set mass of galaxy 2 core
	float t_mass = G_m[0]+G_m[1]; // total mass of two galaxy cores

	float com_x = (G_m[0]*G_r[0][0]+G_m[1]*G_r[1][0])/t_mass; // x component of C.O.M
	float com_y = (G_m[0]*G_r[0][1]+G_m[1]*G_r[1][1])/t_mass; // y component of C.O.M
	float com_z = (G_m[0]*G_r[0][2]+G_m[1]*G_r[1][2])/t_mass; // z component of C.O.M
	// store x,y, and z components
	COM[0].push_back(com_x);   
	COM[1].push_back(com_y);
	COM[2].push_back(com_z);

	tm1.push_back(tt);         // store initial time

	// compute initial kinetic energies of galaxy 1 and 2 cores
	float k1=0, k2=0;          // initialize kinetic energy values
	k1 = 0.5*G_m[0]*(G_v[0][0]*G_v[0][0]+G_v[0][1]*G_v[0][1]+G_v[0][2]*G_v[0][2]); // K.E 1
	k2 = 0.5*G_m[1]*(G_v[1][0]*G_v[1][0]+G_v[1][1]*G_v[1][1]+G_v[1][2]*G_v[1][2]); // K.E 2
	KE_G1.push_back(k1);       // store k.e of galaxy 1
	KE_G2.push_back(k2);       // store k.e of galaxy 2

	// compute initial potential energy between the two galaxy cores
	float G_rji[3];
	for (int k=0; k<3; ++k)
		G_rji[k] = G_r[1][k] - G_r[0][k];  // distance between the two cores
	float G_r2 = 0;
	for (int k=0; k<3; ++k)
		G_r2 += G_rji[k]*G_rji[k];         // distance squared
	float p1= -G_m[0]*G_m[1]/sqrt(G_r2); // P.E
	PE_G1G2.push_back(p1);               // store P.E
}

// evolve galaxy cores and their C.O.M
// numerical integration to compute position, velocity, acceleration and jerk of
// two galaxies due to the influence of each other's gravitational attraction
// forces. It uses fourth order predictor-corrector Hermite scheme.
void evolveGalaxies()
{
	for (int i=0; i<M; ++i)
	{
		for (int j=i+1; j<M; ++j)
		{
			float G_rji[3], G_vji[3];  // position and velocity vectors
			for (int k=0; k<3; ++k)
			{
				G_rji[k] = G_r[j][k]-G_r[i][k]; // position vector form glx 1 to glx 2 cores
				G_vji[k] = G_v[j][k]-G_v[i][k]; // velocity vector
			}
			float G_r2 = 0;                   // initialize distance betn glx 1 and 2 cores
			for (int k=0; k<3; ++k)
				G_r2 += G_rji[k]*G_rji[k];      // compute distance square
			float G_r3 = (G_r2+ep)*sqrt((G_r2+ep)); // compute distance cube
			float G_rv = 0;                   // intialize dot product of position and
			for (int k=0; k<3; ++k)           // velocity vector
				G_rv += G_rji[k]*G_vji[k];      // compute dot product
			// dot product of position and velocity vectors divided by distance squared.
			// Notice the addition of threshold distance 'ep' in the denominator. This
      // is very important because we are treating the galaxies cores and stars as
			// point particles. So, G_r2, which is the distance between two galaxies
			// cores can go to zero (giving G_rv=infinity), but is physically not 
      // posible. So to avoid this situation we have to add some threshold 
			// separation distance
			G_rv /= (G_r2+ep);                
			for (int k=0; k<3; ++k)
			{
				G_a[i][k] += G_m[j] * G_rji[k]/G_r3; // compute acceleration for glx 1 core
				G_a[j][k] -= G_m[i] * G_rji[k]/G_r3; // compute acceleration for glx 2 core
				G_jk[i][k] += G_m[j]*(G_vji[k]-3*G_rv*G_rji[k])/G_r3; // jerk for glx 1 core
				G_jk[j][k] -= G_m[i]*(G_vji[k]-3*G_rv*G_rji[k])/G_r3; // jerk for glx 2 core
			}
		}
	}
	
	// arrays to store previous values of pos, vel, acc, and jerk
	float old_Gr[2][3], old_Gv[2][3], old_Ga[2][3], old_Gjk[2][3];
	
		for (int i=0; i<M; ++i)
		{
			for (int k=0; k<3; ++k)
			{
				// save previous values
				old_Gr[i][k] = G_r[i][k];
				old_Gv[i][k] = G_v[i][k];
				old_Ga[i][k] = G_a[i][k];
				old_Gjk[i][k] = G_jk[i][k];
				// predict new position and velocities (need to be corrected later)
				// this is the predictor step in the Hermite predictor-corrector step
				G_r[i][k] += G_v[i][k]*dt + G_a[i][k]*dt*dt/2 + G_jk[i][k]*dt*dt*dt/6;
				G_v[i][k] += G_a[i][k]*dt + G_jk[i][k]*dt*dt/2;
			}
		}
	
	// following code is for the corrector step
		for (int i=0; i<M; ++i)
			for (int k=0; k<3; ++k)
				G_a[i][k] = G_jk[i][k] = 0.0; // reinitialize acc and jerk to zero

		for (int i=0; i<M; ++i)
		{
			for (int j=i+1; j<M; ++j)
			{
				float G_rji[3], G_vji[3];
				for (int k=0; k<3; ++k)
				{
					G_rji[k] = G_r[j][k]-G_r[i][k]; // position vector
					G_vji[k] = G_v[j][k]-G_v[i][k]; // velocity vector
				}
				float G_r2 = 0;                   // distance squared
				for (int k=0; k<3; ++k)
					G_r2 += G_rji[k]*G_rji[k];      // compute distance squared
				float G_r3 = (G_r2+ep)*sqrt((G_r2+ep)); // distance cubed
				float G_rv = 0;                   // dot product or pos and vel
				for (int k=0; k<3; ++k)
					G_rv += G_rji[k]*G_vji[k];      // compute dot product
				G_rv /= (G_r2+ep);                // dot product divided by dis squared
				for (int k=0; k<3; ++k)
				{
					G_a[i][k] += G_m[j]*G_rji[k]/G_r3;   // new acc of glx 1 core
					G_a[j][k] -= G_m[i] * G_rji[k]/G_r3; // new acc of glx 2 core
					G_jk[i][k] += G_m[j]*(G_vji[k]-3*G_rv*G_rji[k])/G_r3; // new jerk glx 1 core
					G_jk[j][k] -= G_m[i]*(G_vji[k]-3*G_rv*G_rji[k])/G_r3; // new jerk glx 2 core
				}
			}
		}
		for (int i=0; i<M; ++i)
		{
			for (int k=0; k<3; ++k)
			{
				// position and velocity update
				// this is the corrector step
				G_v[i][k] = old_Gv[i][k] + (old_Ga[i][k] + G_a[i][k])*dt/2
																 + (old_Gjk[i][k] - G_jk[i][k])*dt*dt/12;	
				G_r[i][k] = old_Gr[i][k] + (old_Gv[i][k] + G_v[i][k])*dt/2
																 + (old_Ga[i][k] - G_a[i][k])*dt*dt/12;	
			}
		} 

		
		float t_mass = G_m[0]+G_m[1]; // total mass
		// coordinates of C.O.M. of the two galaxy cores
		float com_x = (G_m[0]*G_r[0][0]+G_m[1]*G_r[1][0])/t_mass;
		float com_y = (G_m[0]*G_r[0][1]+G_m[1]*G_r[1][1])/t_mass;
		float com_z = (G_m[0]*G_r[0][2]+G_m[1]*G_r[1][2])/t_mass;
		// store new C.O.M.
		COM[0].push_back(com_x);
		COM[1].push_back(com_y);
		COM[2].push_back(com_z);
		// store new galaxy core positions
		G1_x.push_back(G_r[0][0]);
		G1_y.push_back(G_r[0][1]);
		G1_z.push_back(G_r[0][2]);
		G2_x.push_back(G_r[1][0]);
		G2_y.push_back(G_r[1][1]);
		G2_z.push_back(G_r[1][2]);

		tm1.push_back(tt);
		// compute new kinetic energies of galaxy cores
		float k1=0, k2=0;
		k1 = 0.5*G_m[0]*(G_v[0][0]*G_v[0][0]+G_v[0][1]*G_v[0][1]+G_v[0][2]*G_v[0][2]);	
		k2 = 0.5*G_m[1]*(G_v[1][0]*G_v[1][0]+G_v[1][1]*G_v[1][1]+G_v[1][2]*G_v[1][2]);	
		tt += dt;             // update time
		// store new k.e and time
		KE_G1.push_back(k1);  
		KE_G2.push_back(k2);

		// compute new potential energies of galaxy cores
		float G_rji[3];
		for (int k=0; k<3; ++k)
			G_rji[k] = G_r[1][k] - G_r[0][k]; // distance between galaxy cores
		float G_r2 = 0;
		for (int k=0; k<3; ++k)
			G_r2 += G_rji[k]*G_rji[k];         // distance squared
		float p1= -G_m[0]*G_m[1]/sqrt(G_r2); // new P.E
		PE_G1G2.push_back(p1);               // store new P.E
}

// evolve stars
// this function updates the position, velocity, acceleration, jerk of all the 
// stars due to the gravitational influence of the two galaxy cores.
// all the code in this function is basically similar to the above function.
void evolveStars()
{
	float S_m[N];                              // star-mass array
	for (int i=0; i<N; ++i)
		S_m[i] = 1;                              // set all star-mass to one
	for (int i=0; i<M; ++i)
	{
		for (int j=0; j<N; ++j)
		{
			float S_rji[3], S_vji[3];              // position and velocity arrays
			for (int k=0; k<3; ++k)
			{
				S_rji[k] = S_r[j][k]-G_r[i][k];      // position vector
				S_vji[k] = S_v[j][k]-G_v[i][k];      // velocity vector
			}
			float S_r2 = 0;
			for (int k=0; k<3; ++k)
				S_r2 += S_rji[k]*S_rji[k];            // distance squared
			float S_r3 = (S_r2+ep)*sqrt((S_r2+ep)); // distance cubes
			float S_rv = 0;
			for (int k=0; k<3; ++k)
				S_rv += S_rji[k]*S_vji[k];            // dot product
			S_rv /= (S_r2+ep);                      // dot product divided by dis sq.
			for (int k=0; k<3; ++k)
			{
				G_a[i][k] += S_m[j] * S_rji[k]/S_r3;  // acc of galaxies
				S_a[j][k] -= G_m[i] * S_rji[k]/S_r3;  // acc of stars
				G_jk[i][k] += S_m[j]*(S_vji[k]-3*S_rv*S_rji[k])/S_r3; // jerk of galaxies
				S_jk[j][k] -= G_m[i]*(S_vji[k]-3*S_rv*S_rji[k])/S_r3; // jerk of stars
			}
		}
	}
	// arrays to store previous values
	float old_Sr[N][3], old_Sv[N][3], old_Sa[N][3], old_Sjk[N][3];
		for (int i=0; i<N; ++i)
		{
			for (int k=0; k<3; ++k)
			{
				// save previous time-step data
				old_Sr[i][k] = S_r[i][k];
				old_Sv[i][k] = S_v[i][k];
				old_Sa[i][k] = S_a[i][k];
				old_Sjk[i][k] = S_jk[i][k];
				// predict new position and velocities of stars (predictor)
				S_r[i][k] += S_v[i][k]*dt + S_a[i][k]*dt*dt/2 + S_jk[i][k]*dt*dt*dt/6;
				S_v[i][k] += S_a[i][k]*dt + S_jk[i][k]*dt*dt/2;
			}
		}

		// the following part is for the corrector step
		for (int i=0; i<N; ++i)
			for (int k=0; k<3; ++k)
				S_a[i][k] = S_jk[i][k] = 0.0;           // reset acc and jert to zero

		for (int i=0; i<M; ++i)
		{
			for (int j=0; j<N; ++j)
			{
				float S_rji[3], S_vji[3];
				for (int k=0; k<3; ++k)
				{
					S_rji[k] = S_r[j][k]-G_r[i][k];       // position vector
					S_vji[k] = S_v[j][k]-G_v[i][k];       // velocity vector
				}
				float S_r2 = 0;
				for (int k=0; k<3; ++k)
					S_r2 += S_rji[k]*S_rji[k];            // distance squared
				float S_r3 = (S_r2+ep)*sqrt((S_r2+ep)); // distance cubed
				float S_rv = 0;
				for (int k=0; k<3; ++k)
					S_rv += S_rji[k]*S_vji[k];            // dot product
				S_rv /= (S_r2+ep);                      // dot product divided by dis sq.
				for (int k=0; k<3; ++k)
				{
					G_a[i][k] += S_m[j]*S_rji[k]/S_r3;    // new acc of glx cores
					S_a[j][k] -= G_m[i] * S_rji[k]/S_r3;  // new acc of stars
					G_jk[i][k] += S_m[j]*(S_vji[k]-3*S_rv*S_rji[k])/S_r3; // new jerk glx cores
					S_jk[j][k] -= G_m[i]*(S_vji[k]-3*S_rv*S_rji[k])/S_r3; // new jerk stars
				}
			}
		}
		for (int i=0; i<N; ++i)
		{
			for (int k=0; k<3; ++k)
			{
				// postion and velocity update of stars (corrector step)
				S_v[i][k] = old_Sv[i][k] + (old_Sa[i][k] + S_a[i][k])*dt/2
																 + (old_Sjk[i][k] - S_jk[i][k])*dt*dt/12;	
				S_r[i][k] = old_Sr[i][k] + (old_Sv[i][k] + S_v[i][k])*dt/2
																 + (old_Sa[i][k] - S_a[i][k])*dt*dt/12;	
			}
		} 
}

// draw galaxies
void drawGalaxy()
{
	// velocity varaibles for determing the color of each star
	float vel_dis;  // velocity
	float norm_vel; // normalized velocity
	float v_max;    // maximum velocity
	float v_min;    // minimum velocity
	float v_next;   // velocity of next star in the list

	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);	

	drawAxis();           // draw 3-d axis at the bottom left corner

	// texts to display on the screen
	glColor3f(1,1,1);
	drawBitmapText("Press `h' to view options",p_fact*34,-300,0,2,'h');
	drawBitmapText("Velocity scale factor : ",-p_fact*55,p_fact*38,0,2,'h');
	sprintf(str_v_fact,"%2.2f",v_fact);
	drawBitmapText(str_v_fact,-p_fact*38,p_fact*38,0,2,'h');

	if (t==false)
	{
		drawBitmapText(val_max,p_fact*54,p_fact*37.5,0,2,'h');
		drawBitmapText(val_min,p_fact*54,-p_fact*39,0,2,'h');
		drawColorBar();
	}
	else
	{
		glColor3f(0,1,1);                          // set color for k.e of glx 1
		drawBitmapText("Galaxy 1",p_fact*47,p_fact*41,0,2,'h');
		glColor3f(1,1,0);                          // set color for k.e of glx 1
		drawBitmapText("Galaxy 2",p_fact*47,p_fact*38,0,2,'h');
		glColor3f(1,1,1);                          // set color for k.e of glx 1
		drawBitmapText("C.O.M.",p_fact*47,p_fact*35,0,2,'h');
	}

	//////////////////////// start drawing stars /////////////////////////
	glPushMatrix();
		glScaled(s_f,s_f,s_f); // for zooming effect without zooming the 3d coordinate
		glTranslated(xT,yT,0); // for translation effect
		glRotated(spin_x,0.0f,1.0f,0.0f); // rotate along x-axis
		glRotated(spin_y,1.0f,0.0f,0.0f); // rotate along y-axis
	glColor3f(1,1,1);

	if (o)
		drawBoundary();  // draw 3-d cube

	if (t==false)
	{
		//////// find maximum velocity for the current time step ///////////
		for (int k=0; k<3; ++k)
			v_max += S_v[0][k]*S_v[0][k]/(v_fact*v_fact);
		v_max = sqrt(v_max);
		v_min = v_max;
		for (int i=1; i<N; ++i)
		{
			v_next = 0;
			for (int k=0; k<3; ++k)
				v_next += S_v[i][k]*S_v[i][k]/(v_fact*v_fact);
			v_next = sqrt(v_next);
			if (v_next > v_max)
				v_max = v_next;
			if (v_next < v_min)
				v_min = v_next;
		}	
		sprintf(val_min,"%2.2f",v_min);
		sprintf(val_max,"%2.2f",v_max);
		/////////////////////////////////////////////////////////////////////
		glPointSize(1);                        // set point size

		loadColor();                           // load color data

		glBegin(GL_POINTS);                    // begin drawing stars
			for (int i=0; i<N; ++i)
			{
				vel_dis = 0;
				for (int k=0; k<3; ++k)
					vel_dis += S_v[i][k]*S_v[i][k]/(v_fact*v_fact);	  // velocity squared
				vel_dis = sqrt(vel_dis);
				norm_vel = vel_dis/v_max;           // normalize velocity
				// this norm_vel is used for determing the colors for stars
				for(int j=0; j<257; ++j)          // set star colors
					if (norm_vel>cVal[j].fval && norm_vel<=cVal[j+1].fval)
						glColor3f(cVal[j].rgb[0],cVal[j].rgb[1],cVal[j].rgb[2]);
				glVertex3f(S_r[i][0],S_r[i][1],S_r[i][2]); // draw stars
			}
		glEnd();
		}
		else   // trace the tracks of galaxy cores and their C.O.M
		{
			glColor3f(0,1,1);                      // set color for the path of glx 1 core
			glBegin(GL_LINE_STRIP);  					    
			for (int i=0; i<G1_x.size(); ++i)
				glVertex3f(G1_x[i],G1_y[i],G1_z[i]); // trace of galaxy 1 core
			glEnd();
			glColor3f(1,1,0);                      // set color for the path of glx 2 core
			glBegin(GL_LINE_STRIP);
			for (int i=0; i<G2_x.size(); ++i)
				glVertex3f(G2_x[i],G2_y[i],G2_z[i]); // trace of galaxy 2 core
			glEnd();
			glColor3f(1,1,1);                      // set color for the path of the C.O.M
			glBegin(GL_LINE_STRIP);
			for (int i=0; i<G2_x.size(); ++i)
				glVertex3f(COM[0][i],COM[1][i],COM[2][i]); // trace of the C.O.M
			glEnd();
		}
	glPopMatrix();

	glutSwapBuffers();
}

// draw kinetic and potential energy curves for the galaxy cores
void drawEnergy()
{
	char t_1[5][10];  // time string
	char k_1[5][10];  // k.e. of glx 1 core 
	char k_2[5][10];  // k.e. of glx 2 core
	char p_12[5][10]; // p.e. between glx 1 and 2
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);	

	glColor3f(0,1,1);  
	drawBitmapText("Galaxy 1",p_fact*47,p_fact*41,0,2,'h');
	glColor3f(1,1,0);  
	drawBitmapText("Galaxy 2",p_fact*47,p_fact*38,0,2,'h');

	glPushMatrix();
		glScaled(s_f,s_f,s_f);     // scaling effect
		glTranslated(xT,yT,0);     // translation effect
		glTranslated(-275,50,0);
		glColor3f(1,1,1);

		// upper plotting box
		glBegin(GL_LINES);
				glColor3f(1,1,0);
				glVertex2f(0,205);
				glVertex2f(0,-25);
				glColor3f(0,1,1);
				glVertex2f(550,205);
				glVertex2f(550,-25);
		glEnd();
		glBegin(GL_LINES);
				glColor3f(1,1,1);
				glVertex2f(0,-25);
				glVertex2f(550,-25);
				glVertex2f(0,205);
				glVertex2f(550,205);
		glEnd();
		// lower plotting box
		glBegin(GL_LINES);
				glVertex2f(0,-70);
				glVertex2f(0,-300);
				glVertex2f(550,-300);
				glVertex2f(550,-70);
		glEnd();
		glBegin(GL_LINES);
				glVertex2f(0,-300);
				glVertex2f(550,-300);
				glVertex2f(0,-70);
				glVertex2f(550,-70);
		glEnd();
	
		// draw tick marks on the plotting boxes
		glBegin(GL_LINES);
		for (int i=0; i<4; ++i)
		{
			// upper box -> tick marks for top and bottom boundaries
			glColor3f(1,1,1);
			glVertex2f((i+1)*(550/5),-27.5);
			glVertex2f((i+1)*(550/5),-22.5);
			glVertex2f((i+1)*(550/5),207.5);
			glVertex2f((i+1)*(550/5),202.5);

			// upper box -> tick marks for left and right boundaries
			glColor3f(1,1,0);
			glVertex2f(-2.5,(i+1)*(225/5)-25);
			glVertex2f(2.5,(i+1)*(225/5)-25);
			glColor3f(0,1,1);
			glVertex2f(547.5,(i+1)*(225/5)-25);
			glVertex2f(552.5,(i+1)*(225/5)-25);

			// lower box -> tick marks for top and bottom boundaries
			glColor3f(1,1,1);
			glVertex2f((i+1)*(550/5),-302.5);
			glVertex2f((i+1)*(550/5),-297.5);
			glVertex2f((i+1)*(550/5),-67.5);
			glVertex2f((i+1)*(550/5),-72.5);

			// lower box -> tick marks for left and right boundaries
			glColor3f(1,1,1);
			glVertex2f(547.5,(i+1)*(225/5)-300);
			glVertex2f(552.5,(i+1)*(225/5)-300);
			glColor3f(1,1,1);
			glVertex2f(-2.5,(i+1)*(225/5)-300);
			glVertex2f(2.5,(i+1)*(225/5)-300);
		}
		glEnd();

		// find maximum k.e. values and minimum p.e. value for scaling the plotting
		// axis values appropriately at runtime
		float max_KE_G1 = *max_element(KE_G1.begin(),KE_G1.end());
		float max_KE_G2 = *max_element(KE_G2.begin(),KE_G2.end());
		float min_PE_G1G2 = *min_element(PE_G1G2.begin(),PE_G1G2.end());
		
		// set and level origins for plotting boxes
		glColor3f(1,1,1);
		drawBitmapText("0.00",-5,-38,0,2,'h');
		drawBitmapText("0.00",-5,-313,0,2,'h');

		// set and level axis values dynamically at runtime
		for (int i=0; i<5; ++i)
		{
			sprintf(t_1[i],"%2.2f",((i+1.0)/5)*tm1.back());
			sprintf(k_1[i],"%2.2f",((i+1.0)/5)*max_KE_G1/1e+5);    // Note scaling by 1e+5 to
			sprintf(k_2[i],"%2.2f",((i+1.0)/5)*max_KE_G2/1e+5);    // avoid big numbers on
			sprintf(p_12[i],"%2.2f",((i+1.0)/5)*min_PE_G1G2/1e+5); // axes
			// position the axes value at proper places
			glColor3f(1,1,1);
			drawBitmapText(t_1[i],((i+1.0)/5)*550-10,-38,0,2,'h');
			glColor3f(1,1,0);
			drawBitmapText(k_1[i],-40,((i+1.0)/5)*225-28,0,2,'h');
			glColor3f(0,1,1);
			drawBitmapText(k_2[i],560,((i+1.0)/5)*225-28,0,2,'h');
			glColor3f(1,1,1);
			drawBitmapText(t_1[i],((i+1.0)/5)*550-10,-313,0,2,'h');
			glColor3f(1,1,1);
			drawBitmapText(p_12[i],-45,((i+1.0)/5)*225-303,0,2,'h');
		}

		// set axis titles
		glColor3f(1,1,1);
		drawBitmapText("Time",250,-330,0,2,'h');            // x-axis label
		drawBitmapText("Kinetic Energy",-60,200,0,2,'v');   // y-axis label upper box
		drawBitmapText("Potential Energy",-60,-70,0,2,'v'); // y-axis label lower box

		// draw kinetic energy curve for galaxy 1 core on the upper box
		glColor3f(0,1,1);                         // set color for k.e of glx 1 core
		glBegin(GL_LINE_STRIP);
		for (int i=0; i<KE_G1.size(); ++i)
			glVertex2f(tm1[i]*550/tm1.back(),225*KE_G1[i]/max_KE_G1-25); // k.e of glx 1 core
		glEnd();
		// draw kinetic energy curve for galaxy 2 core on the upper box
		glColor3f(1,1,0);                         // set color for k.e of glx 2 core
		glBegin(GL_LINE_STRIP);
		for (int i=0; i<KE_G2.size(); ++i)
			glVertex2f(tm1[i]*550/tm1.back(),225*KE_G2[i]/max_KE_G2-25); // k.e of glx 2 core
		glEnd();
		// draw potetial energy curve on the lower box
		glColor3f(1,0,0);                         // set color for p.e between the glx cores
		glBegin(GL_LINE_STRIP);
		for (int i=0; i<KE_G2.size(); ++i)
			glVertex2f(tm1[i]*550/tm1.back(),225*PE_G1G2[i]/min_PE_G1G2-300); // p.e of glx cores
		glEnd();
	glPopMatrix();
	
	glutSwapBuffers();
}

// wrapper for drawGalaxies() and drawEnergy()
void draw()
{
	// triggered in menu items
	if (display==GLX) 
		drawGalaxy();
	else if (display==KE)
		drawEnergy();
}

// process menu items
void processMenuEvents(int option)
{
	switch (option)
	{
		case GLX: display=GLX; break; // display glx collision simulation
		case KE : display=KE;  break; // display energy curves
	}
}	
// keyboard functionalities
void myKeyboard(unsigned char theKey, int x, int y)
{
	switch(theKey)
	{
		default: 
			if (theKey == 27) 
			{
				freeMemory();
				exit(0); // exit the program
			} break;
		case 'q':
			freeMemory();
		  exit(0);break;       // quit the program
		case '+':                // zooming in
			s_f += 0.05; break;
		case '-':                // zooming out
			s_f -= 0.05; break;
		case 'c':                // change color
			(changeColor==true) ? changeColor=false : changeColor=true; break;
		case 'f':                 // fullscreen on/off
			if (fullScreen==true)
			{
				glutFullScreen();
				fullScreen = false;
			}
			else
			{
				glutPositionWindow(0,0);
				glutReshapeWindow(WIDTH,HEIGHT);
				fullScreen = true;
			}
			break;
		case 'h':                  // diplay info
			displayInfo(); break;
		case 'i':                  // set everything back to the original setting
			s_f     = 1; xT     = 0; yT     = 0; 
			spin_x = 0; spin_y = 0;
			old_x  = 0; old_y  = 0; break;
		case 's':                  // animation on/off
			(animate==false) ? animate = true : animate = false; break;
		case 'o':                  // boundary on/off
			(o==true) ? o = false : o = true; break;
		case 'r':                  // reset the simulation with different init vel.
		{
			// clear data from previous animation run
			freeMemory();
			tt = 0;
			float range = 21.0;
			// set new scaling factor for initial velocity. This will change the initial
			// velocity for galaxy collision simulation
			v_fact = 5.0+float(range*rand()/(RAND_MAX+1.0));
			initCondition(); // initialize galaxies and stars with new init cond.
			animate=false;
			break;
		}
		case 't':          // trace the path of C.O.M.s of galaxies - on/off
			(t==false) ? t = true : t = false; break;
	}

	MyInit();
	glutPostRedisplay();
}

// deallocate memory of vector datatypes
void freeMemory()
{
	KE_G1.clear();   
	KE_G2.clear();
	PE_G1G2.clear();
	tm1.clear();
	G1_x.clear();G1_y.clear();G1_z.clear();
	G2_x.clear();G2_y.clear();G2_z.clear();
	COM[3].clear();	
}

// mouse functions
void mouse(int button, int state, int x, int y)
{
	old_y = y-spin_y;
	old_x = x-spin_x;

	glutPostRedisplay();
}

// handle mouse motion for rotation
void motion(int x,int y)
{
	spin_x = x-old_x;
	spin_y = y-old_y;

	glutPostRedisplay();
}

// timer function for animation - rotation along y-axis
void TimerFunc(int value)
{
	if (animate == true)
	{	
		evolveStars();    // evolve stars
		evolveGalaxies(); // evolve galaxies
	}
	glutPostRedisplay(); // redisplay
	glutTimerFunc(50,TimerFunc,1); // recursive call to TimerFunc
}

// special keyboard functions
void specialKeys(int key, int x, int y)
{
	if (key == GLUT_KEY_LEFT) // move left
		xT -= 5;
	if (key == GLUT_KEY_RIGHT) // move right
		xT += 5; 
	if (key == GLUT_KEY_DOWN) // move down
		yT -= 5;
	if (key == GLUT_KEY_UP) // move up
		yT += 5;

	MyInit();
	glutPostRedisplay();
}

// load color data
void loadColor()
{
	ifstream myCl, myCl1;
	// cool color data "http://www.sandia.gov/~kmorel/documents/ColorMaps"
	myCl.open("Jet257.dat");   
	assert(myCl);

	myCl1.open("CoolWarmFloat257.dat");
	assert(myCl1);

	if (changeColor==true)
	{
		for (int i=0; i<257; ++i)
		{
			// read warm color data
			myCl1 >> cVal[i].fval >> cVal[i].rgb[0] >> cVal[i].rgb[1] >> cVal[i].rgb[2];
		}
	}
	else if (changeColor==false)
	{
		// read cool color data
		for (int i=0; i<257; ++i)
			myCl >> cVal[i].fval >> cVal[i].rgb[0] >> cVal[i].rgb[1] >> cVal[i].rgb[2];
	}
	myCl1.close();
	myCl.close();	
}

// draw 3-d axis
void drawAxis()
{
	glPushMatrix();
		glTranslated(-p_fact*35,-p_fact*25,p_fact*25); // set position
		glScaled(25,25,25);                            // scale
		glRotated(spin_x,0.0f,1.0f,0.0f); // rotate along x-axis
		glRotated(spin_y,1.0f,0.0f,0.0f); // rotate along y-axis

		glColor3f(1,1,1);
		drawBitmapText("x",1.1,0,0,2,'v'); // x-label
		drawBitmapText("y",0,1.1,0,2,'v'); // y-label
		drawBitmapText("z",0,0,1.1,2,'v'); // z-label
		glBegin(GL_LINES);
		glColor3f(1,0,0);                  // red for x-axis
		glVertex3f(0,0,0);
		glVertex3f(1,0,0);
		glColor3f(0,1,0);                  // green for y-axis
		glVertex3f(0,0,0);
		glVertex3f(0,1,0);
		glColor3f(0,0,1);                  // blue for z-axis
		glVertex3f(0,0,0);
		glVertex3f(0,0,1);
		glEnd();
	glPopMatrix();	
}

// draw velocity colorbar
void drawColorBar()
{
	glPushMatrix();
		glTranslated(p_fact*55,-p_fact*36.5,0); // start from lower right portion of the screen
		// use quads to draw color rectagles
		glBegin(GL_QUAD_STRIP);
			for (int i=0; i<257; ++i)
				{
					glColor3f(cVal[i].rgb[0],cVal[i].rgb[1],cVal[i].rgb[2]); // color of each quad
					glVertex2f(0,2*i);  
					glVertex2f(10,2*i);
				}
		glEnd();
	glPopMatrix();
}

// wrapper for glutBitmapCharacter function
void drawBitmapText(const char* str,float x,float y,float z,int n,char c)
{
  int len = strlen(str);   // string length
	if (c=='h')              // for horizontal display
  	glRasterPos3f(x, y,z); // display position
  if (n==1)
  {
    for (int i=0; i<len; ++i)
		{
			if (c=='v')          // for vertical display
  			glRasterPos3f(x, y-i*5,z); // display position
      glutBitmapCharacter(GLUT_BITMAP_HELVETICA_18, str[i]); // font 1
		}
  }
  else if (n==2)
  {
    for (int i=0; i<len; ++i)
		{
			if (c=='v')          // for vetical display
  			glRasterPos3f(x, y-i*15,z); // display position
      glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12, str[i]); // font 2
		}
  }
}

// display info with all the available options
void displayInfo()
{
	cout<<"------------------Andromeda-MilkyWay Collision Simulation--------------------"<<endl;
	cout<<"Programmer  : Vijay Koju"<<endl;
	cout<<"Class       : CSCI 7300 (Scientific Visualization and Databases)"<<endl;
	cout<<"Instructor  : Dr. Li"<<endl;
	cout<<"Institution : Middle Tennessee State University (MTSU)"<<endl;
	cout<<"Department  : Computational Science Program"<<endl;
	cout<<"-----------------------------------------------------------------------------"<<endl;
	cout<<"This program simulates collision dynamics of Andromeda and MilkyWay galaxies. This kind of problem is often called N-Body problem, and in general the complexity of such a problem is O(n^2) for each simulation step because it requires the computation of gravitational force on each body by all other bodies. There are many other algorithm such as Barnes-Huts Tree algorithm (O(nlogn)), and Fast Multipole Method (O(n)) which use some approximations to reduce the time complexity. This program, however, uses none of such algorithms. It approximates the two galaxies as two huge mass stars located at the center of mass of the galaxies with some huge mass. These two stars interact according to the laws of gravitation. All the stars in both the galaxies are only influenced by these two core stars. The numerical integration in this program is done using the fourth order Hermite Scheme. Detaails of this numerical scheme can be found in 'Movings Stars Around' by Piet Hut and Jun Makino, 2007."<<endl;
	cout<<"-----------------------------------------------------------------------------"<<endl;
	cout<<"--------------------Keyboard options:----------------------"<<endl;
	cout<<"## h                   --> view this page                ##"<<endl;
	cout<<"## s                   --> start/stop animation          ##"<<endl;
	cout<<"## f                   --> fullscreen mode on/off        ##"<<endl;
	cout<<"## c                   --> change color warm/cool        ##"<<endl;
	cout<<"## o                   --> boundary (box) on/off         ##"<<endl;
	cout<<"## i                   --> go back to initial setting    ##"<<endl;
	cout<<"## r                   --> reset with different initial  ##"<<endl;
	cout<<"##                         velocity                      ##"<<endl;
	cout<<"## t                   --> trace the path of the galaxy  ##"<<endl;
	cout<<"##                         and the \"Center of mass\" of   ##"<<endl;
	cout<<"##                         the two galaxies              ##"<<endl;
	cout<<"## q/esc               --> quit the program              ##"<<endl;
	cout<<"## Ctrl +              --> zoom in                       ##"<<endl;
	cout<<"## - (minus sign)      --> zoom out                      ##"<<endl;
	cout<<"## Right arrow         --> move to the right             ##"<<endl;
	cout<<"## Left arrow          --> move to the left              ##"<<endl;
	cout<<"## Up arrow            --> move up                       ##"<<endl;
	cout<<"## Down arrow          --> move down                     ##"<<endl;
	cout<<"-----------------------Mouse options:----------------------"<<endl;
	cout<<"## Left click and drag --> rotate                        ##"<<endl;
	cout<<"## Right click         --> pop menu to switch between    ##"<<endl;
	cout<<"##                     --> galaxy collision mode and     ##"<<endl;
	cout<<"##                         energy mode                   ##"<<endl;
	cout<<"-----------------------------------------------------------------------------"<<endl;
}
