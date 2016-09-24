#include "mesh.h"
//////////////////////////////////////
//The Camera Class
//////////////////////////////////////
class Camera
{
    private:
	Vector3 u, v, n;   // camera coordinates
        Point3 eye;      // camera position
        Point3 look;     // camera view(target)
        Vector3 up;       // camera upvector(tilt)
	float viewAngle, aspect, nearDist, farDist; // view volume shape
	void setModelViewMatrix();  // tell OpenGL where the camera is

    public:
	Camera() {} ;    // constructor

	void set(float ex,  float ey,  float ez,
             float lx, float ly, float lz,
             float ux, float uy, float uz); // like gluLookAt()
	void set(Point3 eye, Point3 look, Vector3 up); // like gluLookAt()

    void setShape(float vAng, float asp, float nearD, float farD);
	void getShape(float &vAng, float &asp, float &nearD, float &farD);

	void roll(float angle);  // roll it

	void pitch(float angle);   // increase pitch

	void yaw(float angle);  // yaw it

	void slide(float delU, float delV, float delN); //slide it
};
