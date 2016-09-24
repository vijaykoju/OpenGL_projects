
#include "Camera.h"
#define GL_PI 2*asin(1) 

void Camera::setModelViewMatrix(void)
{
    // load modelview matrix with existing camera values
    float m[16];
    Vector3 eVec(eye.x, eye.y, eye.z);  // constructor of a vector verion of eye
   
    m[0] = u.x; m[4]=u.y; m[8]=u.z; m[12] = -eVec.dot(u);
    m[1] = v.x; m[5]=v.y; m[9]=v.z; m[13] = -eVec.dot(v);
    m[2] = n.x; m[6]=n.y; m[10]=n.z; m[14]= -eVec.dot(n);
    m[3] = 0;   m[7]=0;   m[11]=0;   m[15]=1;

    glMatrixMode(GL_MODELVIEW);
    glLoadMatrixf(m);
}

void Camera::set(float ex,  float ey,  float ez,                 
				float lx, float ly, float lz,
                float ux, float uy, float uz)
{
    eye.set(ex, ey, ez);
    look.set(lx, ly, lz);
    up.set(ux, uy, uz);

    // create a modelview matrix and send it to OpenGL
    n.set(eye.x - look.x, eye.y - look.y, eye.z - look.z);  // make n
    u.set(up.cross(n));  // make u = up X n

    n.normalize();
    u.normalize();
    v.set(n.cross(u));
    v.normalize();
    setModelViewMatrix();
}

void Camera::set(Point3 Eye, Point3 Look, Vector3 Up)
{
    // create a modelview matrix and send it to OpenGL
    eye.set(Eye); 
    look.set(Look);
    up.set(Up);

    n.set(eye.x - look.x, eye.y - look.y, eye.z - look.z);  // make n
    u.set(up.cross(n));  // make u = up X n
    n.normalize();
    u.normalize();
    v.set(n.cross(u));
    v.normalize();
    setModelViewMatrix();
}

void Camera::setShape(float vAng, float asp, float nearD, float farD)
{
    viewAngle = vAng;
    aspect = asp;
    nearDist = nearD;
    farDist = farD;

    GLdouble fW, fH;
    fH = tan((vAng / 360 * GL_PI) * nearD);
    fW = fH * aspect;

    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
 
    glFrustum(-fW, fW, -fH, fH, nearDist, farDist);
}

void Camera::getShape(float& vAng, float& asp, float& nearD, float& farD)
{
    vAng = viewAngle;
    asp = aspect;
    nearD = nearDist;
    farD = farDist;
}

void Camera::slide(float delU, float delV, float delN)
{
    eye.x += delU * u.x + delV * v.x + delN * n.x;
    eye.y += delU * u.y + delV * v.y + delN * n.y;
    eye.z += delU * u.z + delV * v.z + delN * n.z;

/*
    look.x += delU * u.x + delV * v.x + delN * n.x;
    look.y += delU * u.y + delV * v.y + delN * n.y;
    look.z += delU * u.z + delV * v.z + delN * n.z;
*/

    //setModelViewMatrix();
    set(eye, look, up);
}

//rotating around n
void Camera::roll(float angle)
{
    // roll the camera through angle degrees
    float cs = cos(GL_PI/180 * angle);
    float sn = sin(GL_PI/180 * angle);

    Vector3 t(u);   // remember old u
    u.set(cs*t.x - sn*v.x, cs*t.y-sn*v.y, cs*t.z - sn*v.z);
    v.set(sn*t.x + cs*v.x, sn*t.y+cs*v.y, sn*t.z + cs*v.z);

    setModelViewMatrix();
}

// rotating around v
void Camera::yaw(float angle)
{
    // yaw the camera through angle degrees
    float cs = cos(GL_PI/180 * angle);
    float sn = sin(GL_PI/180 * angle);

    Vector3 t(u);   // remember old u
    u.set(cs*t.x - sn*n.x, cs*t.y-sn*n.y, cs*t.z - sn*n.z);
    n.set(sn*t.x + cs*n.x, sn*t.y+cs*n.y, sn*t.z + cs*n.z);

    setModelViewMatrix();
}

// rotating around u
void Camera::pitch(float angle)
{
    // roll the camera through angle degrees
    float cs = cos(GL_PI/180 * angle);
    float sn = sin(GL_PI/180 * angle);

    Vector3 t(v);   // remember old v
    v.set(cs*t.x - sn*n.x, cs*t.y-sn*n.y, cs*t.z - sn*n.z);
    n.set(sn*t.x + cs*n.x, sn*t.y+cs*n.y, sn*t.z + cs*n.z);

    setModelViewMatrix();
}

