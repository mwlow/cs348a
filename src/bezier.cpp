#include <Eigen/Core>
#include <Eigen/Eigenvalues>
#include <iostream>
#include "bezier.h"
#include <GL/glut.h>
using namespace OpenMesh;
using namespace Eigen;
using namespace std;

float tol = 0.01;

void cubicBezier(Vec3f p[4], vector<Vec3f>& points)
{
    if(((p[0]-p[1]).length()+(p[1]-p[2]).length() + (p[2]-p[3]).length())<tol){
        points.push_back(p[0]);
        //points.push_back(p[3]);
        return;
    }
    //calculate the first new set of control points
    //push into vector and return if the two distance between two ending points are small enough
    Vec3f q[4],r[4];
    q[0] = p[0];
    q[1] = (p[0] + p[1])/2.0;
    q[2] = (p[0] + p[2])/4.0 + p[1]/2.0;
    q[3] = (p[0] + p[3])/8.0 + (p[1] + p[2])*3.0/8.0;
    r[0] = q[3];
    r[1] = (p[1] + p[3])/4.0 + p[2]/2.0;
    r[2] = (p[2] + p[3])/2.0;
    r[3] = p[3];
    //first half
    cubicBezier(q, points);
    //calculate the seconde half
    cubicBezier(r, points);
    return;

}

void bezier(Vec3f p[3], vector<Vec3f>& points){//suppose p[0] p[1] p[2] in clockwise
    //construct points
    Vec3f a[4];
    a[0] = (p[0]*4.0 + p[1] + p[2])/6.0;
    a[1] = (p[0]*2.0 + p[1])/3.0;
    a[2] = (p[0] + p[1]*2.0)/3.0;
    a[3] = (p[0] + p[1]*4.0 + p[2])/6.0;
    cubicBezier(a, points);
    Vec3f b[4];
    b[0] = a[3];
    b[1] = (p[1]*2.0 + p[2])/3.0;
    b[2] = (p[1] + p[2]*2.0)/3.0;
    b[3] = (p[0] + p[1] + p[2]*4.0)/6.0;
    cubicBezier(b, points);
    Vec3f c[4];
    c[0] = b[3];
    c[1] = (p[2]*2.0 + p[0])/3.0;
    c[2] = (p[2] + p[0]*2.0)/3.0;
    c[3] = a[0];
    cubicBezier(c, points);
}

void holy(Vec3f p[3], vector<Vec3f>& points, Vec3f n) {
    //assume control points are in clockwise order

    //Dummy triangle corners, in clockwise order

    // must start at p0 for this code to work
    Vec3f *curr_vertex = &p[0];
    Vec3f *next_vertex = &p[1];

    float min_distance = FLT_MAX;
    vector<Vec3f>::iterator min = points.end();
    vector<Vec3f>::iterator next = points.end();
    for (vector<Vec3f>::iterator it = points.begin(); it != points.end(); ++it) {
        float distance = (*it - *curr_vertex).length();
        if (distance < min_distance) {
            min_distance = distance;
            min = it;
        }
    }

    vector<Vec3f>::iterator curr = min;
    do {
        next = (curr + 1 == points.end()) ? points.begin() : curr + 1;

        if (((*curr_vertex - *next) % (*curr - *next))[2] >0) {
            glBegin(GL_TRIANGLES);
            //anticlockwise ordering here
            //TODO: specify normals per vertex
            glColor3f(1, 0, 0);
            glNormal3f(n[0], n[1], n[2]);
            glVertex3f((*curr_vertex)[0], (*curr_vertex)[1], (*curr_vertex)[2]);
            glVertex3f((*curr)[0], (*curr)[1], (*curr)[2]);
            glVertex3f((*next)[0], (*next)[1], (*next)[2]);
            glEnd();
        }
        else {
            glBegin(GL_TRIANGLES);
            //anticlockwise ordering here
            //TODO: specify normals per vertex
            glColor3f(1, 0, 0);
            glNormal3f(n[0], n[1], n[2]);
            glVertex3f((*curr_vertex)[0], (*curr_vertex)[1], (*curr_vertex)[2]);
            glVertex3f((*curr)[0], (*curr)[1], (*curr)[2]);
            glVertex3f((*next_vertex)[0], (*next_vertex)[1], (*next_vertex)[2]);
            glEnd();

            if (curr_vertex == &p[0]) {
                curr_vertex = &p[1];
                next_vertex = &p[2];
                continue;
            }
            else if (curr_vertex == &p[1]) {
                curr_vertex = &p[2];
                next_vertex = &p[0];
                continue;
            }
            else {
                curr_vertex = &p[0];
                next_vertex = &p[1];
                continue;
            }
        }

        curr = (curr + 1 == points.end()) ? points.begin() : curr + 1;
    }
    while (curr != min);
}

