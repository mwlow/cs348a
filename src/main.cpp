#include <OpenMesh/Core/IO/Options.hh>
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <iostream>
#include <cmath>
#include <cstdio>
#include <GL/glut.h>
#include <Eigen/Core>
#include <Eigen/Eigenvalues>
#include "curvature.h"
#include "mesh_features.h"
#include "image_generation.h"
#include "decimate.h"
#include "bezier.h"
using namespace std;
using namespace OpenMesh;
using namespace Eigen;

VPropHandleT<double> viewCurvature;
FPropHandleT<Vec3f> viewCurvatureDerivative;
VPropHandleT<CurvatureInfo> curvature;
Mesh mesh;

double DWKW_THRESHOLD = 30;
double ANGLE_THRESHOLD = 0.1;
bool showBezier = true, light = false, showSuggestiveContour = false;
bool leftDown = false, rightDown = false, middleDown = false;
int lastPos[2];
float cameraPos[4] = {0,0,4,1};
Vec3f up, pan;
int windowWidth = 640, windowHeight = 480;
bool showSurface = true, showAxes = true, showCurvature = false, showNormals = false;

float specular[] = { 1.0, 1.0, 1.0, 1.0 };
float shininess[] = { 50.0 };

void renderSuggestiveContours(Vec3f actualCamPos) { // use this camera position to account for panning etc.
    glColor3f(1,0,0);

    // RENDER SUGGESTIVE CONTOURS HERE -----------------------------------------------------------------------------
    // -------------------------------------------------------------------------------------------------------------
    for (Mesh::FaceIter it = mesh.faces_begin(); it != mesh.faces_end(); ++it) {

        double kw[3];
        Vec3f p[3];
        Vec3f end_point[2];

        Mesh::ConstFaceVertexIter fvIt = mesh.cfv_iter(it);
        for (int i=0; i<3; i++) {
            p[i] = mesh.point(fvIt.handle());
            kw[i] = mesh.property(viewCurvature, fvIt.handle());
            ++fvIt;
        }

        if(kw[0]*kw[1]>=0 && kw[1]*kw[2]>=0)
            continue;

        //Calculate a representative w
        Vec3f c = (p[0] + p[1] + p[2])/3.0;
        Vec3f n = mesh.calc_face_normal(it.handle());

        Vec3f v = actualCamPos - c;
        Vec3f w = (v - (v | n) * n).normalize();

        //Calculate DwKw
        Vec3f Dw = mesh.property(viewCurvatureDerivative, it);
        double DwKw = Dw | w;

        //Calculate the angle between surface normal and view vector
        double theta = acos(n | v);
        if(DwKw < DWKW_THRESHOLD || (theta < ANGLE_THRESHOLD))
            continue;

        end_point[0] = (kw[0] * kw[1] < 0) ?
            p[0] + (p[1] - p[0]) * -kw[0]/(kw[1] - kw[0]) :
            p[0] + (p[2] - p[0]) * -kw[0]/(kw[2] - kw[0]) ;

        end_point[1] = (kw[1] * kw[2] < 0) ?
            p[1] + (p[2] - p[1]) * -kw[1]/(kw[2] - kw[1]) :
            p[0] + (p[2] - p[0]) * -kw[0]/(kw[2] - kw[0]) ;


        glBegin(GL_LINES);
        glVertex3f(end_point[0][0], end_point[0][1], end_point[0][2]);
        glVertex3f(end_point[1][0], end_point[1][1], end_point[1][2]);
        glEnd();
    }
}

void renderMesh() {
    if (!showSurface) glColorMask(GL_FALSE,GL_FALSE,GL_FALSE,GL_FALSE); // render regardless to remove hidden lines

    glEnable(GL_LIGHTING);
    glLightfv(GL_LIGHT0, GL_POSITION, cameraPos);

    glDepthRange(0.001,1);
    glEnable(GL_NORMALIZE);

    if(!light){
        glDisable(GL_LIGHTING);
        glDepthRange(0, 0.999);
    }
    // WRITE CODE HERE TO RENDER THE TRIANGLES OF THE MESH ---------------------------------------------------------
    for (Mesh::FaceIter f_it=mesh.faces_begin(); f_it!=mesh.faces_end(); ++f_it){
        OpenMesh::Vec3f point[3], normal[3];
        Mesh::ConstFaceVertexIter cfv_it;
        cfv_it = mesh.cfv_iter(f_it.handle());
        Vec3f p0 = mesh.point(cfv_it.handle());
        Vec3f n0 = mesh.normal(cfv_it.handle());
        Vec3f p1 = mesh.point((++cfv_it).handle());
        Vec3f n1 = mesh.normal(cfv_it.handle());
        Vec3f p2 = mesh.point((++cfv_it).handle());
        Vec3f n2 = mesh.normal(cfv_it.handle());

        if (((p1 - p0) % (p2 - p0))[2] < 0) {
            point[0] = p0;
            point[1] = p1;
            point[2] = p2;
            normal[0] = n0;
            normal[1] = n1;
            normal[2] = n2;
        }
        else {
            point[2] = p0;
            point[1] = p1;
            point[0] = p2;
            normal[2] = n0;
            normal[1] = n1;
            normal[0] = n2;
        }

        if(showBezier){
            //glPolygonMode( GL_FRONT_AND_BACK, GL_LINE );
            glColor3f(0,1,0);
            glBegin(GL_LINE_STRIP);
            glVertex3f(point[0][0], point[0][1], point[0][2]);
            glVertex3f(point[1][0], point[1][1], point[1][2]);
            glVertex3f(point[2][0], point[2][1], point[2][2]);
            glEnd();

            std::vector<Vec3f> points;
            bezier(point, points);
            holy(point, points, mesh.calc_face_normal(f_it.handle()));
            glColor3f(0,0,1);
            glBegin(GL_LINE_STRIP);
            for (std::vector<Vec3f>::iterator it = points.begin(); it != points.end(); ++it) {
                glVertex3f((*it)[0], (*it)[1], (*it)[2]);
            }
            glEnd();
            glColor3f(0.2, 0.2, 1);

        }/*
            glBegin(GL_TRIANGLES);
            glNormal3f(normal[0][0], normal[0][1], normal[0][2]);
            glVertex3f(point[0][0], point[0][1], point[0][2]);
            glNormal3f(normal[1][0], normal[1][1], normal[1][2]);
            glVertex3f(point[1][0], point[1][1], point[1][2]);
            glNormal3f(normal[2][0], normal[2][1], normal[2][2]);
            glVertex3f(point[2][0], point[2][1], point[2][2]);
            glEnd();
            */

    }
    // -------------------------------------------------------------------------------------------------------------

    if (!showSurface) glColorMask(GL_TRUE,GL_TRUE,GL_TRUE,GL_TRUE);

    glDisable(GL_LIGHTING);
    glDepthRange(0,0.999);

    Vec3f actualCamPos(cameraPos[0]+pan[0],cameraPos[1]+pan[1],cameraPos[2]+pan[2]);
    if(showSuggestiveContour){
        renderSuggestiveContours(actualCamPos);
    }
    // We'll be nice and provide you with code to render feature edges below
    glBegin(GL_LINES);
    glColor3f(1,0,0);
    glLineWidth(2.0f);
    for (Mesh::ConstEdgeIter it = mesh.edges_begin(); it != mesh.edges_end(); ++it)
        if (isFeatureEdge(mesh,*it,actualCamPos)) {
            Mesh::HalfedgeHandle h0 = mesh.halfedge_handle(it,0);
            Mesh::HalfedgeHandle h1 = mesh.halfedge_handle(it,1);
            Vec3f source(mesh.point(mesh.from_vertex_handle(h0)));
            Vec3f target(mesh.point(mesh.from_vertex_handle(h1)));
            glVertex3f(source[0],source[1],source[2]);
            glVertex3f(target[0],target[1],target[2]);
        }
    glEnd();

    if (showCurvature) {
        // WRITE CODE HERE TO RENDER THE PRINCIPAL DIRECTIONS YOU COMPUTED ---------------------------------------------
        glBegin(GL_LINES);
        for (Mesh::ConstVertexIter it = mesh.vertices_begin(); it != mesh.vertices_end(); ++it) {
            CurvatureInfo info = mesh.property(curvature, it.handle());
            Vec3f p = mesh.point(it.handle());
            Vec3f p0 = p - info.directions[0]*.01;
            Vec3f p1 = p - info.directions[1]*.01;
            Vec3f d0 = p + info.directions[0]*.01;
            Vec3f d1 = p + info.directions[1]*.01;
            glColor3f(1,0,0);
            glVertex3f(p0[0],p0[1],p0[2]);
            glVertex3f(d0[0],d0[1],d0[2]);
            glColor3f(0,0,1);
            glVertex3f(p1[0],p1[1],p1[2]);
            glVertex3f(d1[0],d1[1],d1[2]);
        }

        glEnd();

        // -------------------------------------------------------------------------------------------------------------
    }

    if (showNormals) {
        glBegin(GL_LINES);
        glColor3f(0,1,0);
        for (Mesh::ConstVertexIter it = mesh.vertices_begin(); it != mesh.vertices_end(); ++it) {
            Vec3f n = mesh.normal(it.handle());
            Vec3f p = mesh.point(it.handle());
            Vec3f d = p + n*.01;
            glVertex3f(p[0],p[1],p[2]);
            glVertex3f(d[0],d[1],d[2]);
        }
        glEnd();
    }

    glDepthRange(0,1);
}

void display() {
    glClearColor(1,1,1,1);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glEnable(GL_LINE_SMOOTH);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);

    glEnable(GL_DEPTH_TEST);
    glEnable(GL_LIGHTING);
    glShadeModel(GL_SMOOTH);
    glMaterialfv(GL_FRONT, GL_SPECULAR, specular);
    glMaterialfv(GL_FRONT, GL_SHININESS, shininess);
    glEnable(GL_LIGHT0);

    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glViewport(0,0,windowWidth,windowHeight);

    float ratio = (float)windowWidth / (float)windowHeight;
    gluPerspective(50, ratio, 1, 1000); // 50 degree vertical viewing angle, zNear = 1, zFar = 1000

    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    gluLookAt(cameraPos[0]+pan[0], cameraPos[1]+pan[1], cameraPos[2]+pan[2], pan[0], pan[1], pan[2], up[0], up[1], up[2]);

    // Draw mesh
    renderMesh();

    // Draw axes
    if (showAxes) {
        glDisable(GL_LIGHTING);
        glBegin(GL_LINES);
        glLineWidth(1);
        glColor3f(1,0,0); glVertex3f(0,0,0); glVertex3f(1,0,0); // x axis
        glColor3f(0,1,0); glVertex3f(0,0,0); glVertex3f(0,1,0); // y axis
        glColor3f(0,0,1); glVertex3f(0,0,0); glVertex3f(0,0,1); // z axis
        glEnd(/*GL_LINES*/);
    }

    glutSwapBuffers();
}

void mouse(int button, int state, int x, int y) {
    if (button == GLUT_LEFT_BUTTON) leftDown = (state == GLUT_DOWN);
    else if (button == GLUT_RIGHT_BUTTON) rightDown = (state == GLUT_DOWN);
    else if (button == GLUT_MIDDLE_BUTTON) middleDown = (state == GLUT_DOWN);

    lastPos[0] = x;
    lastPos[1] = y;
}

void mouseMoved(int x, int y) {
    int dx = x - lastPos[0];
    int dy = y - lastPos[1];
    Vec3f curCamera(cameraPos[0],cameraPos[1],cameraPos[2]);
    Vec3f curCameraNormalized = curCamera.normalized();
    Vec3f right = up % curCameraNormalized;

    if (leftDown) {
        // Assume here that up vector is (0,1,0)
        Vec3f newPos = curCamera - 2*(float)((float)dx/(float)windowWidth) * right + 2*(float)((float)dy/(float)windowHeight) * up;
        newPos = newPos.normalized() * curCamera.length();

        up = up - (up | newPos) * newPos / newPos.sqrnorm();
        up.normalize();

        for (int i = 0; i < 3; i++) cameraPos[i] = newPos[i];
    }
    else if (rightDown) for (int i = 0; i < 3; i++) cameraPos[i] *= pow(1.1,dy*.1);
    else if (middleDown) {
        pan += -2*(float)((float)dx/(float)windowWidth) * right + 2*(float)((float)dy/(float)windowHeight) * up;
    }


    lastPos[0] = x;
    lastPos[1] = y;

    Vec3f actualCamPos(cameraPos[0]+pan[0],cameraPos[1]+pan[1],cameraPos[2]+pan[2]);
    computeViewCurvature(mesh,actualCamPos,curvature,viewCurvature,viewCurvatureDerivative);

    glutPostRedisplay();
}

void keyboard(unsigned char key, int x, int y) {
    Vec3f actualCamPos(cameraPos[0]+pan[0],cameraPos[1]+pan[1],cameraPos[2]+pan[2]);

    if (key == 's' || key == 'S') showSurface = !showSurface;
    else if (key == 'a' || key == 'A') showAxes = !showAxes;
    else if (key == 'c' || key == 'C') showCurvature = !showCurvature;
    else if (key == 'n' || key == 'N') showNormals = !showNormals;
    else if (key == 'w' || key == 'W') writeImage(mesh, windowWidth, windowHeight, "renderedImage.svg", actualCamPos);
    else if (key == 'q' || key == 'Q') exit(0);
    else if (key == 'd' || key == 'D') DWKW_THRESHOLD += 1;
    else if (key == 't' || key == 'T') ANGLE_THRESHOLD += 0.1;
    else if (key == 'l' || key == 'L') light = !light;
    else if (key == 'b' || key == 'B') showBezier = !showBezier;
    else if (key == 'u' || key == 'u') showSuggestiveContour = !showSuggestiveContour;
    glutPostRedisplay();
}

void reshape(int width, int height) {
    windowWidth = width;
    windowHeight = height;
    glutPostRedisplay();
}

int main(int argc, char** argv) {
    if (argc < 2) {
        cout << "Usage: " << argv[0] << " mesh_filename\n";
        exit(0);
    }

    IO::Options opt;
    opt += IO::Options::VertexNormal;
    opt += IO::Options::FaceNormal;

    mesh.request_face_normals();
    mesh.request_vertex_normals();

    cout << "Reading from file " << argv[1] << "...\n";
    if ( !IO::read_mesh(mesh, argv[1], opt )) {
        cout << "Read failed.\n";
        exit(0);
    }

    cout << "Mesh stats:\n";
    cout << '\t' << mesh.n_vertices() << " vertices.\n";
    cout << '\t' << mesh.n_edges() << " edges.\n";
    cout << '\t' << mesh.n_faces() << " faces.\n";

    simplify(mesh,.05f);

    mesh.update_normals();

    /*cout << "Writing to file "<<"homer-05.off"<<"..\n";
      IO::write_mesh(mesh, "homer-05.off", opt);*/
    mesh.add_property(viewCurvature);
    mesh.add_property(viewCurvatureDerivative);
    mesh.add_property(curvature);

    // Move center of mass to origin
    Vec3f center(0,0,0);
    for (Mesh::ConstVertexIter vIt = mesh.vertices_begin(); vIt != mesh.vertices_end(); ++vIt) center += mesh.point(vIt);
    center /= mesh.n_vertices();
    for (Mesh::VertexIter vIt = mesh.vertices_begin(); vIt != mesh.vertices_end(); ++vIt) mesh.point(vIt) -= center;

    // Fit in the unit sphere
    float maxLength = 0;
    for (Mesh::ConstVertexIter vIt = mesh.vertices_begin(); vIt != mesh.vertices_end(); ++vIt) maxLength = max(maxLength, mesh.point(vIt).length());
    for (Mesh::VertexIter vIt = mesh.vertices_begin(); vIt != mesh.vertices_end(); ++vIt) mesh.point(vIt) /= maxLength;

    computeCurvature(mesh,curvature);

    up = Vec3f(0,1,0);
    pan = Vec3f(0,0,0);

    Vec3f actualCamPos(cameraPos[0]+pan[0],cameraPos[1]+pan[1],cameraPos[2]+pan[2]);
    computeViewCurvature(mesh,actualCamPos,curvature,viewCurvature,viewCurvatureDerivative);

    /* for(Mesh::VertexIter v_it = mesh.vertices_begin(); v_it!= mesh.vertices_end(); ++ v_it){
       std::cout<<mesh.property(viewCurvature, v_it)<<std::endl;
       }*/
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
    glutInitWindowSize(windowWidth, windowHeight);
    glutCreateWindow(argv[0]);

    glutDisplayFunc(display);
    glutMotionFunc(mouseMoved);
    glutMouseFunc(mouse);
    glutReshapeFunc(reshape);
    glutKeyboardFunc(keyboard);

    glutMainLoop();

    return 0;
}
