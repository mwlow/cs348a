#include <Eigen/Core>
#include <Eigen/Eigenvalues>
#include <iostream>
#include "curvature.h"
using namespace OpenMesh;
using namespace Eigen;
using namespace std;

void computeCurvature(Mesh &mesh, OpenMesh::VPropHandleT<CurvatureInfo> &curvature) {
    for (Mesh::VertexIter it = mesh.vertices_begin(); it != mesh.vertices_end(); ++it) {
        // WRITE CODE HERE TO COMPUTE THE CURVATURE AT THE CURRENT VERTEX ----------------------------------------------
        Vec3f normal = mesh.normal(it.handle());
        Vector3d N(normal[0],normal[1],normal[2]); // example of converting to Eigen's vector class for easier math

        std::vector<double> areas;
        std::vector<Matrix3d> T;

        Vec3f vertex = mesh.point(it.handle());
        Vector3d v_i(vertex[0], vertex[1], vertex[2]);

        for (Mesh::VertexVertexIter vv_it = mesh.vv_iter(it.handle()); vv_it; ++vv_it) {
            Vec3f neighbor = mesh.point(vv_it.handle());
            Vector3d v_j(neighbor[0], neighbor[1], neighbor[2]);
            Vector3d T_ij = ((Matrix3d::Identity() - N * N.transpose()) * (v_i - v_j)).normalized();
            double K_ij = 2 / (v_j - v_i).squaredNorm() * N.transpose() * (v_j - v_i);
            T.push_back(K_ij * T_ij * T_ij.transpose());

            double area = 0.0;
            for (Mesh::VertexFaceIter vf_it = mesh.vf_iter(it.handle()); vf_it; ++vf_it) {
                Mesh::FaceVertexIter fv_it = mesh.fv_iter(vf_it.handle());
                Vec3f p0 = mesh.point(fv_it.handle());
                Vec3f p1 = mesh.point((++fv_it).handle());
                Vec3f p2 = mesh.point((++fv_it).handle());
                if (p0 == neighbor || p1 == neighbor || p2 == neighbor)
                    area += ((p2 - p0) % (p1 - p0)).length() / 2.0;
            }
            areas.push_back(area);
        }

        double sum = std::accumulate(areas.begin(), areas.end(), 0.0);

        Matrix3d M = Matrix3d::Zero();
        for (int i = 0; i < areas.size(); ++i) {
            M += areas[i] / sum * T[i];
        }

        EigenSolver<Matrix3d> solver(M);
        Vector3d v0 = solver.pseudoEigenvectors().block(0,0,3,1);
        Vector3d v1 = solver.pseudoEigenvectors().block(0,1,3,1);
        Vector3d v2 = solver.pseudoEigenvectors().block(0,2,3,1);

        v0.normalize();
        v1.normalize();
        v2.normalize();

        double e0 = real(solver.eigenvalues()(0));
        double e1 = real(solver.eigenvalues()(1));
        double e2 = real(solver.eigenvalues()(2));
        double m1 = (v0.cross(N).norm() < 0.1) ? e2 : e0;
        double m2 = (v1.cross(N).norm() < 0.1) ? e2 : e1;

        CurvatureInfo info;
        info.curvatures[0] = m2 - 3.0*m1;
        info.curvatures[1] = m1 - 3.0*m2;
        info.directions[0] = (v0.cross(N).norm() < 0.1) ? Vec3f(v2[0], v2[1], v2[2]) : Vec3f(v0[0], v0[1], v0[2]);
        info.directions[1] = (v1.cross(N).norm() < 0.1) ? Vec3f(v2[0], v2[1], v2[2]) : Vec3f(v1[0], v1[1], v1[2]);

        mesh.property(curvature,it) = info;
    }
}
void computeViewCurvature(Mesh &mesh, OpenMesh::Vec3f camPos, OpenMesh::VPropHandleT<CurvatureInfo> &curvature, OpenMesh::VPropHandleT<double> &viewCurvature, OpenMesh::FPropHandleT<OpenMesh::Vec3f> &viewCurvatureDerivative) {
    // WRITE CODE HERE TO COMPUTE CURVATURE IN THE VIEW PROJECTION PROJECTED ON THE TANGENT PLANE ------------------
    // Compute vector to viewer and project onto tangent plane, then use components in principal directions to find curvature
    for(Mesh::VertexIter v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); ++v_it){
        Vec3f p = mesh.point(v_it.handle());
        Vec3f n = mesh.normal(v_it.handle());
        Vec3f v = camPos - p;
        Vector3d _v(v[0], v[1], v[2]);
        CurvatureInfo info = mesh.property(curvature, v_it);
        Vec3f T1 = info.directions[0];
        Vec3f T2 = info.directions[1];
        Vector3d _T1(T1[0], T1[1], T1[2]);
        Vector3d _T2(T2[0], T2[1], T2[2]);
        Vector3d N(n[0], n[1], n[2]);
        Vector3d w = _v - _v.dot(N)*N;
        double cos_theta = w.dot(_T1)/w.norm();
        double kw = info.curvatures[0] * cos_theta * cos_theta + info.curvatures[1] * (1.0 - cos_theta * cos_theta);
        mesh.property(viewCurvature, v_it.handle()) = kw;
    }
    // -------------------------------------------------------------------------------------------------------------

    // We'll use the finite elements piecewise hat method to find per-face gradients of the view curvature
    // CS 348a doesn't cover how to differentiate functions on a mesh (Take CS 468! Spring 2013!) so we provide code here

    for (Mesh::FaceIter it = mesh.faces_begin(); it != mesh.faces_end(); ++it) {
        double c[3];
        Vec3f p[3];

        Mesh::ConstFaceVertexIter fvIt = mesh.cfv_iter(it);
        for (int i = 0; i < 3; i++) {
            p[i] = mesh.point(fvIt.handle());
            c[i] = mesh.property(viewCurvature,fvIt.handle());
            ++fvIt;
        }

        Vec3f N = mesh.normal(it.handle());
        double area = mesh.calc_sector_area(mesh.halfedge_handle(it.handle()));

        mesh.property(viewCurvatureDerivative,it) = (N%(p[0]-p[2]))*(c[1]-c[0])/(2*area) + (N%(p[1]-p[0]))*(c[2]-c[0])/(2*area);
    }
}
