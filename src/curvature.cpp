#include <Eigen/Core>
#include <Eigen/Eigenvalues>
#include <iostream>
#include "curvature.h"
using namespace OpenMesh;
using namespace Eigen;
using namespace std;

void computeCurvature(Mesh &mesh, OpenMesh::VPropHandleT<CurvatureInfo> &curvature) {
 double tol = 1e-8;
        Matrix3d I = Matrix3d::Identity();
	for (Mesh::VertexIter it = mesh.vertices_begin(); it != mesh.vertices_end(); ++it) {
		// WRITE CODE HERE TO COMPUTE THE CURVATURE AT THE CURRENT VERTEX ----------------------------------------------
		Vec3f normal = mesh.normal(it.handle());
		Vector3d N(normal[0],normal[1],normal[2]); // example of converting to Eigen's vector class for easier math
                Vec3f vi = mesh.point(it.handle());
                Vector3d v_i(vi[0],vi[1],vi[2]);
                Matrix3d NN_I = N*N.transpose() - I;
                double w = 0;
                //Declare matrix M
                Matrix3d M = Matrix3d::Zero();
                for(Mesh::VertexOHalfedgeIter voh_it = mesh.voh_iter(it); voh_it; ++voh_it){//iterate all outgoing halfedges
                       Mesh::HalfedgeHandle heh = voh_it.handle();
                       Mesh::HalfedgeHandle opposite_heh = mesh.opposite_halfedge_handle(heh);
                       //vj
                       Mesh::VertexHandle vh = mesh.to_vertex_handle(heh); 
                       //vi-vj
                       Vec3f vj = mesh.point(vh);
                       Vector3d v_j(vj[0], vj[1], vi[2]);
                       Vector3d v_ij = v_j - v_i; 
                       //Tij
                       Vector3d T_ij = NN_I*v_ij;
                       T_ij.normalize();
                       //kij
                       double k_ij = N.dot(v_ij);
                       k_ij = 2*k_ij/v_ij.squaredNorm();
                       //wij
                       double w_ij = mesh.calc_sector_area(heh) +
				     mesh.calc_sector_area(opposite_heh);
                       w += w_ij;
                       //M+=
 		       M += w_ij * k_ij * T_ij * T_ij.transpose();
                 }
                 //normalize weight
                 M = M / w;
                //engivector & engivalue
                 EigenSolver<Matrix3d> es(M, true);
                 Vector3d ev1 = es.pseudoEigenvectors().block(0, 0, 3, 1);
		 Vector3d ev2 = es.pseudoEigenvectors().block(0, 1, 3, 1);
                 Vector3d ev3 = es.pseudoEigenvectors().block(0, 2, 3, 1);
                //remove the one parallel to N
                 Vector3d T1, T2;
                 double t1, t2;
                //eigenvalue and curvature
                  double eig1 = real(es.eigenvalues()(0));
                  double eig2 = real(es.eigenvalues()(1));
                  double eig3 = real(es.eigenvalues()(2)); 
                  //printf("%f %f %f\n", eig1, eig2, eig3);
                  if(abs(eig1)<tol){
 			T1 = ev2.normalized();
 			T2 = ev3.normalized();
                        t1 = eig1;
                        t2 = eig2;
		  }
		  else if(abs(eig2)<tol){
			T1 = ev1.normalized();
			T2 = ev3.normalized();
			t1 = eig1;
			t2 = eig3;
			}
			else if(abs(eig3)<tol){
				T1 = ev1.normalized();
				T2 = ev2.normalized();
				t1 = eig1;
				t2 = eig2;
			     }
                             //else
                               // printf("FAIL\n");                
		// In the end you need to fill in this struct
		double k1 = 3*t1 - t2;
                double k2 = 3*t2 - t1;
		CurvatureInfo info;
		info.curvatures[0] = k1;
		info.curvatures[1] = k2;
		info.directions[0] = Vec3f(T1[0], T1[1], T1[2]);
		info.directions[1] = Vec3f(T2[0], T2[1], T2[2]);

		mesh.property(curvature,it) = info;
       }
}
void computeViewCurvature(Mesh &mesh, OpenMesh::Vec3f camPos, OpenMesh::VPropHandleT<CurvatureInfo> &curvature, OpenMesh::VPropHandleT<double> &viewCurvature, OpenMesh::FPropHandleT<OpenMesh::Vec3f> &viewCurvatureDerivative) {
	// WRITE CODE HERE TO COMPUTE CURVATURE IN THE VIEW PROJECTION PROJECTED ON THE TANGENT PLANE ------------------
	// Compute vector to viewer and project onto tangent plane, then use components in principal directions to find curvature
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
