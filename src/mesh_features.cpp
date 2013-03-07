#include "mesh_features.h"
using namespace OpenMesh;

bool isSilhouette(Mesh &mesh, const Mesh::EdgeHandle &e, Vec3f cameraPos)  {
    // CHECK IF e IS A SILHOUETTE HERE -----------------------------------------------------------------------------
    Mesh::HalfedgeHandle h0 = mesh.halfedge_handle(e, 0);
    Mesh::HalfedgeHandle h1 = mesh.halfedge_handle(e, 1);

    Mesh::FaceHandle f0 = mesh.face_handle(h0);
    Mesh::FaceHandle f1 = mesh.face_handle(h1);

    Vec3f n0 = mesh.calc_face_normal(f0);
    Vec3f n1 = mesh.calc_face_normal(f1);

    Vec3f center0(0.0, 0.0, 0.0);
    Vec3f center1(0.0, 0.0, 0.0);

    for (Mesh::FaceVertexIter fv_it = mesh.fv_iter(f0); fv_it; ++fv_it) {
        center0 += mesh.point(fv_it.handle());
        center0 += mesh.point((++fv_it).handle());
        center0 += mesh.point((++fv_it).handle());
        center0 /= 3.0;
    }

    for (Mesh::FaceVertexIter fv_it = mesh.fv_iter(f1); fv_it; ++fv_it) {
        center1 += mesh.point(fv_it.handle());
        center1 += mesh.point((++fv_it).handle());
        center1 += mesh.point((++fv_it).handle());
        center1 /= 3.0;
    }

    Vec3f dir0 = cameraPos - center0;
    Vec3f dir1 = cameraPos - center1;

    return (((dir0|n0) > 0) ^ ((dir1|n1) > 0));

    // -------------------------------------------------------------------------------------------------------------
}

bool isSharpEdge(Mesh &mesh, const Mesh::EdgeHandle &e) {
    // CHECK IF e IS SHARP HERE ------------------------------------------------------------------------------------

    Mesh::HalfedgeHandle h0 = mesh.halfedge_handle(e, 0);
    Mesh::HalfedgeHandle h1 = mesh.halfedge_handle(e, 1);

    Mesh::FaceHandle f0 = mesh.face_handle(h0);
    Mesh::FaceHandle f1 = mesh.face_handle(h1);

    Vec3f n0 = mesh.calc_face_normal(f0);
    Vec3f n1 = mesh.calc_face_normal(f1);

    return (n0 | n1) < 0.5;

    // -------------------------------------------------------------------------------------------------------------
}

bool isFeatureEdge(Mesh &mesh, const Mesh::EdgeHandle &e, Vec3f cameraPos) {
    return mesh.is_boundary(e) || isSilhouette(mesh,e, cameraPos) || isSharpEdge(mesh,e);
}

