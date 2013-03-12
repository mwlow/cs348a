#ifndef BAZIER_H
#define BAZIER_H

#include "mesh_definitions.h"

void cubicBezier(OpenMesh::Vec3f p[4], std::vector<OpenMesh::Vec3f>&);

void bezier(OpenMesh::Vec3f p[3], std::vector<OpenMesh::Vec3f>&);

void holy(OpenMesh::Vec3f p[3], std::vector<OpenMesh::Vec3f>&);
#endif
