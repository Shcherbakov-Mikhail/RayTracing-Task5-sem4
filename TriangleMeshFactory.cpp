#include "Object.hpp"
#include "MeshTriangle.hpp"
#include "Factory.hpp"
#include "TriangleMeshFactory.hpp"

Object* TriangleMeshFactory::Create(const Vec3f *verts, const int &numTris){
        return new MeshTriangle(verts, numTris);
    }
