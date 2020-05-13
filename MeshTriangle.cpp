#include "Object.hpp"
#include "MeshTriangle.hpp"
#include "functions.hpp"

MeshTriangle::MeshTriangle(const Vec3f *verts, const int &numTris)
    {
        int vertsIndex[12] = {2,1,0,3,2,0,0,1,3,1,2,3}; 
        Vec2f st[4] = {{0, 0}, {1, 0}, {1, 1}, {0, 1}}; 
        int maxIndex = 0;
        for (int i = 0; i < numTris * 3; ++i)
            if (vertsIndex[i] > maxIndex) maxIndex = vertsIndex[i];
        maxIndex += 1;
        vertices = std::unique_ptr<Vec3f[]>(new Vec3f[maxIndex]);
        memcpy(vertices.get(), verts, sizeof(Vec3f) * maxIndex);
        vertexIndex = std::unique_ptr<int[]>(new int[numTris * 3]);
        memcpy(vertexIndex.get(), vertsIndex, sizeof(int) * numTris * 3);
        numTriangles = numTris;
        massCenter[0] = (verts[0][0] + verts[1][0] + verts[2][0] + verts[3][0])/4;
        massCenter[1] = (verts[0][1] + verts[1][1] + verts[2][1] + verts[3][1])/4;
        massCenter[2] = (verts[0][2] + verts[1][2] + verts[2][2] + verts[3][2])/4;
    }

bool MeshTriangle::intersect(const Vec3f &orig, const Vec3f &dir, float &tnear, int &index, Vec2f &uv) const  //BOOST
    {
        bool intersect = false;
        //#pragma omp parallel for
        for (int k = 0; k < numTriangles; ++k) {
            const Vec3f &v0 = vertices[vertexIndex[k * 3]];
            const Vec3f &v1 = vertices[vertexIndex[k * 3 + 1]];
            const Vec3f &v2 = vertices[vertexIndex[k * 3 + 2]];
            float t, u, v;
            if (rayTriangleIntersect(v0, v1, v2, orig, dir, t, u, v) && t < tnear) {
                tnear = t;
                uv.x = u;
                uv.y = v;
                index = k;
                intersect |= true;
            }
        }
        return intersect;
    }

Vec3f MeshTriangle::getMassCenter() const{
        return massCenter;
    }