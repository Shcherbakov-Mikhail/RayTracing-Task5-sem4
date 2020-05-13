class MeshTriangle : public Object
{
    public:
    std::unique_ptr<Vec3f[]> vertices;
    int numTriangles;
    std::unique_ptr<int[]> vertexIndex;
    Vec3f massCenter;

    MeshTriangle(const Vec3f *verts, const int &numTris);
    bool intersect(const Vec3f &orig, const Vec3f &dir, float &tnear, int &index, Vec2f &uv) const;
    Vec3f getMassCenter() const;
};
