class TriangleMeshFactory : public Factory
{
    public:
	Object* Create(const Vec3f *verts, const int &numTris);
};
