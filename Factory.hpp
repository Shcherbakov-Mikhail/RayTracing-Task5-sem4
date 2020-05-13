class Factory
{
    public:
	virtual Object* Create(const Vec3f *verts, const int &numTris) = 0;
};
