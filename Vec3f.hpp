class Vec3f 
{
    public:
    float x, y, z;
    Vec3f();
    Vec3f(float xx);
    Vec3f(float xx, float yy, float zz);
    Vec3f operator * (const float &r) const;
    Vec3f operator * (const Vec3f &v) const;
    Vec3f operator - (const Vec3f &v) const;
    Vec3f operator + (const Vec3f &v) const;
    Vec3f operator - () const;
    Vec3f& operator += (const Vec3f &v);
    const float& operator [] (int i) const;
    float& operator [] (int i);
};