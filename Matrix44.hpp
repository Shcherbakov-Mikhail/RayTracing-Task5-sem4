class Matrix44
{
    public:
    Matrix44();
    Matrix44 (float a, float b, float c, float d, float e, float f, float g, float h, float i, float j, float k, float l, float m, float n, float o, float p);
    const float* operator [] (int i) const;
    float* operator [] (int i);
    static void multiply(const Matrix44 &a, const Matrix44& b, Matrix44 &c);
    Matrix44 transposed() const;
    Matrix44& transpose ();
    void multVecMatrix(const Vec3f& src, Vec3f& dst) const;
    void multDirMatrix(const Vec3f& src, Vec3f& dst) const;
    Matrix44 inverse() const;
    const Matrix44& invert();
};
