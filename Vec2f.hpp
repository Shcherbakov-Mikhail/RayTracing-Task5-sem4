class Vec2f
{
    public:
    float x, y;
    Vec2f();
    Vec2f(float xx);
    Vec2f(float xx, float yy);
    Vec2f operator * (const float &r) const;
    Vec2f operator + (const Vec2f &v) const;
};
