#include "Vec3f.hpp"

Vec3f::Vec3f() : x(0), y(0), z(0) {}
Vec3f::Vec3f(float xx) : x(xx), y(xx), z(xx) {}
Vec3f::Vec3f(float xx, float yy, float zz) : x(xx), y(yy), z(zz) {}
Vec3f Vec3f::operator * (const float &r) const { return Vec3f(x * r, y * r, z * r); }
Vec3f Vec3f::operator * (const Vec3f &v) const { return Vec3f(x * v.x, y * v.y, z * v.z); }
Vec3f Vec3f::operator - (const Vec3f &v) const { return Vec3f(x - v.x, y - v.y, z - v.z); }
Vec3f Vec3f::operator + (const Vec3f &v) const { return Vec3f(x + v.x, y + v.y, z + v.z); }
Vec3f Vec3f::operator - () const { return Vec3f(-x, -y, -z); }
Vec3f& Vec3f::operator += (const Vec3f &v) { x += v.x, y += v.y, z += v.z; return *this; }
const float& Vec3f::operator [] (int i) const { return (&x)[i]; }
float& Vec3f::operator [] (int i) { return (&x)[i]; }