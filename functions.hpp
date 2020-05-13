#define _USE_MATH_DEFINES
#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>
#include <memory>
#include <string>
#include <chrono>
#include <sstream>
#include "omp.h"
#include <algorithm>

Vec3f normalize(const Vec3f &v);

inline
float dotProduct(const Vec3f &a, const Vec3f &b)
{ return a.x * b.x + a.y * b.y + a.z * b.z; }

Vec3f crossProduct(const Vec3f &a, const Vec3f &b);

inline
float clamp(const float &lo, const float &hi, const float &v)
{ return std::max(lo, std::min(hi, v)); }

inline
float deg2rad(const float &deg)
{ return deg * M_PI / 180; }

inline
Vec3f mix(const Vec3f &a, const Vec3f& b, const float &mixValue)
{ return a * (1 - mixValue) + b * mixValue; }

bool rayTriangleIntersect(const Vec3f &v0, const Vec3f &v1, const Vec3f &v2, const Vec3f &orig, const Vec3f &dir, float &tnear, float &u, float &v);

bool trace(const Vec3f &orig, const Vec3f &dir, const std::vector<Object*> &objects, float &tNear, int &index, Vec2f &uv, Object **hitObject, Options options);

Vec3f castRay(const Vec3f &orig, const Vec3f &dir, const std::vector<Object*> &objects, const Options &options);

void render(const Options &options, const std::vector<Object*> &objects);

Options GetSettings(std::string FileIn, std::string FileOut);

std::vector<Object*> GetFigures(std::string FileName);

void autotest1();

void autotest2();

void autotest3();
