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
#include "CImg.h"
#include "Vec2f.hpp"
#include "Vec3f.hpp"
#include "Matrix44.hpp"
#include <map>

class Object
{
    public:
    int pos;
    std::string name = "";
    Object();
    virtual ~Object();
    virtual bool intersect(const Vec3f &, const Vec3f &, float &, int &, Vec2f &) const = 0;
    virtual Vec3f getMassCenter() const = 0;
};

struct Options
{
    int width; 
    int height; 
    float fov;
    Vec3f backgroundColor;
    double dist_to_screen;
    Vec3f vup;
    Vec3f camera;
    Vec3f normal;
    Vec3f right;
    Matrix44 cameraToWorld;
    float max_dist;
    Vec3f target;
    std::string output_file;
    int num_objects;
    std::map <float, Object*> distances;
};