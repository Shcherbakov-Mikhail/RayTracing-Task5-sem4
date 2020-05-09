#include <cstdio> 
#include <cstdlib> 
#include <memory> 
#include <vector> 
#include <utility> 
#include <cstdint> 
#include <iostream> 
#include <fstream> 
#include <cmath> 
#include "omp.h"
using namespace std;

//const double M_PI = 3.14;
const float kInfinity = std::numeric_limits<float>::max(); 

class Vect3 {
    public:

    double X;
    double Y;
    double Z;

    Vect3(double x = 0.0, double y = 0.0, double z = 0.0) : X(x), Y(y), Z(z) {}

    Vect3 operator +(const Vect3& vec) const{
        return Vect3(X+vec.X,Y+vec.Y,Z+vec.Z);
    }

    Vect3 operator -(const Vect3& vec) const {
        return Vect3(X-vec.X,Y-vec.Y,Z-vec.Z);
    }

    void fill(double x, double y, double z) const {
        Vect3(x,y,z);
    }

    Vect3 operator *(double scalar) const{
        return Vect3(scalar*X,scalar*Y,scalar*Z);
    }

	Vect3 crossProduct(const Vect3& vec) const {
		return Vect3(vec.Y*Z - vec.Z*Y, vec.Z*X - vec.X*Z, vec.X*Y - vec.Y*X);
	}

    double dot(const Vect3& vec) const {
        return X*vec.X+Y*vec.Y+Z*vec.Z;
    }

    double norm() const {
        return sqrt(dot(*this));
    }

    Vect3 normalize() {
        return (*this)*(1/norm());
    }

    void print(){
        cout << X << ", " << Y << ", " << Z << endl;
    }

};

double clamp(const double &lo, const double &hi, const double &v) { 
    return max(lo, min(hi, v)); 
} 

double deg2rad(const float &deg){ 
    return deg * M_PI / 180; 
} 

typedef struct Options{
    int width; 
    int height; 
    double fov; 
    double imageAspectRatio; 
    double maxDepth; 
    Vect3 backgroundColor;
    Vect3 cameraPosition;
}Options;

class Ray {
    public:
    Vect3 orig;
    Vect3 dir;
    Ray() {}
    Ray(const Vect3& a, const Vect3& b): orig(a), dir(b) {}
    Vect3 origin() const { return orig; }
    Vect3 direction() const { return dir; }
    Vect3 point_at_parameter(double t) const{ return orig + dir * t; }
};

class Sphere{
    public:
    double radius;
    Vect3 center;
    Sphere(Vect3 cent, double rad){
        center = cent;
        radius = rad;
    }
    bool intersect(const Vect3 &orig, const Vect3 &dir, double &tnear, int &index, Vect3 &uv) const {
        double t0, t1;
        Vect3 L = center - orig; 
        double tca = L.dot(dir); 
        if (tca < 0){
            return false;
        }
        double d2 = L.dot(L) - tca * tca; 
        if (d2 > radius){
            return false;
        } 
        double thc = sqrt(radius - d2); 
        t0 = tca - thc; // first point
        t1 = tca + thc; // second point
        
        if (t0 > t1){
            swap(t0, t1); 
        }
 
        if (t0 < 0){ 
            t0 = t1; // if t0 is negative, let's use t1 instead 
            if (t0 < 0){
                return false;
            }
        } 
 
        return true;
    }
};

class Screen {
    public:
    Vect3 normal;
    Vect3 up;
    double dist;
    double limit;
    double angle;
    double width;
    double height;
    Screen(double dist = 0.0, double limit = 0.0, double angle = 0.0, double width = 0.0, double height = 0.0){
        this->dist = dist;
        this->limit = limit;
        this->angle = angle;
        this->width = width;
        this->height = height;
    }
    void add_normal(Vect3 norm){
        normal = norm;
    }
    void add_up(Vect3 vup){
        up = vup;
    }
};

bool trace( 
    const Vect3 &orig, const Vect3 &dir, 
    const vector<Sphere> &objects, 
    double &t_nearest, int &index, 
    Vect3 &uv, Sphere hitObject) 
{ 
    hitObject = NULL; 
    for (int k = 0; k < objects.size(); ++k) { 
        double tNearK = kInfinity; 
        int indexK; 
        Vec2f uvK; 
        if (objects[k].intersect(orig, dir, tNearK, indexK, uvK) && (tNearK < t_nearest)) { 
            hitObject = objects[k]; 
            t_nearest = tNearK; 
            index = indexK; 
            uv = uvK; 
        } 
    } 
 
    return (*hitObject != nullptr); 
} 

Vect3 castRay(
        const Vect3 &orig, 
        const Vect3 &dir, 
        const vector<Sphere> &objects, 
        const Options &options,
        bool test = false) 
{
    double t_nearest = kInfinity; // distance to the intersected object
    Sphere *hitObject = nullptr; // pointer to the intersected object 
    Vect3 hitColor = options.backgroundColor; // the color of the intersected point

    Vect3 uv; //???
    int index = 0; //???

    if (trace(orig, dir, objects, t_nearest, index, uv, &hitObject)) { 
        //Vect3 hitPoint = orig + dir * t_nearest;
        //Vec2f st; // st coordinates 
        //hitObject->getSurfaceProperties(hitPoint, dir, index, uv, N, st); 
        hitColor = Vect3(1, 1, 1);

    } 
 
    return hitColor; 
} 

void render(const Options &options, const vector<Sphere> &objects) 
{ 
    Vect3 *framebuffer = new Vect3[options.width * options.height]; 
    Vect3 *pix = framebuffer; 
    double scale = tan(deg2rad(options.fov * 0.5)); 
    double imageAspectRatio = options.width / (float)options.height; 
    Vect3 orig = options.cameraPosition; 

    for (int j = 0; j < options.height; ++j) { 
        for (int i = 0; i < options.width; ++i) {
            double x = (2 * (i + 0.5) / (double)options.width - 1) * imageAspectRatio * scale; 
            double y = (1 - 2 * (j + 0.5) / (double)options.height) * scale; 
            Vect3 dir = Vect3(x, y, -1).normalize(); 
            *(pix++) = castRay(orig, dir, objects, options); // returns black or backgroundColor
        } 
    }
 
    ofstream ofs; 
    ofs.open("./render.bmp"); 
    ofs << "P6\n" << options.width << " " << options.height << "\n255\n"; 
    for (uint32_t i = 0; i < options.height * options.width; ++i) { 
        char r = (char)(255 * clamp(0, 1, framebuffer[i].X)); 
        char g = (char)(255 * clamp(0, 1, framebuffer[i].Y)); 
        char b = (char)(255 * clamp(0, 1, framebuffer[i].Z)); 
        ofs << r << g << b; 
    } 
 
    ofs.close(); 
 
    delete[] framebuffer; 
} 

int main()
{

    Screen scr(1.0, 20.0, 30.0, 640, 480); //dist to screen, distance of view, vertical camera angle, screen params
    scr.add_normal(Vect3(1.0, 0.0, 0.0)); //towards object
    scr.add_up(Vect3(0.0, 0.0, 1.0)); // up direction

    vector<Sphere> objects; 

    Sphere sph1(Vect3(5.0, 0.0, 0.0), 1.0);

    objects.push_back(sph1);

    // Mesh Triangle creation:::
    // Vec3f verts[4] = {{-5,-3,-6}, {5,-3,-6}, {5,-3,-16}, {-5,-3,-16}}; 
    // uint32_t vertIndex[6] = {0, 1, 3, 1, 2, 3}; 
    // Vec2f st[4] = {{0, 0}, {1, 0}, {1, 1}, {0, 1}}; 
    // MeshTriangle *mesh = new MeshTriangle(verts, vertIndex, 2, st); 

    Options options; 
    options.width = 640; 
    options.height = 480; 
    options.fov = 60; 
    options.backgroundColor = Vect3(255, 255, 255); 
    options.maxDepth = 20;
    options.cameraPosition = Vect3(0.0, 0.0, 0.0);

    //render(options, objects); 

    return 0;
}