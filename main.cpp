#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <chrono>
#include <sstream>
#include "omp.h"

class Vec3f {
    public:
    float x, y, z;
    Vec3f() : x(0), y(0), z(0) {}
    Vec3f(float xx) : x(xx), y(xx), z(xx) {}
    Vec3f(float xx, float yy, float zz) : x(xx), y(yy), z(zz) {}
    Vec3f operator * (const float &r) const { return Vec3f(x * r, y * r, z * r); }
    Vec3f operator * (const Vec3f &v) const { return Vec3f(x * v.x, y * v.y, z * v.z); }
    Vec3f operator - (const Vec3f &v) const { return Vec3f(x - v.x, y - v.y, z - v.z); }
    Vec3f operator + (const Vec3f &v) const { return Vec3f(x + v.x, y + v.y, z + v.z); }
    Vec3f operator - () const { return Vec3f(-x, -y, -z); }
    Vec3f& operator += (const Vec3f &v) { x += v.x, y += v.y, z += v.z; return *this; }
    const float& operator [] (int i) const { return (&x)[i]; }
    float& operator [] (int i) { return (&x)[i]; }
};

class Vec2f
{
    public:
    float x, y;
    Vec2f() : x(0), y(0) {}
    Vec2f(float xx) : x(xx), y(xx) {}
    Vec2f(float xx, float yy) : x(xx), y(yy) {}
    Vec2f operator * (const float &r) const { return Vec2f(x * r, y * r); }
    Vec2f operator + (const Vec2f &v) const { return Vec2f(x + v.x, y + v.y); }
};

Vec3f normalize(const Vec3f &v)
{
    float mag2 = v.x * v.x + v.y * v.y + v.z * v.z;
    if (mag2 > 0) {
        float invMag = 1 / sqrtf(mag2);
        return Vec3f(v.x * invMag, v.y * invMag, v.z * invMag);
    }
    return v;
}

inline
float dotProduct(const Vec3f &a, const Vec3f &b)
{ return a.x * b.x + a.y * b.y + a.z * b.z; }

Vec3f crossProduct(const Vec3f &a, const Vec3f &b)
{
    return Vec3f(
        a.y * b.z - a.z * b.y,
        a.z * b.x - a.x * b.z,
        a.x * b.y - a.y * b.x
    );
}

class Matrix44
{
    public:
    float x[4][4] = {{1,0,0,0},{0,1,0,0},{0,0,1,0},{0,0,0,1}};

    Matrix44() {}

    Matrix44 (float a, float b, float c, float d, float e, float f, float g, float h,
              float i, float j, float k, float l, float m, float n, float o, float p)
    {
        x[0][0] = a;
        x[0][1] = b;
        x[0][2] = c;
        x[0][3] = d;
        x[1][0] = e;
        x[1][1] = f;
        x[1][2] = g;
        x[1][3] = h;
        x[2][0] = i;
        x[2][1] = j;
        x[2][2] = k;
        x[2][3] = l;
        x[3][0] = m;
        x[3][1] = n;
        x[3][2] = o;
        x[3][3] = p;
    }
    
    const float* operator [] (int i) const { return x[i]; }
    float* operator [] (int i) { return x[i]; }

    static void multiply(const Matrix44 &a, const Matrix44& b, Matrix44 &c) //BOOST
    {
        #pragma omp parallel for
        for (int i = 0; i < 4; ++i) {
            for (int j = 0; j < 4; ++j) {
                c[i][j] = a[i][0] * b[0][j] + a[i][1] * b[1][j] +
                    a[i][2] * b[2][j] + a[i][3] * b[3][j];
            }
        }
    }
    
    Matrix44 transposed() const //BOOST
    {
        #pragma omp parallel for private(t, x)
        Matrix44 t;
        for (int i = 0; i < 4; ++i) {
            for (int j = 0; j < 4; ++j) {
                t[i][j] = x[j][i];
            }
        }
        return t;
    }

    Matrix44& transpose ()
    {
        Matrix44 tmp (x[0][0],
                      x[1][0],
                      x[2][0],
                      x[3][0],
                      x[0][1],
                      x[1][1],
                      x[2][1],
                      x[3][1],
                      x[0][2],
                      x[1][2],
                      x[2][2],
                      x[3][2],
                      x[0][3],
                      x[1][3],
                      x[2][3],
                      x[3][3]);
        *this = tmp;
        return *this;
    }

    
    
    void multVecMatrix(const Vec3f &src, Vec3f &dst) const
    {
        float a, b, c, w;
        
        a = src[0] * x[0][0] + src[1] * x[1][0] + src[2] * x[2][0] + x[3][0];
        b = src[0] * x[0][1] + src[1] * x[1][1] + src[2] * x[2][1] + x[3][1];
        c = src[0] * x[0][2] + src[1] * x[1][2] + src[2] * x[2][2] + x[3][2];
        w = src[0] * x[0][3] + src[1] * x[1][3] + src[2] * x[2][3] + x[3][3];
        
        dst.x = a / w;
        dst.y = b / w;
        dst.z = c / w;
    }
    
    void multDirMatrix(const Vec3f &src, Vec3f &dst) const
    {
        float a, b, c;
        
        a = src[0] * x[0][0] + src[1] * x[1][0] + src[2] * x[2][0];
        b = src[0] * x[0][1] + src[1] * x[1][1] + src[2] * x[2][1];
        c = src[0] * x[0][2] + src[1] * x[1][2] + src[2] * x[2][2];
        
        dst.x = a;
        dst.y = b;
        dst.z = c;
    }

    Matrix44 inverse() const //BOOST
    {
        int i, j, k;
        Matrix44 s;
        Matrix44 t (*this);
        
        for (i = 0; i < 3 ; i++) {
            int pivot = i;
            
            float pivotsize = t[i][i];
            
            if (pivotsize < 0)
                pivotsize = -pivotsize;
                
                for (j = i + 1; j < 4; j++) {
                    float tmp = t[j][i];
                    
                    if (tmp < 0)
                        tmp = -tmp;
                        
                        if (tmp > pivotsize) {
                            pivot = j;
                            pivotsize = tmp;
                        }
                }
            
            if (pivotsize == 0) {
                return Matrix44();
            }
            
            if (pivot != i) {
                for (j = 0; j < 4; j++) {
                    float tmp;
                    
                    tmp = t[i][j];
                    t[i][j] = t[pivot][j];
                    t[pivot][j] = tmp;
                    
                    tmp = s[i][j];
                    s[i][j] = s[pivot][j];
                    s[pivot][j] = tmp;
                }
            }
            
            for (j = i + 1; j < 4; j++) {
                float f = t[j][i] / t[i][i];
                
                for (k = 0; k < 4; k++) {
                    t[j][k] -= f * t[i][k];
                    s[j][k] -= f * s[i][k];
                }
            }
        }
        
        for (i = 3; i >= 0; --i) {
            float f;
            
            if ((f = t[i][i]) == 0) {
                return Matrix44();
            }
            
            for (j = 0; j < 4; j++) {
                t[i][j] /= f;
                s[i][j] /= f;
            }
            
            for (j = 0; j < i; j++) {
                f = t[j][i];
                
                #pragma omp parallel for private(t, s)
                for (k = 0; k < 4; k++) {
                    t[j][k] -= f * t[i][k];
                    s[j][k] -= f * s[i][k];
                }
            }
        }
        return s;
    }

    const Matrix44& invert()
    {
        *this = inverse();
        return *this;
    }
};

inline
float clamp(const float &lo, const float &hi, const float &v)
{ return std::max(lo, std::min(hi, v)); }

inline
float deg2rad(const float &deg)
{ return deg * M_PI / 180; }

inline
Vec3f mix(const Vec3f &a, const Vec3f& b, const float &mixValue)
{ return a * (1 - mixValue) + b * mixValue; }

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
};

class Object
{
    public:
    std::string name = "";
    Object(){}
    virtual ~Object() {}
    virtual bool intersect(const Vec3f &, const Vec3f &, float &, int &, Vec2f &) const = 0;
    virtual Vec3f getMassCenter() const = 0;
};

bool rayTriangleIntersect(
    const Vec3f &v0, const Vec3f &v1, const Vec3f &v2,
    const Vec3f &orig, const Vec3f &dir,
    float &tnear, float &u, float &v)
{
    Vec3f edge1 = v1 - v0;
    Vec3f edge2 = v2 - v0;
    Vec3f pvec = crossProduct(dir, edge2);
    float det = dotProduct(edge1, pvec);
    if (det == 0 || det < 0) return false;

    Vec3f tvec = orig - v0;
    u = dotProduct(tvec, pvec);
    if (u < 0 || u > det) return false;

    Vec3f qvec = crossProduct(tvec, edge1);
    v = dotProduct(dir, qvec);
    if (v < 0 || u + v > det) return false;

    float invDet = 1 / det;
    
    tnear = dotProduct(edge2, qvec) * invDet;
    u *= invDet;
    v *= invDet;

    return true;
}

class MeshTriangle : public Object
{
    public:
    std::unique_ptr<Vec3f[]> vertices;
    int numTriangles;
    std::unique_ptr<int[]> vertexIndex;
    Vec3f massCenter;

    MeshTriangle(
        const Vec3f *verts,
        const int &numTris)
    {
        int vertsIndex[12] = {2,1,0,3,2,0,0,1,3,1,2,3}; 
        Vec2f st[4] = {{0, 0}, {1, 0}, {1, 1}, {0, 1}}; 
        int maxIndex = 0;
        for (int i = 0; i < numTris * 3; ++i)
            if (vertsIndex[i] > maxIndex) maxIndex = vertsIndex[i];
        maxIndex += 1;
        vertices = std::unique_ptr<Vec3f[]>(new Vec3f[maxIndex]);
        memcpy(vertices.get(), verts, sizeof(Vec3f) * maxIndex);
        vertexIndex = std::unique_ptr<int[]>(new int[numTris * 3]);
        memcpy(vertexIndex.get(), vertsIndex, sizeof(int) * numTris * 3);
        numTriangles = numTris;
        massCenter[0] = (verts[0][0] + verts[1][0] + verts[2][0] + verts[3][0])/4;
        massCenter[1] = (verts[0][1] + verts[1][1] + verts[2][1] + verts[3][1])/4;
        massCenter[2] = (verts[0][2] + verts[1][2] + verts[2][2] + verts[3][2])/4;
    }

    bool intersect(const Vec3f &orig, const Vec3f &dir, float &tnear, int &index, Vec2f &uv) const  //BOOST
    {
        bool intersect = false;
        #pragma omp parallel for private(vertices, vertexIndex)
        for (int k = 0; k < numTriangles; ++k) {
            const Vec3f &v0 = vertices[vertexIndex[k * 3]];
            const Vec3f &v1 = vertices[vertexIndex[k * 3 + 1]];
            const Vec3f &v2 = vertices[vertexIndex[k * 3 + 2]];
            float t, u, v;
            if (rayTriangleIntersect(v0, v1, v2, orig, dir, t, u, v) && t < tnear) {
                tnear = t;
                uv.x = u;
                uv.y = v;
                index = k;
                intersect |= true;
            }
        }
        return intersect;
    }

    Vec3f getMassCenter() const{
        return massCenter;
    }
};

class Factory
{
public:
	virtual Object* Create(const Vec3f *verts, const int &numTris) = 0;
	virtual ~Factory() {}
};

class TriangleMeshFactory : public Factory
{
    public:
	Object* Create(const Vec3f *verts, const int &numTris){
        return new MeshTriangle(verts, numTris);
    }
};

bool trace(
    const Vec3f &orig, const Vec3f &dir,
    const std::vector<Object*> &objects,
    float &tNear, int &index, Vec2f &uv, Object **hitObject, Options options)
{
    *hitObject = nullptr;
    #pragma omp parallel for private(objects) //BOOST
    for (int k = 0; k < objects.size(); ++k) {
        float tNearK = options.max_dist;
        int indexK;
        Vec2f uvK;
        if (objects[k]->intersect(orig, dir, tNearK, indexK, uvK) && tNearK < tNear) {
            *hitObject = objects[k];
            tNear = tNearK;
            index = indexK;
            uv = uvK;
        }
    }
    return (*hitObject != nullptr);
}

Vec3f castRay(
    const Vec3f &orig, const Vec3f &dir,
    const std::vector<Object*> &objects,
    const Options &options)
{
    Vec3f hitColor = options.backgroundColor;
    float tnear = options.max_dist;
    Vec2f uv;
    int index = 0;
    Object *hitObject = nullptr;
    if (trace(orig, dir, objects, tnear, index, uv, &hitObject, options)){
        hitColor = Vec3f(64, 64, 64) + 10*dotProduct(hitObject->getMassCenter(), hitObject->getMassCenter());
        if (hitColor.x > 191){
            hitColor = Vec3f(191, 191, 191);
        }
    }
    return hitColor;
}

void render(
    const Options &options,
    const std::vector<Object*> &objects)
{
    Vec3f *framebuffer = new Vec3f[options.width * options.height];
    Vec3f *pix = framebuffer;
    float scale = tan(deg2rad(options.fov * 0.5))*options.dist_to_screen;
    float imageAspectRatio = options.width / (float)options.height;
    Vec3f orig;
    options.cameraToWorld.multVecMatrix(Vec3f(0), orig);
    #pragma opm parallel for private(pix, x, y)
    for (int j = 0; j < options.height; ++j) {
        for (int i = 0; i < options.width; ++i) {
            float x = (2 * (i + 0.5) / (float)options.width - 1) * imageAspectRatio * scale;
            float y = (1 - 2 * (j + 0.5) / (float)options.height) * scale;
            Vec3f dir;
            options.cameraToWorld.multDirMatrix(Vec3f(x, y, -1), dir);
            normalize(dir);
            *(pix++) = castRay(orig, dir, objects, options);
        }
    }

    std::ofstream ofs;
    ofs.open(options.output_file, std::ios::out | std::ios::binary);
    ofs << "P6\n" << options.width << " " << options.height << "\n255\n";
    for (int i = 0; i < options.height * options.width; ++i) {
        char r = (char)(framebuffer[i].x);
        char g = (char)(framebuffer[i].y);
        char b = (char)(framebuffer[i].z);
        ofs << r << g << b;
    }

    ofs.close();
    delete[] framebuffer;
}

Options GetSettings(std::string FileIn, std::string FileOut){

    Options options;

    std::ifstream in(FileIn);
    if (!in)
	{
		std::cout << "Settings file doesn't exist" << std::endl;
		return options;
	}

    for(std::string line; std::getline(in, line); )
    {
        std::istringstream in(line);
        std::string option;
        double x, y, z, val;
        in >> option;

        if(option == "cam") {
            in >> x >> y >> z;
            options.camera = Vec3f(x, y, z);
        }
        else if(option == "target") {
            in >> x >> y >> z;
            options.target = Vec3f(x, y, z);
        }
        else if(option == "up") {
            in >> x >> y >> z;
            options.vup = normalize(Vec3f(x, y, z));
        }
        else if(option == "distance") {
            in >> val;
            options.dist_to_screen = val;
        }
        else if(option == "limit") {
            in >> val;
            options.max_dist = val;
        }
        else if(option == "alpha") {
            in >> val;
            options.fov = val;
        }
        else if(option == "width") {
            in >> val;
            options.width = val;
        }
        else if(option == "height") {
            in >> val;
            options.height = val;
        }
        else{
            std::cout << "Smth wrong in 'settings.txt'" << std::endl;
        }
    }

    options.output_file = FileOut;
    options.backgroundColor = Vec3f(255, 255, 255);
    options.normal = normalize((options.camera - options.target));
    options.right = normalize(crossProduct(options.normal, options.vup));
    options.vup = normalize(crossProduct(options.right,options.normal));
    options.cameraToWorld = Matrix44(options.right[0], options.right[1], options.right[2], 0,
                                     options.vup[0], options.vup[1], options.vup[2], 0, 
                                     options.normal[0], options.normal[1], options.normal[2], 0, 
                                     options.camera[0], options.camera[1], options.camera[2], 1); 
    in.close();
    return options;
}

std::vector<Object*> GetFigures(std::string FileName){

    std::vector<Object*> objects;
    std::ifstream in(FileName);
    //in.open(R"(полный путь к файлу)");
    if (!in)
	{
		std::cout << "Figures file doesn't exist" << std::endl;
		return objects;
	}

    TriangleMeshFactory* TriangleMesh_Factory = new TriangleMeshFactory;
    float tmp[12];

    for(std::string line; std::getline(in, line); )
    {
        std::istringstream in(line);
        std::string type;
        in >> type;

        if(type == "tetra") {
            in >> tmp[0] >> tmp[1] >> tmp[2] >> 
            tmp[3] >> tmp[4] >> tmp[5] >> tmp[6] >> 
            tmp[7] >> tmp[8] >> tmp[9] >> tmp[10] 
            >> tmp[11];
            Vec3f verts[4] = {{tmp[0],tmp[1],tmp[2]}, {tmp[3], tmp[4],tmp[5]}, {tmp[6],tmp[7],tmp[8]}, {tmp[9],tmp[10],tmp[11]}};
            objects.push_back((TriangleMesh_Factory->Create(verts, 4)));
        }
        else{
            std::cout << "Wrong info in 'figures.txt'" << std::endl;
		    return objects;
        }
    }

    in.close();
    return objects;
}

void autotest1(){ // Front view

    Options options = GetSettings("settings1.txt", "front.bmp");
    std::vector<Object*> objects = GetFigures("figures.txt");

    render(options, objects);

    for (int i = 0; i < objects.size(); i++){
        delete objects[i];
    }
    
    std::cout << "Test1: OK (Front view)" << std::endl;
}

void autotest2(){ // Side view

    Options options = GetSettings("settings2.txt", "side.bmp");
    std::vector<Object*> objects = GetFigures("figures.txt");

    render(options, objects);

    for (int i = 0; i < objects.size(); i++){
        delete objects[i];
    }

    std::cout << "Test2: OK (Side view)" << std::endl;
}

void autotest3(){ // Back view

    Options options = GetSettings("settings3.txt", "back.bmp");
    std::vector<Object*> objects = GetFigures("figures.txt");

    render(options, objects);

    for (int i = 0; i < objects.size(); i++){
        delete objects[i];
    }

    std::cout << "Test3: OK (Back view)" << std::endl;
}

int main()
{
    auto start = std::chrono::system_clock::now();

    autotest1();
    autotest2();
    autotest3();

    auto end = std::chrono::system_clock::now(); 
    int elapsed_ms = static_cast<int>(std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count()); 
    std::cout << "Program runtime is " << elapsed_ms << " ms\n";

    return 0;
}