#include "Object.hpp"
#include "Factory.hpp"
#include "TriangleMeshFactory.hpp"
#include "functions.hpp"

Vec3f normalize(const Vec3f &v)
{
    float mag2 = v.x * v.x + v.y * v.y + v.z * v.z;
    if (mag2 > 0) {
        float invMag = 1 / sqrtf(mag2);
        return Vec3f(v.x * invMag, v.y * invMag, v.z * invMag);
    }
    return v;
}

Vec3f crossProduct(const Vec3f &a, const Vec3f &b)
{
    return Vec3f(
        a.y * b.z - a.z * b.y,
        a.z * b.x - a.x * b.z,
        a.x * b.y - a.y * b.x
    );
}

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

bool trace(
    const Vec3f &orig, const Vec3f &dir,
    const std::vector<Object*> &objects,
    float &tNear, int &index, Vec2f &uv, Object **hitObject, Options options)
{
    *hitObject = nullptr;
    #pragma omp parallel for
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
        if (options.num_objects == 1){
            hitColor = Vec3f(64, 64, 64);
        }
        else{
            hitColor = Vec3f(64, 64, 64) + Vec3f(static_cast<int>(127/(options.num_objects-1))) * (hitObject->pos-1);
        }
    }
    return hitColor;
}

void render(
    const Options &options,
    const std::vector<Object*> &objects)
{
    float scale = tan(deg2rad(options.fov * 0.5))*options.dist_to_screen;
    float imageAspectRatio = options.width / (float)options.height;
    Vec3f orig, tmp;
    float color[3];
    options.cameraToWorld.multVecMatrix(Vec3f(0), orig);
    cimg_library::CImg<float> img(options.width, options.height, 1, 3); 
#pragma omp parallel for private(color,tmp)
    for (int j = 0; j < options.height; ++j) {
        for (int i = 0; i < options.width; ++i) {
            float x = (2 * (i + 0.5) / (float)options.width - 1) * imageAspectRatio * scale;
            float y = (1 - 2 * (j + 0.5) / (float)options.height) * scale;
            Vec3f dir;
            options.cameraToWorld.multDirMatrix(Vec3f(x, y, -1), dir);
            normalize(dir);
            tmp = castRay(orig, dir, objects, options);
            color[0] = tmp[0];
			color[1] = tmp[1];
			color[2] = tmp[2];
            img.draw_point(i, j, color);
        }
    }
    img.save(options.output_file.c_str());
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
    options.backgroundColor = Vec3f(252, 255, 255);
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

std::vector<Object*> GetFigures(std::string FileName, Options &options){

    std::vector<Object*> objects;
    std::ifstream in(FileName);
    //in.open(R"(полный путь к файлу)");
    if (!in)
	{
		std::cout << "Figures file doesn't exist!" << std::endl;
		return objects;
	}

    TriangleMeshFactory* TriangleMesh_Factory = new TriangleMeshFactory;
    float tmp[12];
    options.num_objects = 0;

    for(std::string line; std::getline(in, line); )
    {
        std::istringstream in(line);
        std::string type;
        in >> type;

        if(type == "tetra") {
            options.num_objects++;
            in >> tmp[0] >> tmp[1] >> tmp[2] >> 
            tmp[3] >> tmp[4] >> tmp[5] >> tmp[6] >> 
            tmp[7] >> tmp[8] >> tmp[9] >> tmp[10] 
            >> tmp[11];
            Vec3f verts[4] = {{tmp[0],tmp[1],tmp[2]}, {tmp[3], tmp[4],tmp[5]}, {tmp[6],tmp[7],tmp[8]}, {tmp[9],tmp[10],tmp[11]}};
            objects.push_back((TriangleMesh_Factory->Create(verts, 4)));
            float dist = dotProduct(objects.back()->getMassCenter(), objects.back()->getMassCenter());
            options.distances.insert(std::make_pair(dist, objects.back()));
        }
        else{
            std::cout << "Wrong info in 'figures.txt'" << std::endl;
		    return objects;
        }
    }

    std::map<int, Object*> tmp1;
    std::map<Object*, int> poss;
    int i = 1;
    for(std::map<float, Object*>::const_iterator it = options.distances.begin();
    it != options.distances.end(); ++it){
        tmp1.insert(std::make_pair(i, it->second));
        i++;
    }

    for(std::map<int, Object*>::const_iterator it = tmp1.begin();
    it != tmp1.end(); ++it){
        poss.insert(std::make_pair(it->second, it->first));
    }

    for (std::vector<Object*>::iterator it = objects.begin() ; it != objects.end(); ++it){
        (*it)->pos = poss.find(*it)->second;
    }

    in.close();
    return objects;
}

void autotest1(){ // Front view

    Options options = GetSettings("settings1.txt", "front.bmp");
    std::vector<Object*> objects = GetFigures("figures.txt", options);
    if (objects.size() == 0){
        return;
    }
    render(options, objects);

    for (int i = 0; i < objects.size(); i++){
        delete objects[i];
    }
    
    std::cout << "Test1: OK (Front view)" << std::endl;
}

void autotest2(){ // Side view

    Options options = GetSettings("settings2.txt", "side.bmp");
    std::vector<Object*> objects = GetFigures("figures.txt", options);
    if (objects.size() == 0){
        return;
    }
    render(options, objects);

    for (int i = 0; i < objects.size(); i++){
        delete objects[i];
    }

    std::cout << "Test2: OK (Side view)" << std::endl;
}

void autotest3(){ // Back view

    Options options = GetSettings("settings3.txt", "back.bmp");
    std::vector<Object*> objects = GetFigures("figures.txt", options);
    if (objects.size() == 0){
        return;
    }
    render(options, objects);

    for (int i = 0; i < objects.size(); i++){
        delete objects[i];
    }

    std::cout << "Test3: OK (Back view)" << std::endl;
}