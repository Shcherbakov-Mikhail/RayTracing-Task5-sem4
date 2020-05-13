cmake_minimum_required(VERSION 3.15)


project(RayTracer)

set(CMAKE_CXX_STANDARD 14)



add_executable(RayTracer main.cpp Vec3f.cpp Vec3f.hpp Vec2f.cpp Vec2f.hpp functions.cpp functions.hpp TriangleMeshFactory.cpp TriangleMeshFactory.hpp Factory.hpp MeshTriangle.cpp MeshTriangle.hpp Object.cpp Object.hpp Matrix44.hpp Matrix44.cpp CImg.h)