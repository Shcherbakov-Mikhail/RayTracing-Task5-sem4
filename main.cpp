#include "Object.hpp"
#include "MeshTriangle.hpp"
#include "functions.hpp"

int main(void)
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