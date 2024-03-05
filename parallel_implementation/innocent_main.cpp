#include "include/access.hpp"
#include "include/simulation.hpp"
#include "include/utils.hpp"
#include "include/two_lattice_sequential.hpp"
#include <hpx/hpx_init.hpp>

int hpx_main()
{
    std::cout << "Hi there! I am just an innocent main. I just want to be executed. Please. Have mercy with me HPX. "<< std::endl;
    return hpx::finalize();
}

int main()
{
    // Initialize HPX, run hpx_main as the first HPX thread, and
    // wait for hpx::finalize being called.
    return hpx::init();
}