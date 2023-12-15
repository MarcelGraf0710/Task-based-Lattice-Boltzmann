#include "../include/update.hpp"
#include "../include/simulation_in_cool.hpp"

#include "../include/access.hpp"
#include "../include/boundaries.hpp"
#include "../include/macroscopic.hpp"

#include "../include/utils.hpp"

void print_happy_life()
{
    std::vector<double> ways_to_die{1,2,3,4,5,6,7,8,9};
    std::vector<double> more_ways_to_die{0,0,0,0,0,0,0,0,0};
    access_function f = access::collision;
    std::cout << "I think I'm gonna kill myself, just a little suicide like it's 19" << matrix_access(3,3,3) << std::endl;
    std::cout << "I think I'm gonna kill myself, just a little suicide like it's 19" << access::collision(3,3) << std::endl;
    boundaries::lower_inlet_boundary_stream(ways_to_die, more_ways_to_die, f);
    std::cout << "This is fine." << std::endl;
    std::cout << "I think I'm gonna kill myself, just a little suicide like it's 19" << more_ways_to_die[3] << std::endl;
    std::cout << "I think I'm gonna kill myself, just a little suicide like it's 19" << macroscopic::density(more_ways_to_die) << std::endl;
    ways_to_die = collision::collide_bgk(more_ways_to_die, {.5, .5}, 2);
    std::cout << "I think I'm gonna kill myself, just a little suicide like it's 19" << ways_to_die[1] << std::endl;
}