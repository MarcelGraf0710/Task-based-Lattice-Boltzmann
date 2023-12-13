#include "boundaries.hpp"
//#include "access.hpp"
//#include "defines.hpp"

    /**
     * @brief This map contains the directions that each kind of border node is facing a neighbor node,
     *        i.e. it contains all directions that perform actual streams.
     */
    std::map<boundaries::boundary_tuple, std::list<int>> boundaries::neighbor_directions
    {
        {boundary_scenarios::inlet, {1,2,5,7,8}},
        {boundary_scenarios::outlet, {0,1,3,6,7}},
        {boundary_scenarios::lower_inlet, {5,7,8}},
        {boundary_scenarios::upper_inlet, {1,2,5}},
        {boundary_scenarios::lower_outlet, {3,6,7}},
        {boundary_scenarios::upper_outlet, {0,1,3}},
        {boundary_scenarios::wall_up, {0,1,2,3,5}},
        {boundary_scenarios::wall_down, {3,5,6,7,8}},
    };

    void boundaries::inlet_boundary_stream(
        all_distributions &source, 
        all_distributions &destination, 
        access_function &access_function, 
        int y,
        double velocity_x, 
        double density)
    {
        vec_of_dist_val result;
        //unsigned int node_index = access::get_node_index(0, y);
        //result = access::get_all_distribution_values(source, node_index, access_function);
        double rho_times_u_x = 1.0/6 * density * velocity_x;

        std::vector<unsigned int> changes = {5, 8, 2};
        result[5] = result[3] + 2.0/3 * density * velocity_x;
        result[8] = result[0] - 0.5 * (result[7] - result[1]) + rho_times_u_x;
        result[2] = result[6] + 0.5 * (result[7] - result[1]) + rho_times_u_x;
        
        //access::set_all_distribution_values(result, destination, node_index, access_function);
    }

    void boundaries::outlet_boundary_stream(
        all_distributions &source, 
        all_distributions &destination, 
        access_function &access_function, 
        int y, 
        double velocity_x, 
        double density)
    {
        vec_of_dist_val result;
        //unsigned int node_index = access::get_node_index(HORIZONTAL_NODES - 1, y);
        //result = access::get_all_distribution_values(source, node_index, access_function);
        double rho_times_u_x = 1.0/6 * density * velocity_x;

        std::vector<unsigned int> changes = {5, 8, 2};
        result[3] = result[5] - 2.0/3 * density * velocity_x;
        result[0] = result[8] + 0.5 * (result[7] - result[1]) - rho_times_u_x;
        result[6] = result[2] - 0.5 * (result[7] - result[1]) - rho_times_u_x;

        //access::set_all_distribution_values(result, destination, node_index, access_function);
    }

    double boundaries::lower_wall_boundary_stream(  
        all_distributions &source,      
        all_distributions &destination, 
        access_function &access_function, 
        int x, 
        velocity wall_velocity)
    {
        vec_of_dist_val result;
        //unsigned int node_index = access::get_node_index(x, 0);
        //result = access::get_all_distribution_values(source, node_index, access_function);
        double density = 1 / (1 - wall_velocity[1]) * (result[4] + result[5] + result[3] + 2 * (result[1] + result[0] + result[2]));

        double rho_times_u_x = 0.5 * density * wall_velocity[0];
        double rho_times_u_y = 1.0/6 * density * wall_velocity[1];

        std::vector<unsigned int> changes = {6,7,8};
        result[7] = result[1] + 2.0/3 * density * wall_velocity[1];
        result[8] = result[0] + 0.5 * (result[5] - result[3]) + rho_times_u_x + rho_times_u_y;
        result[6] = result[2] + 0.5 * (result[5] - result[3]) - rho_times_u_x + rho_times_u_y;

        //access::set_all_distribution_values(result, destination, node_index, access_function);
        return density;
    }

    double boundaries::upper_wall_boundary_stream(   
        all_distributions &source,      
        all_distributions &destination, 
        access_function &access_function, 
        int x, 
        velocity wall_velocity)
    {
        vec_of_dist_val result;
        //unsigned int node_index = access::get_node_index(x, VERTICAL_NODES - 1);
        //result = access::get_all_distribution_values(source, node_index, access_function);
        double density = 1 / (1 - wall_velocity[1]) * (result[4] + result[5] + result[3] + 2 * (result[1] + result[0] + result[2]));

        double rho_times_u_x = 0.5 * density * wall_velocity[0];
        double rho_times_u_y = 1.0/6 * density * wall_velocity[1];

        std::vector<unsigned int> changes = {0,1,2};
        result[1] = result[7] - 2.0/3 * density * wall_velocity[1];
        result[0] = result[8] - 0.5 * (result[3] - result[5]) - rho_times_u_x - rho_times_u_y;
        result[2] = result[6] + 0.5 * (result[3] - result[5]) + rho_times_u_x - rho_times_u_y;

        //access::set_all_distribution_values(result, destination, node_index, access_function);
        return density;
    }

    void boundaries::lower_inlet_boundary_stream(
        all_distributions &source, 
        all_distributions &destination,
        access_function access_function,
        double density)
    {
        vec_of_dist_val result;
        //unsigned int lower_inlet_node = access::get_node_index(0, 0);
        //result = access::get_all_distribution_values(source, lower_inlet_node, access_function);

        std::vector<unsigned int> changes = {5, 7, 8, 2, 6};
        result[5] = result[3];
        result[7] = result[1];
        result[8] = result[0];
        result[2] = 0.5 * (density - (result[4] + 2 * result[0] + 2 * result[1] + 2 * result[3]));
        result[6] = result[2];

        //access::set_all_distribution_values(result, destination, lower_inlet_node, access_function);
    }

    void boundaries::upper_inlet_boundary_stream(
        all_distributions &source, 
        all_distributions &destination, 
        access_function access_function,
        double density)
    {
        vec_of_dist_val result;
        //unsigned int upper_inlet_node = access::get_node_index(0, HORIZONTAL_NODES - 1);
        //result = access::get_all_distribution_values(source, upper_inlet_node, access_function);

        std::vector<unsigned int> changes = {5, 1, 2, 0, 8};
        result[5] = result[3];
        result[1] = result[7];
        result[2] = result[6];
        result[0] = 0.5 * (density - (result[4] + 2 * result[3] + 2 * result[6] + 2 * result[7]));
        result[8] = result[0];

        //access::set_all_distribution_values(result, destination, upper_inlet_node, access_function);
    }

    void boundaries::lower_outlet_boundary_stream(
        all_distributions &source,        
        all_distributions &destination,        
        access_function access_function,
        double density)
    {
        vec_of_dist_val result;
        //unsigned int lower_outlet_node = access::get_node_index(VERTICAL_NODES - 1, 0);
        //result = access::get_all_distribution_values(source, lower_outlet_node, access_function);
        
        std::vector<unsigned int> changes = {3,7,6,0,8};
        result[3] = result[5];
        result[7] = result[1];
        result[6] = result[2];
        result[0] = 0.5 * (density - (result[4] + 2 * result[1] + 2 * result[2] + 2 * result[5]));
        result[8] = result[0];

        //access::set_all_distribution_values(result, destination, lower_outlet_node, access_function);
    }

    void boundaries::upper_outlet_boundary_stream(   
        all_distributions &source, 
        all_distributions &destination,         
        access_function access_function,
        double density)
    {
        vec_of_dist_val result;
        //unsigned int lower_outlet_node = access::get_node_index(VERTICAL_NODES - 1, HORIZONTAL_NODES - 1);
        //result = access::get_all_distribution_values(source, lower_outlet_node, access_function);
        
        std::vector<unsigned int> changes = {3,1,0,2,6};
        result[3] = result[5];
        result[1] = result[7];
        result[0] = result[8];
        result[2] = 0.5 * (density - (result[4] + 2 * result[5] + 2 * result[7] + 2 * result[8]));
        result[6] = result[2];

        //access::set_all_distribution_values(result, destination, lower_outlet_node, access_function);
    }

