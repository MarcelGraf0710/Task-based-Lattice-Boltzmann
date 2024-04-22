#include <vector>
#include <string>
#include <cstring>
#include <iostream>
#include <fstream>
#include "./include/defines.hpp"
#include "./include/file_interaction.hpp"

char* algorithm_picker(const std::string &algorithm, bool parallel, unsigned int number_of_cores)
{
    std::string command_line_instruction = "";

    // Select if algorithm is parallel or not
    if(parallel)
    {
        command_line_instruction = "./parallel_";
    }
    else
    {
        command_line_instruction = "./sequential_";
    }

    // Select algorithm that is to be executed
    if(algorithm != "two_lattice" && algorithm != "two_step" && algorithm != "shift" && algorithm != "swap")
    {
        std::cout << "Unknown algorithm: " << algorithm << std::endl;
    }
    else
    {
        command_line_instruction.append(algorithm);
    }

    // Specify HPX settings
    command_line_instruction.append(" -t");
    command_line_instruction.append(std::to_string(number_of_cores));
    command_line_instruction.append(" --hpx:bind=thread:0-");
    command_line_instruction.append(std::to_string(number_of_cores - 1));
    command_line_instruction.append("=core:0-");
    command_line_instruction.append(std::to_string(number_of_cores - 1));
    command_line_instruction.append(".pu:0 --hpx:dump-config");

    // Convert to char array
    const int length = command_line_instruction.length(); 
    char* result = new char[length + 1]; 
    strcpy(result, command_line_instruction.c_str()); 

    return result;
}



int main(int argc, char* argv[])
{
    /* Selections that actually vary */
    std::vector<bool> select_parallel_algorithm{false, true};
    std::vector<std::string> algorithms{"two_lattice", "two_step", "swap", "shift"};
    std::vector<std::string> access_patterns{"collision", "stream", "bundle"};
    std::vector<unsigned int> pow_2_vertical_nodes{4,5,6}; // Will be: {4,5,6,7,8,9,10,11,12,13,14}
    std::vector<unsigned int> horizontal_nodes{10,20,50,100,200,500}; // Will be: {10,20,50,100,200,500,1000,2000,5000,10000,20000,50000}

    /* Selections assumed static */
    const double relaxation_time = 1.4;
    const unsigned int time_steps = 20;
    const int write_csv = 0;
    const int debug_mode = 0;

    /* Examples */

    // system("./parallel_two_lattice_framework -t6 --hpx:bind=thread:0-5=core:0-5.pu:0 --hpx:dump-config");
    // system("./parallel_two_lattice_framework -t4 --hpx:bind=thread:0-3=core:0-3.pu:0 --hpx:dump-config");
    // system("./parallel_two_lattice_framework -t2 --hpx:bind=thread:0-1=core:0-1.pu:0 --hpx:dump-config");

    /* Example execution */
    std::cout << algorithm_picker("two_lattice", true, 4) << std::endl;
    system(algorithm_picker("two_lattice", true, 4));

    for(bool parallel : select_parallel_algorithm)
    {
        for(std::string access_pattern : access_patterns)
        {
            for(std::string algorithm : algorithms)
            {
                for(unsigned int current_pow_2_vertical_nodes : pow_2_vertical_nodes)
                {  
                   for(unsigned int current_horizontal_nodes : horizontal_nodes) 
                   {
                        // Write options file
                        // execute algorithm
                   }
                }
            }
        }
    }


}

    // Settings example1;
    // example1.is_parallel = true;
    // example1.access_pattern = "stream";
    // example1.relaxation_time = 1.0;
    // example1.vertical_nodes_excluding_buffers = 32;
    // example1.subdomain_count = 4;
    // example1.horizontal_nodes = 10;
    // example1.inlet_velocity = {0.2,0.0};
    // example1.outlet_velocity = {0.0, 0.0};
    // example1.inlet_density = 1.0;
    // example1.outlet_density = 1.0;

    // // write_csv_config_file(example1);


    // // Settings example2;
    // // example2.is_parallel = false;
    // // example2.access_pattern = "stream";
    // // example2.relaxation_time = 1.6;
    // // example2.vertical_nodes_excluding_buffers = 21;
    // // example2.horizontal_nodes = 11;
    // // example2.inlet_velocity = {0.25,0.0};
    // // example2.outlet_velocity = {0.0, 0.0};
    // // example2.inlet_density = 1.0;
    // // example2.outlet_density = 1.0;

    // Settings example3 = retrieve_settings_from_csv("retrieve_me.csv");
    // write_csv_config_file(example3);