#include <string>
#include <algorithm>
#include <fstream>
#include <iostream>
#include <string>
#include "../include/file_interaction.hpp"


void load_settings_from_file(std::string &filename)
{
    /* Open file */
    std::ifstream settings_file{filename};

    if(!settings_file.is_open())
    {
        std::cout << "Could not open file " << filename << std::endl;
    }
    else
    {
        /* Line-wise reading of settings file */
        std::string current_line{};
		std::string parameter_name{};
		std::string parameter_value_string{};
        size_t position_of_equality_sign{};

        while (std::getline(settings_file, current_line))
        {
            // Only lines that are NOT a comment or empty or without an equality sign require treatment
            if(((current_line[0] != '#') && (current_line.find_first_of("= \t") != std::string::npos)))
            {
                // Set parameter name to all signs before the equality sign
                parameter_name = current_line.substr(0, current_line.find_first_of('='));

                // Delete any white spaces
                parameter_name.erase(std::remove(parameter_name.begin(), parameter_name.end(), ' '), parameter_name.end());

                position_of_equality_sign = current_line.find_first_of('=');
				parameter_value_string = current_line.substr(position_of_equality_sign + 1,
											current_line.find_first_of("#\t\n") - position_of_equality_sign - 1);
            }
        }
        // TODO: FINISH
    }
}

void sim_data_to_csv(std::vector<sim_data_tuple> &data, const std::string &filename)
{
    std::ofstream file;
    unsigned int current_node = 0;
    file.open("results.csv");
    file << "iteration,x,y,vx,vy,density\n"; 
    for(auto time = 0; time < data.size(); ++time)
    {
        for(auto y = 1; y < VERTICAL_NODES - 1; ++y)
        {
            for(auto x = 1; x < HORIZONTAL_NODES - 1; ++x)
            {
                current_node = lbm_access::get_node_index(x,y);
                file << time << ',' << x << ',' << y << ',' 
                << std::get<0>(data[time])[current_node][0] << ',' 
                << std::get<0>(data[time])[current_node][1] << ',' << std::get<1>(data[time])[current_node] << '\n';
            }
        }
    }
    file.close();
}