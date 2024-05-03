#include <vector>
#include <string>
#include <cstring>
#include <iostream>
#include <fstream>
#include <cmath>
#include "./include/defines.hpp"
#include "./include/file_interaction.hpp"

#include <hpx/hpx_init.hpp>

#include<sys/sysinfo.h>

char* algorithm_picker(unsigned int number_of_cores)
{
    std::string command_line_instruction = "./lattice_boltzmann";

    // Specify HPX settings
    command_line_instruction.append(" -t");
    command_line_instruction.append(std::to_string(number_of_cores));
    command_line_instruction.append(" --hpx:bind=thread:0-");
    command_line_instruction.append(std::to_string(number_of_cores - 1));
    command_line_instruction.append("=core:0-");
    command_line_instruction.append(std::to_string(number_of_cores - 1));
    command_line_instruction.append(".pu:0");

    // Convert to char array
    const int length = command_line_instruction.length(); 
    char* result = new char[length + 1]; 
    strcpy(result, command_line_instruction.c_str()); 

    return result;
}

Settings csv_test(const std::string& algorithm)
{
    Settings result;
    result.debug_mode = 0;
    result.results_to_csv = 1;
    result.horizontal_nodes = 200;
    result.algorithm = algorithm;
    result.vertical_nodes_excluding_buffers = 50;
    result.subdomain_count = 5;
    result.access_pattern = "bundle";
    result.time_steps = 800;
    return result;
}

void strong_scaling_new
(
    const std::vector<std::string> &sequential_algorithms,
    const std::vector<std::string> &parallel_algorithms,
    const std::vector<std::string> &access_patterns,
    const std::vector<unsigned int> &multi_core_counts,
    double relaxation_time,
    unsigned int time_steps
)
{
    hpx::chrono::high_resolution_timer timer;
    unsigned int test_runs = 20;

    std::vector<std::string> result_lines;
    std::string current_line{"algorithm,access_pattern,cores"};
    for(auto i = 0; i < test_runs; ++i)
    {
        current_line += ", runtime test " + std::to_string(i) + " [s]";
    }
    current_line += "\n";

    result_lines.push_back(current_line);
    current_line = {};

    Settings settings;
    settings.debug_mode = 0;
    settings.results_to_csv = 0;
    settings.horizontal_nodes = 512;
    settings.vertical_nodes_excluding_buffers = 128;
    settings.time_steps = time_steps;

    std::vector<double> runtimes{};

    double total = ((parallel_algorithms.size()) * multi_core_counts.size() + sequential_algorithms.size()) * access_patterns.size();
    double progress = 0; 

    std::cout << "Starting strong scaling test." << std::endl;
    std::cout << "Progress: " << floor(100* progress / total) << " %\r" << std::flush;

    /* Sequential tests */ 

    settings.subdomain_count = 0;

    for(std::string algorithm : sequential_algorithms)
    {
        settings.algorithm = algorithm;
        
        for(const std::string& access_pattern : access_patterns) 
        {
            settings.access_pattern = access_pattern;
        
            // Write options file
            write_csv_config_file(settings);

            // Execute algorithm
            for(auto i = 0; i < test_runs; ++i)
            {
                timer.restart();
                system("./lattice_boltzmann");     
                runtimes.push_back(timer.elapsed());
            }

            // Evaluate data
            current_line = algorithm + "," + access_pattern + "," + std::to_string(1);
            
            
            for(double runtime : runtimes)
            {
                current_line += "," + std::to_string(runtime);
            }

            current_line += "\n";

            result_lines.push_back(current_line);
            current_line = {};
            runtimes = {};

            progress++;
            std::cout << "Progress: " << floor(100* progress / total) << " %\r" << std::flush;
        }
    }

    /* Parallel tests */ 
    for(std::string algorithm : parallel_algorithms)
    {
        settings.algorithm = algorithm;

        for(const std::string& access_pattern : access_patterns) 
        {
            settings.access_pattern = access_pattern;

            for(const auto current_cores : multi_core_counts)
            {
                settings.subdomain_count = current_cores;

                // Write options file
                write_csv_config_file(settings);

                // Execute algorithm
                for(auto i = 0; i < test_runs; ++i)
                {
                    timer.restart();
                    system(algorithm_picker(current_cores));       
                    runtimes.push_back(timer.elapsed());
                }

                // Evaluate data
                current_line = algorithm + "," + access_pattern + "," + std::to_string(current_cores);

                for(double runtime : runtimes)
                {
                    current_line += "," + std::to_string(runtime);
                }

                current_line += "\n";

                result_lines.push_back(current_line);
                current_line = {};
                runtimes = {};

                progress++;
                std::cout << "Progress: " << floor(100* progress / total) << " %\r" << std::flush;
            }
        }
    } 
    std::cout << "Storing results to 'strong_scaling_results.csv'..." << std::endl;

    std::ofstream file;
    file.open("strong_scaling_results.csv");

    for(const auto& line : result_lines)
    {
        file << line;
    }
    file.close();
    
    std::cout << "Done." << std::endl;
    std::cout << std::endl;
}

void weak_scaling_new
(
    const std::vector<std::string> &sequential_algorithms,
    const std::vector<std::string> &parallel_algorithms,
    const std::vector<std::string> &access_patterns,
    const std::vector<unsigned int> &multi_core_counts,
    double relaxation_time,
    unsigned int time_steps  
)
{
    hpx::chrono::high_resolution_timer timer;
    unsigned int test_runs = 20;

    std::vector<std::string> result_lines;
    std::string current_line{"algorithm,access_pattern,cores"};
    for(auto i = 0; i < test_runs; ++i)
    {
        current_line += ", runtime test " + std::to_string(i) + " [s]";
    }
    current_line += "\n";

    result_lines.push_back(current_line);
    current_line = {};

    unsigned int base_subdomain_height = 128;
    unsigned int horizontal_nodes = 512;

    Settings settings;
    settings.debug_mode = 0;
    settings.results_to_csv = 0;
    settings.horizontal_nodes = horizontal_nodes;
    settings.time_steps = time_steps;

    std::vector<double> runtimes{};

    double total = ((parallel_algorithms.size()) * multi_core_counts.size() + sequential_algorithms.size()) * access_patterns.size();
    double progress = 0; 

    std::cout << "Starting weak scaling test." << std::endl;
    std::cout << "Progress: " << floor(100* progress / total) << " %\r" << std::flush;

    /* Sequential tests */ 

    settings.subdomain_count = 0;
    settings.vertical_nodes_excluding_buffers = base_subdomain_height * 1;

    for(std::string algorithm : sequential_algorithms)
    {
        settings.algorithm = algorithm;

        for(const std::string& access_pattern : access_patterns) 
        {
            settings.access_pattern = access_pattern;
        
            // Write options file
            write_csv_config_file(settings);

            // Execute algorithm
            for(auto i = 0; i < test_runs; ++i)
            {
                timer.restart();
                system("./lattice_boltzmann");     
                runtimes.push_back(timer.elapsed());
            }

            // Evaluate data
            current_line = algorithm + "," + access_pattern + "," + std::to_string(1);
            
            for(double runtime : runtimes)
            {
                current_line += "," + std::to_string(runtime);
            }

            current_line += "\n";

            result_lines.push_back(current_line);
            current_line = {};
            runtimes = {};

            progress++;
            std::cout << "Progress: " << floor(100* progress / total) << " %\r" << std::flush;
        }
    }

    /* Parallel tests */
    
    for(std::string algorithm : parallel_algorithms)
    {
        settings.algorithm = algorithm;

        for(const std::string& access_pattern : access_patterns) 
        {
            settings.access_pattern = access_pattern;

            for(unsigned int core_count : multi_core_counts)
            {
                settings.vertical_nodes_excluding_buffers = base_subdomain_height * core_count;
                settings.subdomain_count = core_count;
            
                // Write options file
                write_csv_config_file(settings);

                // Execute algorithm
                for(auto i = 0; i < test_runs; ++i)
                {
                    timer.restart();
                    system(algorithm_picker(core_count));     
                    runtimes.push_back(timer.elapsed());
                }

                // Evaluate data
                current_line = algorithm + "," + access_pattern + "," + std::to_string(core_count);

                for(double runtime : runtimes)
                {
                    current_line += "," + std::to_string(runtime);
                }

                current_line += "\n";

                result_lines.push_back(current_line);
                current_line = {};
                runtimes = {};

                progress++;
                std::cout << "Progress: " << floor(100* progress / total) << " %\r" << std::flush;
            }
        }
    }

    std::cout << "Storing results to 'weak_scaling_results.csv'..." << std::endl;

    std::ofstream file;
    file.open("weak_scaling_results.csv");

    for(const auto& line : result_lines)
    {
        file << line;
    }
    file.close();
    
    std::cout << "Done." << std::endl;
    std::cout << std::endl;
}

int main(int argc, char* argv[])
{
    /* Selections that actually vary */
    std::vector<std::string> sequential_algorithms{"sequential_two_lattice", "sequential_two_step", "sequential_swap", "sequential_shift"};
    std::vector<std::string> parallel_algorithms{"parallel_two_lattice", "parallel_two_lattice_framework", "parallel_two_step", "parallel_swap", "parallel_shift"};
    std::vector<std::string> access_patterns{"collision", "stream", "bundle"};

    /* Selections assumed static */
    const double relaxation_time = 1.4;
    const unsigned int time_steps = 20;
    const int write_csv = 0;
    const int debug_mode = 0;

    unsigned int available_cores = std::thread::hardware_concurrency() / 2;
    std::cout << "Up to " << available_cores << " concurrent threads are supported.\n";

    std::vector<unsigned int> multicore_setups;
    unsigned int current_max_core_count = 2;

    while(current_max_core_count < available_cores)
    {
        multicore_setups.push_back(current_max_core_count);
        current_max_core_count *= 2;
    }

    strong_scaling_new(sequential_algorithms, parallel_algorithms, access_patterns, multicore_setups, relaxation_time, time_steps);
    weak_scaling_new(sequential_algorithms, parallel_algorithms, access_patterns, multicore_setups, relaxation_time, time_steps);

    std::cout << "Benchmark finished." << std::endl;
}

    /* Examples */

    // system("./parallel_two_lattice_framework -t6 --hpx:bind=thread:0-5=core:0-5.pu:0 --hpx:dump-config");
    // system("./parallel_two_lattice_framework -t4 --hpx:bind=thread:0-3=core:0-3.pu:0 --hpx:dump-config");
    // system("./parallel_two_lattice_framework -t2 --hpx:bind=thread:0-1=core:0-1.pu:0 --hpx:dump-config");

    /* Example execution */
    // std::cout << algorithm_picker("two_lattice", true, 4) << std::endl;
    // system(algorithm_picker("two_lattice", true, 4));

    // for(bool parallel : select_parallel_algorithm)
    // {
    //     for(std::string access_pattern : access_patterns)
    //     {
    //         for(std::string algorithm : algorithms)
    //         {
    //             for(unsigned int current_pow_2_vertical_nodes : pow_2_vertical_nodes)
    //             {  
    //                for(unsigned int current_horizontal_nodes : horizontal_nodes) 
    //                {
    //                     // Write options file
    //                     // execute algorithm
    //                }
    //             }
    //         }
    //     }
    // }

    // Settings example = csv_test("parallel_two_lattice_framework");
    // write_csv_config_file(example);
