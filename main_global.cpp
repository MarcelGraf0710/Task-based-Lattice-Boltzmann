#include <hpx/hpx_init.hpp>
#include <hpx/execution.hpp>
#include <string>

int hpx_main(hpx::program_options::variables_map& vm)
{
    std::string algorithm = vm["algorithm"].as<std::string>();
    std::string data_layout = vm["data_layout"].as<std::string>();
    double relaxation_time = vm["relaxation_time"].as<double>();
    std::uint64_t time_steps = vm["time_steps"].as<std::uint64_t>();
    std::uint64_t horizontal_nodes = vm["horizontal_nodes"].as<std::uint64_t>();
    std::uint64_t vertical_nodes = vm["vertical_nodes"].as<std::uint64_t>();

    std::cout << "Trying to launch sequential " << algorithm << " algorithm for " << time_steps << " time steps with the following parameters: " << std::endl;
    std::cout << "Data layout: " << data_layout << ", relaxation_time: " << relaxation_time << ", horizontal nodes: " << horizontal_nodes << ", vertical nodes: " << vertical_nodes << std::endl;
    // std::cout << "Lol you wish! Here, have some gedit instead..." << std::endl;
    // system("gedit");

    if(algorithm == "two-lattice")
    {
        system("./parallel_two_lattice_framework --hpx:bind=compact --hpx:cores=3 --hpx:dump-config");
    }
    else if(algorithm == "two-step")
    {
        system("./sequential_two_step");
    }
    else if(algorithm == "swap")
    {
        system("./sequential_swap");
    }
    else if(algorithm == "shift")
    {
        system("./sequential_shift");
    }
    else
    {
        std::cout << "Unknown algorithm: " << algorithm << std::endl;
    }
    return hpx::local::finalize();
}

int main(int argc, char* argv[])
{
    hpx::program_options::options_description desc_commandline("Usage: " HPX_APPLICATION_STRING " [options]");
    
    /* General definitions */
    desc_commandline.add_options()("algorithm",
    hpx::program_options::value<std::string>()->default_value("two-lattice"),
    "This string is used to determine the algorithm that is to be used");

    desc_commandline.add_options()("data_layout",
    hpx::program_options::value<std::string>()->default_value("collision"),
    "This string is used to determine the data layout that is to be used");

    desc_commandline.add_options()("relaxation_time",
    hpx::program_options::value<double>()->default_value(1.4),
    "Relaxation time within the simulation");

    desc_commandline.add_options()("time_steps",
    hpx::program_options::value<std::uint64_t>()->default_value(50),
    "This many iterations of the algorithm will be executed");

    desc_commandline.add_options()("horizontal_nodes",
    hpx::program_options::value<std::uint64_t>()->default_value(7),
    "The amount of horizontal nodes including ghost nodes");

    desc_commandline.add_options()("vertical_nodes",
    hpx::program_options::value<std::uint64_t>()->default_value(24),
    "The amount of vertical nodes including ghost nodes");

    hpx::local::init_params init_args;
    init_args.desc_cmdline = desc_commandline;

    return hpx::local::init(hpx_main, argc, argv, init_args);

    // system("./parallel_two_lattice_framework -t6 --hpx:bind=thread:0-5=core:0-5.pu:0 --hpx:dump-config");
    // system("./parallel_two_lattice_framework -t4 --hpx:bind=thread:0-3=core:0-3.pu:0 --hpx:dump-config");
    // system("./parallel_two_lattice_framework -t2 --hpx:bind=thread:0-1=core:0-1.pu:0 --hpx:dump-config");
}