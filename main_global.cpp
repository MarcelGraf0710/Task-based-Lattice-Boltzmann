#include <hpx/hpx_init.hpp>
#include <hpx/execution.hpp>
#include <string>
#include "./include/defines.hpp"
#include "./include/file_interaction.hpp"
#include "include/parallel_two_lattice_framework.hpp"
#include "include/sequential_shift.hpp"
#include "include/parallel_shift_framework.hpp"
#include "include/lbm_execution.hpp"

int hpx_main(hpx::program_options::variables_map& vm)
{
    Settings settings = retrieve_settings_from_csv("config.csv");
    setup_global_variables(settings);
    select_and_execute(settings.algorithm);
    return hpx::local::finalize();
}


int main(int argc, char* argv[])
{
    hpx::program_options::options_description desc_commandline("Usage: " HPX_APPLICATION_STRING " [options]");

    hpx::local::init_params init_args;
    init_args.desc_cmdline = desc_commandline;

    return hpx::local::init(hpx_main, argc, argv, init_args);

    // system("./parallel_two_lattice_framework -t6 --hpx:bind=thread:0-5=core:0-5.pu:0 --hpx:dump-config");
    // system("./parallel_two_lattice_framework -t4 --hpx:bind=thread:0-3=core:0-3.pu:0 --hpx:dump-config");
    // system("./parallel_two_lattice_framework -t2 --hpx:bind=thread:0-1=core:0-1.pu:0 --hpx:dump-config");
}