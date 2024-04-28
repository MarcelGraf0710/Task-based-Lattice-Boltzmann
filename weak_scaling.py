# Imports
import csv
import numpy
import matplotlib.pyplot as plot

# Initializations
core_ticks = []

read_content = []

access_patterns = ["collision", "stream", "bundle"]

# Store content of file with weak scaling results
with open("./build/weak_scaling_results.csv", "r") as f:
    reader = csv.reader(f, delimiter=",")
    for i, line in enumerate(reader):
        read_content.append(line)

# Determine core ticks
max_core_count = 0
for i in range(1, len(read_content)):
    if int(read_content[i][2]) > max_core_count:
        core_ticks.append(int(read_content[i][2]))
        max_core_count = int(read_content[i][2])
    elif int(read_content[i][2]) == max_core_count:
        max_core_count = int(read_content[i][2])
    else:
        break

# One plot for each access pattern
for access_pattern in access_patterns:

    # Initializations
    parallel_two_lattice_results = []
    parallel_two_lattice_framework_results = []
    parallel_two_step_results = []
    parallel_swap_results = []
    parallel_shift_results = []

    parallel_two_lattice_weak_scaling = []
    parallel_two_lattice_framework_weak_scaling = []
    parallel_two_step_weak_scaling = []
    parallel_swap_weak_scaling = []
    parallel_shift_weak_scaling = []

    # Assign times to corresponding list
    for line in read_content:
        if line[1] == access_pattern:
            if line[0] == "parallel_two_lattice":
                parallel_two_lattice_results.append(line[-1])
            elif line[0] == "parallel_two_lattice_framework":
                parallel_two_lattice_framework_results.append(line[-1])
            elif line[0] == "parallel_two_step" or line[0] == "sequential_two_step":
                parallel_two_step_results.append(line[-1])
            elif line[0] == "parallel_swap" or line[0] == "sequential_swap":
                parallel_swap_results.append(line[-1])
            elif line[0] == "parallel_shift" or line[0] == "sequential_shift":
                parallel_shift_results.append(line[-1])
            elif line[0] == "sequential_two_lattice":
                parallel_two_lattice_results.append(line[-1])
                parallel_two_lattice_framework_results.append(line[-1])

    # Determine weak scaling factor based on results
    for i in range(0, len(core_ticks)):
        core_count = core_ticks[i]
        parallel_two_lattice_weak_scaling.append((core_count * float(parallel_two_lattice_results[0])) / float(parallel_two_lattice_results[i]))
        parallel_two_lattice_framework_weak_scaling.append((core_count * float(parallel_two_lattice_framework_results[0])) / float(parallel_two_lattice_framework_results[i]))
        parallel_two_step_weak_scaling.append((core_count * float(parallel_two_step_results[0])) / float(parallel_two_step_results[i]))
        parallel_swap_weak_scaling.append((core_count * float(parallel_swap_results[0])) / float(parallel_swap_results[i]))
        parallel_shift_weak_scaling.append((core_count * float(parallel_shift_results[0])) / float(parallel_shift_results[i]))

    # Regression with polynomial
    model_tl = numpy.poly1d(numpy.polyfit(core_ticks, parallel_two_lattice_weak_scaling, len(core_ticks)-1))
    model_tlf = numpy.poly1d(numpy.polyfit(core_ticks, parallel_two_lattice_framework_weak_scaling, len(core_ticks)-1))
    model_tf = numpy.poly1d(numpy.polyfit(core_ticks, parallel_two_lattice_framework_weak_scaling, len(core_ticks)-1))
    model_ts = numpy.poly1d(numpy.polyfit(core_ticks, parallel_two_step_weak_scaling, len(core_ticks)-1))
    model_shift = numpy.poly1d(numpy.polyfit(core_ticks, parallel_shift_weak_scaling, len(core_ticks)-1))
    model_swap = numpy.poly1d(numpy.polyfit(core_ticks, parallel_swap_weak_scaling, len(core_ticks)-1))
    line_linspace = numpy.linspace(1, core_ticks[-1], 100)

    ### Plot
    plot.figure(dpi=200)

    plot.plot(core_ticks, parallel_two_lattice_weak_scaling, 'o', color='k', markersize=5)
    plot.plot(core_ticks, parallel_two_lattice_framework_weak_scaling, '*', color='y', markersize=5)
    plot.plot(core_ticks, parallel_two_step_weak_scaling, 's', color='r', markersize=5)
    plot.plot(core_ticks, parallel_swap_weak_scaling, 'D', color='g', markersize=5)
    plot.plot(core_ticks, parallel_shift_weak_scaling, 'p', color='b', markersize=5)

    plot.plot(core_ticks, core_ticks, color='k', linestyle='--')

    plot.plot(line_linspace, model_tl(line_linspace), color='k')
    plot.plot(line_linspace, model_tlf(line_linspace), color='y')
    plot.plot(line_linspace, model_ts(line_linspace), color='r')
    plot.plot(line_linspace, model_swap(line_linspace), color='g')
    plot.plot(line_linspace, model_shift(line_linspace), color='b')

    plot.grid(color='#808080', linestyle='--', linewidth=0.5)
    plot.legend(['Non-framework two-lattice', 'Framework two-lattice', 'Two-step', 'Swap', 'Shift', 'Ideal'], loc='best')
    plot.yticks(numpy.linspace(0, core_ticks[-1], core_ticks[-1]*2 + 1))
    plot.xticks(core_ticks)
    plot.title("Weak scaling: 32x128 nodes per core, " + access_pattern + " layout")
    plot.xlabel('cores')
    plot.ylabel('speedup')
    plot.savefig("./images/weak_scaling_" + access_pattern + ".pdf", format="pdf", bbox_inches="tight")
    plot.close()