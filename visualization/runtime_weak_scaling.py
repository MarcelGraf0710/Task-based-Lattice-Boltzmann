# Imports
import csv
import numpy
import matplotlib.pyplot as plot

# Initializations
core_ticks = []

read_content = []

access_patterns = ["collision", "stream", "bundle"]

# Store content of file with weak scaling results
with open("../build/weak_scaling_results.csv", "r") as f:
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

    # Assign times to corresponding list
    for line in read_content:
        if line[1] == access_pattern:
            if line[0] == "parallel_two_lattice":
                parallel_two_lattice_results.append(float(line[-1]))
            elif line[0] == "parallel_two_lattice_framework":
                parallel_two_lattice_framework_results.append(float(line[-1]))
            elif line[0] == "parallel_two_step" or line[0] == "sequential_two_step":
                parallel_two_step_results.append(float(line[-1]))
            elif line[0] == "parallel_swap" or line[0] == "sequential_swap":
                parallel_swap_results.append(float(line[-1]))
            elif line[0] == "parallel_shift" or line[0] == "sequential_shift":
                parallel_shift_results.append(float(line[-1]))
            elif line[0] == "sequential_two_lattice":
                parallel_two_lattice_results.append(float(line[-1]))
                parallel_two_lattice_framework_results.append(float(line[-1]))

    # Regression with polynomial
    model_tl = numpy.poly1d(numpy.polyfit(core_ticks, parallel_two_lattice_results, len(core_ticks)-1))
    model_tlf = numpy.poly1d(numpy.polyfit(core_ticks, parallel_two_lattice_framework_results, len(core_ticks)-1))
    model_ts = numpy.poly1d(numpy.polyfit(core_ticks, parallel_two_step_results, len(core_ticks)-1))
    model_shift = numpy.poly1d(numpy.polyfit(core_ticks, parallel_shift_results, len(core_ticks)-1))
    model_swap = numpy.poly1d(numpy.polyfit(core_ticks, parallel_swap_results, len(core_ticks)-1))
    line_linspace = numpy.linspace(1, core_ticks[-1], 100)

    ### Plot
    plot.figure(dpi=200)

    plot.plot(core_ticks, parallel_two_lattice_results, 'o', color='k', markersize=5)
    plot.plot(core_ticks, parallel_two_lattice_framework_results, '*', color='y', markersize=5)
    plot.plot(core_ticks, parallel_two_step_results, 's', color='r', markersize=5)
    plot.plot(core_ticks, parallel_swap_results, 'D', color='g', markersize=5)
    plot.plot(core_ticks, parallel_shift_results, 'p', color='b', markersize=5)

    plot.plot(line_linspace, model_tl(line_linspace), color='k')
    plot.plot(line_linspace, model_tlf(line_linspace), color='y')
    plot.plot(line_linspace, model_ts(line_linspace), color='r')
    plot.plot(line_linspace, model_swap(line_linspace), color='g')
    plot.plot(line_linspace, model_shift(line_linspace), color='b')

    plot.grid(color='#808080', linestyle='--', linewidth=0.5)
    plot.legend(['Non-framework two-lattice', 'Framework two-lattice', 'Two-step', 'Swap', 'Shift'], loc='best')
    plot.xticks(core_ticks)
    plot.title("Weak scaling: 32x128 nodes per core, " + access_pattern + " layout")
    plot.xlabel('cores')
    plot.ylabel('time[s]')
    plot.savefig("../images/runtimes_weak_scaling_" + access_pattern + ".pdf", format="pdf", bbox_inches="tight")
    plot.close()