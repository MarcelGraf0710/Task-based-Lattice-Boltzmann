# Imports
import csv
import numpy
import matplotlib.pyplot as plot

# Initializations
core_ticks = []

read_content_weak = []
read_content_strong = []

access_patterns = ["collision", "stream", "bundle"]

# Store content of file with weak scaling results
with open("../build/weak_scaling_results.csv", "r") as f:
    reader = csv.reader(f, delimiter=",")
    for i, line in enumerate(reader):
        read_content_weak.append(line)

with open("../build/strong_scaling_results.csv", "r") as f:
    reader = csv.reader(f, delimiter=",")
    for i, line in enumerate(reader):
        read_content_strong.append(line)

# Determine core ticks
max_core_count = 0
for i in range(1, len(read_content_weak)):
    if int(read_content_weak[i][2]) > max_core_count:
        core_ticks.append(int(read_content_weak[i][2]))
        max_core_count = int(read_content_weak[i][2])
    elif int(read_content_weak[i][2]) == max_core_count:
        max_core_count = int(read_content_weak[i][2])
    else:
        break

### Weak scaling
print("Creating weak scaling plots...")

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
    for line in read_content_weak:
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

    # Determine weak scaling factor based on results
    for i in range(0, len(core_ticks)):
        core_count = core_ticks[i]
        parallel_two_lattice_weak_scaling.append((core_count * float(parallel_two_lattice_results[0])) / float(parallel_two_lattice_results[i]))
        parallel_two_lattice_framework_weak_scaling.append((core_count * float(parallel_two_lattice_framework_results[0])) / float(parallel_two_lattice_framework_results[i]))
        parallel_two_step_weak_scaling.append((core_count * float(parallel_two_step_results[0])) / float(parallel_two_step_results[i]))
        parallel_swap_weak_scaling.append((core_count * float(parallel_swap_results[0])) / float(parallel_swap_results[i]))
        parallel_shift_weak_scaling.append((core_count * float(parallel_shift_results[0])) / float(parallel_shift_results[i]))

    ## Plot runtimes
    plot.figure(dpi=200)

    plot.plot(core_ticks, parallel_two_lattice_results, 'o-', color='k', markersize=5)
    plot.plot(core_ticks, parallel_two_lattice_framework_results, '*-', color='y', markersize=5)
    plot.plot(core_ticks, parallel_two_step_results, 's-', color='r', markersize=5)
    plot.plot(core_ticks, parallel_swap_results, 'D-', color='g', markersize=5)
    plot.plot(core_ticks, parallel_shift_results, 'p-', color='b', markersize=5)

    plot.grid(color='#808080', linestyle='--', linewidth=0.5)
    plot.legend(['Non-framework two-lattice', 'Framework two-lattice', 'Two-step', 'Swap', 'Shift'], loc='best')
    plot.xticks(core_ticks)
    plot.title("Weak scaling: 32x128 nodes per core, " + access_pattern + " layout")
    plot.xlabel('cores')
    plot.ylabel('time[s]')
    plot.savefig("../images/runtimes_weak_scaling_" + access_pattern + ".pdf", format="pdf", bbox_inches="tight")
    plot.close()

    ## Plot scaling results
    plot.figure(dpi=200)

    plot.plot(core_ticks, parallel_two_lattice_weak_scaling, 'o-', color='k', markersize=5)
    plot.plot(core_ticks, parallel_two_lattice_framework_weak_scaling, '*-', color='y', markersize=5)
    plot.plot(core_ticks, parallel_two_step_weak_scaling, 's-', color='r', markersize=5)
    plot.plot(core_ticks, parallel_swap_weak_scaling, 'D-', color='g', markersize=5)
    plot.plot(core_ticks, parallel_shift_weak_scaling, 'p-', color='b', markersize=5)

    plot.plot(core_ticks, core_ticks, color='k', linestyle='--')

    plot.grid(color='#808080', linestyle='--', linewidth=0.5)
    plot.legend(['Non-framework two-lattice', 'Framework two-lattice', 'Two-step', 'Swap', 'Shift', 'Ideal'], loc='best')
    plot.yticks(numpy.linspace(0, core_ticks[-1], core_ticks[-1]*2 + 1))
    plot.xticks(core_ticks)
    plot.title("Weak scaling: 32x128 nodes per core, " + access_pattern + " layout")
    plot.xlabel('cores')
    plot.ylabel('speedup')
    plot.savefig("../images/weak_scaling_" + access_pattern + ".pdf", format="pdf", bbox_inches="tight")
    plot.close()

### Strong scaling
print("Creating strong scaling plots...")

# One plot for each access pattern
for access_pattern in access_patterns:

    # Initializations
    parallel_two_lattice_results = []
    parallel_two_lattice_framework_results = []
    parallel_two_step_results = []
    parallel_swap_results = []
    parallel_shift_results = []

    parallel_two_lattice_strong_scaling = []
    parallel_two_lattice_framework_strong_scaling = []
    parallel_two_step_strong_scaling = []
    parallel_swap_strong_scaling = []
    parallel_shift_strong_scaling = []

    # Assign times to corresponding list
    for line in read_content_strong:
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

    # Determine strong scaling factor based on results
    for i in range(0, len(core_ticks)):
        parallel_two_lattice_strong_scaling.append((float(parallel_two_lattice_results[0])) / float(parallel_two_lattice_results[i]))
        parallel_two_lattice_framework_strong_scaling.append((float(parallel_two_lattice_framework_results[0])) / float(parallel_two_lattice_framework_results[i]))
        parallel_two_step_strong_scaling.append((float(parallel_two_step_results[0])) / float(parallel_two_step_results[i]))
        parallel_swap_strong_scaling.append((float(parallel_swap_results[0])) / float(parallel_swap_results[i]))
        parallel_shift_strong_scaling.append((float(parallel_shift_results[0])) / float(parallel_shift_results[i]))

    ## Plot runtimes
    plot.figure(dpi=200)

    plot.plot(core_ticks, parallel_two_lattice_results, 'o-', color='k', markersize=5)
    plot.plot(core_ticks, parallel_two_lattice_framework_results, '*-', color='y', markersize=5)
    plot.plot(core_ticks, parallel_two_step_results, 's-', color='r', markersize=5)
    plot.plot(core_ticks, parallel_swap_results, 'D-', color='g', markersize=5)
    plot.plot(core_ticks, parallel_shift_results, 'p-', color='b', markersize=5)

    plot.grid(color='#808080', linestyle='--', linewidth=0.5)
    plot.legend(['Non-framework two-lattice', 'Framework two-lattice', 'Two-step', 'Swap', 'Shift'], loc='best')
    plot.xticks(core_ticks)
    plot.title("Strong scaling: 64x256 nodes, " + access_pattern + " layout")
    plot.xlabel('cores')
    plot.ylabel('time[s]')
    plot.yscale('log')
    plot.savefig("../images/runtimes_strong_scaling_" + access_pattern + ".pdf", format="pdf", bbox_inches="tight")
    plot.close()

    ## Plot scaling results
    plot.figure(dpi=200)

    plot.plot(core_ticks, parallel_two_lattice_strong_scaling, 'o-', color='k', markersize=5)
    plot.plot(core_ticks, parallel_two_lattice_framework_strong_scaling, '*-', color='y', markersize=5)
    plot.plot(core_ticks, parallel_two_step_strong_scaling, 's-', color='r', markersize=5)
    plot.plot(core_ticks, parallel_swap_strong_scaling, 'D-', color='g', markersize=5)
    plot.plot(core_ticks, parallel_shift_strong_scaling, 'p-', color='b', markersize=5)

    plot.plot(core_ticks, core_ticks, color='k', linestyle='--')

    plot.grid(color='#808080', linestyle='--', linewidth=0.5)
    plot.legend(['Non-framework two-lattice', 'Framework two-lattice', 'Two-step', 'Swap', 'Shift', 'Ideal'], loc='best')
    plot.yticks(numpy.linspace(0, core_ticks[-1], core_ticks[-1]*2 + 1))
    plot.xticks(core_ticks)
    plot.title("Strong scaling: 64x256 nodes, " + access_pattern + " layout")
    plot.xlabel('cores')
    plot.ylabel('speedup')
    plot.savefig("../images/strong_scaling_" + access_pattern + ".pdf", format="pdf", bbox_inches="tight")
    plot.close()