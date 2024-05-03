# Imports
import csv
import math
import matplotlib
import numpy
import matplotlib.pyplot as plot
from scipy.stats import chi2

# Matplotlib configurations
matplotlib.rcParams['text.usetex'] = True
matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'
matplotlib.rcParams.update({'font.size': 12})
matplotlib.rcParams.update({'errorbar.capsize': 5})
matplotlib.rcParams['figure.dpi'] = 300

# Global definitions
ACCESS_PATTERNS = ["collision", "stream", "bundle"]
ALGORITHMS = ["two_lattice", "two_lattice_framework", "two_step", "swap", "shift"]
CONTENTS = ["averages", "lower_bounds", "upper_bounds"]
ACCESS_PATTERNS = ["collision", "stream", "bundle"]
ALGORITHMS = ["two_lattice", "two_lattice_framework", "two_step", "swap", "shift"]
COLORS = ['k', 'y', 'r', 'g', 'b']
MARKERS = ['o', '*', 's', 'D', 'p']


def confidence_error(standard_deviation: float, n: int, confidence:float=0.95) -> list[float]:
    """
    Determines the confidence error for the specified confidence, standard deviation and number of samples.

    Args:
        standard_deviation (float): Standard deviation of the test
        n (int): number of samples
        confidence (float, optional): The confidence used. Defaults to 0.95.

    Returns:
        list[float]: A list containing the lower error bound and the upper error bound in this order.
    """

    alpha = 1.0 - confidence
    lower_limit = standard_deviation * numpy.sqrt(n / chi2.ppf(1 - (alpha / 2), df=n))
    upper_limit = standard_deviation * numpy.sqrt(n / chi2.ppf(alpha / 2, df=n))
    return [lower_limit, upper_limit]



def determine_algorithm(algorithm_name: str) -> list[str]:
    """
    Determines to which data set the algorithm with the given name belongs.
    This is necessary to group the results of the corresponding sequential and parallel algorithms together.
    Notice that the sequential two-lattice algorithm contributes to two data sets, 
    namely "two_lattice" and "two_lattice_framework".

    Returns:
        list[str]: A list containing the string representations of the data sets this algorithm contributes to.
    """

    if algorithm_name == "parallel_two_lattice":
        return ["two_lattice"]
    elif line[0] == "parallel_two_lattice_framework":
        return ["two_lattice_framework"]
    elif line[0] == "parallel_two_step" or line[0] == "sequential_two_step":
        return ["two_step"]
    elif line[0] == "parallel_swap" or line[0] == "sequential_swap":
        return ["swap"]
    elif line[0] == "parallel_shift" or line[0] == "sequential_shift":
        return ["shift"]
    elif line[0] == "sequential_two_lattice":
        return ["two_lattice", "two_lattice_framework"]
    else:
        print("Unknown algorithm name:", algorithm_name)
    

def extract_runtimes_from_line(line: list[str]) -> list[float]:
    """
    Extracts all runtimes in a given line and returns them.

    Args:
        line (list[str]): The line in question

    Returns:
        list[float]: A list containing all runtimes extracted from the given line
    """
    result = []
    for i in range(3, len(line)):
        result.append(float(line[i]))
    return result


def set_result(
        results: numpy.ndarray, 
        algorithm: str, 
        access_pattern: str, 
        content: str, 
        core_tick: list[int], 
        value: float
    ) -> None:
    """
    Sets a value at the specified result data structure.
    The indices are determined using the remaining parameters.

    Args:
        results (numpy.ndarray): the data structure storing the runtimes
        algorithm (str): string representation of an algorithm family
        access_pattern (str): string representation of the access pattern used
        content (str): either "averages" or "lower_bounds" or "upper_bounds"
        core_tick (list[int]): a list containing the surveyed amount of cores in ascending order
        value (float): the value that is to be set
    """

    access_pattern_index = 0
    algorithm_index = 0
    content_index = 0

    for i in range(0, len(ALGORITHMS)):
        if algorithm == ALGORITHMS[i]:
            algorithm_index = i
            break

    for i in range(0, len(ACCESS_PATTERNS)):
        if access_pattern == ACCESS_PATTERNS[i]:
            access_pattern_index = i
            break
            
    for i in range(0, len(CONTENTS)):
        if content == CONTENTS[i]:
            content_index = i
            break

    results[algorithm_index][access_pattern_index][content_index][core_tick] = value


def get_result(
        results: numpy.ndarray, 
        algorithm: str, 
        access_pattern: str, 
        content: str, 
        core_tick: list[int], 
    ) -> float:
    """
    Returns a value read from the specified result data structure.
    The indices are determined using the remaining parameters.

    Args:
        results (numpy.ndarray): the data structure storing the runtimes
        algorithm (str): string representation of an algorithm family
        access_pattern (str): string representation of the access pattern used
        content (str): either "averages" or "lower_bounds" or "upper_bounds"
        core_tick (list[int]): a list containing the surveyed amount of cores in ascending order

     Returns:
        float: The content of the respective test run
    """

    access_pattern_index = 0
    algorithm_index = 0
    content_index = 0

    for i in range(0, len(ALGORITHMS)):
        if algorithm == ALGORITHMS[i]:
            algorithm_index = i
            break

    for i in range(0, len(ACCESS_PATTERNS)):
        if access_pattern == ACCESS_PATTERNS[i]:
            access_pattern_index = i
            break
            
    for i in range(0, len(CONTENTS)):
        if content == CONTENTS[i]:
            content_index = i
            break

    return results[algorithm_index][access_pattern_index][content_index][core_tick]


def set_scaling_result(scaling_results: numpy.ndarray, algorithm: str, core_tick: list[int], value: float) -> None:
    """
    Sets a value at the specified data structure containing the scaling results.
    The indices are determined using the remaining parameters.

    Args:
        scaling_results (numpy.ndarray): the data structure storing the scaling factors
        algorithm (str): string representation of an algorithm family
        core_tick (list[int]): a list containing the surveyed amount of cores in ascending order
        value (float): the value that is to be set
    """

    algorithm_index = 0

    for i in range(0, len(ALGORITHMS)):
        if algorithm == ALGORITHMS[i]:
            algorithm_index = i
            break
    
    scaling_results[algorithm_index][core_tick] = value


def get_scaling_result(scaling_results: numpy.ndarray, algorithm: str, core_tick: list[int]) -> float:
    """
    Returns a value read from the specified result data structure.
    The indices are determined using the remaining parameters.

    Args:
        scaling_results (numpy.ndarray): the data structure storing the scaling factors
        algorithm (str): string representation of an algorithm family
        core_tick (list[int]): a list containing the surveyed amount of cores in ascending order

     Returns:
        float: The content of the respective test run
    """

    algorithm_index = 0

    for i in range(0, len(ALGORITHMS)):
        if algorithm == ALGORITHMS[i]:
            algorithm_index = i
    
    return scaling_results[algorithm_index][core_tick]

'''
main
'''
if __name__ == "__main__":

    # Initializations
    core_ticks = []
    read_content_weak = []
    read_content_strong = []

    # Store content of file with weak scaling results
    with open("../build/weak_scaling_results.csv", "r") as f:
        reader = csv.reader(f, delimiter=",")
        for i, line in enumerate(reader):
            read_content_weak.append(line)

    # Store content of file with strong scaling results   
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
    print("Core ticks are: ", core_ticks)

    ### Weak scaling
    print("Creating weak scaling plots...")

    runtimes = []
    lower_bound = 0
    upper_bound = 0
    average = 0
    standard_deviation = 0
    algorithm = ""
    results = numpy.zeros(shape=(5,3,3,len(core_ticks)), dtype=float)

    ## Assign times to corresponding data set
    for i in range(1, len(read_content_weak)):
        line = read_content_weak[i]
        algorithm = determine_algorithm(line[0])
        length = len(line)
        runtimes = extract_runtimes_from_line(line)
        average = numpy.average(runtimes)

        for current_algorithm in algorithm:
            set_result(results, current_algorithm, line[1], "averages", int(math.log2(int(line[2]))), average)

        standard_deviation = numpy.std(runtimes)
        lower, upper = confidence_error(standard_deviation, len(runtimes))

        for current_algorithm in algorithm:
            set_result(results, current_algorithm, line[1], "lower_bounds", core_ticks.index(int(line[2])), lower)
            set_result(results, current_algorithm, line[1], "upper_bounds", core_ticks.index(int(line[2])), upper)

    ## Determine scaling factors and plot results
    for access_pattern in ACCESS_PATTERNS:

        # Initializations
        scaling_results = numpy.zeros(shape=(5,len(core_ticks)), dtype=float)

        # Determine weak scaling factor based on results
        for i in range(0, len(core_ticks)):
            core_count = core_ticks[i]

            for current_algo in ALGORITHMS: 
                set_scaling_result(scaling_results, current_algo, i, 
                                   (core_count * get_result(results, current_algo, access_pattern, "averages", 0)) / 
                                   get_result(results, current_algo, access_pattern, "averages", i)) 

        # Plot runtime results
        fig, ax = plot.subplots()

        for i in range(0,5):
            y_values = [get_result(results, ALGORITHMS[i], access_pattern, "average", j) for j in range(0, len(core_ticks))]
            error_bounds = [[get_result(results, ALGORITHMS[i], access_pattern, content, j) 
                             for j in range(0, len(core_ticks))] 
                             for content in ["lower_bounds", "upper_bounds"]]
            ax.errorbar(core_ticks, y_values, yerr=error_bounds, capsize=5, elinewidth=2, markeredgewidth=2, color=COLORS[i])

        plot.grid(color='#808080', linestyle='--', linewidth=0.5)
        plot.legend(['Non-framework two-lattice', 'Framework two-lattice', 'Two-step', 'Swap', 'Shift'], loc='best')
        plot.title(r"Weak scaling: $\displaystyle 32 \times 128$ nodes per core, runtimes with $\displaystyle\gamma=95\%$, " + access_pattern + " layout")
        plot.xlabel('cores')
        plot.ylabel('time[s]')
        ax.set_xscale('log', base=2)
        ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
        plot.xticks(core_ticks)
        plot.savefig("../images/runtimes_weak_scaling_" + access_pattern + ".pdf", format="pdf", bbox_inches="tight")
        plot.close()

        # Plot scaling results
        fig, ax = plot.subplots()

        for i in range(0,5):
            y_values = [get_scaling_result(scaling_results, ALGORITHMS[i], j) for j in range(0, len(core_ticks))]
            plot.plot(core_ticks, y_values, color=COLORS[i], linestyle='-', marker=MARKERS[i])

        plot.plot(core_ticks, core_ticks, color='k', linestyle='--')

        plot.grid(color='#808080', linestyle='--', linewidth=0.5)
        plot.legend(['Non-framework two-lattice', 'Framework two-lattice', 'Two-step', 'Swap', 'Shift', 'Ideal'], loc='best')
        plot.title(r"Weak scaling: $\displaystyle 32 \times 128$ nodes per core, " + access_pattern + " layout")
        plot.xlabel('cores')
        plot.ylabel('speedup')
        ax.set_xscale('log', base=2)
        ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
        plot.xticks(core_ticks)
        plot.savefig("../images/weak_scaling_" + access_pattern + ".pdf", format="pdf", bbox_inches="tight")
        plot.close()

    print("Creating strong scaling plots...")

    runtimes = []
    lower_bound = 0
    upper_bound = 0
    average = 0
    standard_deviation = 0
    algorithm = ""
    results = numpy.zeros(shape=(5,3,3,len(core_ticks)), dtype=float)

    # Assign times to corresponding data set
    for i in range(1, len(read_content_strong)):
        line = read_content_strong[i]
        algorithm = determine_algorithm(line[0])
        length = len(line)
        runtimes = extract_runtimes_from_line(line)
        average = numpy.average(runtimes)

        for current_algorithm in algorithm:
            set_result(results, current_algorithm, line[1], "averages", core_ticks.index(int(line[2])), average)

        standard_deviation = numpy.std(runtimes)
        lower, upper = confidence_error(standard_deviation, len(runtimes))

        for current_algorithm in algorithm:
            set_result(results, current_algorithm, line[1], "lower_bounds", core_ticks.index(int(line[2])), lower)
            set_result(results, current_algorithm, line[1], "upper_bounds", core_ticks.index(int(line[2])), upper)

   ## Determine scaling factors and plot results
    for access_pattern in ACCESS_PATTERNS:

        # Initializations
        scaling_results = numpy.zeros(shape=(5,len(core_ticks)), dtype=float)

        # Determine strong scaling factor based on results
        for i in range(0, len(core_ticks)):
            core_count = core_ticks[i]

            for current_algo in ALGORITHMS: 
                set_scaling_result(scaling_results, current_algo, i, 
                                   (get_result(results, current_algo, access_pattern, "averages", 0)) / 
                                   get_result(results, current_algo, access_pattern, "averages", i)) 

        # Plot runtime results
        fig, ax = plot.subplots()

        for i in range(0,5):
            y_values = [get_result(results, ALGORITHMS[i], access_pattern, "average", j) for j in range(0, len(core_ticks))]
            error_bounds = [[get_result(results, ALGORITHMS[i], access_pattern, content, j) 
                             for j in range(0, len(core_ticks))] 
                             for content in ["lower_bounds", "upper_bounds"]]
            ax.errorbar(core_ticks, y_values, yerr=error_bounds, capsize=5, elinewidth=2, markeredgewidth=2, color=COLORS[i])

        plot.grid(color='#808080', linestyle='--', linewidth=0.5)
        plot.legend(['Non-framework two-lattice', 'Framework two-lattice', 'Two-step', 'Swap', 'Shift'], loc='best')
        plot.title(r"Strong scaling: $\displaystyle 64 \times 256$ nodes, runtimes with $\displaystyle\gamma=95\%$, " + access_pattern + " layout")
        plot.xlabel('cores')
        plot.ylabel('time[s]')
        ax.set_xscale('log', base=2)
        ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
        ax.set_xticks(core_ticks)
        plot.yscale('log')
        plot.savefig("../images/runtimes_strong_scaling_" + access_pattern + ".pdf", format="pdf", bbox_inches="tight")
        plot.close()

        # Plot scaling results
        fig, ax = plot.subplots()

        for i in range(0,5):
            y_values = [get_scaling_result(scaling_results, ALGORITHMS[i], j) for j in range(0, len(core_ticks))]
            plot.plot(core_ticks, y_values, color=COLORS[i], linestyle='-', marker=MARKERS[i])

        plot.plot(core_ticks, core_ticks, color='k', linestyle='--')

        plot.grid(color='#808080', linestyle='--', linewidth=0.5)
        plot.legend(['Non-framework two-lattice', 'Framework two-lattice', 'Two-step', 'Swap', 'Shift', 'Ideal'], loc='best')
        plot.title(r"Strong scaling: $\displaystyle 64 \times 256$ nodes, " + access_pattern + " layout")
        plot.xlabel('cores')
        plot.ylabel('speedup')
        ax.set_xscale('log', base=2)
        ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
        plot.xticks(core_ticks)
        plot.savefig("../images/strong_scaling_" + access_pattern + ".pdf", format="pdf", bbox_inches="tight")
        plot.close()

    print("All done.")

'''
/main
'''