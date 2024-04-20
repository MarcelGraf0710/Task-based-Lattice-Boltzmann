import csv

import numpy

import matplotlib.pyplot as plot

min_x = 0
min_y = 0
min_time_step = 0

max_x = 0
max_y = 0
max_time_step = 0

read_content=[]

with open("results.csv", "r") as f:
    reader = csv.reader(f, delimiter=",")
    for i, line in enumerate(reader):
        #print("line[{}] = {}".format(i, line))
        read_content.append(line)

min_time_step = int(read_content[1][0])
min_x = int(read_content[1][1])
min_y = int(read_content[1][2])


max_time_step = int(read_content[-1][0])
max_x = int(read_content[-1][1])
max_y = int(read_content[-1][2])

v_x = [numpy.array([max_x - min_x, max_y - min_y]) for i in range(0, max_time_step+1)]
v_y = [numpy.array([max_x - min_x, max_y - min_y]) for i in range(0, max_time_step+1)]
densities = [numpy.array([max_x - min_x, max_y - min_y]) for i in range(0, max_time_step+1)]

print("Got min_x = {}, max_x = {}, min_y = {}, max_y = {}, min_time_step = {}, max_time_step = {}".format(min_x, max_x, min_y, max_y, min_time_step, max_time_step))

# for line in read_content:
#     densities[int(line[0])]

#densities = [numpy.array([line[-1]]) for line in read_content]

densities = [numpy.ndarray([max_y, max_x]) for time in range(0,max_time_step+1)]
velocities_x = [numpy.ndarray([max_y, max_x]) for time in range(0,max_time_step+1)]

for i in range(1, len(read_content)):
    line = read_content[i]
    densities[int(line[0])-min_time_step][int(line[2])-min_y][int(line[1])-min_x] = line[-1]
    velocities_x[int(line[0])-min_time_step][int(line[2])-min_y][int(line[1])-min_x] = numpy.sqrt(pow(float(line[3]),2) ) #+ pow(float(line[4]),2)

print(densities[0])

plot.figure(dpi=200)
plot.title('Rosenbrock function as 2d heat map')
p = plot.imshow(velocities_x[-1], extent=[min_x, max_x, min_y, max_y], cmap='rainbow', vmin=0, vmax=0.1)
plot.colorbar(p)
plot.xticks(numpy.arange(min_x, max_x+1, step=10))
plot.yticks(numpy.arange(min_y, max_y+1, step=10))
plot.xlabel('x')
plot.ylabel('y')
plot.show()