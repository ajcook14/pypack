import cPickle
from pack import Polyhedron
from pack import create_circles
from pack import real
import systems
import matplotlib.pyplot as plt
import matplotlib.lines as lines
import sys

plot_spinors = False

p = sys.argv[1]

f = open('./polyhedra/' + p, 'r')

string = f.read()
polyhedron = cPickle.loads(string)
f.close()

graph = polyhedron.graph
faces = polyhedron.faces
checkerboarding = polyhedron.checkerboarding
solution = polyhedron.solution
print('Volume = ' + str(polyhedron.volume))

blue_faces = []
red_faces = []
for i in range(len(faces)):
    if 0 in faces[i]: # don't include faces adjacent to the vertex at infinity. Note that it is guaranteed that the spinor with index zero will get taken to infinity.
        continue

    if checkerboarding[i] == 0:
        blue_faces = blue_faces + [faces[i]] # faces must be clockwise oriented and not be adjacent to 0
    else:
        red_faces = red_faces + [faces[i]] # faces must be clockwise oriented and not be adjacent to 0

# find centres and radii
blue_centres = []
blue_radii   = []

red_centres = []
red_radii   = []

create_circles(solution, blue_centres, blue_radii, red_centres, red_radii, blue_faces, red_faces)

# plot circles
circles = []

for i in range(len(blue_centres)):
    circles += [plt.Circle(blue_centres[i], blue_radii[i], color='b', fill=False)]

for i in range(len(red_centres)):
    circles += [plt.Circle(red_centres[i], red_radii[i], color='r', fill=False)]

fig, ax = plt.subplots()

for i in range(len(circles)):
    ax.add_artist(circles[i])

# warning: this method of finding the opposite vertex for printing is unjustified.
vertex = real(solution[graph[0][systems.search(graph[0],1) - 2]][0])
ax.axhline(y=0, color='b', linewidth=1.1)
ax.axhline(y=vertex[1], color='b', linewidth=1.1)
ax.axvline(x=0, color='r', linewidth=1.1)
ax.axvline(x=vertex[0], color='r', linewidth=1.1)

if plot_spinors:
    for i in range(1,len(solution)):
        position = real(solution[i][0])
        ax.plot(position[0], position[1], 'ro')
        ax.text(position[0], position[1], str(i), fontsize=12)

ax.set_xlim(-vertex[0]/2, 1.5 * vertex[0])
ax.set_ylim((vertex[1]/2) - abs(vertex[1]),(vertex[1]/2) + abs(vertex[1]))

#fig.savefig('circlepacking.png')
plt.show()

