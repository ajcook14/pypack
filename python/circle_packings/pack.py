from scipy import optimize
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.lines as lines
import copy
import systems
from scipy import integrate
import cPickle
import sys

# Care has been taken here not to confuse the namespaces between cmath and math.
# Take care when modifying these headings.
import cmath
from math import pi, cos, asin, sin, log

transformation = [] # global variable

width = 3   # range of guesses for starting points in root finding

def compress(spinors):
    '''
    spinors: list/ndarray of real numbers, where every 4 real numbers specify a spinor.
    Returns: a list of pairs of complex numbers, where every pair corresponds to 4 numbers in spinors.
    '''
    n = len(spinors)/4  # number of spinors

    compressed = [[]]*n

    for i in range(n):
        k = 4 * i

        compressed[i] = [spinors[k] + 1j*spinors[k+1],spinors[k+2] + 1j * spinors[k+3]]

    return(compressed)
    
def convert(compressed):
    '''
    compressed: list/ndarray of pairs of complex numbers, where every pair corresponds to a spinor.
    Returns: a list of pairs [z, w], where z is a position and w is a direction.
    '''
    for i in range(len(compressed)):
        compressed[i][0] = compressed[i][0]/compressed[i][1]
        compressed[i][1] = 1/(compressed[i][1]**2)

    return(compressed)

def real(z):
    return([z.real, z.imag])

def circle(pair):
    '''
    pair: list of the form [[x + iy, u + iv], [a + ib, c + id]], where x + iy and a + ib are positions,
    u + iv and c + id are directions.
    Returns: centre and radius of the unique circle, in the form [(x, y), r].
    '''
    p = pair[0][0]
    w = pair[0][1]
    q = pair[1][0]
    v = pair[1][1]

    alpha = cmath.phase((p - q)/v) - pi/2

    r_sgn = abs(p - q)/(2 * cos(alpha)) # signed radius

    to_c = r_sgn * (1j * v) / abs(v)

    c = tuple(real(q + to_c))

    r = abs(r_sgn)

    return([c, r])

def f(x):
    global transformation

    # tack on the pre-determined spinor coordinates from fixing three points
    x = list([-x[0], -x[1], x[0], x[1], -x[3], x[2], x[2], x[3], x[4], x[5], x[4], x[5]]) + list(x[6:])

    kappa = compress(x)
    
    return(systems.system(kappa, transformation))

def create_circles(solution, blue_centres, blue_radii, red_centres, red_radii, blue_faces, red_faces):
    # find centres and radii of circles
    # blue circles
    for face in blue_faces:
        circ = circle([solution[face[0]], solution[face[1]]])

        blue_centres += [circ[0]]

        blue_radii += [circ[1]]

    # red circles
    for face in red_faces:
        pair = copy.deepcopy([solution[face[0]], solution[face[1]]])

        pair[0][1] *= 1j # rotate each spinor in pair by pi/2 first
        pair[1][1] *= 1j

        circ = circle(pair)

        red_centres += [circ[0]]

        red_radii += [circ[1]]

def lobachevsky(theta):
    result = integrate.quad(lambda u: log(abs(2*sin(u))), 0, theta)
    return(-result[0])

class Polyhedron():
    def __init__(self, graph, checkerboarding, faces, solution, volume):
        self.graph = graph
        self.checkerboarding = checkerboarding
        self.faces = faces
        self.solution = solution
        self.volume = volume
        
def main():
    global transformation

    g = open('graphs', 'rb')

    pol = int(sys.argv[1]) # number of polyhedra to enumerate

    for iteration in range(pol):
        print('finding polyhedron = ' + str(iteration))

        polyhedron = open('./polyhedra/' + str(iteration), 'w')

        graph = systems.read_file(g)

        transformation, faces, checkerboarding = systems.generate(graph)

        equns = len(transformation)

        blue_faces = []
        red_faces = []
        for i in range(len(faces)):
            if 0 in faces[i]: # don't include faces adjacent to the vertex at infinity. Note that it is guaranteed that the spinor with index zero will get taken to infinity.
                continue

            if checkerboarding[i] == 0:
                blue_faces = blue_faces + [faces[i]] # faces must be clockwise oriented and not be adjacent to 0
            else:
                red_faces = red_faces + [faces[i]] # faces must be clockwise oriented and not be adjacent to 0

        # solve system of equations with random guesses for initial vectors
        i = 0 ############## diagnostic
        #print("Number of initial guesses made:")
        while True:
            guess = 2 * width * np.random.rand(equns) - width * np.ones(equns)
            result = optimize.root(f,guess)
            if result.success:
                break
            #if i % 50 == 0:
                #print(i) ############## diagnostic
            i = i + 1 ############## diagnostic
        #print(i)

        zeros = result.x

        zeros = list([-zeros[0], -zeros[1], zeros[0], zeros[1], -zeros[3], zeros[2], zeros[2], zeros[3], zeros[4], zeros[5], zeros[4], zeros[5]]) + list(zeros[6:])

        spinors = compress(zeros)

        # NOTE: Mobius transformation must have determinant 1!
        phi = [[1, -1j],[0.5 - 0.5j, 0.5 - 0.5j]]
        for i in range(1,len(spinors)): # it is unnecessary and computationally infeasible to actually take the 0th spinor to infinity
            spinors[i] = np.dot(phi, spinors[i])

        solution = convert(spinors)

        # find centres and radii
        blue_centres = []
        blue_radii   = []

        red_centres = []
        red_radii   = []

        create_circles(solution, blue_centres, blue_radii, red_centres, red_radii, blue_faces, red_faces)

        # compute volume
        zetas = []

        for i in range(len(blue_faces)):
            face = blue_faces[i]

            centre = blue_centres[i][0] + 1j * blue_centres[i][1]

            for j in range(len(face)):
                zetas += [0.5 * cmath.phase( (solution[face[j - 1]][0] - centre) / (solution[face[j]][0] - centre) )]

        for i in range(len(red_faces)):
            face = red_faces[i]

            centre = red_centres[i][0] + 1j * red_centres[i][1]

            for j in range(len(face)):
                zetas += [0.5 * cmath.phase( (solution[face[j - 1]][0] - centre) / (solution[face[j]][0] - centre) )]


        volume = 0

        for zeta in zetas:
            volume += lobachevsky(zeta) # * 2 / 2

        volume = abs(volume) # volume is negative in the case of finding the solution with reverse orientation

        # write the full combinatorial and geometric data to a file
        p = Polyhedron(graph, checkerboarding, faces, solution, volume)

        string = cPickle.dumps(p)

        polyhedron.write(string)

        polyhedron.close()

    g.close()


if __name__ == '__main__':
    main()
