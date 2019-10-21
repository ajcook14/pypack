import cmath
from math import pi, cos
import copy
import numpy as np

def product(kappa,omega):
    return(kappa[0]*omega[1] - kappa[1]*omega[0])

def reg(z, comp):
    '''
    z: complex variable to be projected and regularized
    comp: component of z that must be zero - comp = 'i' for imaginary, comp = 'r' for real.
    Returns: projection of z back with an added fourth order 'bump' function of the other component at zero.
    '''
    eps = 1e-1  # regularization epsilon (width of bump function)

    bump = lambda x: ((x + eps)**2) * ((x - eps)**2) / eps**4

    if comp == 'i':
        if -eps < z.real and z.real < eps:
            return(z.imag**2 + bump(z.real))
        else:
            return(z.imag**2)
    elif comp == 'r':
        if -eps < z.imag and z.imag < eps:
            return(z.real**2 + bump(z.imag))
        else:
            return(z.real**2)
    else:
        raise AttributeError("attribute comp must be 'i' or 'r'")

def rot(spinor):
    spinor = copy.deepcopy(spinor)
    z = cmath.exp(pi * 1j / 4) 

    return([z * spinor[0], z * spinor[1]])

def read_file(f):
    '''
    f: binary planar code file containing a series of graphs output from the command
    "plantri -qc4m3od <verts>" where <verts> is the number of vertices of each graph on output
    Returns: list of lists of integers, where each inner list L is a clockwise-oriented adjacency
    list of the vertex given by the index of L in the outermost list, starting at the adjacent
    vertex with the lowest index.
    '''
    while True:
        string = f.read(15)

        if string == '>>planar_code<<':
            pass # seek to the next existent graph
        elif string == '':
            return([])
        else:
            break

    f.seek(-15,1)

    V = ord(f.read(1)) # number of vertices

    faces = list(range(V + 2))

    vertices = list(range(V))

    graph = [[]] * V

    for i in range(V):
        while True:
            edge = ord(f.read(1))

            if edge != 0:
                graph[i] = graph[i] + [edge]
            else:
                break

    for i in range(V):  # set indices of vertices to start at 0 and not 1
        for j in range(len(graph[i])):
            graph[i][j] -= 1

    return(graph)

def search(a_list, target):
    for i in range(len(a_list)):
        if a_list[i] == target:
            return(i)

    return(-1)

def generate(graph):
    '''
    Usage:
    faces, checkerboarding = generate(graph)

    graph: see output of the read_file function.

    faces: A list of lists of integers representing vertices. Each inner list corresponds to
    exactly one face of graph, its elements are the adjacent vertices in clockwise order.

    checkerboarding: list whose length is the face count of graph. This list specifies a
    proper 2-colouring of the faces of graph, where the elements of checkerboarding are elements
    of the binary colour set {0, 1}, and where the index of each colour is the index of the
    corresponding face in the list faces.
    '''
    # Find faces

    V = len(graph)

    vertices = V * [[]] # each entry is a list of adjacent faces to the corresponding vertex in clockwise order

    faces = []

    adjacency = -np.ones([V,V],dtype=int) # n = -1 ==> 'not traversed', n > -1 ==> 'traversed, face = n'

    for vertex in range(V):
        for neighbor in graph[vertex]:
            if adjacency[vertex, neighbor] >= 0: # has this face been done already?
                vertices[vertex] = vertices[vertex] + [adjacency[vertex, neighbor]]

                continue

            f = len(faces) # index of face being added

            adjacency[vertex, neighbor] = f

            face = [vertex]

            prev = vertex
            current = neighbor
            while current != vertex:
                face = face + [current]

                index = search(graph[current], prev)

                prev = current

                current = graph[current][index - 1]

                assert adjacency[prev, current] == -1, "Entered a face that has already been counted. This should not be happening."

                adjacency[prev, current] = f

            faces = faces + [face]

            vertices[vertex] = vertices[vertex] + [f] # append index of face


    # Checkerboard colour the faces

    F = V + 2

    checkerboarding = F * [-1]

    stack = [0]

    covered = [0] # vertices of the graph currently covered

    failed = False # flag to indicate whether the last step was a push or a pop

    curr = 0 # top of stack.

    prev = graph[curr][0] # second from top of stack. This is an initialization trick to begin colouring without tampering with the stack.

    checkerboarding[adjacency[curr, prev]] = 0

    while len(stack) > 0:

        if not failed: # checkerboard around the current vertex
            index = search(vertices[curr], adjacency[curr, prev])

            assert index >= 0, "Failed to find starting face for checkerboarding."

            colour = checkerboarding[adjacency[curr, prev]]

            assert colour >= 0, "Face should already be coloured."

            for i in range(4):
                if checkerboarding[vertices[curr][index - i]] >= 0:
                    assert checkerboarding[vertices[curr][index - i]] == colour, "Overlapping colours do not match."
                else:
                    checkerboarding[vertices[curr][index - i]] = colour

                colour = (colour + 1) % 2

        failed = True
        for neighbor in graph[stack[-1]]:

            if search(covered, neighbor) == -1:
                covered.append(neighbor)

                prev = stack[-1]

                curr = neighbor

                stack.append(neighbor)

                failed = False

                break

        if failed:
            stack.pop()

    assert len(covered) >= len(graph), "Checkerboarding did not cover graph."
    assert len(covered) <= len(graph), "Checkerboarding covered too many vertices."


    # Choose directions of spinors

    directions = np.zeros([V,V],dtype=int) # directions of spinors. Orientations of circles induce orientations of faces in graph. Each spinor direction is represented by two 1's on the two directed edges in agreeing orientation with the two corresponding oriented circles that are tangent at the spinor.

    for i in range(V):
        if checkerboarding[adjacency[i, graph[i][0]]] == 0:
            directions[i, graph[i][1]] = 1
            directions[i, graph[i][2]] = 1
        else:
            directions[i, graph[i][0]] = 1
            directions[i, graph[i][1]] = 1


    # Generate equations

    transformation = []

    for i in range(len(faces)):
        if checkerboarding[i] == 0:
            for j in range(len(faces[i])):
                if directions[faces[i][j - 2], faces[i][j - 1]] == directions[faces[i][j - 1], faces[i][j]]:
                    transformation = transformation + [lambda kappa, k=faces[i][j - 2], l=faces[i][j - 1]: reg(product(kappa[k], kappa[l]), 'i')]
                else:
                    transformation = transformation + [lambda kappa, k=faces[i][j - 2], l=faces[i][j - 1]: reg(product(kappa[k], kappa[l]), 'r')]

            for j in range(len(faces[i]) - 3):
                if directions[faces[i][0], faces[i][1]] == directions[faces[i][j + 2], faces[i][j + 3]]:
                    transformation = transformation + [lambda kappa, k=faces[i][0], l=faces[i][j + 2]: reg(product(kappa[k], kappa[l]), 'i')]
                else:
                    transformation = transformation + [lambda kappa, k=faces[i][0], l=faces[i][j + 2]: reg(product(kappa[k], kappa[l]), 'r')]

        elif checkerboarding[i] == 1:
            for j in range(len(faces[i]) - 3):
                if directions[faces[i][0], faces[i][1]] == directions[faces[i][j + 2], faces[i][j + 3]]:
                    transformation = transformation + [lambda kappa, k=faces[i][0], l=faces[i][j + 2]: reg(product(rot(kappa[k]), rot(kappa[l])), 'i')]
                else:
                    transformation = transformation + [lambda kappa, k=faces[i][0], l=faces[i][j + 2]: reg(product(rot(kappa[k]), rot(kappa[l])), 'r')]
        else:
            assert False, "checkerboarding list has unfilled entries."

    for i in range(V):
        transformation = transformation + [lambda kappa, k=i: abs(kappa[k][1]) - 1]

    assert len(transformation) == 4 * V - 6, "Poorly posed system: Q != 4V - 6."

    return(transformation, faces, checkerboarding)

def system(kappa, transformation):
    vector = len(transformation) * [0] ### creates a new list each iteration - potential bottleneck in time complexity

    for i in range(len(transformation)):
        function = transformation[i]

        vector[i] = function(kappa)

    return(vector)

    f.close()
