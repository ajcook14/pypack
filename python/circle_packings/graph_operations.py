import numpy as np



def search(a_list, target):
    for i in range(len(a_list)):
        if a_list[i] == target:
            return(i)

    return(-1)



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


def find_faces(graph):
    """
    Inputs:
    graph: see output of the read_file function.
    
    
    Outputs:
    faces: A list of lists of integers representing vertices. Each inner list corresponds to
    exactly one face of graph, its elements are the adjacent vertices in clockwise order.
    
    vertices: A list of lists, each inner list contains the indices of the faces adjacent to the
    vertex corresponding to the index in the outer list, in clockwise order.
    
    adjacency: A numpy 2-dimensional ndarray whose [i, j] entry is -1 if there is no edge from
    vertex i to vertex j, otherwise adjacency[i, j] = f, where f is the index of the face to the
    right of the oriented edge (i, j).
    """
    
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
            
    return(faces, vertices, adjacency)


def checkerboard(graph, faces, vertices, adjacency):
    """
    Inputs:
    faces, vertices, adjacency: see outputs of find_faces function.
    
    
    Outputs:
    checkerboarding: list whose length is the face count of graph. This list specifies a
    proper 2-colouring of the faces of graph, where the elements of checkerboarding are elements
    of the binary colour set {0, 1}, and where the index of each colour is the index of the
    corresponding face in the list faces.
    """
    
    V = len(graph)
    
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
    
    return(checkerboarding)