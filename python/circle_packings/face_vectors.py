import numpy as np
import sys
from scipy import integrate
from math import sin, log, pi

pol = int(sys.argv[1]) # number of polyhedra to enumerate

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

        if string.decode("utf-8") == '>>planar_code<<':
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

def find_faces(graph):

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

    return(faces)

def face_vector(faces):
    """
    return: face_vector, a list such that face_vector[i] is the number of faces of valence i in the graph represented by faces
    """
    
    F = len(faces)
    face_vector = [0] * ((F // 2 - 1) + 1)    # add 1 to make F // 2 - 1 the maximum index
    
    for face in faces:
        face_vector[len(face)] += 1
    
    return(face_vector)

def lobachevsky(theta):
    result = integrate.quad(lambda u: log(abs(2*sin(u))), 0, theta)
    return(-result[0])

def bipyramid_vol(n):
    tetrahedron = lobachevsky(2*pi/n) + 2 * lobachevsky((n - 2) * pi / (2 * n))
    return(n * tetrahedron)

def vbound(fv, pyramid_vols):
    vol = 0
    
    for i in range(len(fv)):
        vol += pyramid_vols[i] * fv[i]
        
    return(vol)

def main():
    f = open('graphs', 'rb')

    pyramid_vols = [0, 0, 0]
    fvectors = [[]]
    bounds = []
    F = 0

    for i in range(pol):
        graph = read_file(f)
        faces = find_faces(graph)

        while True:
            if len(faces) > F: # reached the end of this number of faces
                if F > 2:
                    pyramid_vols += [bipyramid_vol(F) / 2.0]
                    
                maxi = 0
                for fv in fvectors[F]:
                    if F // 2 - 1 >= 5:
                        if sum(fv[5:]) == 0:    # only count polyhedra with large faces
                            continue
                    else:
                        continue
                            
                    vol = vbound(fv, pyramid_vols)
                    
                    if vol > maxi:
                        maxi = vol
                        
                bounds += [maxi]
                fvectors += [[]]
                F += 1
            else:
                break

        fv = face_vector(faces)
        fvectors[-1] += [fv]
        
    maxi = 0
    for fv in fvectors[F]:
        if F // 2 - 1 >= 5:
            if sum(fv[5:]) == 0:    # only count polyhedra with large faces
                continue
        else:
            continue

        vol = vbound(fv, pyramid_vols)

        if vol > maxi:
            maxi = vol

    bounds += [maxi]   
    f.close()

    print(fvectors)
    print(bounds)




    
if __name__ == '__main__':
    main()



