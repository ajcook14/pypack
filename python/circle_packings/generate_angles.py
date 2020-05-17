#!/usr/bin/sage -python

import sage.all

#from sage.matrix.matrix_space import MatrixSpace
from sage.matrix.constructor import Matrix
from sage.modules.free_module_element import vector
#from sage.geometry.polyhedron.parent import Polyhedra
from sage.geometry.polyhedron.all import Polyhedron #base import Polyhedron_base
from sage.rings.all import QQ, RDF

from lobachevsky import lobachevsky as lb

from scipy.optimize import minimize

import copy
import numpy as np
import time
import sys

try:
    import cPickle
except ModuleNotFoundError:
    import pickle as cPickle

from angle_structure import AngleStructure


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

def generate(graph):
    '''
    graph: see output of the read_file function.

    faces: A list of lists of integers representing vertices. Each inner list corresponds to
    exactly one face of graph, its elements are the adjacent vertices in clockwise order.
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

    # Find angle structure
    
    A = 0 # number of angles
    
    tetrahedra = [] # triangulation of P. Elements are faces. Nested (inner) elements are the triangles of that face.
    
    adjacency_triangulated = [[[-1,-1,-1] for _ in range(V)] for _ in range(V)] # edges are represented by triples: [face index, tetrahedron index, already and edge? 1 else 0]
    
    for i in range(len(faces)): # triangulate P into tetrahedra
        face = faces[i]
        
        if 0 in face: # don't count vertical faces
            tetrahedra = tetrahedra + [[]]
            
            continue
        
        assert len(face) >= 3, "Face is a bigon or point! Not allowed..."
        
        triangulated = [[face[0], face[1], face[2]]]
        
        adjacency_triangulated[face[0]][face[1]] = [i,0,1]
        adjacency_triangulated[face[1]][face[2]] = [i,0,1]
        adjacency_triangulated[face[2]][face[0]] = [i,0,0]

        for j in range(2, len(face) - 1):
            triangulated = triangulated + [[face[0], face[j], face[j + 1]]]

            adjacency_triangulated[face[0]][face[j]] = [i,j-1,0]
            adjacency_triangulated[face[j]][face[j + 1]] = [i,j-1,1]
            adjacency_triangulated[face[j + 1]][face[0]] = [i,j-1,0]
            
        adjacency_triangulated[face[-1]][face[0]] = [i,len(face) - 3,1]
        
        tetrahedra = tetrahedra + [triangulated]
        
        A = A + 3 * len(triangulated)
        
    e = np.eye(A) # standard basis vectors
    
    vert_degree_excess = 0
    for face in vertices[0]:
        vert_degree_excess += len(faces[face])
    
    equations = np.zeros([(A // 3) + (len(graph) - 1) + (3 * len(faces) - vert_degree_excess - 4), A]) # type (2) + type (3) + (types (5) and (6))
    # see notes for how to derive this expression
    
    constants = np.zeros((A // 3) + (len(graph) - 1) + (3 * len(faces) - vert_degree_excess - 4)) # type (2) + type (3) + (types (5) and (6))
    
    for i in range(A // 3): # add type (2) conditions
        equations[i] = np.array(e[3*i] + e[3*i + 1] + e[3*i + 2])
        
        constants[i] = 1
            
    
    count = 0
    
    for i in range(len(tetrahedra)): # add type (3) conditions
        face = tetrahedra[i]
        
        for j in range(len(face)):
            tetrahedron = face[j]
            
            for k in range(len(tetrahedron)):
                vertex = tetrahedron[k]
                
                equations[vertex + A // 3 - 1][count + k] = 1
            
            count += 3
    
    for i in range(1,len(vertices)): # add type (3) conditions
        vertex = vertices[i]
        
        if 0 in graph[i]: # vertex on the corner of boundary projection?
            constants[i + A // 3 - 1] = 0.5
            
        else:
            constants[i + A // 3 - 1] = 2 # assume not on corner of edge of projection
            
            for face in vertex:
                if 0 in faces[face]: # vertex on edge of boundary projection?
                    constants[i + A // 3 - 1] = 1
                    
                    break
                    
    
    indices = [] # angle index of each tetrahedron in equations matrix
    
    count = 0
    
    for i in range(len(tetrahedra)):
        face = tetrahedra[i]
        
        indices = indices + [[]]
        
        for j in range(len(face)):
            indices[i] = indices[i] + [[count, count + 1, count + 2]]
            
            count += 3
    
    row = (A // 3) + (len(graph) - 1)
    
    for i in range(1, V): # add type (5) and (6) conditions
        for j in range(i + 1, V):
            triple1 = adjacency_triangulated[i][j]
            triple2 = adjacency_triangulated[j][i]
            
            if triple1[0] >= 0: # exclude vertical faces
                tetrahedron = tetrahedra[triple1[0]][triple1[1]]

                for k in range(len(tetrahedron)):
                    vertex = tetrahedron[k]

                    if vertex != i and vertex != j:
                        equations[row][indices[triple1[0]][triple1[1]][k]] = 1
                        
                        break
            
            if triple2[0] >= 0: # exclude vertical faces
                tetrahedron = tetrahedra[triple2[0]][triple2[1]]

                for k in range(len(tetrahedron)):
                    vertex = tetrahedron[k]

                    if vertex != i and vertex != j:
                        equations[row][indices[triple2[0]][triple2[1]][k]] = 1
                        
                        break
            
            if triple1[2] == 0 or triple2[2] == 0: # type (5) condition
                assert triple1[2] == triple2[2], "Edge types in adjacency_triangulated do not agree for opposing orientations."
                
                constants[row] = 1
                
                row += 1
                
            elif triple1[2] == 1 or triple2[2] == 1: # type (6) condition
                constants[row] = 0.5
                
                row += 1
            
            
    
    equations = Matrix(QQ, equations.tolist())
    
    constants = vector(QQ, constants.tolist())
    
    M = equations.augment(constants)
    
    free = M.nonpivots()
    
    freedom = len(free) - 1 # degrees of freedom
    
    M.echelonize()
    
    np_ech = M.numpy()
    
    I = np.eye(A + 1)
    
    for i in range(freedom):
        var = free[i]
        
        np_ech = np.insert(np_ech,var,-I[var,:],axis=0) # include free variables. There may be a faster way of doing this.
    
    b = np.pi * np_ech[0:A,-1]
    
    coeff = np.zeros((A, freedom), dtype='float64')
    
    for i in range(freedom):
        coeff[:,i] = -np_ech[0:A,free[i]]
    
    return(coeff, b)

def objective(t,coeff,b):
    angles = np.dot(coeff,t) + b
    
    volume = 0
    
    for a in angles:
        volume += lb(a)
        
    return(-volume)


def main():
    g = open('./graphs', 'rb')
    
    f = open('./angle_structures/%d'%(int(sys.argv[1])), 'ab')
    
    for i in range(int(sys.argv[2])):
        graph = read_file(g)
    
    try:
        iteration = int(sys.argv[2]) + 1
        
        while True:
            start = time.time()

            graph = read_file(g)

            if len(graph) == 0:
                break

            read_time = time.time() - start
            start = time.time()

            coeff, b = generate(graph)

            constraints1 = np.insert(coeff, 0, b, axis=1).tolist()
            constraints2 = np.insert(-coeff, 0, np.pi - b, axis=1)

            constraints = np.concatenate((constraints1, constraints2), axis=0).tolist() # angle constraints for free variables

            for i in range(len(constraints)):
                constraints[i] = tuple(constraints[i])

            pol = Polyhedron(ieqs=constraints)

            pol_vertices = pol.Vrepresentation()

            total = np.zeros(len(list(pol_vertices[0])))

            for i in range(len(pol_vertices)):
                total = total + np.array(list(pol_vertices[i]))

            x0 = total/float(len(pol_vertices)) # centroid of the polytope vertices

            gen_time = time.time() - start
            start = time.time()

            result = minimize(objective, x0, args=(coeff, b))

            opt_time = time.time() - start

            result.x = angles = np.dot(coeff,result.x) + b
            
            angle_structure = AngleStructure(graph, result, read_time, gen_time, opt_time)

            #########################################################################################
            ########################## Should make this transaction atomic ########################## 
            
            cPickle.dump(angle_structure, f)
            
            print('\rnumber of vertices = %d\tnumber in sequence (by Plantri output) = %d of %d'%(int(sys.argv[1]),iteration,int(sys.argv[3]))), # python 3:, end='')
            
            sys.stdout.flush()
            
            iteration += 1
            
            ########################## Should make this transaction atomic ########################## 
            #########################################################################################
            
            
    except KeyboardInterrupt:
        g.close()
        
        sys.exit(1)





if __name__ == '__main__':
    main()