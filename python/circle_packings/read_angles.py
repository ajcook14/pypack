import numpy as np
import sys

from math import sin

from graph_operations import find_faces, checkerboard, search

try:
    import cPickle
except ModuleNotFoundError:
    import pickle as cPickle

from angle_structure import AngleStructure
from queue import Queue

f = open('./angle_structures/%d'%(int(sys.argv[1])), 'rb')

for i in range(int(sys.argv[2])):
    try:
        angle_structure = cPickle.load(f)
        
    except EOFError:
        print('angle structure not enumerated or not that many in polyhedra with %s vertices'%(sys.argv[1]))
        
        sys.exit(1)

f.close()

faces, vertices, adjacency = find_faces(angle_structure.graph)

checkerboarding = checkerboard(angle_structure.graph, faces, vertices, adjacency)

V = len(angle_structure.graph)

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




def find_geometry(graph, faces, adjacency, tetrahedra, adjacency_triangulated, angles):
    
    centres = [''] * len(faces)  # complex numbers
    
    radii = [-1] * len(faces)
    
    edge_lengths = [[-1 for _ in range(V)] for _ in range(V)]
    
    positions = [''] * V  # positions of vertices in boundary projection (complex numbers)
    
    l = 0
    
    angled_tetrahedra = []
    
    for i in range(len(tetrahedra)):
        
        angled_tetrahedra = angled_tetrahedra + [[]]
        
        face = tetrahedra[i]
        
        for j in range(len(face)):
            
            angled_tetrahedra[-1] = angled_tetrahedra[-1] + [[]]
            
            tetrahedron = face[j]
            
            for k in range(len(tetrahedron)):
            
                angled_tetrahedra[-1][-1] = angled_tetrahedra[-1][-1] + [angles[l]]
            
                l += 1
    
    positions[0] = 'inf'
    
    u = graph[0][0]
    
    index = search(graph[u], 0)
    
    v = graph[u][index - 1]
    
    positions[u] = 0
    
    positions[v] = 1
    
    edge_lengths[u][v] = 1.0
    
    edge_lengths[v][u] = 1.0
    
    queue = Queue()
    
    queue.append(adjacency[v][u])
    
    while len(queue) > 0:  # breadth first traversal of faces
        
        #queue.print_queue()
        
        face = queue.serve()
        
        # find all lengths of edges in current face
        
        n = len(faces[face]) - 1
        
        if edge_lengths[faces[face][n]][faces[face][0]] >= 0:
            
            known_edge = [n, 0]
        
        for i in range(len(faces[face]) - 1):
            
            if edge_lengths[faces[face][i]][faces[face][i + 1]] >= 0:
                
                known_edge = [i, i + 1]
                
                break
        
        tetrahedron_index = adjacency_triangulated[faces[face][known_edge[0]]][faces[face][known_edge[1]]][1]
        
        tetrahedron = tetrahedra[face][tetrahedron_index]
        

        if tetrahedron_index == 0:
            
            # make sure the first edge is known
            
            if known_edge[0] == 1:

                if edge_lengths[faces[face][0]][faces[face][1]] < 0:

                    index_beta = search(tetrahedron, faces[face][2])

                    beta = angled_tetrahedra[face][tetrahedron_index][index_beta]

                    index_alpha = search(tetrahedron, faces[face][0])

                    alpha = angled_tetrahedra[face][tetrahedron_index][index_alpha]

                    length = edge_lengths[faces[face][known_edge[0]]][faces[face][known_edge[1]]] * sin(beta) / sin(alpha)
                    
                    edge_lengths[faces[face][0]][faces[face][1]] = length
                    
                    edge_lengths[faces[face][1]][faces[face][0]] = length
             
            if known_edge[0] == 2:  # face is triangular

                if edge_lengths[faces[face][0]][faces[face][1]] < 0:

                    index_beta = search(tetrahedron, faces[face][2])

                    beta = angled_tetrahedra[face][tetrahedron_index][index_beta]

                    index_alpha = search(tetrahedron, faces[face][1])

                    alpha = angled_tetrahedra[face][tetrahedron_index][index_alpha]

                    length = edge_lengths[faces[face][known_edge[0]]][faces[face][known_edge[1]]] * sin(beta) / sin(alpha)
                    
                    edge_lengths[faces[face][0]][faces[face][1]] = length
                    
                    edge_lengths[faces[face][1]][faces[face][0]] = length
                           
            # go clockwise only
        
            first = True
            
            while adjacency_triangulated[faces[face][tetrahedron_index + 1]][faces[face][0]][2] == 0 or first:
                
                tetrahedron = tetrahedra[face][tetrahedron_index]
                
                index_alpha = search(tetrahedron, faces[face][tetrahedron_index + 2])
                
                alpha = angled_tetrahedra[face][tetrahedron_index][index_alpha]
                
                if edge_lengths[faces[face][tetrahedron_index + 1]][faces[face][tetrahedron_index + 2]] < 0:
                    
                    index_beta = search(tetrahedron, faces[face][0])
                    
                    beta = angled_tetrahedra[face][tetrahedron_index][index_beta]
                    
                    length = edge_lengths[faces[face][0]][faces[face][tetrahedron_index + 1]] * sin(beta) / sin(alpha)
                    
                    edge_lengths[faces[face][tetrahedron_index + 1]][faces[face][tetrahedron_index + 2]] = length
                    
                    edge_lengths[faces[face][tetrahedron_index + 2]][faces[face][tetrahedron_index + 1]] = length
                    
                if edge_lengths[faces[face][0]][faces[face][tetrahedron_index + 2]] < 0:
                    
                    index_beta = search(tetrahedron, faces[face][tetrahedron_index + 1])

                    beta = angled_tetrahedra[face][tetrahedron_index][index_beta]

                    length = edge_lengths[faces[face][0]][faces[face][tetrahedron_index + 1]] * sin(beta) / sin(alpha)

                    edge_lengths[faces[face][0]][faces[face][tetrahedron_index + 2]] = length

                    edge_lengths[faces[face][tetrahedron_index + 2]][faces[face][0]] = length
                
                
                tetrahedron_index += 1
                
                first = False
                
        elif tetrahedron_index == len(tetrahedra[face]) - 1:
            
            n = len(faces[face])
        
            if known_edge[0] == n - 2:
                
                if edge_lengths[faces[face][-1]][faces[face][0]] < 0:  # make sure last edge is known
                    
                    index_alpha = search(tetrahedron, faces[face][0])
                    
                    alpha = angled_tetrahedra[face][tetrahedron_index][index_alpha]
                    
                    index_beta = search(tetrahedron, faces[face][n - 2])
                    
                    beta = angled_tetrahedra[face][tetrahedron_index][index_beta]
                    
                    length = edge_lengths[faces[face][n - 2]][faces[face][n - 1]] * sin(beta) / sin(alpha)
                    
                    edge_lengths[faces[face][n - 1]][faces[face][0]] = length
                    
                    edge_lengths[faces[face][0]][faces[face][n - 1]] = length
            
            # go anticlockwise only

            first = True

            while adjacency_triangulated[faces[face][0]][faces[face][tetrahedron_index + 2]][2] == 0 or first:
                
                tetrahedron = tetrahedra[face][tetrahedron_index]
                
                index_alpha = search(tetrahedron, faces[face][tetrahedron_index + 1])
                
                alpha = angled_tetrahedra[face][tetrahedron_index][index_alpha]
                
                if edge_lengths[faces[face][tetrahedron_index + 1]][faces[face][tetrahedron_index + 2]] < 0:
                    
                    index_beta = search(tetrahedron, faces[face][0])
                    
                    beta = angled_tetrahedra[face][tetrahedron_index][index_beta]
                    
                    length = edge_lengths[faces[face][0]][faces[face][tetrahedron_index + 2]] * sin(beta) / sin(alpha)
                    
                    edge_lengths[faces[face][tetrahedron_index + 1]][faces[face][tetrahedron_index + 2]] = length
                    
                    edge_lengths[faces[face][tetrahedron_index + 2]][faces[face][tetrahedron_index + 1]] = length
                    
                if edge_lengths[faces[face][0]][faces[face][tetrahedron_index + 1]] < 0:
                    
                    index_beta = search(tetrahedron, faces[face][tetrahedron_index + 2])

                    beta = angled_tetrahedra[face][tetrahedron_index][index_beta]

                    length = edge_lengths[faces[face][0]][faces[face][tetrahedron_index + 2]] * sin(beta) / sin(alpha)

                    edge_lengths[faces[face][0]][faces[face][tetrahedron_index + 1]] = length

                    edge_lengths[faces[face][tetrahedron_index + 1]][faces[face][0]] = length
                
                
                tetrahedron_index -= 1
                
                first = False
                
                
        else:
            
            # go both clockwise and anticlockwise
            
            initial = tetrahedron_index
            
            # go clockwise
            
            if edge_lengths[faces[face][0]][faces[face][tetrahedron_index + 1]] < 0:
                
                index_alpha = search(tetrahedron, faces[face][0])
                
                alpha = angled_tetrahedra[face][tetrahedron_index][index_alpha]
                
                index_beta = search(tetrahedron, faces[face][tetrahedron_index + 2])
                
                beta = angled_tetrahedra[face][tetrahedron_index][index_beta]
                
                length = edge_lengths[faces[face][tetrahedron_index + 2]][faces[face][tetrahedron_index + 1]] * sin(beta) / sin(alpha)

                edge_lengths[faces[face][0]][faces[face][tetrahedron_index + 1]] = length

                edge_lengths[faces[face][tetrahedron_index + 1]][faces[face][0]] = length
                
        
            while adjacency_triangulated[faces[face][tetrahedron_index + 1]][faces[face][0]][2] == 0:
                
                tetrahedron = tetrahedra[face][tetrahedron_index]
                
                index_alpha = search(tetrahedron, faces[face][tetrahedron_index + 2])
                
                alpha = angled_tetrahedra[face][tetrahedron_index][index_alpha]
                
                if edge_lengths[faces[face][tetrahedron_index + 1]][faces[face][tetrahedron_index + 2]] < 0:
                    
                    index_beta = search(tetrahedron, faces[face][0])
                    
                    beta = angled_tetrahedra[face][tetrahedron_index][index_beta]
                    
                    length = edge_lengths[faces[face][0]][faces[face][tetrahedron_index + 1]] * sin(beta) / sin(alpha)
                    
                    edge_lengths[faces[face][tetrahedron_index + 1]][faces[face][tetrahedron_index + 2]] = length
                    
                    edge_lengths[faces[face][tetrahedron_index + 2]][faces[face][tetrahedron_index + 1]] = length
                    
                if edge_lengths[faces[face][0]][faces[face][tetrahedron_index + 2]] < 0:
                    
                    index_beta = search(tetrahedron, faces[face][tetrahedron_index + 1])

                    beta = angled_tetrahedra[face][tetrahedron_index][index_beta]

                    length = edge_lengths[faces[face][0]][faces[face][tetrahedron_index + 1]] * sin(beta) / sin(alpha)

                    edge_lengths[faces[face][0]][faces[face][tetrahedron_index + 2]] = length

                    edge_lengths[faces[face][tetrahedron_index + 2]][faces[face][0]] = length
                
                
                tetrahedron_index += 1
                
            # go anticlockwise
            
            tetrahedron_index = initial
                
            if edge_lengths[faces[face][0]][faces[face][tetrahedron_index + 2]] < 0:
                
                index_alpha = search(tetrahedron, faces[face][0])
                
                alpha = angled_tetrahedra[face][tetrahedron_index][index_alpha]
                
                index_beta = search(tetrahedron, faces[face][tetrahedron_index + 1])
                
                beta = angled_tetrahedra[face][tetrahedron_index][index_beta]
                
                length = edge_lengths[faces[face][tetrahedron_index + 2]][faces[face][tetrahedron_index + 1]] * sin(beta) / sin(alpha)

                edge_lengths[faces[face][0]][faces[face][tetrahedron_index + 2]] = length

                edge_lengths[faces[face][tetrahedron_index + 2]][faces[face][0]] = length
                

            while adjacency_triangulated[faces[face][0]][faces[face][tetrahedron_index + 2]][2] == 0:
                
                tetrahedron = tetrahedra[face][tetrahedron_index]
        
                index_alpha = search(tetrahedron, faces[face][tetrahedron_index + 1])
                
                alpha = angled_tetrahedra[face][tetrahedron_index][index_alpha]
                
                if edge_lengths[faces[face][tetrahedron_index + 1]][faces[face][tetrahedron_index + 2]] < 0:
                    
                    index_beta = search(tetrahedron, faces[face][0])
                    
                    beta = angled_tetrahedra[face][tetrahedron_index][index_beta]
                    
                    length = edge_lengths[faces[face][0]][faces[face][tetrahedron_index + 2]] * sin(beta) / sin(alpha)
                    
                    edge_lengths[faces[face][tetrahedron_index + 1]][faces[face][tetrahedron_index + 2]] = length
                    
                    edge_lengths[faces[face][tetrahedron_index + 2]][faces[face][tetrahedron_index + 1]] = length
                    
                if edge_lengths[faces[face][0]][faces[face][tetrahedron_index + 1]] < 0:
                    
                    index_beta = search(tetrahedron, faces[face][tetrahedron_index + 2])

                    beta = angled_tetrahedra[face][tetrahedron_index][index_beta]

                    length = edge_lengths[faces[face][0]][faces[face][tetrahedron_index + 2]] * sin(beta) / sin(alpha)

                    edge_lengths[faces[face][0]][faces[face][tetrahedron_index + 1]] = length

                    edge_lengths[faces[face][tetrahedron_index + 1]][faces[face][0]] = length
                
                
                tetrahedron_index -= 1
                
        
        # find centre and radius of the corresponding circle
        
        radii[face] = 1
        
        # find the positions of the vertices around the current face
        
        for i in range(len(faces[face]) - 1):
            
            neighbor = adjacency[faces[face][i + 1]][faces[face][i]]
            
            if 0 in faces[neighbor]:  # neighbor is adjacent to vertex at infinity
                
                pass
            
            elif queue.contains(neighbor):  # face appended, but not yet processed
                
                pass
            
            elif radii[neighbor] < 0:  # neighbor has not yet been traversed
                
                queue.append(neighbor)

    return(edge_lengths)

adjacency = adjacency.tolist()

print('graph = %s'%(angle_structure.graph))
print('faces = %s'%(faces))
print('tetrahedra = %s'%(tetrahedra))

edge_lengths = find_geometry(angle_structure.graph, faces, adjacency, tetrahedra, adjacency_triangulated, angle_structure.angles)

print('lengths = %s'%(edge_lengths))

print('success = %r'%(angle_structure.success))
print('volume = %f'%(angle_structure.volume))
print('angles = %s'%(str(angle_structure.angles)))
print('read time = %f seconds'%(angle_structure.read_time))
print('generation time = %f seconds'%(angle_structure.gen_time))
print('optimization time = %f seconds'%(angle_structure.opt_time))














