import numpy as np

import sys

try:
    import cPickle
except ModuleNotFoundError:
    import pickle as cPickle

from angle_structure import AngleStructure

f = open('./angle_structures/%d'%(int(sys.argv[1])), 'rb')

sort_list = []

while True:
    try:
        angle_structure = cPickle.load(f)
        
    except EOFError:
        break
    
    sort_list = sort_list + [angle_structure]

f.close()

def search(a_list, target):
    for i in range(len(a_list)):
        if a_list[i] == target:
            return(i)

    return(-1)

def extract(angle_structure):
    return(angle_structure.volume)

sort_list.sort(key=extract)

count_3_4_facial = 0

for angle_structure in sort_list:
    
    graph = angle_structure.graph

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
            
    max_deg = max([len(f) for f in faces])
    
    if max_deg <= 4:
        
        count_3_4_facial += 1
        
    
    #print("number of vertices = %d, \tvolume = %s,\tgraph = %s"%(V, str(angle_structure.volume), str(graph)))
    #print("number of vertices = %d, \tvolume = %f,\tgraph = %s"%(V, angle_structure.volume, str(graph)))
    print("max_deg = %d, volume/v_10 = %f, graph = %s"%(max_deg, angle_structure.volume/6.02304600917, str(graph)))
    #print("max_deg = %d, volume = %f, graph = %s"%(max_deg, angle_structure.volume, str(graph)))

#print('number of 3-4 facial polyhedra = %d, number of polyhedra = %d, proportion = %f'%(count_3_4_facial, len(sort_list), float(count_3_4_facial) / float(len(sort_list))))

#print(sort_list[0].volume)
#print(sort_list[-1].volume)





