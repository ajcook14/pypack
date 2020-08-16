import sys

try:
    
    import cPickle
    
except ModuleNotFoundError:
    
    import pickle as cPickle

from angle_structure import AngleStructure

from graph_operations import find_faces



f = open('./angle_structures/%d'%(int(sys.argv[1])), 'rb')

sort_list = []

while True:
    
    try:
        
        angle_structure = cPickle.load(f)
        
    except EOFError:
        
        break
    
    sort_list = sort_list + [angle_structure]

f.close()

def extract(angle_structure):
    
    return(angle_structure.volume)

sort_list.sort(key=extract)

count_3_4_facial = 0

f = open('./sorted_volumes/%d'%(int(sys.argv[1])), 'wb')

for i in range(len(sort_list)):
    
    angle_structure = sort_list[i]
    
    cPickle.dump((i, angle_structure), f)
    
    graph = angle_structure.graph
    
    faces, vertices, adjacency = find_faces(graph)

    """
    max_deg = max([len(f) for f in faces])
    
    if max_deg <= 4:
        
        count_3_4_facial += 1
    """
        
    
    #print("number of vertices = %d, \tvolume = %s,\tgraph = %s"%(V, str(angle_structure.volume), str(graph)))
    #print("number of vertices = %d, \tvolume = %f,\tgraph = %s"%(V, angle_structure.volume, str(graph)))
    #print("max_deg = %d, volume/v_10 = %f, graph = %s"%(max_deg, angle_structure.volume/6.02304600917, str(graph)))
    #print("max_deg = %d, volume = %f, graph = %s"%(max_deg, angle_structure.volume, str(graph)))

f.close()

#print('number of 3-4 facial polyhedra = %d, number of polyhedra = %d, proportion = %f'%(count_3_4_facial, len(sort_list), float(count_3_4_facial) / float(len(sort_list))))

#print(sort_list[0].volume)
#print(sort_list[-1].volume)





