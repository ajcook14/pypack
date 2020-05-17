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

def extract(angle_structure):
    return(angle_structure.volume)

sort_list.sort(key=extract)

print(sort_list[0].volume)
print(sort_list[-1].volume)