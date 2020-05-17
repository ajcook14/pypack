import sys

try:
    import cPickle
except ModuleNotFoundError:
    import pickle as cPickle

from angle_structure import AngleStructure

f = open('./angle_structures/%d'%(int(sys.argv[1])), 'rb')

for i in range(int(sys.argv[2])):
    try:
        angle_structure = cPickle.load(f)
        
    except EOFError:
        print('angle structure not enumerated or not that many in polyhedra with %s vertices'%(sys.argv[1]))
        
        sys.exit(1)

f.close()

print('graph = %s'%(angle_structure.graph))
print('success = %r'%(angle_structure.success))
print('volume = %f'%(angle_structure.volume))
print('angles = %s'%(str(angle_structure.angles)))
print('read time = %f seconds'%(angle_structure.read_time))
print('generation time = %f seconds'%(angle_structure.gen_time))
print('optimization time = %f seconds'%(angle_structure.opt_time))