
class AngleStructure():
    def __init__(self, graph, result, read_time, gen_time, opt_time):
        self.graph = graph
        self.success = result.success
        self.volume = -result.fun
        self.angles = result.x
        self.read_time = read_time
        self.gen_time = gen_time
        self.opt_time = opt_time
