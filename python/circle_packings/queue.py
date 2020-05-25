class Queue():
    
    def __init__(self):
        
        self.queue = []
    
    def __len__(self):
        
        return(len(self.queue))
    
    def print_queue(self):
        
        print(self.queue)
        
    def append(self, item):
        
        assert type(item) == type(0), "Item must be a non-negative integer."
        
        assert item >= 0, "Item must be a non-negative integer."
        
        self.queue.append(item)
        
    def serve(self):
        
        if len(self.queue) > 0:
            
            item = self.queue[0]
            
            del self.queue[0]
            
            return(item)
        
        else:
            
            return(-1)
    
    def contains(self, item):
        
        if item in self.queue:
            
            return(True)
        
        else:
            
            return(False)