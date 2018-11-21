class NiiObject:
    def __init__(self):
        self.file = None
        self.reader = None
        self.extent = ()
        self.labels = []
        self.image_mapper = None
        self.scalar_range = None
