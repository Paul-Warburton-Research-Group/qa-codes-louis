import pickle

# Init a working directory?

# Save object to binary file
def saveBin(obj,fname):
    fd = open(fname,"wb")
    pickle.dump(obj,fd)
    fd.close()

# Load binary file
def loadBin(fname):
    fd = open(fname,"rb")
    obj = pickle.load(fd)
    fd.close()
    return obj
