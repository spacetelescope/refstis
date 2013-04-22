import os

def RemoveIfThere(item):
    if os.path.exists(item):
        os.remove(item)
