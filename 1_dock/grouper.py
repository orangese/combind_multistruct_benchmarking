import itertools

def grouper(n, iterable, fillvalue=None):
    args = [iter(iterable)] * n 
    return itertools.izip_longest(*args, fillvalue=fillvalue)

