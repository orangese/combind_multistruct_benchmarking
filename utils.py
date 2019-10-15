import re
import os

class StringFunction:
    """
    Wraps function such that it's string representation
    can be read with the eval() function.
    """
    def __init__(self, func):
        self.str = func
        self.func = eval(func)
    def __call__(self, *args):
        return self.func(*args)
    def __str__(self):
        return self.str
    def __repr__(self):
        return self.str

def grouper(n, iterable):
	iterable = list(iterable)
	return [iterable[i: i+n] for i in range(0, len(iterable), n)]

def resolve(paths):
    keys = {k: v for k, v in paths.items()}
    for path in paths.values():
        for wild in re.findall('{(.*?)}', path):
            if wild.lower() == wild:
                keys[wild] = '{' + wild + '}'

    resolved = {}
    for name, path in paths.items():
        resolved[name] = path.format(**keys)
    if resolved == paths:
        return resolved
    return resolve(resolved)

def get_proteins(paths, exclude):
	return [p for p in os.listdir(paths['DATA'])
	        if p[0] != '.' and p not in exclude]
