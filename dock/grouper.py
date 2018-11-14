def grouper(n, iterable, fillvalue=None):
    iterable = list(iterable)
    out = []
    for i in range(0, len(iterable), n):
        out += [iterable[i: i+n]]
    return out
