#Container for information about a Pose
class Pose:
    def __init__(self, rmsd, fp, num, gscore):
        self.rmsd = rmsd
        self.fp = fp
        self.num = num
        self.gscore = gscore

    def __str__(self):
        return ';'.join(map(str, [self.num, self.rmsd, self.fp, self.gscore]))

