SCHRODINGER_RUN = "/share/PI/rondror/software/schrodinger2017-1/run"
import subprocess

class RMSD:
    def __init__(self, data):
        self.data = data

    @classmethod
    def read(cls, line):
        return cls(map(float, line.strip().split(',')))
    
    @classmethod
    def create(cls, reference, pv):
        temp = reference.split('/')[-1].split('.')[0] + pv.split('/')[-1].split('.')[0]
        calc = subprocess.Popen("{} rmsd.py -pv second -d {} {} {}".format(SCHRODINGER_RUN, temp, reference, pv), shell = True)
        calc.wait()
        data = map(float, open(temp).read().strip().split())
        subprocess.Popen("rm {}".format(temp), shell = True)
        return cls(data)

    def get_rmsd(self, pose):
        return self.data[pose]
    
    def getData(self):
        return self.data

    def getMinPose(self):
        return self.data.index(min(self.data))

    def best(self):
        return min(self.data)

    def __str__(self):
        return ','.join(map(str, self.data))
        
if __name__ == '__main__':
    import sys
    rmsd = RMSD(sys.argv[1], sys.argv[2])
    print rmsd
