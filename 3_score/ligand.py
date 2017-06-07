#Container for information about a ligand
class Ligand:
    def __init__(self, crystal):
        self.crystal = crystal
        self.poses = {}

    def add_pose(self, pose, pose_num):
        self.poses[pose_num] = pose

    def __str__(self):
        return str(self.crystal)

