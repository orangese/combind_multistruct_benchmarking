import os
import sys

def redo(struct):
    subdirs = ['structures','ifp','docking']
    for subdir in subdirs:
        for ssdir in os.listdir(subdir):
            if ssdir[0] == '.' or ssdir in ['downloads','raw_files','processed_files']: continue
            for item in os.listdir('{}/{}'.format(subdir, ssdir)):
                if struct in item:
                    print '{}/{}/{}'.format(subdir, ssdir, item)
                    os.system('rm -rf {}/{}/{}'.format(subdir, ssdir, item))


