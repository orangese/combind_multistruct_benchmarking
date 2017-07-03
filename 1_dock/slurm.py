import os

#Process to submit a job that will wait until it terminates (salloc)
#Important Notes:
#  -By default, submits to the rondror partition
#  -Inherits the environment of the parent process (includes path of the parent's working directory)
#  -Asks for all of the processors on one node
def salloc(jobString, numProcessors, timeLimit):
    print(os.getcwd())
    os.system("salloc -n {} -p rondror -t {} /usr/bin/srun --ntasks=1 --nodes=1 --preserve-env {}".format(numProcessors, timeLimit, jobString))
