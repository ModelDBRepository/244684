# This code was used in: Masquelier & Kheradpisheh (2018) Optimal localist and distributed coding of spatiotemporal spike patterns through STDP and coincidence detection. Frontiers in Computational Neuroscience.
# Aug 2018
# timothee.masquelier@cnrs.fr
#
# This is a Python script to launch several threads of main.m (with different random seeds). This is useful if multiple cores are available.
# The number of threads can be specfied (currently 16, see below)
# For example batch.py -i 0 -f 8 -s 1 will launch 8 threads, each thread will save its results in ../data/
#
# A priori, you want to remove all mat files in ../data/ before launching this batch.

import subprocess
from numpy import *
import sys
import getopt

try:
    opts, args = getopt.getopt(sys.argv[1:], 'i:f:s:', ['initial=', 'final=','step='])
except getopt.GetoptError, err:
    # print help information and exit:
    print str(err) # will print something like "option -a not recognized"
#    usage()
    sys.exit(2)
#opts, extraparams = getopt.getopt(sys.argv[1:],'r:w:',['randomSeed='])
# starts at the second element of argv since the first one is the script name
# extraparms are extra arguments passed after all option/keywords are assigned
# opts is a list containing the pair "option"/"value"
for o,p in opts:
  if o in ['-i','--initial']:
     initial = double(p)
  elif o in ['-f','--final']:
     final = double(p)
  elif o in ['-s','--step']:
     step = double(p)

        
import multiprocessing
import subprocess

def calculate(alpha):
    
    subprocess.call('matlab -singleCompThread -nojvm -nosplash -nodesktop -r  "seed='+ str(alpha) +';main "; ',shell=True)
    

if __name__ == '__main__':
    # pool = multiprocessing.Pool(None) # Use this if you want to use all the CPUs
    pool = multiprocessing.Pool(16) # Use this to specifify the number of threads
    tasks = arange(initial,final,step)
    results = []
    r = pool.map_async(calculate, tasks, callback=results.append)
    r.wait() # Wait on the results
    print results
