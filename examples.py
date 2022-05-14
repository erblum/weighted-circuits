import circuitgraph as cg
import mysubcircuits as sc
import myutils as mu
import eval as ev
import subprocess, math, os
from pathlib import Path
import csv

def exampleCircuit1():
    c = cg.Circuit()
    #add gates
    c.add('in_0','input')
    c.add('in_1','input')
    c.add('and_0','and', fanin=['in_0','in_1'])
    c.add('and_1','and', fanin=['in_0','in_1'])
    c.add('out','or',fanin=['and_0','and_1'], output=True)
    return c

def exampleCircuit2():
    c = cg.Circuit()
    c.add('in_1','input')
    c.add('in_2','input')
    c.add('in_3','input')
    c.add('and_1','and', fanin=['in_1','in_2'])
    c.add('and_2','and', fanin=['in_2','in_3'])
    c.add('or_1','or', fanin=['and_1','and_2'])
    c.add('and_3','and', fanin=['and_1','and_2'])
    c.add('out','or', fanin=['and_3','or_1'],output=True)
    return c

def comparatorExample():
    c = sc.comparator()
    return c

def sorterExample(n):
    c = mu.removeBuffers(sc.sorter(n))
    return c

def mergeHalveExample(n):
    c = mu.removeBuffers(sc.mergeHalve(n))
    return c

def BWMergeHalveExample(n):
    c = mu.removeBuffers(sc.BWMergeHalve(n))
    return c

def u_nExample(n,original=False):
    #if original==True, builds bw version; otherwise, builds opt1 version
    tau = math.ceil(n*math.log2(n))
    c = sc.universalThresholdCircuit(n,tau,original=False)
    return c


bw = u_nExample(3,original=True)
mu.plainVisualize(bw,"./u_3_bw.jpg")

# build opt1, opt2 versions of a 4-party weighted threshold circuit with weights [1,1,2,3]
opt1 = sc.weightedThresholdCircuit(4,[1,1,2,3],3,original=False,prebuilt=False,prebuilt_circuit=None)
cg.to_file(opt1,"./weighted_example_opt1.v")
mu.yosysOpt("./weighted_example_opt1.v","./weighted_example_opt2.v")
opt2 = cg.from_file("./weighted_example_opt2.v")

# visualize circuits and export as .jpg
mu.plainVisualize(opt1,"./weighted_example_opt1.jpg")
mu.optVisualize(opt2,"./weighted_example_opt2.jpg")

# build u_{n}_opt2.v for n=4,5,6 in basedir:
ev.batchBuildU_n([4,5,6],"./basedir/")


#Generates tests according to n=8, t="half", distr="uniform", params=[1,64], num_samples = 20,
#   then calls batchBuildFn on n=8,basedir="./basedir",tests,"my_batch","my_csv"
ev.batch(8,"half","uniform",[1,64],20,"./basedir/","my_batch","my_csv")
