from pathlib import Path
from tempfile import NamedTemporaryFile
import subprocess, math, random, csv, copy
import numpy as np

from circuitgraph.io import circuit_to_verilog, from_file, to_file

def ceildiv(a, b):
    #from https://stackoverflow.com/questions/14822184/is-there-a-ceiling-equivalent-of-operator-in-python
    return -(a // -b)

def removeBuffers(c):
    # removes redundant buffers from circuit (leaving outputs alone)
    for buf in c.filter_type("buf"):
        if buf not in c.outputs():
            fanin = c.fanin(buf)
            fanout = c.fanout(buf)
            #disconnect buf
            c.disconnect(fanin,buf)
            c.disconnect(buf,fanout)
            #connect fanin of buf to fanout of buf
            c.connect(fanin,fanout)
            #delete buf
            c.remove(buf)
    return c

def plainVisualize(c,output_file):
    """
    Visualize a circuit using Yosys. (adapted from circuitgraph's visualize)

    Parameters
    ----------
    c: Circuit
            Circuit to visualize.
    output_file: str
            Where to write the image to.
    """
    verilog = circuit_to_verilog(c)
    output_file = Path(output_file)
    fmt = output_file.suffix[1:]
    prefix = output_file.with_suffix("")
    with NamedTemporaryFile(
        prefix="circuitgraph_synthesis_input", suffix=".v"
    ) as tmp_in:
        tmp_in.write(bytes(verilog, "ascii"))
        tmp_in.flush()
        cmd = [
            "yosys",
            "-p",
            f"read_verilog {tmp_in.name}; rename -hide; show -format {fmt} -prefix {prefix}",
        ]
        subprocess.run(cmd)

    if fmt != "dot":
        prefix.with_suffix(".dot").unlink()

def optVisualize(c, output_file):
    """
    Optimize and then visualize a circuit using Yosys. (adapted from circuitgraph's visualize)

    Parameters
    ----------
    c: Circuit
            Circuit to visualize.
    output_file: str
            Where to write the image to.
    """
    verilog = circuit_to_verilog(c)
    output_file = Path(output_file)
    fmt = output_file.suffix[1:]
    prefix = output_file.with_suffix("")
    with NamedTemporaryFile(
        prefix="circuitgraph_synthesis_input", suffix=".v"
    ) as tmp_in:
        tmp_in.write(bytes(verilog, "ascii"))
        tmp_in.flush()
        cmd = [
            "yosys",
            "-p",
            f"read_verilog {tmp_in.name}; " f"rename -hide; opt; " f"show -format {fmt} -prefix {prefix}",
        ]
        subprocess.run(cmd)

    if fmt != "dot":
        prefix.with_suffix(".dot").unlink()

def debugVisualize(c, output_file):
    """
    Optimize and then visualize a circuit using Yosys with debug messages on. (adapted from circuitgraph's visualize)

    Parameters
    ----------
    c: Circuit
            Circuit to visualize.
    output_file: str
            Where to write the image to.
    """
    verilog = circuit_to_verilog(c)
    output_file = Path(output_file)
    fmt = output_file.suffix[1:]
    prefix = output_file.with_suffix("")
    with NamedTemporaryFile(
        prefix="circuitgraph_synthesis_input", suffix=".v"
    ) as tmp_in:
        tmp_in.write(bytes(verilog, "ascii"))
        tmp_in.flush()
        cmd = [
            "yosys",
            "-p",
            f"read_verilog {tmp_in.name}; " f"rename -hide; debug opt; " f"show -format {fmt} -prefix {prefix}",
        ]
        subprocess.run(cmd)

    if fmt != "dot":
        prefix.with_suffix(".dot").unlink()

def greatestPowerOfTwoLessThan(n):
    # helper function for bitonic merge
    #    input: int n
    #    output: int k
    # src: https://www.inf.hs-flensburg.de/lang/algorithmen/sortieren/bitonic/oddn.htm

    k = 1

    while(k>0 and k<n):
        k = k<<1
    return k>>1

def countAnds(c):
    #lazy alias for number of ands
    return len(c.filter_type('and'))

def countOrs(c):
    #lazy alias for number of ors
    return len(c.filter_type('or'))

def countInputs(c):
    #lazy alias for number of inputs
    return len(c.inputs())

def countOutputs(c):
    #lazy alias for number of outputs
    return len(c.outputs())

def maxShareSize(c):
    #maximum number of shares given to any party
    # (if c is a circuit, share size for Yao CSS; if c is a formula, share size for BL ITSS)
    max = 0
    for input in c.inputs():
        fanout = len(c.fanout(input))
        if fanout > max:
            max = fanout
    return max

def countWires(c):
    #total fan-in of ANDs, ORs, and outputs
    wires = 0
    for g in c.nodes():
        if ((c.type(g)=='or') or (c.type(g)=='and') or (c.is_output(g)==True)):
            wires += len(c.fanin(g))
    return wires

def getDepth(c):
    #returns maximum depth over all outputs
    #usually want to call removeBuffers on c before calling this function

    depth = 0
    for out in c.outputs():
        d = c.fanin_depth(out)
        depth = max(depth,d)
    return depth

def expandToFormula(c):
    #returns a formula (i.e., a circuit where all non-input nodes have fanout<=1)
    f = c.copy()
    while (not isFormula(f)):
        for node in f.topo_sort():
            #get an "unfinished" node (node with fanout >1)
            if ((f.type(node) != 'input') and (len(f.fanout(node))>1)):
                ins = f.fanin(node)
                outs = f.fanout(node)
                ancestors = f.transitive_fanin(node)
                ancestors.add(node)
                f.disconnect(node,outs)
                # make dictionary mapping original ancestors to list of copies
                copies = {}
                for ancestor in ancestors:
                    if (f.type(ancestor) != 'input'):
                        #make list of copies of the ancestor and store in dictionary
                        copylist = []
                        for i in range(len(outs)-1):
                            copy = f.add(f'copy_{i}',f.type(ancestor),uid=True)
                            copylist.append(copy)
                        copies[ancestor] = copylist
                # after making all the copies, connect them in the same way as originals
                for original in copies:
                    # for each node in the fanin of original node:
                    # if fanin is an input: for each i, connect i^th copy to input directly
                    # else (node is not an input):
                    #   for each i, connect i^th copy of original to i^th copy of fanin
                    for fanin in f.fanin(original):
                        if (f.type(fanin)=='input'):
                            f.connect(fanin,copies[original])
                        else:
                            for i in range(len(outs)-1):
                                f.connect(copies[fanin][i],copies[original][i])
                #finally connect i^th copy of node to outs[i]
                i=0
                for out in outs:
                    if i==0:
                        f.connect(node,out)
                    else:
                        f.connect(copies[node][i-1],out)
                    i+=1

    #sanity check
    if (isFormula(f)):
        return f
    else:
        return "error"

def isFormula(c):
    for n in c.nodes():
        fanout = c.fanout(n)
        if ((c.type(n) != 'input') and (len(fanout)>1)):
            return False
    return True

def yosysOpt(in_path,out_path):
    """
    Read, opt, write

    Parameters
    ----------
    in_path : path
            Verilog file containing circuit to optimize
    Output
    ----------
    out_path : path
            Path to destination for optimized circuit

    """
    cmd = [
        "yosys",
        "-p",
        f"read_verilog {in_path}; rename -hide; opt; proc; write_verilog -noattr {out_path} ",
    ]
    subprocess.run(cmd)
    return 0

def writeRowsToCSV(header,rows,filename):
    with open(filename, 'w', encoding='UTF8', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(header)
        writer.writerows(rows)

def writeRowToCSV(header,row,filename):
    with open(filename, 'w', encoding='UTF8', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(header)
        writer.writerow(row)

def appendRowToCSV(row,filename):
    with open(filename, 'a', encoding='UTF8', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(row)

def appendRowsToCSV(rows,filename):
    with open(filename, 'a', encoding='UTF8', newline='') as f:
        writer = csv.writer(f)
        writer.writerows(rows)

def generateTests(n,t,distr_type,distr_params,num_samples):
    #returns a dictionary with test names as keys and parameters [n,W,T,distr_type,params,num_samples] as values
    #   ns : int (number of parties to generate tests for)
    #   t: str (recognizes "half" and "twothirds")
    #   distr : str (recognizes )
    #       if type == geometric, expects params [lo, hi, p]
    #       if type == uniform, expects params [lo,hi]
    #       if type == zipf, expects [s]

    known_thrs = ["half", "twothirds"]
    known_types = ["uniform","geometric","zipf"]

    tests = {}

    if (not (distr_type in known_types)):
        print(f"distribution type {distr_type} not recognized")
    elif (not (t in known_thrs)):
        print(f"threshold type {t} not recognized")
    else:
        if (distr_type == "zipf"):
            s = distr_params[0]
            i = 0
            fail = 0
            while i < num_samples:
                L = np.random.default_rng().zipf(s, size=n)
                W = L.tolist()
                totalW = sum(W)

                if (t == "half"):
                    thr = ceildiv(totalW - 2, 2)
                elif (t == "twothirds"):
                    thr = ceildiv(totalW - 3, 3)

                if max(W) >= thr:
                    fail+=1
                else:
                    tests[f'sample_{i}'] = [n,W,thr,distr_type,t,distr_params]
                    i+=1
            print(f"threw out {fail} tests, kept {i} tests")
        else:
            if (distr_type == "uniform"):
                lo = distr_params[0]
                hi = distr_params[1]
                values = range(lo,hi)
                #using random.choices with weights = None chooses uniformly from values
                weights = None
            elif (distr_type =="geometric"):
                lo = distr_params[0]
                hi = distr_params[1]
                p = distr_params[2]
                values = range(lo,hi)
                weights = []
                for j in range(lo,hi):
                    #create list of weights for random.choices (same length as values)
                    weights.append(math.pow(1-p,j-lo)*p)
            for i in range(num_samples):
                W = random.choices(values,weights,k=n)
                totalW = sum(W)
                if (t == "half"):
                    thr = ceildiv(totalW - 2, 2)
                elif (t == "twothirds"):
                    thr = ceildiv(totalW - 3, 3)
                tests[f'sample_{i}'] = [n,W,thr,distr_type,t,distr_params]

    return tests
