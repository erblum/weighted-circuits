import csv
import circuitgraph as cg
import mysubcircuits as sc
import myutils as mu
import subprocess, math
from pathlib import Path
import os

def getStats(c):
    #given circuit c, return [header,row]
    header = ['num_gates','num_wires','depth','max_input_fanout', 'total_input_fanout']
    c = mu.removeBuffers(c)
    numAnds = mu.countAnds(c)
    numOrs = mu.countOrs(c)
    size = numAnds + numOrs
    numWires = mu.countWires(c)
    fanouts = []
    for node in c.inputs():
        fanouts.append(len(c.fanout(node)))
    maxfanout = max(fanouts)
    totalfanout = sum(fanouts)
    depth = mu.getDepth(c)
    row = [size,numWires,depth,maxfanout,totalfanout]
    if not len(row)==len(header):
        print("check that number of values = number of columns")
    return [header,row]

def buildToFile(n,W,T,dir,filename):
    """
    Builds and optimizes circuit with parameters n, W, T
    Circuit is written to {dir}/{filename}_opt2.v (overwriting if file already exists)
    Parameters n, W, T are saved to {filename}_args.txt
    """
    if not os.path.exists(dir):
        os.mkdir(dir)
    with open(dir+f'{filename}_args.txt', 'w') as f:
        f.write(f'n: {n}\nW: {W}\nT: {T}\ndistribution: {filename}')
    #build circuit
    opt1 = sc.weightedThresholdCircuit(n,W,T,original=False)
    #write to file
    cg.to_file(opt1,dir+f"{filename}_opt1.v")
    #yosysOpt
    mu.yosysOpt(dir+f"{filename}_opt1.v",dir+f"{filename}_opt2.v")
    #delete intermediate file
    os.remove(dir+f"{filename}_opt1.v")

def batchEval(batchdir,csvname):
    """
    for each Verilog file in batchdir:
        - load circuit from file
        - compute stats
        - write stats to csv
    """
    header = ['name','num_gates','num_wires','depth','max_input_fanout', 'total_input_fanout']
    csvpath = batchdir + f'{csvname}.csv'
    if (not (os.path.exists(csvpath))):
        with open(csvpath, 'w', encoding='UTF8', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(header)
    for filename in os.scandir(batchdir):
        filename = Path(filename)
        fmt = filename.suffix
        if (fmt == ".v"):
            name = filename.stem
            c = cg.from_file(filename)
            data = getStats(c)[1]
            print(f"data: {data}")
            row = [name] + data
            print(f"row: {row}")
            if (len(row)!=len(header)):
                print("check that number of values = number of columns")
            else:
                with open(csvpath, 'a', encoding='UTF8', newline='') as f:
                    writer = csv.writer(f)
                    writer.writerow(row)

def batchCompareAll(ns,basedir,csvname):
    """
    for each n in ns, build bw, opt1, and opt2 versions of u_n (or load
        from file if they already exist in basedir) and write stats to basedir/csvname.csv
    """
    #set up csv file
    header = ['circuit type', 'n','num_gates',
                        'num_wires','depth','max_input_fanout', 'total_input_fanout']

    if (not (os.path.exists(basedir+f'{csvname}.csv'))):
        with open(basedir+f'{csvname}.csv', 'w', encoding='UTF8', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(header)

    for n in ns:

        tau = math.ceil(n*math.log2(n))
        if (not os.path.exists(basedir+f"bw_{n}.v")):
            bw = sc.universalThresholdCircuit(n,tau,original=True)
            cg.to_file(bw,basedir+f"u_{n}_bw.v")
        else:
            bw = cg.from_file(basedir+f"u_{n}_bw.v")

        if (not os.path.exists(basedir+f"u_{n}_opt1.v")):
            opt1 = sc.universalThresholdCircuit(n,tau,original=False)
            cg.to_file(opt1,basedir+f"u_{n}_opt1.v")
        else:
            opt1 = cg.from_file(basedir+f"u_{n}_opt1.v")

        if (not os.path.exists(basedir+f"u_{n}_opt2.v")):
            mu.yosysOpt(basedir+f"u_{n}_opt1.v",basedir+f"u_{n}_opt2.v")
            opt2 = cg.from_file(basedir+f"u_{n}_opt2.v")
        else:
            opt2 = cg.from_file(basedir+f"u_{n}_opt2.v")

        cs = [[bw,"bw"],[opt1,"opt1"],[opt2,"opt2"]]

        for c in cs:
            res = getStats(c[0])
            data = res[1]
            # row: type,n, #gates,#wires,max_input_fanout,total_input_fanout,depth
            row = [c[1]]+[n]+data
            if (len(row)!=len(header)):
                print("check that number of values = number of columns")
            else:
                with open(basedir+f'{csvname}.csv', 'a', encoding='UTF8', newline='') as f:
                    writer = csv.writer(f)
                    writer.writerow(row)

def batchBuildU_n(ns,basedir):
    #build opt2 version of U_n for each n in ns, located in basedir
    for n in ns:
        if (not os.path.exists(basedir+f"u_{n}_opt2.v")):
            tau = math.ceil(n*math.log2(n))
            u = sc.universalThresholdCircuit(n,tau,original=False)
            cg.to_file(u,basedir+f"u_{n}_opt1.v")
            res = mu.yosysOpt(basedir+f"u_{n}_opt1.v",basedir+f"u_{n}_opt2.v")
            os.remove(basedir+f"u_{n}_opt1.v")

def batchBuildFn(n,basedir,tests,batchname,csvname):
    """
    optimizes u_{n+1} then builds f_n
       test circuits will be put in directory batchdir = {basedir}+{batchname}
       if opt_u_{n+1}.v exists in {basedir}, will load from file instead of building from scratch
       writes stats to {csvname}.csv in batchdir
    """
    #set up directory and csv file
    batchdir = basedir+f"{batchname}/"
    header = ['test_name', 'n','W','T','distr_type', 'thr_type', 'distr_params',
            'num_gates','num_wires','depth','max_input_fanout', 'total_input_fanout']
    if (not os.path.exists(batchdir)):
        os.mkdir(batchdir)

    if (not (os.path.exists(batchdir+f'{csvname}.csv'))):
        with open(batchdir+f'{csvname}.csv', 'w', encoding='UTF8', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(header)

    #build base circuit u_n (if it's not already there)
    if (not os.path.exists(basedir+f"u_{n+1}_opt1.v")):
        tau = math.ceil(n*math.log2(n))
        u = sc.universalThresholdCircuit(n+1,tau)
        cg.to_file(u,basedir+f"u_{n+1}_opt1.v")
        mu.yosysOpt(basedir+f"u_{n+1}_opt1.v",basedir+f"u_{n+1}_opt2.v")

    u_opt = cg.from_file(basedir+f"u_{n+1}_opt2.v")

    for t in tests:
        if ( len(tests[t]) != 6):
            print("check that tests have 6 values: n,W,thr,distr_type,thr_type,distr_params")
        else:
            #parse test parameters
            n = tests[t][0]
            W = tests[t][1]
            T = tests[t][2]
            distr_type = tests[t][3]
            thr_type = tests[t][4]
            distr_params = tests[t][5]

            #build circuit using prebuilt base circuit
            opt1 = sc.weightedThresholdCircuit(n,W,T,original=False,prebuilt=True,prebuilt_circuit=u_opt)
            cg.to_file(opt1,batchdir+f"{t}_opt1.v")
            mu.yosysOpt(batchdir+f"{t}_opt1.v",batchdir+f"{t}_opt2.v")
            opt2 = cg.from_file(batchdir+f"{t}_opt2.v")
            #delete intermediate file
            os.remove(batchdir+f"{t}_opt1.v")

            cs = [[opt1,f"{t}_opt1"],[opt2,f"{t}_opt2"]]
            for c in cs:
                data = getStats(c[0])
                # row: type,n, #gates,#wires,max_input_fanout,total_input_fanout,depth
                row = [c[1]]+[n]+[W]+[T]+[distr_type] + [thr_type] + [distr_params] + data[1]
                if (len(row)!=len(header)):
                    print("error: length of row != length of header")
                    print(f"row: {row}")
                    print(f"header: {header}")
                else:
                    with open(batchdir+f'{csvname}.csv', 'a', encoding='UTF8', newline='') as f:
                        writer = csv.writer(f)
                        writer.writerow(row)

def batch(n,t,distr,params,num_samples,basedir,batchname,csvname):
    """
        Generates tests according to n, t, distr, params, k
            and then calls batchBuildFn on n,basedir,tests,batchname

        n: int (number of parties)
        t: "twothirds" or "half"
        distr: one of the types recognized by generateTests
        params: params recognized by generateTests
        k: number of sample vectors
        basedir: location of opt_n, if it exists
        batchname: where to put (or look for) circuits
        csvname: where to put stats
    """
    tests = mu.generateTests(n,t,distr,params,num_samples)
    batchBuildFn(n,basedir,tests,batchname,csvname)
