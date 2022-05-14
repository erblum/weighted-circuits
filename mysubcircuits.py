from circuitgraph import Circuit, BlackBox
import circuitgraph.tx as ct
import circuitgraph.io as io
from myutils import greatestPowerOfTwoLessThan
import math

#extends circuitgraph library, found at circuitgraph.github.io/circuitgraph/index.html
def mergeHalve(n):
    """
    Create a (lazy) merge-halve circuit.

    Parameters
    ----------
    n : int
            Number of outputs
    Returns
    -------
    Circuit
            Merge-halve circuit
             - 2n inputs, n outputs
             - If inputs are pre-sorted (high to low), output is sorted
             - (sum of output bits) = floor(sum of inputs bits/2)
    """
    c = Circuit(name='merge_halve')

    # create inputs and outputs
    for i in range(1, n+1):
        c.add(f'in_L_{i}', 'input')
        c.add(f'in_R_{i}', 'input')

    for l in range(2, 2*n+1, 2):
        if (l <= n):
            #case 1: out_{l \\ 2) = OR of {x_l, y_l, {(x_i AND y_j) for each i,j s.t. i+j=l}}
            c.add(f'out_{l//2}',
                'or',
                fanin=[c.add(f'AND_{i}_{l-i}', 'and', fanin=[f'in_L_{i}',f'in_R_{l-i}']) for i in range(1,l)]
                    + [f'in_L_{l}']
                    + [f'in_R_{l}'],
                output=True)
        elif (l < 2*n):
            #case 2: out_{l \\ 2) = OR of {(x_i AND y_j) for each i,j s.t. i+j=l}
            c.add(f'out_{l//2}',
                'or',
                fanin=[c.add(f'AND_{i}_{l-i}', 'and', fanin=[f'in_L_{i}',f'in_R_{l-i}']) for i in range(l-n,n+1)]
                    + [f'in_L_{2*n-l}']
                    + [f'in_R_{2*n-l}'],
                output=True)
        else:
            #case 3: out_n = x_n AND y_n
            c.add(f'out_{l//2}',
            'and',
            fanin = [f'in_L_{n}',f'in_R_{n}'],
            output=True)

    return c

def BWMergeHalve(n):
    """
    Create a merge-halve circuit as specified in BW.

    Parameters
    ----------
    n : int
            Number of outputs
    Returns
    -------
    Circuit
            Merge-halve circuit
             - 2n inputs, n outputs
             - If inputs are pre-sorted (high to low), output is sorted
             - (sum of output bits) = floor(sum of inputs bits/2)
    """
    c = Circuit(name='bw_merge_halve')

    # create inputs and outputs
    for i in range(1, n+1):
        c.add(f'in_L_{i}', 'input')
        c.add(f'in_R_{i}', 'input')
        c.add(f'out_{i}', 'buf',output=True)

    for l in range(1, 2*n+1):
        if (l == 1):
            c.add(f'OR_{l}',
            'or',
            fanin=[f'in_L_{l}',f'in_R_{l}'])
        elif (l <= n):
            #case 1: out_{l \\ 2) = OR of {x_l, y_l, {(x_i AND y_j) for each i,j s.t. i+j=l}}
            if (l % 2 == 0):
                c.add(f'OR_{l}',
                    'or',
                    fanin=[c.add(f'AND_{i}_{l-i}', 'and', fanin=[f'in_L_{i}',f'in_R_{l-i}']) for i in range(1,l)]
                        + [f'in_L_{l}']
                        + [f'in_R_{l}'],
                    fanout=[f'out_{l//2}'])
            else:
                    c.add(f'OR_{l}',
                        'or',
                        fanin=[c.add(f'AND_{i}_{l-i}', 'and', fanin=[f'in_L_{i}',f'in_R_{l-i}']) for i in range(1,l)]
                            + [f'in_L_{l}']
                            + [f'in_R_{l}'])

        elif (l < 2*n):
            #case 2: out_{l \\ 2) = OR of {(x_i AND y_j) for each i,j s.t. i+j=l}
            if (l % 2 == 0):
                c.add(f'OR_{l}',
                    'or',
                    fanin=[c.add(f'AND_{i}_{l-i}', 'and', fanin=[f'in_L_{i}',f'in_R_{l-i}']) for i in range(l-n,n+1)]
                        + [f'in_L_{2*n-l}']
                        + [f'in_R_{2*n-l}'],
                    fanout=[f'out_{l//2}'])
            else:
                c.add(f'OR_{l}',
                    'or',
                    fanin=[c.add(f'AND_{i}_{l-i}', 'and', fanin=[f'in_L_{i}',f'in_R_{l-i}']) for i in range(l-n,n+1)]
                        + [f'in_L_{2*n-l}']
                        + [f'in_R_{2*n-l}'])

        else:
            #case 3: out_n = x_n AND y_n
            c.add(f'AND_{n}_{n}',
            'and',
            fanin = [f'in_L_{n}',f'in_R_{n}'],
            fanout=[f'out_{l//2}'])

    return c

def mergeHalveBB(n):
    """
    Create a merge-halve blackbox.

    Parameters
    ----------
    n : int
            Number of outputs
    Returns
    -------
    BlackBox
            Merge-halve circuit
             - 2n inputs, n outputs
    """
    inputs = [f"in_{i}" for i in range(1,2*n+1)]
    outputs = [f"out_{i}" for i in range(1,n+1)]
    return BlackBox("merge_and_halve",inputs,outputs)

    return c

def comparatorBB(dir=True):
    """
    Comparator blackbox

    Parameters
    ----------
    dir : Boolean
            sorting direction; default is descending (high to low)

    Returns
    -------
    c : Blackbox
            Comparator
            - 2 inputs in_0, in_1
            - 2 outputs out_0, out_1
    """

    inputs = ["in_0","in_1"]
    outputs = ["out_0","out_1"]
    return BlackBox("comparator",inputs,outputs)

def comparator(dir=True):
    """
    Create a comparison subcircuit (implements comparison gate with AND/OR)
    with direction dir

    Parameters
    ----------
    dir : Boolean
            sorting direction; default is descending (high to low)

    Returns
    -------
    c : Circuit
            Comparator
            - 2 inputs in_0, in_1
            - 2 outputs out_0, out_1
            - if dir=True (default) outputs (in_0 OR in_1, in_0 AND in_1)
            - otherwise outputs (in_0 AND in_1, in_0 OR in_1)
    """

    c = Circuit()
    c.add('in_0','input')
    c.add('in_1','input')
    c.add('g_OR','or',fanin=["in_0","in_1"],output=True)
    c.add('g_AND','and',fanin=["in_0","in_1"],output=True)
    if (dir==True):
        c.relabel({'g_OR':'out_0'})
        c.relabel({'g_AND':'out_1'})
    else:
        c.relabel({'g_AND':'out_0'})
        c.relabel({'g_OR':'out_1'})

    return c


def sorterBB(n,dir=True):
    """
    Sorter blackbox

    Parameters
    ----------
    n : int
            number of inputs
    dir : Boolean
            sorting direction; default is descending (high to low)

    Returns
    -------
    c : BlackBox
            sorter bb
            - inputs: in_1,...,in_{n}
            - outputs: out_1,...,out_{n}
    """

    inputs = [f"in_{i}" for i in range(1,n+1)]
    outputs = [f"out_{i}" for i in range(1,n+1)]
    return BlackBox("sort",inputs,outputs)

def sorter(n,dir=True):
    """
    Wrapper for bitonicSorter

    Parameters
    ----------
    n : int
            number of inputs
    dir : Boolean
            sorting direction; default is descending (high to low)

    Returns
    -------
    c : Circuit
            sorter
            - inputs: in_1,...,in_{n}
            - outputs: out_1,...,out_{n}
            - sum(inputs)=sum(outputs)
            - outputs are sorted
    """
    c = Circuit()
    s = bitonicSorter(0,n)
    c.add_subcircuit(s,"bsort")
    for i in range(n):
        c.relabel({f"bsort_in_{i}" : f"in_{i+1}"})
        c.relabel({f"bsort_out_{i}" : f"out_{i+1}"})
        c.set_type(f"in_{i+1}","input")
        c.set_output(f"out_{i+1}",True)
    return c

def bitonicSorter(lo,n,dir=True):
    """
    Create sorting subcircuit for n inputs

    based on https://www.inf.hs-flensburg.de/lang/algorithmen/sortieren/bitonic/oddn.htm
    further discussion: https://jix.one/proving-50-year-old-sorting-networks-optimal-part-1/

    Parameters
    ----------
    lo : int
            starting index
    n : int
            number of inputs
    dir : Boolean
            sorting direction; default is descending (high to low)

    Returns
    -------
    c : Circuit
            bitonic sorter
            - inputs: in_{lo},...,in_{lo+n}
            - outputs: out_{lo},...,out_{lo+n}
            - sum(inputs)=sum(outputs)
            - outputs are sorted
    """
    c = Circuit()
    if (n>1):
        mid = n // 2

        #build subcircuits
        c.add_subcircuit(bitonicSorter(lo, mid, not dir), "L0")
        c.add_subcircuit(bitonicSorter(lo+mid, n-mid, dir), "L1")
        c.add_subcircuit(bitonicMerger(lo, n, dir),"R")

        #relabel inputs/outputs

        for i in range(lo,lo+n):
            c.relabel({f"R_out_{i}" : f"out_{i}"})
            c.set_output(f"out_{i}",True)
            if (i<lo+mid):
                c.relabel({f"L0_in_{i}" : f"in_{i}"})
            else:
                c.relabel({f"L1_in_{i}" : f"in_{i}"})

        # connect outputs of bsorts to inputs of merge
        for i in range(lo,lo+mid):
            c.connect(f"L0_out_{i}", f"R_in_{i}")
        for i in range(lo+mid,lo+n):
            c.connect(f"L1_out_{i}", f"R_in_{i}")

    else:
        # otherwise return trivial circuit
        c.add(f"in_{lo}","input")
        c.add(f"out_{lo}","buf",fanin=[f"in_{lo}"],output=True)

    return c

def mergeHelper(lo,n,dir=True):
    """
    Creates the comparators for bitonicMerger

    Parameters
    ----------
    lo : int
            starting index
    n : int
            number of inputs
    dir : Boolean
            sorting direction; default is descending (high to low)

    Returns
    -------
    c : Circuit
    """
    c = Circuit()
    if (n>1):
        pow = greatestPowerOfTwoLessThan(n)
        #create comparators and relabel inputs/outputs
        for i in range(lo,lo+n-pow):
            c.add_subcircuit(comparator(dir),f"C{i}")
            c.relabel({f"C{i}_in_0" : f"in_{i}"})
            c.relabel({f"C{i}_in_1" : f"in_{i+pow}"})
            c.relabel({f"C{i}_out_0" : f"out_{i}"})
            c.relabel({f"C{i}_out_1" : f"out_{i+pow}"})
            c.set_type([f"in_{i}",f"in_{i+pow}"],"input")
            c.set_output([f"out_{i}",f"out_{i+pow}"],True)
        #create other wires
        for i in range(lo+n-pow,lo+pow):
            c.add(f"in_{i}","input")
            c.add(f"out_{i}","buf",fanin=[f"in_{i}"],output=True)
    else:
        c.add(f"in_{lo}","input")
        c.add(f"out_{lo}","buf",fanin=[f"in_{lo}"],output=True)
    return c

def bitonicMerger(lo,n,dir=True):
    """
    Merges two bitonic sequences into a single bitonic sequence

    Parameters
    ----------
    lo : int
            starting index
    n : int
            number of inputs
    dir : Boolean
            sorting direction; default is descending (high to low)

    Returns
    -------
    c : Circuit
            merger
            - inputs: in_lo,...,in_{lo+n}
            - outputs: out_lo,...,out_{lo+n}
            - outputs are a bitonic sequence
    """
    c = Circuit()
    if (n>1):
        pow = greatestPowerOfTwoLessThan(n)

        #build subcircuits
        c.add_subcircuit(mergeHelper(lo,n,dir),"L")
        c.add_subcircuit(bitonicMerger(lo,pow,dir),"R0")
        c.add_subcircuit(bitonicMerger(lo+pow,n-pow,dir),"R1")

        #relabel inputs/outputs
        for i in range(lo,lo+n):
            c.relabel({f"L_in_{i}" : f"in_{i}"})
            if (i<lo+pow):
                c.relabel({f"R0_out_{i}" : f"out_{i}"})
            else:
                c.relabel({f"R1_out_{i}" : f"out_{i}"})
            c.set_type(f"in_{i}","input")
            c.set_output(f"out_{i}",True)

        #connect outputs of helper to inputs of merges
        for i in range(lo,lo+n):
            if (i<lo+pow):
                c.connect(f"L_out_{i}",f"R0_in_{i}")
            else:
                c.connect(f"L_out_{i}",f"R1_in_{i}")

    else:
        # otherwise return trivial circuit
        c.add(f"in_{lo}","input")
        c.add(f"out_{lo}","buf",fanin=[f"in_{lo}"],output=True)
    return c

#def recursiveBlock(lo,hi,n): **
def recursiveBlock(lo,hi,n,original=False):
    """
    Implements Beimel and Weinreb's shallow circuit recursive building block.

    Parameters
    ----------
    lo : int
            starting index ('i' in original paper)
    hi : int
            ending index ('j' in original paper)
    n : int
            number of inputs of the top-level circuit
    dir : Boolean
            sorting direction; default is descending (high to low)

    Returns
    -------
    c : Circuit
            recursive block D_{lo,hi}
            - inputs: presorted input variables for levels lo to hi (inclusive)
                indexed as in_i_j = jth variable of level i
                i=lo,...,hi; j=1,...n
            - outputs: n+1 lists of n output variables,
                indexed as out_i_j = jth output of ith list
                i=0,...,n; j=1,...n

    """
    c = Circuit()

    #recursive case (lo<hi)
    if lo<hi:
        mid = (hi+lo) // 2
        #build D_lo,mid
        #DL = recursiveBlock(lo,mid,n) **
        DL = recursiveBlock(lo,mid,n,original)
        c.add_subcircuit(DL,'DL',strip_io=False)
        #relabel inputs DL_in_{i}_{j} as in_{i}_{j}
        new_DL_labels = {}
        for i in range(lo,mid+1):
            for j in range(1,n+1):
                new_DL_labels[f'DL_in_{i}_{j}']=f'in_{i}_{j}'
        c.relabel(new_DL_labels)

        #build D_mid+1,hi
        #DR = recursiveBlock(mid+1,hi,n) **
        DR = recursiveBlock(mid+1,hi,n,original)
        c.add_subcircuit(DR,'DR',strip_io=False)
        #relabel inputs DL_in_{i}_{j} as in_{i}_{j}
        new_DR_labels = {}
        for i in range(mid+1,hi+1):
            for j in range(1,n+1):
                new_DR_labels[f'DR_in_{i}_{j}']=f'in_{i}_{j}'
        c.relabel(new_DR_labels)

        #create n selectors indexed k=0,...,n
        #   the input to selector k is n+2 lists of n variables:
        #       DL_out_k_1,...,DL_out_k_n, DR_out_0_1,...,DR_out_0_n, ... DR_out_n_1,..., DR_out_n_n
        for k in range(0,n+1):
            for i in range(1,n+1):
                for j in range(1,n+1):
                    #connect and_k_i_j to DL_out_k_i and DR_out_i_j
                    c.add(f'and_{k}_{i}_{j}', 'and',
                        fanin = [f'DL_out_{k}_{i}',f'DR_out_{i}_{j}'])
            for i in range(1,n+1):
                ins = [f'DR_out_0_{i}']
                #connect out_k_i (the 'or' gate) to DR_out_0_i, and_k_1_i, ..., and_k_n_i
                for j in range(1,n+1):
                    ins.append(f'and_{k}_{j}_{i}')
                c.add(f'out_{k}_{i}','or', fanin=ins)

    #base case (lo=hi)
    # create constants 1,0
    # for i = 0,...,n:
    #  make a mergeHalve subcircuit with 2n inputs, n outputs
    #  connect level lo inputs to first n inputs
    #  connect 1 to next i inputs, then connect 0 to remaining n-i inputs

    else: #lo == hi
        c.add('const0','0')
        c.add('const1','1')
        #create inputs 'in_{lo}_1',...,'in_{lo}_n'
        for i in range(1,n+1):
            c.add(f'in_{lo}_{i}','input')

        for i in range(0,n+1):
            #create and hook up M0,...,Mn
            cnx = {}
            for j in range(1,n+1):
                #input in_L_j of merge connects to in_lo_j
                cnx[f'in_L_{j}']=f'in_{lo}_{j}'
                if j<=i:
                    #inputs in_R_j up to i of merge connects to const1
                    cnx[f'in_R_{j}']='const1'
                else:
                    #remaining inputs in_R_j of merge connect to const0
                    cnx[f'in_R_{j}']='const0'
                    #output out_j of merge becomes out_i_j

            if (original==True):
                c.add_subcircuit(BWMergeHalve(n),f'M{i}',connections=cnx)
            else:
                c.add_subcircuit(mergeHalve(n),f'M{i}',connections=cnx)

            new_labels = {}
            for j in range(1,n+1):
                new_labels[f'M{i}_out_{j}']=f'out_{i}_{j}'
            c.relabel(new_labels)

    return c

def universalThresholdCircuit(n,tau,original=False):
    """
    Implements Beimel and Weinreb's universal threshold circuit u_n.

    Parameters
    ----------
    n

    Returns
    -------
    c : Circuit
        universal threshold circuit u_n
        - inputs: in_0_1,...,in_0_n,...,in_{tau-1}_1,in_{tau-1}_n
        - output: out

    """
    c = Circuit()

    #D = recursiveBlock(0,tau-1,n)
    D = recursiveBlock(0,tau-1,n,original)

    c.add_subcircuit(D,"D")

    #sort inputs and connect to inputs of D
    #D inputs: D_in_0_1,...,D_in_0_n,...,D_in_{tau-1}_1,...,D_in_{tau-1}_n

    sort = sorter(n)
    for i in range(0,tau):
        c.add_subcircuit(sort,f"sort_{i}")
        for j in range(1,n+1):
            c.add(f"in_{i}_{j}",'input')
            c.connect(f"in_{i}_{j}",f"sort_{i}_in_{j}")
            c.connect(f"sort_{i}_out_{j}",f"D_in_{i}_{j}")

    # connect first n outputs of D to single OR
    c.add('out','or',fanin=[f"D_out_0_{i}" for i in range(1,n+1)],output=True)

    return c

#def weightedThresholdCircuit(n,W_f,T_f):
def weightedThresholdCircuit(n,W_f,T_f,original=False,prebuilt=False,prebuilt_circuit=None):
    """
    Implements Beimel and Weinreb's universal threshold circuit f_n for
    weights W=[w_1,...w_n] and threshold T_f.

    Parameters
    ----------
    n : int
        number of parties
    W_f : list of ints
        W[i] = weight of party i
    T_f : int
        threshold

    Returns
    -------
    c : Circuit
        threshold circuit f_n
        - inputs: in_1,...,in_n
        - outputs: out

    """
    c = Circuit()
    tau = math.ceil(n*math.log2(n))

    if (not (len(W_f)==n)):
        print("Error: number of weights not equal to n")
        return c
    else:
        if (prebuilt==False):
        #u_n_plus_1 = universalThresholdCircuit(n+1,tau)
            u_n_plus_1 = universalThresholdCircuit(n+1,tau,original)
        else:
            if (prebuilt_circuit == None):
                printf(f"prebuilt == True but no prebuilt circuit passed as argument")
                return c
            else:
                u_n_plus_1 = prebuilt_circuit

        c.add('const0', '0')
        c.add('const1', '1')
        c.add('out','buf',output=True)
        for i in range(1,n+1):
            c.add(f'in_{i}','input')


        c.add_subcircuit(u_n_plus_1,'U')
        #get bit representation of each weight and use to connect variables:
        #       if binstr(W[j])[::-1][i] = 1, connect x_j to x_i_j.
        #       otherwise, connect const0 to x_i_j
        c.connect('U_out','out')
        for j in range(0,n):
            wt = W_f[j]
            bin = f'{wt:b}'.zfill(tau)
            for k in range(0,tau):
                if (bin[::-1][k] == '1'):
                    c.connect(f'in_{j+1}',f'U_in_{k}_{j+1}')
                else:
                    c.connect('const0',f'U_in_{k}_{j+1}')

        #now connect padding variables to in_i_{n+1} to make thresholds match
        dif = pow(2,tau) - T_f
        bin = f'{dif:b}'.zfill(tau)
        for k in range(0,tau):
            if (bin[::-1][k] == '1'):
                c.connect(f'const1',f'U_in_{k}_{n+1}')
            else:
                c.connect('const0',f'U_in_{k}_{n+1}')
    return c
