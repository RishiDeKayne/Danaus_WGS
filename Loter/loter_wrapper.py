import gzip, sys, argparse
import numpy as np
import loter.locanc.local_ancestry as lc

def extend_good_calls(anc, score, minimum):
    n = anc.shape[0]
    good = np.zeros(n, dtype=int)
    current = -1
    #forward search
    for i in range(n):
        if score[i] >= minimum:
            good[i] = 1
            current = anc[i]
        elif anc[i] == current:
            good[i] = 1
    #backwwad search
    current = -1
    for i in reversed(range(n)):
        if score[i] >= minimum:
            good[i] = 1
            current = anc[i]
        elif anc[i] == current:
            good[i] = 1
        
        return good

def reduce_to_tracts(anc):
    start = 0
    starts = [0]
    ends = []
    for i in range(1, anc.shape[1]):
        if not np.all(anc[:,i] == anc[:,start]):
            ends.append(i-1)
            start = i
            starts.append(i)
    
    ends.append(i)
    n_tracts = len(starts)
    
    tract_anc = np.zeros(shape=(anc.shape[0], n_tracts), dtype=int)
    
    for i in range(n_tracts):
        tract_anc[:,i] = anc[:,starts[i]]
    
    return (tract_anc, starts, ends)


### parse arguments

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input_file", help="Input file giving derived allele counts per individual", action = "store", required=True)
parser.add_argument("-o", "--output", help="Output file containing ancestry calls", action = "store", required=True)

#reference populations
parser.add_argument("-g", "--group", help="Group name and optionally sample names (separated by commas)",
                    required = True, action='append', nargs="+", metavar=("group_name","[samples]"))
parser.add_argument("--groups_file", help="Optional file of sample names and the group they belong to", action = "store", required = False)
parser.add_argument("--haploidify_IDs", help="Split each reference ID provided into '_A' and '_B' versions", action='store_true')

parser.add_argument("--min_score", help="Minimum score to set ancestry calls as indecisive (0). Maximum is 160.", type=int, action = "store", default = 0)

parser.add_argument("--extend", help="Extend outward from decisive runs if ancestry call is the same", action='store_true')

parser.add_argument("-t", "--threads", help="Number of processors to use", type=int, action = "store", default = 1)


args = parser.parse_args()

#input_file = "test.haplo.gz"

input_file = args.input_file

with gzip.open(input_file, "rt") as f:
    IDs = f.readline().split()[2:]

n = len(IDs)


assert len(args.group) >= 2, "You must specify at least two reference groups."

ref_IDs = {}
for g in args.group:
    ref_IDs[g[0]] = g[1].split(",") if len(g) > 1 else []

if args.groups_file:
    with open(args.groups_file, "rt") as gf: groupDict = dict([ln.split() for ln in gf])
    for ID in groupDict.keys():
        try: ref_IDs[groupDict[ID]].append(ID)
        except: pass

for ref in ref_IDs:
    assert len(ref_IDs[ref]) >= 1, f"Reference group {ref} is empty."

#ref_IDs = {"orientis":["SM21TS03","SM21TS05","SM21TS16","SM21TS33","SM21TS34","SM21TS35","SM21TS36"],
        #"klugii":["SM15W61","SM15W66","SM15W69","SM15W72","SM15W74","SM18W01"],
        #"chrysippus":["SM16N01","SM16N04","SM16N05","SM16N06","SM16N20","SM16N37"]}

if args.haploidify_IDs:
    ref_hap_IDs = dict([(name, [label + suffix for label in ref_IDs[name] for suffix in ["_A", "_B"]],) for name in ref_IDs])
else:
    ref_hap_IDs = ref_IDs

for ref in ref_hap_IDs:
    for ID in ref_hap_IDs[ref]:
        assert ID in IDs, f"{ID} not in file header line. Did you mistype, or do you need to add the --haploidify_IDs option?"


#load allele counts
print("\nLoading complete haplotype file into memory...", file=sys.stderr)
alleles = np.loadtxt(input_file, skiprows=1, usecols = list(range(2, n+2)), dtype=int)



ref_hap_indices = dict([(name, [IDs.index(ID) for ID in ref_hap_IDs[name]],) for name in ref_hap_IDs])

ref_haps_list = [alleles[:,ref_hap_indices[name]].transpose() for name in ref_hap_IDs]

all_haps = alleles.transpose()

del alleles

## Loter with bagging only (phase correction not possible when there are >2 sources)
print(f"\nRunning loter inference using {args.threads} threads...", file=sys.stderr)
threads = args.threads
#threads = 8
res_loter = lc.loter_local_ancestry(l_H=ref_haps_list, h_adm=all_haps, num_threads=threads, rate_vote=0, nb_bagging=20)

del ref_haps_list

del all_haps

ancestry = res_loter[0] + 1

if args.min_score > 0:
    assert res_loter[1].max() >= args.min_score, "None of the ancestry scores ae above the specified minimum. Try a lower value."
    if args.extend:
        print("\nCleaning and extending ancestry calls...", file=sys.stderr)
        
        onePercent = int(round(n/100))
        
        for i in range(n):
            if i % onePercent == 0: print(".", end="", file=sys.stderr, flush=True)
            good = extend_good_calls(ancestry[i], res_loter[1][i], minimum=args.min_score)
            not_good = np.where(np.logical_not(good))[0]
            ancestry[i, not_good] = 0
    else:
        ancestry[res_loter[1] < args.min_score] = 0


print("\nReducing to tracts of identical ancestry calls...", file=sys.stderr)

tract_anc, tract_start, tract_end = reduce_to_tracts(ancestry)

position = np.loadtxt(input_file, skiprows=1, usecols = 1, dtype=int)

print("\nWriting output...", file=sys.stderr)

out_array = np.column_stack([position[tract_start], position[tract_end], tract_anc.transpose()])

out_array = np.row_stack([np.concatenate([["start","end"], IDs]), out_array.astype(str)])

np.savetxt(args.output, out_array, delimiter = "\t", fmt='%s')
