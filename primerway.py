#!/usr/bin/env python3
import configparser
import getopt
import gzip
import json
import math
import os
import re
import subprocess
import sys
import primer3
import pysam
import vcf


def get_pseudo_qualities_from_VCF(VCF_file, target_region):
    pseudo_qualities = []
    if VCF_file != "":
        vcf_reader = vcf.Reader(filename=VCF_file)
        for i in target_region[2]:
            try:
                variants = vcf_reader.fetch(target_region[1], i, i + 1)
            except:
                try:
                    variants = vcf_reader.fetch(target_region[1].replace("chrM", "chrMT", 1).replace("chr", "", 1),
                                                i, i + 1)
                except:
                    variants = []
            qual = 0.0001
            ref = "-"
            for variant in variants:
                try:
                    ref = variant.REF
                    q = 1.0 - float(variant.INFO["CAF"][0])
                    if q > qual:
                        qual = q
                except:
                    pass
            pseudo_qualities.append(int(-10.0 * math.log10(qual)))
            if verbose:
                sys.stdout.write(ref[0])
        if verbose:
            print("")
            print(len(pseudo_qualities))
            print(pseudo_qualities)
    return pseudo_qualities


def get_candidate_pairs(sequence, left_edge, right_edge, step):
    candidate_pairs = {}
    p3seq = {
        'SEQUENCE_ID': 'example',
        'SEQUENCE_TEMPLATE': ''
    }
    p3prim = dict([(x[0].upper(), json.loads(x[1])) for x in config.items("Primer3")])

    for i in range(left_edge, right_edge, step):
        target_length = right_edge - i
        if target_length > step:
            target_length = step

        p3seq['SEQUENCE_TEMPLATE'] = sequence

        if variant_file != "":
            p3seq['SEQUENCE_QUALITY'] = l_quality
        else:
            p3prim['PRIMER_WT_SEQ_QUAL'] = 0
            p3prim['PRIMER_WT_END_QUAL'] = 0

        if as_coord == "":
            if deletion == "":
                p3seq['SEQUENCE_TARGET'] = [i, target_length]
            else:
                p3seq['SEQUENCE_TARGET'] = [config.getint("PrimerWay", "flanking") + 1,
                                            config.getint("PrimerWay", "min_overlap") * 2 + 1]
        else:
            if exon < 2:
                p3seq['SEQUENCE_FORCE_LEFT_END'] = 50
            else:
                p3seq['SEQUENCE_FORCE_RIGHT_END'] = 500

        result = primer3.bindings.designPrimers(p3seq, p3prim)

        for i in range(p3prim['PRIMER_NUM_RETURN']):
            si = str(i)
            spair = (
                result.get('PRIMER_LEFT_' + si + '_SEQUENCE', ''), result.get('PRIMER_RIGHT_' + si + '_SEQUENCE', ''))
            if spair[0] != '':
                candidate_pairs[spair] = [
                    result['PRIMER_LEFT_' + si][0] + result['PRIMER_LEFT_' + si][1] - 1,
                    result['PRIMER_RIGHT_' + si][0] - result['PRIMER_RIGHT_' + si][1] + 1,
                    result['PRIMER_PAIR_' + si + '_PENALTY'],
                    0
                ]

        for pair in candidate_pairs:
            if (pair[0][-1] in set("ATat")) | (pair[1][-1] in set("ATat")):
                candidate_pairs[pair][2] += 1
    return candidate_pairs


def add_penalty_for_non_specific_pairs(candidate_pairs):
    s = ""
    i = 0
    complete = 0
    for pair in candidate_pairs:
        i += 1
        complete += 1
        s += "_".join(pair) + "\t" + pair[0] + "\t" + pair[1] + "\n"
        if (i >= config.getint("tntBLAST", "pairs_per_run")) | (pair == list(candidate_pairs.keys())[-1]):
            if verbose: print(s)
            print("Tntblast: " + str(complete) + " of " + str(len(candidate_pairs)))
            proc = subprocess.Popen(
                (config.get("tntBLAST", "run_command") + " -i /dev/stdin -d " + reference_file).split(),
                stdout=subprocess.PIPE, stdin=subprocess.PIPE
            )
            proc.stdin.write(s.encode())
            s = ""
            i = 0
            proc.stdin.close()
            result = proc.stdout.read().decode()
            proc.wait()
            pcrs = re.findall("name = (\w+)", result)
            for pcr in pcrs:
                candidate_pairs[tuple(pcr.split('_'))][3] += 1
                candidate_pairs[tuple(pcr.split('_'))][2] += config.getfloat("PrimerWay", "penalty_PCR_product")
    return candidate_pairs


def get_the_best_way(penalised_pairs, left_edge, right_edge):
    graph = {"START": []}
    for pair1 in penalised_pairs:
        if penalised_pairs[pair1][0] <= left_edge: graph["START"].append(pair1)
        if penalised_pairs[pair1][1] >= right_edge:
            graph[pair1] = graph.get(pair1, [])
            graph[pair1].append("END")
        for pair2 in penalised_pairs:
            if penalised_pairs[pair2][0] < (penalised_pairs[pair1][1] - config.getint("PrimerWay", "min_overlap")):
                graph[pair1] = graph.get(pair1, [])
                graph[pair1].append(pair2)

    penalised_pairs["END"] = [0, 0, 0, 0]
    visited = {}
    to_visit = {"START": 0}
    paths = {"START": ["START"]}
    while to_visit:
        v = min(to_visit, key=lambda x: to_visit[x])
        visited[v] = to_visit[v]
        del to_visit[v]
        for w in graph.get(v, []):
            if w not in visited:
                vwLength = (visited[v] + ((penalised_pairs[w][2] - (
                        (deletion == "") * config.getfloat("PrimerWay", "penalty_PCR_product"))) ** 2) +
                            config.getfloat("PrimerWay", "base_penalty_for_pair"))
                if (w not in to_visit) or (vwLength < to_visit[w]):
                    to_visit[w] = vwLength
                    paths[w] = paths[v] + [w]
    try:
        the_best_way = paths["END"][1:-1]
    except:
        the_best_way = []
    return the_best_way


try:
    opts, args = getopt.getopt(sys.argv[1:], 'hvG:R:V:i:o:p:r:a:n:s:d:c:')
except:
    opts, args = [], []
if len(opts) == 0: print("Type -h for help")
text_width = 60
verbose = False
job_name = "Exon"
protein_id = ""
user_region = ""
deletion = ""
as_coord = ""
output_d = ""
input_d = ""
config_file = "primerway.cfg"
reference_file = ""
reference_gff_file = ""
variant_file = ""
target_regions = []
searches_coord = []
start_exon = 1

for o, a in opts:
    if o == "-v":
        verbose = True
    elif o in ("-h", "--help"):
        print("-R : reference unzipped fasta file")
        print("-G : reference gzipped GFF description file")
        print("-V : variants gzipped VCF dbSNP file")
        print("-v : verbose mode for debug")
        print("-i : input directory")
        print("-o : output directory")
        print("-p : protein ID for automatic exon CDS extraction (example NP_000442.1)")
        print("-d : deletion region coordinates (example: 11:121212-454545)")
        print("-r : target region coordinates (example: 11:121212-454545)")
        print("-a : target nucleotide coordinate for allele specific (example: 11:121212)")
        print("-n : prefix for primer names")
        print("-s : start from exon number (default 1)")
        print("-c : set configuration file (default primerway.cfg)")
        sys.exit()

    elif o in "-s":
        try:
            start_exon = int(a)
        except:
            pass
    elif o in "-p":
        protein_id = a
    elif o in "-r":
        user_region = a
    elif o in "-d":
        deletion = a
    elif o in "-a":
        as_coord = a
    elif o in "-o":
        output_d = a
    elif o in "-c":
        config_file = a
    elif o in "-i":
        input_d = a
    elif o in "-R":
        reference_file = a
    elif o in "-G":
        reference_gff_file = a
    elif o in "-V":
        variant_file = a
    elif o in "-n":
        job_name = a

config = configparser.ConfigParser()
config.read(config_file)

if output_d == "":
    print("Error: Output directory is needed")
    sys.exit()
try:
    os.stat(output_d)
except:
    try:
        os.mkdir(output_d)
    except:
        print("Error: Can not create output directory")

        sys.exit()

if (as_coord != "") & (reference_file != ""):
    exons = []
    regreg = re.findall("(.+):(\d+)", as_coord)
    with pysam.FastaFile(reference_file) as fh:
        try:
            start = int(regreg[0][1])
            end = int(regreg[0][1])
        except:
            start = 0
            end = 0
        target_regions.append(
            [fh.fetch(reference=regreg[0][0], start=(start - 50), end=(end + 500)), regreg[0][0],
             range(start - 50, end + 500)])
        target_regions.append(
            [fh.fetch(reference=regreg[0][0], start=(start - 500), end=(end + 50)), regreg[0][0],
             range(start - 500, end + 50)])
    protein_id = ""
    user_region = ""

if (deletion != "") & (reference_file != ""):
    exons = []
    regreg = re.findall("(.+):(\d+)-(\d+)", deletion)
    with pysam.FastaFile(reference_file) as fh:
        try:
            start = int(regreg[0][1]) - 1
            end = int(regreg[0][2])
        except:
            start = 0
            end = 0
        target_regions.append([
            fh.fetch(reference=regreg[0][0], start=(
                    start - config.getint("PrimerWay", "flanking") - config.getint("PrimerWay", "min_overlap") - 0),
                     end=(start + 1)) +
            fh.fetch(reference=regreg[0][0], start=end - 1,
                     end=(end + config.getint("PrimerWay", "flanking") + config.getint("PrimerWay", "min_overlap"))),
            regreg[0][0],
            list(range(start - config.getint("PrimerWay", "flanking") - config.getint("PrimerWay", "min_overlap") - 0,
                       start + 1)) +
            list(range(end -1,
                       end + config.getint("PrimerWay", "flanking") + config.getint("PrimerWay", "min_overlap")))
        ])
    protein_id = ""
    user_region = ""

if (user_region != "") & (reference_file != ""):
    exons = []
    regreg = re.findall("(.+):(\d+)-(\d+)", user_region)
    with pysam.FastaFile(reference_file) as fh:
        try:
            start = int(regreg[0][1]) - config.getint("PrimerWay", "flanking") - config.getint("PrimerWay",
                                                                                               "min_overlap") - 1
            end = int(regreg[0][2]) + config.getint("PrimerWay", "flanking") + config.getint("PrimerWay", "min_overlap")
        except:
            start = 0
            end = 0
        target_regions.append(
            [fh.fetch(reference=regreg[0][0], start=(start), end=(end)), regreg[0][0], range(start, end)])
    protein_id = ""

if (protein_id != "") & (reference_gff_file != ""):
    with gzip.open(reference_gff_file, "r") as csvfile:
        print("Searching exons for transcript ID = " + protein_id + " ...")
        exons = re.findall("(chr[\dXY]+)\t.+\tCDS\t(\d+)\t(\d+)\t.+protein_id=" + protein_id, csvfile.read().decode())
    print("Found " + str(len(exons)) + " exons.")
    fex = open(output_d + "/Exon_list.txt", 'w')
    fex.write("Coord without flanking " + str(config.getint("PrimerWay", "flanking") +
                                              config.getint("PrimerWay", "min_overlap")) + "\n")
    for i in exons:
        fex.write("\t".join(i) + "\n")
        with pysam.FastaFile(reference_file) as fh:
            try:
                start = int(i[1]) - config.getint("PrimerWay", "flanking") - config.getint("PrimerWay",
                                                                                           "min_overlap") - 1
                end = int(i[2]) + config.getint("PrimerWay", "flanking") + config.getint("PrimerWay", "min_overlap")
            except:
                start = 0
                end = 0
            if verbose:
                print(i[0])
            target_regions.append(
                [fh.fetch(reference=i[0], start=(start), end=(end)), i[0], range(start, end)])
    fex.flush()
    fex.write("\n")

fsint = open(output_d + "/resultforsintes.txt", "w")
exon = 0
for target_region in target_regions:
    exon += 1
    if (user_region == "") & (as_coord == ""):
        print("Exon " + str(exon) + " of " + str(len(exons)) + " is assaying...")
    if exon < start_exon:
        print("Skipping...")
        continue
    pairs = {}
    sequence = re.sub('[^ACGTNacgtn]+', '', target_region[0])
    if verbose:
        print(len(sequence))
        print("===")
        print(sequence)
        print("===")
    l_quality = get_pseudo_qualities_from_VCF(variant_file, target_region)

    right_edge = len(sequence) - config.getint("PrimerWay", "flanking") - 1
    left_edge = config.getint("PrimerWay", "flanking")
    if as_coord != "":
        right_edge = left_edge + 1

    pairs = get_candidate_pairs(sequence, left_edge, right_edge, config.getint("PrimerWay", "Primer3_step"))

    print("primers: " + str(len(pairs)))

    pairs = add_penalty_for_non_specific_pairs(pairs)

    print("Computing best way...")
    f = open(output_d + "/Exon_" + str(exon) + '_specprimers.txt', 'w')
    f.write("left edge " + str(left_edge) + ", right edge " + str(right_edge) + "\n")

    for pair in pairs:
        if (pairs[pair][3] == 0) & (deletion == ""):
            pairs[pair][2] += 2 * config.getfloat("PrimerWay", "penalty_PCR_product")

    sorted_pairs = sorted(pairs, key=lambda d: pairs[d][2])
    for pair in sorted_pairs:
        f.write(pair[0] + "\t" + pair[1] + "\t" + "\t" + str(pairs[pair][0]) + "\t" + str(pairs[pair][1]) + "\t" + str(
            pairs[pair][2]) + "\n")
    f.close()

    best_way = get_the_best_way(pairs, left_edge, right_edge)

    f = open(output_d + "/Exon_" + str(exon) + '_resultprimers.txt', 'w')
    f.write(job_name + "_" + str(exon) + "\n")
    f.write("left edge " + str(left_edge) + ", right edge " + str(right_edge) + "\n")
    f.write("The best way consists " + str(len(best_way)) + " pairs.\n\n")
    print("The best way consists " + str(len(best_way)) + " pairs.")
    penalty = 0
    prpair = 0

    if verbose:
        print(left_edge, right_edge, config.getint("PrimerWay", "min_overlap"))
    specchar = dict()
    specchar[left_edge + config.getint("PrimerWay", "min_overlap")] = "{"
    specchar[right_edge - config.getint("PrimerWay", "min_overlap") + 1] = "}"
    for pair in best_way:
        prpair += 1
        if deletion != "":
            penalty += (pairs[pair][2] - 0) ** 2
        else:
            penalty += (pairs[pair][2] - config.getfloat("PrimerWay", "penalty_PCR_product")) ** 2
        f.write(pair[0] + "\t" + pair[1] + "\t" + "\t" + str(pairs[pair][0]) + "\t" + str(pairs[pair][1]) + "\t" + str(
            pairs[pair][2]) + "\n")
        f.flush()
        prname = job_name
        if (pairs[pair][2] >= config.getfloat("PrimerWay", "penalty_PCR_product") * 2) or (
                deletion != "" and pairs[pair][2] >= config.getfloat("PrimerWay", "penalty_PCR_product")):
            prname = "_NS_" + prname
        if user_region == "":
            prname += "_" + str(exon)
        if len(best_way) > 1:
            prname += "_" + str(prpair)
        fsint.write((prname + "f," + pair[0] + "\n").lower())
        fsint.write((prname + "r," + pair[1] + "\n").lower())
        fsint.flush()
        if (protein_id != "") & (reference_gff_file != ""):
            fex.write(prname + "\t" + str(pairs[pair][2]) + "\n")
        fex.flush()
        specchar[pairs[pair][0] - len(pair[0]) + 1] = specchar.get(pairs[pair][0] - len(pair[0]) + 1, "") + "["
        specchar[pairs[pair][0] + 1] = ">" + specchar.get(pairs[pair][0] + 1, "")
        specchar[pairs[pair][1]] = specchar.get(pairs[pair][1], "") + "<"
        specchar[pairs[pair][1] + len(pair[1])] = "]" + specchar.get(pairs[pair][1] + len(pair[1]), "")
        f.write("\n")
        print("Penalty summ sqr = " + str(penalty))
        f.write("Penalty saumm sqr = " + str(penalty) + "\n\n")
        n = 0
        if verbose:
            print(specchar)
        for i in range(0, len(sequence)):
            f.write(specchar.get(i, "") + sequence[i])
        if n > text_width:
            f.write("\n")
        n = 0

        n += 1
        f.write("\n\n{...} - target DNA sequence\n[...> - forward primer\n<...] - reverse primer")
        f.close()
        fsint.close()
        if (protein_id != "") & (reference_gff_file != ""):
            fex.close()
