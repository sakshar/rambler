from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import csv
import numpy as np
import sys
import os
import time

def create_directory(path):
    if not os.path.exists(path):
        os.makedirs(path)

def get_unikmer_map(unikmer_file):
    unikmer_map = dict()
    counter = 0
    while True:
        line = unikmer_file.readline()
        if not line:
            break
        unikmer_map[line.strip()] = counter
        counter += 1
    return unikmer_map


def write_unikmer_map(out_file, unikmer_map):
    output = open(out_file, "w")
    for unikmer in unikmer_map.keys():
        output.write(str(unikmer) + "," + str(unikmer_map[unikmer]) + "\n")
    output.close()


def get_read_to_unikmer_and_unikmer_to_read_maps(unikmer_map, read_file, k):
    read_to_unikmer_map = dict()
    unikmer_to_read_map = dict()
    unikmers = unikmer_map.keys()
    for record in read_file:
        read = record.seq
        id = record.id
        for i in range(len(read) - k + 1):
            kmer = read[i:i + k]
            str_kmer = str(kmer)
            str_rev_kmer = str(kmer.reverse_complement())
            if str_kmer in unikmers:
                if id not in read_to_unikmer_map.keys():
                    read_to_unikmer_map[id] = []
                read_to_unikmer_map[id] += [(unikmer_map[str_kmer], i)]  # (unikmer_id, pos)
                if unikmer_map[str_kmer] not in unikmer_to_read_map.keys():
                    unikmer_to_read_map[unikmer_map[str_kmer]] = []
                unikmer_to_read_map[unikmer_map[str_kmer]] += [(id, i)]  # (read_id, pos)
            elif str_rev_kmer in unikmers:
                if id not in read_to_unikmer_map.keys():
                    read_to_unikmer_map[id] = []
                read_to_unikmer_map[id] += [(unikmer_map[str_rev_kmer], i)]  # (unikmer_id, pos)
                if unikmer_map[str_rev_kmer] not in unikmer_to_read_map.keys():
                    unikmer_to_read_map[unikmer_map[str_rev_kmer]] = []
                unikmer_to_read_map[unikmer_map[str_rev_kmer]] += [(id, i)]  # (read_id, pos)
    return read_to_unikmer_map, unikmer_to_read_map


def write_read_to_unikmer_and_unikmer_to_read_maps(out_r2u, out_u2r, read_to_unikmer_map, unikmer_to_read_map):
    out_r2u_file = open(out_r2u, "w")
    out_u2r_file = open(out_u2r, "w")
    for read in read_to_unikmer_map.keys():
        out_r2u_file.write(str(read) + ">")
        for tup in read_to_unikmer_map[read]:
            out_r2u_file.write(str(tup) + ";")
        out_r2u_file.write("\n")
    out_r2u_file.close()
    for unikmer in unikmer_to_read_map.keys():
        out_u2r_file.write(str(unikmer) + ">")
        for tup in unikmer_to_read_map[unikmer]:
            out_u2r_file.write(str(tup) + ";")
        out_u2r_file.write("\n")
    out_u2r_file.close()


def get_read_to_unikmer_map(infile):
    read_to_unikmer_map = dict()
    read_to_unikmer_list = dict()
    read_to_unikmerpos_map = dict()
    while True:
        line = infile.readline()
        if not line:
            break
        words = line.strip().split(">")
        read = int(words[0].split("_")[1])
        read_to_unikmer_map[read] = dict()
        read_to_unikmer_list[read] = []
        read_to_unikmerpos_map[read] = dict()
        pairs = words[1].split(";")
        pos = 0
        for pair in pairs[:-1]:
            values = pair[1:-1].split(",")
            read_to_unikmer_map[read][int(values[0])] = int(values[1])
            read_to_unikmer_list[read] += [(int(values[0]), int(values[1]))]
            read_to_unikmerpos_map[read][int(values[1])] = pos
            pos += 1
    return read_to_unikmer_map, read_to_unikmer_list, read_to_unikmerpos_map


def get_unikmer_to_read_map(infile):
    unikmer_to_read_map = dict()
    while True:
        line = infile.readline()
        if not line:
            break
        words = line.strip().split(">")
        unikmer_to_read_map[int(words[0])] = dict()
        pairs = words[1].split(";")
        for pair in pairs[:-1]:
            values = pair[1:-1].split(",")
            unikmer_to_read_map[int(words[0])][int(values[0].split("_")[1][:-1])] = int(values[1])
    return unikmer_to_read_map


# older version of get_pairwise_profiles method
# def get_pairwise_profiles(unikmer_to_read_map, read_to_unikmer_map, read_to_unikmer_list, read_to_unikmerpos_map, tolerance):
#     pairwise_profiles = dict()
#     unikmers = unikmer_to_read_map.keys()
#     count = 0
#     for unikmer in unikmers:
#         reads = list(unikmer_to_read_map[unikmer].keys())
#         for i in range(len(reads)):
#             r1 = reads[i]
#             for j in range(i+1, len(reads)):
#                 r2 = reads[j]
#                 pos1, pos2 = read_to_unikmer_map[r1][unikmer], read_to_unikmer_map[r2][unikmer]
#                 pairwise_profiles[(r1, r2)] = [(unikmer, pos1, pos2)]
#                 unikmer_index1, unikmer_index2 = read_to_unikmerpos_map[r1][pos1], read_to_unikmerpos_map[r2][pos2]
#                 r1_list, r2_list = read_to_unikmer_list[r1], read_to_unikmer_list[r2]
#                 while unikmer_index1 + 1 < len(r1_list) and unikmer_index2 + 1 < len(r2_list):
#                     u1, u2 = r1_list[unikmer_index1+1][0], r2_list[unikmer_index2+1][0]
#                     delta1, delta2 = r1_list[unikmer_index1+1][1] - pos1, r2_list[unikmer_index2+1][1] - pos2
#                     if u1 == u2 and delta1 - tolerance <= delta2 <= delta1 + tolerance:
#                         unikmer_index1 += 1
#                         unikmer_index2 += 1
#                         pos1, pos2 = r1_list[unikmer_index1][1], r2_list[unikmer_index2][1]
#                         pairwise_profiles[(r1, r2)] += [(u1, pos1, pos2)]
#                     else:
#                         break
#         count += 1
#     return pairwise_profiles


# newer version of get_pairwise_profiles method
def get_pairwise_profiles(unikmer_to_read_map, read_to_unikmer_map, read_to_unikmer_list, read_to_unikmerpos_map, tolerance):
    pairwise_profiles = dict()
    reads = list(read_to_unikmer_map.keys())
    count = 0
    for i in range(len(reads)):
        r1 = reads[i]
        r1_unikmers = set(read_to_unikmer_map[r1].keys())
        for j in range(i+1, len(reads)):
            r2 = reads[j]
            r2_unikmers = set(read_to_unikmer_map[r2].keys())
            intersection = r1_unikmers & r2_unikmers
            for unikmer in intersection:
                pos1, pos2 = read_to_unikmer_map[r1][unikmer], read_to_unikmer_map[r2][unikmer]
                pairwise_profiles[(r1, r2)] = [(unikmer, pos1, pos2)]
                unikmer_index1, unikmer_index2 = read_to_unikmerpos_map[r1][pos1], read_to_unikmerpos_map[r2][pos2]
                r1_list, r2_list = read_to_unikmer_list[r1], read_to_unikmer_list[r2]
                while unikmer_index1 + 1 < len(r1_list) and unikmer_index2 + 1 < len(r2_list):
                    u1, u2 = r1_list[unikmer_index1+1][0], r2_list[unikmer_index2+1][0]
                    delta1, delta2 = r1_list[unikmer_index1+1][1] - pos1, r2_list[unikmer_index2+1][1] - pos2
                    if u1 == u2 and delta1 - tolerance <= delta2 <= delta1 + tolerance:
                        unikmer_index1 += 1
                        unikmer_index2 += 1
                        pos1, pos2 = r1_list[unikmer_index1][1], r2_list[unikmer_index2][1]
                        pairwise_profiles[(r1, r2)] += [(u1, pos1, pos2)]
                    else:
                        break
                pos1, pos2 = read_to_unikmer_map[r1][unikmer], read_to_unikmer_map[r2][unikmer]
                unikmer_index1, unikmer_index2 = read_to_unikmerpos_map[r1][pos1], read_to_unikmerpos_map[r2][pos2]
                while unikmer_index1 - 1 >= 0 and unikmer_index2 - 1 >= 0:
                    u1, u2 = r1_list[unikmer_index1-1][0], r2_list[unikmer_index2-1][0]
                    delta1, delta2 = r1_list[unikmer_index1-1][1] - pos1, r2_list[unikmer_index2-1][1] - pos2
                    if u1 == u2 and delta1 - tolerance <= delta2 <= delta1 + tolerance:
                        unikmer_index1 -= 1
                        unikmer_index2 -= 1
                        pos1, pos2 = r1_list[unikmer_index1][1], r2_list[unikmer_index2][1]
                        pairwise_profiles[(r1, r2)] = [(u1, pos1, pos2)] + pairwise_profiles[(r1, r2)]
                    else:
                        break
                count += 1
                break
    return pairwise_profiles


def write_pairwise_profiles(pairwise_profiles, file_path):
    output_pairwise_profiles = open(file_path, "w")
    for pairs in pairwise_profiles.keys():
        output_pairwise_profiles.write("(" + str(pairs[0]) + "," + str(pairs[1]) + ")>")
        for profile in pairwise_profiles[pairs]:
            output_pairwise_profiles.write(
                "(" + str(profile[0]) + "," + str(profile[1]) + "," + str(profile[2]) + ");")
        output_pairwise_profiles.write("\n")
    output_pairwise_profiles.close()


def read_pairwise_profiles(infile):
    pairwise_profiles = dict()
    while True:
        line = infile.readline()
        if not line:
            break
        key_values = line.strip().split(">")
        read_pair = key_values[0][1:-1].split(",")
        pairwise_profiles[(int(read_pair[0]), int(read_pair[1]))] = []
        triplets = key_values[1].split(";")
        for triplet in triplets[:-1]:
            values = triplet[1:-1].split(",")
            pairwise_profiles[(int(read_pair[0]), int(read_pair[1]))] += [(int(values[0]), int(values[1]), int(values[2]))]
    return pairwise_profiles


# unikmer_delta = unikmerPos_r1 - unikmerPos_r2, true_delta = start_r1 - start_r2
def get_pairwise_location_data(pw_profiles, threshold):
    overlap_locations = [["r1", "r2", "unikmer_delta_mean", "unikmer_delta_sd", "#_of_unikmers"]]
    read_set = set()
    for key in pw_profiles.keys():
        r1, r2 = key[0], key[1]
        shared_unikmers = pw_profiles[(r1, r2)]
        unikmer_delta_mean = 0
        no_shared_unikmer = len(shared_unikmers)
        unikmer_deltas = np.zeros(no_shared_unikmer)
        for i in range(no_shared_unikmer):
            unikmer_delta = shared_unikmers[i][1] - shared_unikmers[i][2]
            unikmer_delta_mean += unikmer_delta/no_shared_unikmer
            unikmer_deltas[i] = unikmer_delta
        unikmer_delta_sd = np.sqrt(np.sum((unikmer_deltas - unikmer_delta_mean)**2)/no_shared_unikmer)
        if len(shared_unikmers) > threshold:
            read_set.add(r1)
            read_set.add(r2)
            overlap_locations.append([r1, r2, np.round(unikmer_delta_mean, 2), np.round(unikmer_delta_sd, 6), no_shared_unikmer])
    return overlap_locations


def pairwise_location_data_writer(outfile, overlap_locations):
    with open(outfile, "w") as csvfile:
        writer = csv.writer(csvfile, delimiter="\t")
        [writer.writerow(r) for r in overlap_locations]


class DisjointSet:
    def __init__(self, vertices, parent):
        self.vertices = vertices
        self.parent = parent

    def find(self, item):
        if self.parent[item] == item:
            return item
        else:
            res = self.find(self.parent[item])
            self.parent[item] = res
            return res

    def union(self, set1, set2):
        root1 = self.find(set1)
        root2 = self.find(set2)
        self.parent[root1] = root2


def csv_reader(infile):
    data = []
    with open(infile, newline='') as csvfile:
        reader = csv.reader(csvfile, delimiter='\t')
        for row in reader:
            data.append(row)
    return data


def get_read_set(data):
    reads = set()
    for i in range(1, len(data)):
        r1, r2 = int(data[i][0]), int(data[i][1])
        reads.add(r1)
        reads.add(r2)
    return reads


def convert_data_to_float(data):
    new_data = []
    size = len(data)
    for i in range(1, size):
        new_row = [float(r) for r in data[i]]
        new_data.append(new_row)
    return new_data


def get_disjoint_set(reads, data):
    data = convert_data_to_float(data)
    parent = {}
    for r in reads:
        parent[r] = r
    ds = DisjointSet(reads, parent)
    size = len(data)
    for i in range(size):
        r1, r2 = int(data[i][0]), int(data[i][1])
        ds.union(r1, r2)
    return ds


def get_clusters(reads, ds):
    clusters = dict()
    for r in reads:
        parent = ds.find(r)
        if parent not in clusters:
            clusters[parent] = [r]
        else:
            clusters[parent] += [r]
    return clusters


def write_clusters(outfile, clusters):
    outfile.write(str(len(clusters.keys()))+"\n")
    for key in clusters.keys():
        outfile.write(','.join(str(i) for i in clusters[key]))
        outfile.write("\n")


def read_clusters(cluster_file):
    no_of_clusters = int(cluster_file.readline().strip())
    clusters = []
    for i in range(no_of_clusters):
        reads = cluster_file.readline().strip().split(",")
        clusters.append([int(j) for j in reads])
    return clusters, no_of_clusters


def write_read_clusters_to_fasta(cluster_file, read_file, path_to_output):
    clusters, no_of_clusters = read_clusters(cluster_file)
    read_dict = dict()
    for record in read_file:
        id = int(str(record.id).split("_")[1])
        read_dict[id] = record.seq
    for i in range(no_of_clusters):
        seqs = []
        for j in range(len(clusters[i])):
            record = SeqRecord(Seq(read_dict[clusters[i][j]]), id=str(clusters[i][j]),
                               description="cluster_"+str(i)+"_read_"+str(clusters[i][j]))
            seqs.append(record)
        SeqIO.write(seqs, path_to_output + str(i) + ".fasta", "fasta")


def get_clusters_from_unikmers_and_reads(unikmer_file, read_file1, read_file2, k, out_dir, tolerance, threshold):
    unikmer_map = get_unikmer_map(unikmer_file)
    out_unikmer_map = out_dir+ "/intermediates/unikmer_map.txt"
    write_unikmer_map(out_unikmer_map, unikmer_map)
    r2u, u2r = get_read_to_unikmer_and_unikmer_to_read_maps(unikmer_map, read_file1, k)
    out_r2u = out_dir+ "/intermediates/read_to_unikmer_map.txt"
    out_u2r = out_dir+ "/intermediates/unikmer_to_read_map.txt"
    write_read_to_unikmer_and_unikmer_to_read_maps(out_r2u, out_u2r, r2u, u2r)
    print("status: --- unikmer barcoding done ---")
    in_r2u = open(out_r2u, "r")
    in_u2r = open(out_u2r, "r")
    r2u_map, r2u_list, r2up_map = get_read_to_unikmer_map(in_r2u)
    u2r_map = get_unikmer_to_read_map(in_u2r)
    pw_profiles = get_pairwise_profiles(u2r_map, r2u_map, r2u_list, r2up_map, tolerance)
    out_pw_profiles = out_dir+ "/intermediates/pairwise_profiles.txt"
    write_pairwise_profiles(pw_profiles, out_pw_profiles)
    print("status: --- pairwise profiling done ---")
    overlap_data = get_pairwise_location_data(pw_profiles, threshold)
    out_overlap = out_dir+ "/intermediates/overlap.csv"
    pairwise_location_data_writer(out_overlap, overlap_data)
    reads = list(get_read_set(overlap_data))
    out_path = out_dir + "/clusters/"
    ds = get_disjoint_set(reads, overlap_data)
    clusters = get_clusters(reads, ds)
    out_path = out_path + "cluster_log.txt"
    outfile = open(out_path, "w")
    write_clusters(outfile, clusters)
    outfile.close()
    cluster_file = open(out_path, "r")
    out_clusters = out_dir + "/clusters/"
    write_read_clusters_to_fasta(cluster_file, read_file2, out_clusters)
    print("status: --- clustering done ---")


def rambler_run():

    print("----------------- RAmbler initializing -----------------")
    print("read file:", sys.argv[1])
    print("unikmer file:", sys.argv[2])
    print("unikmer size:", sys.argv[5])
    print("estimated reference size:", sys.argv[4])
    print("tolerance:", sys.argv[6])
    print("threshold:", sys.argv[7])
    print("output directory:", sys.argv[3])
    create_directory(sys.argv[3]+"/assembly/hifiasm")
    start_time = time.time()
    unikmer_file = open(sys.argv[2], "r")
    read_file1 = SeqIO.parse(sys.argv[3]+"/intermediates/reads.fasta", "fasta")
    read_file2 = SeqIO.parse(sys.argv[3]+"/intermediates/reads.fasta", "fasta")
    get_clusters_from_unikmers_and_reads(unikmer_file, read_file1, read_file2, int(sys.argv[5]), sys.argv[3], int(sys.argv[6]), int(sys.argv[7]))
    end_time = time.time()
    print("Runtime:", end_time - start_time, "seconds")


if __name__ == "__main__":
    rambler_run()