import numpy as np
import time
import math
 

green_exon_pos = np.array([
    [149288166, 149288277], # G1
    [149293258, 149293554], # G2
    [149295542, 149295710], # G3
    [149297178, 149297343], # G4
    [149298898, 149299137], # G5
    [149301420, 149301530]  # G6
    ])


red_exon_pos = np.array([
    [149249757, 149249868], # R1
    [149256127, 149256423], # R2
    [149258412, 149258580], # R3
    [149260048, 149260213], # R4
    [149261768, 149262007], # R5
    [149264290, 149264400]  # R6
    ])


def load_seq_ref(filename):
     
    seq_ref = ''.join(np.loadtxt(filename, dtype=str)[1:])
    return seq_ref  


def load_last_cols(filename):
     
    last_column = ''.join(np.loadtxt(filename, dtype=str))
    return last_column 


def load_map_seq(filename):
     
    map_seq_ref = np.loadtxt(filename, dtype=int)
    return map_seq_ref 


def load_reads(filename):

    _reads = np.loadtxt(filename, dtype=str)
    return _reads   


def srting_nparray(string):    
    return np.fromiter(string, (str,1))



def combinations(n, r):
    r = min(r, n - r)
    n = math.prod(range(n, n - r, -1))
    d = math.prod(range(1, r + 1))
    return n/d



def matching_read_locations(read):
    positions = []

    read_length = len(read)
    parts = [read[:read_length // 3], read[read_length // 3: 2 * (read_length // 3)], read[2 * (read_length // 3):]]

    for part in parts:
        exact_match_indices = calc_matching(part)
        positions.extend(exact_match_indices)

    positions = list(set(positions))

    valid_positions = []
    for position in positions:
        if position + read_length <= len(seq_ref):
            reference_part = seq_ref[position: position + read_length]
            mismatch_count = np.sum(srting_nparray(reference_part) != srting_nparray(read))
            if mismatch_count <= 2:
                valid_positions.append(position)

    return valid_positions



def calc_matching(read):
    read = read.replace("N", "A")   
    batch = first_cols

    try:
        batch_start = batch.index(read[-1])
        batch_end = batch.rindex(read[-1])
    except ValueError:
        return []

    for i in range(2, len(read) + 1):
        last = last_column[batch_start: batch_end + 1]

        try:
            left_char = read[-i]
            sub_batch_start = Rank[batch_start + last.index(left_char)]
            sub_batch_end = Rank[batch_start + last.rindex(left_char)]
        except (ValueError, IndexError):
            return []

        batch_start = sub_batch_start 
        batch_end = sub_batch_end
        batch = first_cols[batch_start: batch_end + 1]

    return [_map[i] for i in range(batch_start, batch_end + 1)]



def compute_probabilities(exon_match_counts):
    red_exon_counts = [int(c) for c in exon_match_counts[:6]]
    green_exon_counts = [int(c) for c in exon_match_counts[6:]]
    total_counts = [r + g for r, g in zip(red_exon_counts, green_exon_counts)]

    configurations = [
        ([1/3, 1/3, 1/3, 1/3], [2/3, 2/3, 2/3, 2/3]),
        ([1/2, 1/2, 0, 0], [1/2, 1/2, 1.0, 1.0]),
        ([1/4, 1/4, 1/2, 1/2], [3/4, 3/4, 1/2, 1/2]),
        ([1/4, 1/4, 1/4, 1/2], [3/4, 3/4, 3/4, 1/2])
    ]

    probabilities = []
    for config_red, config_green in configurations:
        probability = 1.0
        for i in range(len(config_red)):
            probability *= combinations(total_counts[i + 1], red_exon_counts[i + 1]) * \
                pow(config_red[i], red_exon_counts[i + 1]) * \
                pow(config_green[i], green_exon_counts[i + 1])
        probabilities.append(probability)

    return probabilities



def particular_exon(positions):
    r_counts = [0.0] * 6
    g_counts = [0.0] * 6
    red_gene_positions = red_exon_pos
    green_gene_positions = green_exon_pos

    for j in positions:
        for i in range(6):
            red_overlap = red_gene_positions[i][0] <= j <= red_gene_positions[i][1]
            green_overlap = green_gene_positions[i][0] <= j <= green_gene_positions[i][1]
            
            if red_overlap and green_overlap:
                r_counts[i] += 0.5
                g_counts[i] += 0.5
            elif red_overlap:
                r_counts[i] += 1.0
            elif green_overlap:
                g_counts[i] += 1.0

    return np.array(r_counts + g_counts)



def best_match(list_prob):
     
    likely_config = np.asarray(list_prob).argmax()
    return likely_config  



if __name__ == "__main__":
    start = time.time()
    
    seq_ref = load_seq_ref("../data/chrX.fa") 
    last_column = load_last_cols("../data/chrX_last_col.txt")   
    _map = load_map_seq("../data/chrX_map.txt")  
    _reads = load_reads("../data/reads") 
    print("LOADED THE DATA")

    a = last_column.count('A')
    c = last_column.count('C')
    g = last_column.count('G')
    t = last_column.count('T')
    first_cols = 'A'*a + 'C'*c + 'G'*g + 'T'*t +'$'  
    
    A = 0 
    C =  a 
    G = a+c
    T = a+c+g
    Rank = []
    for i in last_column:
        if i in {'A', 'C', 'G', 'T'}:
            Rank.append({'A': A, 'C': C, 'G': G, 'T': T}[i])
            A += (i == 'A')
            C += (i == 'C')
            G += (i == 'G')
            T += (i == 'T')
        elif i == '$':
            Rank.append(len(first_cols) - 1)

    print("CALCULATED THE RANK")

    exon_count = np.zeros(12)  
    for read in _reads[2936000:2948000]:  
        positions = matching_read_locations(read)  
        exon_count += particular_exon(positions)  
    print("CALCULATED THE EXONS")
    print("EXONS MATCH COUNT : ", exon_count)

    list_prob = compute_probabilities(exon_count)  
    print("PROBABILITIES : ", list_prob)

    MostLikely = best_match(list_prob)  
    print("CONFIGURATION %d IS THE BEST MATCH" %MostLikely)

    end = time.time()
    print(f"TOTAL TIME ELAPSED = {end-start:0.2f} s.\n")