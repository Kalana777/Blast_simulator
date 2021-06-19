from Bio import SeqIO
import numpy as np
from Bio.SubsMat import MatrixInfo

QUERY_FILE = 'query.fasta'
DB_FILE  = 'target.fasta'
THRESHOLD = 10  #When expanding hits, the current score is allowed to fall to (MAX_SCORE-THRESHOLD)


def parse_db(db_seq):
    genes=[]
    header_dict = {}

    file = SeqIO.parse(db_seq, "fasta")
    for i in file :
        genes.append(split(str(i.seq)))
        header_dict[str(i.id)] = 0
    return genes, header_dict

def split(word):
    return [char for char in word]

def parse_query(input_seq):
    query_seq=[]
    file = SeqIO.parse(input_seq, "fasta")
    for i in file :
        query_seq = split(str(i.seq))

    return query_seq


def hsp(query_words_seq, target_words_seq):

    score_matrix =[]
    for i in range(0, len(query_words_seq)):
        # for each matching word, in both the query and databse sequence, we record the offset and the positions (indexes) of the word
        for j in range(0, len(target_words_seq)):
            if query_words_seq[i] == target_words_seq[j]:
                offset = j - i
                query_index = i
                db_seq_index = j
                score_matrix.append([offset,i,j])

    return score_matrix

def score_match(pair):
    matrix = MatrixInfo.blosum62
    if pair not in matrix:
        return matrix[(tuple(reversed(pair)))]
    else:
        return matrix[pair]

def expand_hit(max_diag, score_matrix, query_seq, gene, threshold):
    max_hits = []
    for hit in score_matrix:
        if hit[0]==max_diag:
            curr_score = 0
            max_score = 0

            # expand to the right
            for j in range(hit[1], len(query_seq)):
                # blossom
                if hit[2]+j-hit[1] >= len(gene):
                    break
                pair = (query_seq[j],gene[hit[2]+j-hit[1]] )
                curr_score += score_match(pair)
                if curr_score>= max_score:
                    max_score = max_score +  score_match(pair)
                elif max_score - curr_score > threshold:
                    break

            # expand to the left
            curr_score=max_score
            for k in range(hit[1]-1, -1,-1):
                # blossom
                if hit[2]+k-hit[1] < 0:
                    break
                pair = (query_seq[k],gene[hit[2]+k-hit[1]] )
                curr_score += score_match(pair)
                if curr_score>= max_score:
                    max_score = max_score +  score_match(pair)
                elif max_score - curr_score > threshold:
                    break
            max_hits.append(max_score)
    return max(max_hits)



def run_BLAST(input_file, db_file, threshold):

    # parsing the data in the sequence database
    genes , header_dict = parse_db(db_file)
    # parsing the query sequence in the query file
    query_seq = parse_query(input_file)

    query_words_seq=[]

    # making the kmer (3mers) for the query sequence
    for i in range(0, len(query_seq)):
        if i + 3 <= len(query_seq):
            query_words_seq.append(query_seq[i:i + 3])



    max_scores_of_each_db_seq = []
    for db_seq in genes:
        target_words_seq = []

        # making the kmer (3mers) for each of the sequences in the database
        for i in range(0, len(db_seq)):
            if i + 3 <= len(db_seq):
                target_words_seq.append(db_seq[i:i + 3])

        # THis function executes a tuple join
        score_matrix = hsp(query_words_seq, target_words_seq)

        arr = np.array(score_matrix, dtype=int)

        # finding the diagonal with largest count of matching kmers
        cnts = {}
        for n in arr:
            if n[0] in cnts:
                cnts[n[0]] = cnts.get(n[0], 0) + 1
            else:
                cnts[n[0]] = 1

        max_diag = max(cnts.items(), key=lambda t: t[1])[0]

        # expanding to the right and left to the given threshold.
        # the current score is allowed to fall to a max_score-threshold
        # returns the maximum score for the particular sequence in the database
        max_hits = expand_hit(max_diag, arr,query_seq, db_seq,threshold)


        max_scores_of_each_db_seq.append(max_hits)

    # presenting the data
    keys_list = list(header_dict)
    for i in range(0,len(keys_list)):
        header_dict[keys_list[i]] = max_scores_of_each_db_seq[i]

    return header_dict



print("Scores for each sequence in DB:")

print(run_BLAST(QUERY_FILE,DB_FILE, THRESHOLD))
