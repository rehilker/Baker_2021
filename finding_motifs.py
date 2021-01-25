""" This code takes a text file of upper or lower case DNA sequences in
fasta format and returns the starting positions (starting at 0) of all consensus
sequences with a specified number of mismatches. Sequences must ONLY contain
A/G/C/T characters """

#Put your consensus motif here; all NT codes listed in NT_code are accepted.
#Letters can be upper or lower case
consensus = "TANA"
#Put the tolerated number of mismatches here, should equal the number of "N"s
#in the consensus sequence
mis = 1

#The Single Letter NT code for reference
NT_code = {"W": ["A", "T"], "R": ["G", "A"], "S": ["G", "C"], 
"N": ["A", "G", "C", "T"], "Y": ["C", "T"], "K": ["G", "T"], 
"M": ["A", "C"], "B": ["G", "T", "C"], "D": ["G", "A", "T"], 
"H": ["A", "C", "T"], "V": ["G", "C", "A"]}

def fasta2dict(fil):
    """
    Read fasta-format file fil, return dict of form scaffold:sequence.
    Note: Uses only the unique identifier of each sequence, rather than the 
    entire header, for dict keys. 
    """
    dic = {}
    cur_scaf = ''
    cur_seq = []
    for line in open(fil):
        if line.startswith(">") and cur_scaf == '':
            cur_scaf = line.split(' ')[0]
        elif line.startswith(">") and cur_scaf != '':
            dic[cur_scaf] = ''.join(cur_seq)
            cur_scaf = line.split(' ')[0]
            cur_seq = []
        else:
            cur_seq.append(line.rstrip())
    dic[cur_scaf] = ''.join(cur_seq)
    return dic

#Replace the underline with the name of your text file in fasta format
new = fasta2dict("_.txt")

def HammingDistance(seq,con):
    count = 0
    seq = seq.upper()
    con = con.upper()
    ls = len(seq)
    for i in range(0,ls):
        if con[i] == seq[i]:
            count += 0
        elif con[i] != "A" or con[i] != "G" or con[i] != "C" or con[i] != "T" or con[i] != "N":
            if seq[i] == "A":
                if con[i] == "W" or con[i] == "R" or con[i] == "M" or con[i] == "V" or con[i] == "H" or con[i] == "D":
                    count += 0
                else:
                    count += 1
            elif seq[i] == "T":
                if con[i] == "W" or con[i] == "Y" or con[i] == "K" or con[i] == "B" or con[i] == "H" or con[i] == "D":
                    count += 0
                else:
                    count += 1
            elif seq[i] == "C":
                if con[i] == "S" or con[i] == "Y" or con[i] == "M" or con[i] == "V" or con[i] == "H" or con[i] == "B":
                    count += 0
                else:
                    count += 1
            elif seq[i] == "G":
                if con[i] == "S" or con[i] == "R" or con[i] == "K" or con[i] == "V" or con[i] == "B" or con[i] == "D":
                    count += 0
                else:
                    count += 1
        elif con[i] != seq[i]:
            count += 1
    return count

def HammingCount(Text, Pattern):
    HD = {}
    count = 0
    n = len(Text)-len(Pattern)+1
    for i in range(0,n):
        HD[count] = 0
        count += 1
    for i in range(n):
        x = HammingDistance(Text[i:i+len(Pattern)],Pattern)
        HD[len(Text[0:i])] = x
    return HD
    
def ApproximatePatternMatching_1(Text, Pattern, d):
    positions = []
    dic = HammingCount(Text, Pattern)
    count = 0
    for key in dic:
        if dic[count] <= d:
            positions.append(count)
            count += 1
        else:
            count += 1
    return positions
    

for key in new:
    seq = new[key]
    x = ApproximatePatternMatching_1(seq, consensus, mis)
    print(x, key)

