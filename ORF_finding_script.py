fasta_file = open('/Users/ritikaaggarwal/Bioinfo/dna2.fasta','r')
#checking if the file to open exist
try:
    file = open('/Users/ritikaaggarwal/Bioinfo/dna2.fasta')
except IOError:
    print('oops!, file does not exist')

#make a dictionary of the sequences having identifier as keys and seq as values
seqs = {}
for line in fasta_file:
    line = line.rstrip()
    if line[0]=='>':
        words = line.split()
        name = words[0][1:]
        seqs[name] = ''
    else:
        seqs[name] = seqs[name]+line

#q1 How many records are in the file? 
print(len(seqs.keys()))

#q2.1 What are the lengths of the sequences in the file? 
#q2.2What is the longest sequence and what is the shortest sequence? 
#q2.3Is there more than one longest or shortest sequence? 
#q2.4 What are their identifiers? 
seq_length = {}
for k,v in seqs.items():
    seq_length[k] = len(v)
    
print(seq_length)
print(sorted(seq_length.values()))
maxi = max(seq_length.values())
print('max length identifier is',[key for key in seq_length if seq_length[key] == maxi],'and maximum value is',maxi)
mini = min(seq_length.values())
print('min length identifier is',[key for key in seq_length if seq_length[key] == mini],'and minimumm value is',mini)

#define a function to find the start,end of ORFs in a sequence
def find_orfs1(sequence,frame):
    'find orfs in a given reading frame'
    orfs = [] #create a empty list to append the start,end of the orfs found
    start_codon = ['ATG'] #orf starts only at ATG
    stop_codons = ['TAG','TAA','TGA'] #orf can stop when encountering any one of them in the same frame
    seq_to_be_screened = sequence[frame-1:] #sequence to be screened must start from the base of the frame
    start = 0 #assume the start position is 0
    while True:
        start = seq_to_be_screened.find('ATG',start) #we will first find the start site in the frame
        if start == -1: #if none found, break the loop
            break
        for i in range(start,len(seq_to_be_screened),3): #if found, we will screen the sequence for stop codon
            if seq_to_be_screened[i:i+3] in stop_codons:
                orfs.append((start,i+3,seq_to_be_screened[start:(i+3)]))
                break
        start += 1 #after this seq gets added to the orfs list, we will start from the next base for screening of next ORF
    return orfs

#define a function to find the longest ORF found, with length and seqID
def longest_orf_all(seqs,frame):
    longest_length = 0
    longest_start = 0
    seq_ID = ''
    for k,v in seqs.items():
        orfs = find_orfs1(v,frame)
        for start,end,sequence in orfs:
            orf_length = end-start
            if orf_length > longest_length:
                longest_length = orf_length
                longest_start = start
                seq_ID = k
    return longest_length, longest_start, seq_ID

longest_length1,longest_start1,seq_ID1 = longest_orf_all(seqs,1)
print('longest length ORF foun in reading frame 1 is {}, its start position is {}, and its ID is {}'.format(longest_length1,longest_start1,seq_ID1))
longest_length2,longest_start2,seq_ID2 = longest_orf_all(seqs,2)
print('longest length ORF found in reading frame 2 is {}, its start position is {}, and its ID is {}'.format(longest_length2,longest_start2,seq_ID2))
longest_length3,longest_start3,seq_ID3 = longest_orf_all(seqs,3)
print('longest length ORF found in reading frame 3 is {}, its start position is {}, and its ID is {}'.format(longest_length3,longest_start3,seq_ID3))

def finding_repeats(fasta_seq_dictionary, length):
    dict_repeats = {}
    seqs = fasta_seq_dictionary
    for k,v in seqs.items():
        dict_repeats[k] = {}
        for i in range(0,len(v)-length,1):
            subseq = v[i:i+length]
            if subseq in dict_repeats[k].keys():
                dict_repeats[k][subseq] += 1
            else:
                dict_repeats[k][subseq] = 1

    dict_repeats = {geneID: {subseq: number_of_repeats for subseq, number_of_repeats in dict_repeats[geneID].items() if number_of_repeats > 3} for geneID in dict_repeats}

    return dict_repeats

dictionary_of_repeats_6 = finding_repeats(seqs,6)


'''
        if dict_repeats[k][subseq] == 3:
            dict_repeats[k].pop(subseq)

        if dict_repeats[k][subseq] == 2:
            del dict_repeats[k][subseq]

        if dict_repeats[k][subseq] == 1:
            del dict_repeats[k][subseq] 

'''
sum_of_max_value = 0
most_repetitive_ID = {}
list_of_max_value_sequences = []
for ID,seq_repeats in dictionary_of_repeats_6.items():
    if dictionary_of_repeats_6[ID]:
        list_of_v = []
        for k,v in dictionary_of_repeats_6[ID].items():
            list_of_v.append(v)
    max_value_repeat = max(list_of_v)
    for k,v in dictionary_of_repeats_6[ID].items():
        if v == max_value_repeat:
            list_of_max_value_sequences.append(k)
for sequences_repeats in list_of_max_value_sequences:
    most_repetitive_ID[sequences_repeats] = []
    for ID,seq_repeat_data in dictionary_of_repeats_6.items():
        for subseq,repeat_val in dictionary_of_repeats_6[ID].items():
            if subseq == sequences_repeats:
              most_repetitive_ID[subseq].append(repeat_val)
print(most_repetitive_ID) 

for sequence_repeats,value_list in most_repetitive_ID.items():
    sum_of_max_value += sum(value_list)

print(sum_of_max_value)


'''
for ID,repeat_values in dictionary_of_repeats_6.items():
    list_most_repetitive_ID[(repeats for repeats,values in dictionary_of_repeats_6[ID].items() if values in list_of_max_value_sequences)]
        if v == max_value_repeat:
            list_most_repetitive_ID[k] = []
    for gene_ID in list_most_repetitive_ID.keys():
        for key,value in dictionary_of_repeats_6[ID].items():
            if key == gene_ID:
                list_most_repetitive_ID[gene_ID].append(value)

'''




        


