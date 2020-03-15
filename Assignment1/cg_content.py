DNA_sequence = "ACAGTCGACTAGCTTGCACGTAC"


def base_content(sequence, base):
    noOfGivenBase = sequence.count(base)
    return noOfGivenBase


def percentage_of_base(sequence,bases):
    percentage=0
    for i in range(len(bases)):
        percentage=percentage+base_content(sequence,bases[i])

    total_percentage=(percentage/len(sequence))
    return total_percentage

###########Testing############
Cper=base_content(DNA_sequence,'C')
Gper=base_content(DNA_sequence,'G')
#print(Cper, Gper)
CGper=percentage_of_base(DNA_sequence,"CG")
#print(CGper)
