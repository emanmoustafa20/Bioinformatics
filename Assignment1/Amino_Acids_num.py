import DNA_to_RNA

def Amino_Acids_No(RNASequence):
     noOfAminoAcids=int(len(RNASequence)/3)
     Remainder=(len(RNASequence)%3)
     return  noOfAminoAcids,Remainder


#######Testing#######
no,remainder=Amino_Acids_No(DNA_to_RNA.RNA)
print(no,remainder)
