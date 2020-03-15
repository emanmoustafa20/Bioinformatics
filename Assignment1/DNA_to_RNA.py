import cg_content as content

def DNA_to_RNA(DNASequence):
    RNA_Sequence=DNASequence.replace('T','U')
    return RNA_Sequence



######Testing######
RNA=DNA_to_RNA(content.DNA_sequence)

