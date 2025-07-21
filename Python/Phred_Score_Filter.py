# This code selects sequences from a fastq file based on the quality of the reads. 
# It uses the Phred score to filter the sequences and select only those with a high quality score.
# The code uses the Biopython library to read and write fastq files, and to perform pairwise alignment of sequences.
# The input files are the assembled sequences and the R1 and R2 reads, and the output is a fastq file with the selected sequences.

from Bio.SeqIO import write
from Bio import pairwise2
from statistics import mean

def ver_qualidade (R1_tamanho_id, R2_tamanho_id ,R1_seqs, R2_seqs, assemble_ids, i, seq_selecionada):
    id_geral = assemble_ids[i].id
    seq_ref = assemble_ids[i].seq
    seq_len = len(assemble_ids[i].letter_annotations["phred_quality"])
    
    #ALIGNING SEQ_R1

    index_temp = R1_tamanho_id.index(id_geral) #Get the index of the assembly's sequence ID from R1's list of IDs
    alinhamento_R1_teste = pairwise2.align.localxx(R1_seqs[index_temp], seq_ref[0:len(R1_seqs[index_temp])])
    if len(alinhamento_R1_teste) == 0:
        return 0
    alinhamento_R1_assemble_local = alinhamento_R1_teste[0]
    R1_len = len(R1_seqs[index_temp])
    
    #ALIGNING SEQ_R2
    
    index_temp = R2_tamanho_id.index(id_geral) #Get the index of the assembly's sequence ID from R1's list of IDs
    R2_ini = len(seq_ref)-len(R2_seqs[index_temp]) #Index of probable start of R2 in the assembly
    seq_aln_R2 = seq_ref[R2_ini:seq_len]
    alinhamento_R2_teste = pairwise2.align.localxx(R2_seqs[index_temp].reverse_complement(), seq_aln_R2)    
    if len(alinhamento_R2_teste) == 0:
        return 0
    alinhamento_R2_assemble_local = alinhamento_R2_teste[0]
    R2_len = len(R2_seqs[index_temp])

    interreads = range((seq_len-R2_len),R1_len)
    
    '''
    Conditions:
        The region without a gap is GREATER than or equal to the start of the region where the reads intersect:
            The region before the intersection is used as the starting point. If there is a gap at the beginning, discard it.
        
        Region without a gap is LESS than the start of the region where the reads intersect:
            The region without a gap is used as the starting point
        
    '''
    
    R1_no_gap = alinhamento_R1_assemble_local[0].split("-")[0] #region of R1 that aligns with assembly without gaps
    if (len(R1_no_gap) -1) not in interreads:
        return None
    R2_no_gap = alinhamento_R2_assemble_local[0].split("-")[-1] #region of R2 that aligns with assembly without gaps
    if ((seq_len - len(R2_no_gap))-1) not in interreads:
        return None

    #Let's take the regions that don't overlap with the assembly and take the ones with good quality.
    regiao_inicial = assemble_ids[i].letter_annotations["phred_quality"][0:(seq_len - R2_len)] #The start region has to end where R2 begins
    regiao_final = assemble_ids[i].letter_annotations["phred_quality"][R1_len:seq_len] #The final region has to start where R1 ends
    regiao_meio = assemble_ids[i].letter_annotations["phred_quality"][(seq_len - R2_len):R1_len] #R2 gap-free region is the final
    
    try:
        if (min(regiao_inicial) > 20) and (min(regiao_final) > 20) and (mean(regiao_meio) > 30):    
            write(assemble_ids[i], seq_selecionada, "fastq")
    except:
        return None
    
    return None


from Bio.SeqIO import parse

assemble_ids = list(parse("844-NI-CPesado_assemble-pass.fastq", "fastq")) #Place the files of interest

R1_ids = list(parse("844-NI-Cpesada_L001_R1_Output.fastq", "fastq")) #Place the files of interest
R2_ids = list(parse("844-NI-Cpesada_L001_R2_Output.fastq", "fastq")) #Place the files of interest

R1_tamanho_id = [R1_ids[x].id for x in  range (len(R1_ids))] #For each element of R1_ids, write the id with the seq number, it only writes the ID.
R1_seqs = [R1_ids[x].seq for x in range (len(R1_ids))] #For each element of R1_ids, write the id with the seq number, it only writes the ID.

R2_tamanho_id = [R2_ids[x].id for x in  range (len(R2_ids))] #For each element of R2_ids, write the id with the seq number, it only writes the ID.
R2_seqs = [R2_ids[x].seq for x in range (len(R2_ids))] #For each element of R2_ids, write the id with the seq number, it only writes the ID.

seq_selecionada = open("seq_selecionada_844_NI.fastq", "a+") 

for i in range (len(assemble_ids)):
    if (i % 1000 == 0) and (i != 0):
        print(i)
        seq_selecionada.close()
        seq_selecionada = open("seq_selecionada_844_NI.fastq", "a+")
    ver_qualidade(R1_tamanho_id, R2_tamanho_id ,R1_seqs, R2_seqs, assemble_ids, i, seq_selecionada)
seq_selecionada.close()
