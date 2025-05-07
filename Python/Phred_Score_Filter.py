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
    
    #FAZENDO ALINHAMENTO SEQ_R1

    index_temp = R1_tamanho_id.index(id_geral) #Pega o index do ID da sequência do assemble na lista de IDs do R1
    alinhamento_R1_teste = pairwise2.align.localxx(R1_seqs[index_temp], seq_ref[0:len(R1_seqs[index_temp])])
    if len(alinhamento_R1_teste) == 0:
        return 0
    alinhamento_R1_assemble_local = alinhamento_R1_teste[0]
    R1_len = len(R1_seqs[index_temp])
    
    #FAZENDO ALINHAMENTO SEQ_R2
    
    index_temp = R2_tamanho_id.index(id_geral) #Pega o index do ID da sequência do assemble na lista de IDs do R1
    R2_ini = len(seq_ref)-len(R2_seqs[index_temp]) #Índice de início provável do R2 no assemble
    seq_aln_R2 = seq_ref[R2_ini:seq_len]
    alinhamento_R2_teste = pairwise2.align.localxx(R2_seqs[index_temp].reverse_complement(), seq_aln_R2)    
    if len(alinhamento_R2_teste) == 0:
        return 0
    alinhamento_R2_assemble_local = alinhamento_R2_teste[0]
    R2_len = len(R2_seqs[index_temp])

    interreads = range((seq_len-R2_len),R1_len)
    
    '''
    Condições:
        Regiao sem gap é MAIOR ou igual ao início da região de interssecção das reads:
            Utiliza-se a região anterior à intersecção como inicial
            Se tem um gap no começo joga fora.
        
        Regiao sem gap é MENOR que o início da região de interssecção das reads:
            Utiliza-se a região sem gap como inicial
        
    '''
    
    R1_no_gap = alinhamento_R1_assemble_local[0].split("-")[0] #região do R1 que alinha com assemble sem gaps
    if (len(R1_no_gap) -1) not in interreads:
        return None
    R2_no_gap = alinhamento_R2_assemble_local[0].split("-")[-1] #região do R2 que alinha com assemble sem gaps
    if ((seq_len - len(R2_no_gap))-1) not in interreads:
        return None

    #Vamos pegar as regiões que não se sobrepõem com o assemble e pegar as com qualidade boa.
    regiao_inicial = assemble_ids[i].letter_annotations["phred_quality"][0:(seq_len - R2_len)] #A região inicial tem que terminar onde o R2 começa
    regiao_final = assemble_ids[i].letter_annotations["phred_quality"][R1_len:seq_len] #A região final tem que começar onde o R1 termina
    regiao_meio = assemble_ids[i].letter_annotations["phred_quality"][(seq_len - R2_len):R1_len] #regiao sem gap do R2 é a final
    
    try:
        if (min(regiao_inicial) > 20) and (min(regiao_final) > 20) and (mean(regiao_meio) > 30):    
            write(assemble_ids[i], seq_selecionada, "fastq")
    except:
        return None
    
    return None


from Bio.SeqIO import parse

assemble_ids = list(parse("844-NI-CPesado_assemble-pass.fastq", "fastq")) #Colocar os arquivos de interesse

R1_ids = list(parse("844-NI-Cpesada_L001_R1_Output.fastq", "fastq")) #Colocar os arquivos de interesse
R2_ids = list(parse("844-NI-Cpesada_L001_R2_Output.fastq", "fastq")) #Colocar os arquivos de interesse

R1_tamanho_id = [R1_ids[x].id for x in  range (len(R1_ids))] #Pra cada elemento do R1_ids, escreve o id com o número da seq, ele escreve só o ID.
R1_seqs = [R1_ids[x].seq for x in range (len(R1_ids))] #Pra cada elemento do R1_ids, escreve o id com o número da seq, ele escreve só o ID.

R2_tamanho_id = [R2_ids[x].id for x in  range (len(R2_ids))] #Pra cada elemento do R2_ids, escreve o id com o número da seq, ele escreve só o ID.
R2_seqs = [R2_ids[x].seq for x in range (len(R2_ids))] #Pra cada elemento do R2_ids, escreve o id com o número da seq, ele escreve só o ID.

seq_selecionada = open("seq_selecionada_844_NI.fastq", "a+") 

for i in range (len(assemble_ids)):
    if (i % 1000 == 0) and (i != 0):
        print(i)
        seq_selecionada.close()
        seq_selecionada = open("seq_selecionada_844_NI.fastq", "a+")
    ver_qualidade(R1_tamanho_id, R2_tamanho_id ,R1_seqs, R2_seqs, assemble_ids, i, seq_selecionada)
seq_selecionada.close()
