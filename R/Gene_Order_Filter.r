#Código para selecionar sequências com base na ordem do lócus no IMGT.
#O racional é que se o gene tiver orientação direta não escolhe nada que vier depois dele.
#Se tiver orientação indireta não escolhe nada que vier antes, só o que vier depois.

# Primeiro pega o gene funcional, ver se ele está presente no dicionário direto ou oposto.
# E pega o índice de i.
# Se for direto exclui tudo que tiver depois do i
# Se for oposto exclui tudo o que tiver antes de i
# pega o pseudogene e vê se ele está presente no dicionário, se estiver guarda o valor dele.

#Esse código está usando especificamente a tabela de cavalo do imgt.

direto <- c(
  "IGHV4S1"	= 0,
  "IGHV1-41N"	=	1,
  "IGHV3-40N"	=	2,
  "IGHV3-39N"	=	3,
  "IGHV4-38N"	=	4,
  "IGHV4-37N"	=	5,
  "IGHV1-84/ORF" =	6,
  "IGHV4-83" =	7,
  "IGHV4-82"	=	8,
  "IGHV4-81"	=	9,
  "IGHV3-80"	=	10,
  "IGHV1-79"	=	11,
  "IGHV3-78"	=	12,
  "IGHV2-77/ORF"	=	13,
  "IGHV(II)-43N"	=	46,
  "IGHV4-42N/ORF"	=	47,
  "IGHV(II)-43"	=	54,
  "IGHV4-42/ORF"	=	55,
  "IGHV1-41D"	=	56,
  "IGHV3-40D"	=	57,
  "IGHV3-39D"	=	58,
  "IGHV4-38D"	=	59,
  "IGHV4-37D"	=	60,
  "IGHV1-41"	=	61,
  "IGHV3-40"	=	62,
  "IGHV3-39"	=	63,
  "IGHV4-38"	=	64,
  "IGHV4-37"	=	65,
  "IGHV3-36"	=	66,
  "IGHV4-35"	=	67,
  "IGHV7-34"	=	68,
  "IGHV4-33"	=	69,
  "IGHV3-32"	=	70,
  "IGHV1-31"	=	71,
  "IGHV4-30"	=	72,
  "IGHV4-29"	=	73,
  "IGHV4-28"	=	74,
  "IGHV4-20N"	=	75,
  "IGHV3-27"	=	76,
  "IGHV3-26"	=	77,
  "IGHV4-25"	=	78,
  "IGHV3-24"	=	79,
  "IGHV4-23"	=	80,
  "IGHV4-22"	=	81,
  "IGHV4-20"	=	84,
  "IGHV3-19"	=	85,
  "IGHV3-18"	=	86,
  "IGHV4-17"	=	87,
  "IGHV3-16"	=	88,
  "IGHV(II)-15"	=	89,
  "IGHV4-14"	=	90,
  "IGHV3-13"	=	91,
  "IGHV(II)-12"	=	92,
  "IGHV4-11"	=	93,
  "IGHV4-10"	=	94,
  "IGHV3-9"	=	95,
  "IGHV4-8"	=	96,
  "IGHV3-7"	=	97,
  "IGHV4-6"	=	98,
  "IGHV1-5"	=	99,
  "IGHV4-4"	=	100,
  "IGHV4-3"	=	101,
  "IGHV4-2"	=	102,
  "IGHV3-1"	=	103
)

oposto <- c(
  "IGHV4-33D"	=	14,
  "IGHV7-34D"	=	15,
  "IGHV4-35D"	=	16,
  "IGHV3-76"	=	17,
  "IGHV4-75"	=	18,
  "IGHV4-74"	=	19,
  "IGHV3-73"	=	20,
  "IGHV4-72"	=	21,
  "IGHV5-71"	=	22,
  "IGHV1-70"	=	23,
  "IGHV3-69"	=	24,
  "IGHV4-68"	=	25,
  "IGHV5-67/ORF"	=	26,
  "IGHV9-66"	=	27,
  "IGHV4-65"	=	28,
  "IGHV3-64"	=	29,
  "IGHV2-63"	=	30,
  "IGHV3-62/ORF"	=	31,
  "IGHV1-61"	=	32,
  "IGHV3-60/ORF"	=	33,
  "IGHV4-59"	=	34,
  "IGHV3-58"	=	35,
  "IGHV4-57"	=	36,
  "IGHV3-56"	=	37,
  "IGHV4-55"	=	38,
  "IGHV(II)-54"	=	39,
  "IGHV4-53"	=	40,
  "IGHV4-52"	=	41,
  "IGHV4-51"	=	42,
  "IGHV3-50"	=	43,
  "IGHV4-49"	=	44,
  "IGHV4-48"	=	45,
  "IGHV3-47"	=	48,
  "IGHV3-46"	=	49,
  "IGHV1-45"	=	50,
  "IGHV4-42D/ORF"	=	51,
  "IGHV(II)-43D"	=	52,
  "IGHV3-44"	=	53,
  "IGHV4-20D"	=	82,
  "IGHV4-21"	=	83
)

### Importando banco de dados 
data_834 <- read.csv("resultado_final_teste_834.csv", sep = "\t", stringsAsFactors = FALSE, header = TRUE)

### Manter na coluna apenas o 1o anotação e deletar tudo depois da vírgula (nesse caso: ;) 
data_834$gene <- c(sub(";(.*)", "",data_834$gene))

for (i in data_834$V.GENE.and.allele) {
  i = gsub("\\*01(.*)", "", i) #Removendo tudo que vem depois do *01 
  if (i %in% names(direto)) { #Verificando pertencimento no dicionário direto
    i = direto[[i]] #Se pertence o i vira o índice do dicionário, que é a localização no lócus.
    dic_direto <- direto[direto <= i] #Removendo o que for maior que o valor do índice.
    #print(dic_direto)
  } else if (i %in% names(oposto)) { 
    i = oposto[[i]]
    dic_oposto <- oposto[oposto >= i] #Removendo o que for menor que o valor do índice.
    #print(dic_oposto)
  }
  pseudogenes_selecionados = list()
  for (j in data_834$gene) {
    j = gsub("\\*(.*)", "", j) #Removendo tudo depois do *
    #print(j)
    if (j %in% names(dic_direto) | j %in% names(dic_oposto)) {
      pseudogenes_selecionados <- c(pseudogenes_selecionados, j)
    }
  }
  }
