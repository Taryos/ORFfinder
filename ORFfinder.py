#!/usr/bin/env python3
import argparse
import argcomplete

parser = argparse.ArgumentParser(description="ORF detection")
parser.add_argument("-i", "--input_file", required=True, help="Path to the input fasta file")
parser.add_argument("-o", "--output_file", required=True, help="Path to the output fasta file")
parser.add_argument("-t", "--translate", action="store_true", help="Write translated ORF sequences")
parser.add_argument("-m", "--min_length", type=int, default=100, help="Minimum legnth of ORF to consider")
parser.add_argument("-r", "--reverse_frames", action="store_true", help="find reverse strand ORFs")
parser.add_argument("-s", "--strict_start_codon", action="store_true", help="only ATG codon is used as start codon")

argcomplete.autocomplete(parser)
args = parser.parse_args()

print(f"input_path: {args.input_file}")
print(f"output_path: {args.output_file}")
print(f"translate: {args.translate}")
print(f"min_length: {args.min_length}")

def lit_fasta(input_fasta):
    sequences = []
    nom = ""
    sequence = "" 

    with open(input_fasta, 'r') as f:
        for ligne in f:
            if ligne[-1] == '\n':
                ligne = ligne[:-1]

            if ligne[0] == '>': 
                if nom:
                    sequences.append({"nom": nom, "sequence": sequence})

                nom = ligne
                sequence = ""
            else:
                sequence += ligne
                
        if nom:     
            sequences.append({"nom": nom, "sequence": sequence})

    return sequences

def find_ORFs(sequence, cadre_de_lecture, min_len, start_codon_list):
    stop_codons_liste = ["TAA", "TAG", "TGA"]
    ORFs = []

    ORF = ""
    start_index = -1
    id = 0

    for i in range(cadre_de_lecture - 1, len(sequence) - 2, 3):
        codon = sequence[i:i+3]

        if ORF:  # ORF en cours
            ORF += codon

            if codon in stop_codons_liste:
                if len(ORF) // 3 >= min_len:
                    id += 1
                    ORF_info = {
                        "id" : id,
                        "sequence" : ORF,
                        "start" : start_index + 1,
                        "end" : i + 3
                    }
                    ORFs.append(ORF_info)

                ORF = ""
                start_index = -1

        elif codon in start_codon_list:  # début d'ORF
            ORF = codon
            start_index = i

    return ORFs
          
def reverse_complement(sequence):
    reverse = ""
    for n in sequence:
        match n:
            case 'A':
                reverse += 'T'
            case 'T':
                reverse += 'A'
            case 'C':
                reverse += 'G'
            case 'G':
                reverse += 'C'
            case _:
                reverse += n
            
    return reverse[::-1]

code_genetique = {
    'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
    'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
    'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
    'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
    'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
    'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
    'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
    'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
    'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
    'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
    'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
    'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
    'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
    'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
    'TAC':'Y', 'TAT':'Y', 'TAA':'*', 'TAG':'*',
    'TGC':'C', 'TGT':'C', 'TGA':'*', 'TGG':'W',
}

def traduire(sequence):
    prot = ""

    for i in range(0, len(sequence) - 2, 3):
        codon = sequence[i:i+3]

        if codon in code_genetique:
            aa = code_genetique[codon]
            
        else:
            aa = 'X'

        prot += aa
    
    return prot

start_codons_liste = ["ATG", "GTG", "CTG", "TTG"]
strict_start_codons_liste = ["ATG"]


def fonction_principale(fasta_path=args.input_file, output_path=args.output_file, 
                        min_len=args.min_length, strict_start_codon=args.strict_start_codon, 
                        translate=args.translate, reverse=args.reverse_frames):
    
    sequences = lit_fasta(fasta_path)

    if strict_start_codon:
        start_liste = strict_start_codons_liste
    else:
        start_liste = start_codons_liste


    with open(output_path, "w") as f:
        for sequence in sequences:

            brin = "direct"
            seq_a_analyser = sequence["sequence"]

            if reverse:
                seq_a_analyser = reverse_complement(seq_a_analyser)
                brin = "reverse"
            

            for i in range(3):

                ORFs = find_ORFs(seq_a_analyser, i+1, min_len, start_codon_list = start_liste)

                for ORF in ORFs:
                    f.write(">ORF number " + str(ORF["id"]) + " in reading frame " + str((i+1)) + " on the " + brin + " strand extends from base "
                            + str(ORF["start"]) + " to " + str(ORF["end"]) + ".\n" + ORF["sequence"] + "\n\n")
                    if translate:
                        f.write(">Translation of ORF number " + str(ORF["id"]) + " in reading frame " + str((i+1)) + " on the " + brin + " strand." + "\n"
                                + traduire(ORF["sequence"]) + "\n\n")


if __name__ == "__main__":
    fonction_principale()







