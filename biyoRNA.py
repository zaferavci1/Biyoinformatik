from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqUtils import gc_fraction
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
from Bio import AlignIO
from Bio.Align.Applications import ClustalOmegaCommandline

# Dosyaya yazma işlemi için bir dosya adı belirtin
output_file = "/Users/zaferavci/Documents/GitHub/Biyoinformatik/outputRNA.txt"

def read_fasta(fasta_file):
    sequences = []
    with open(fasta_file, "r") as f:
        for record in SeqIO.parse(f, "fasta"):
            sequences.append(str(record.seq))
    return sequences

# Dizileri okuma
sequences = read_fasta("/Users/zaferavci/Documents/GitHub/Biyoinformatik/RNA.fasta")

def find_motifs_in_protein(fasta_file, motif_length=5, specific_motif=None, output_file=None):
    with open(output_file, 'a') as f:
        f.write("#################### Motif Arama Sonuçları ####################\n")
        for record in SeqIO.parse(fasta_file, "fasta"):
            sequence = str(record.seq)
            for i in range(len(sequence) - motif_length + 1):
                motif = sequence[i:i+motif_length]
                if motif == specific_motif:
                    result = f"Dizi ID: {record.id} - Motif: {motif} - Pozisyon: {i + 1}\n"
                    print(result)
                    f.write(result)
        f.write("\n")

# Motif arama
specific_motifs = "TAFMG"
print("\n\n#################### 2.a  #####################\n\n")
find_motifs_in_protein("/Users/zaferavci/Documents/GitHub/Biyoinformatik/RNA.fasta", 5, specific_motifs, output_file)

def global_alignment(seq1, seq2, match_score, mismatch_score, gap_penalty, output_file=None):
    alignment = pairwise2.align.globalms(seq1, seq2, match_score, mismatch_score, gap_penalty, gap_penalty)[0]
    result = f"Global hizalama skoru: {alignment.score}\n"
    result += f"{alignment.seqA}\n{alignment.seqB}\n"
    if output_file:
        with open(output_file, 'a') as f:
            f.write(result)
    return alignment.score, alignment.seqA, alignment.seqB

# Varsayılan skorlar ile hizalama (d şıkkı için b şıkkı)
print("\n\nVarsayılan skorlar d şıkkı için b şıkkının denemeleri. Varsayılan skorlar ile.\n\n")
match_score = 2
mismatch_score = -1
gap_penalty = -2

with open(output_file, 'a') as f:
    f.write("#################### Global Hizalama (Varsayılan Skorlar) ####################\n")
    for i in range(len(sequences)):
        for j in range(i + 1, len(sequences)):
            score, aligned_seq1, aligned_seq2 = global_alignment(sequences[i], sequences[j], match_score, mismatch_score, gap_penalty, output_file)
            f.write(f"Diziler {i+1} ve {j+1} arasındaki global hizalama skoru: {score}\n")
            f.write(f"{aligned_seq1}\n")
            f.write(f"{aligned_seq2}\n")
        f.write("\n")

# Farklı skorlar ile hizalama
print("\n\nVarsayılan skorlar d şıkkı için b şıkkının denemeleri. Farklı skorlar ile.\n\n")
match_score = 3
mismatch_score = -2
gap_penalty = -3

with open(output_file, 'a') as f:
    f.write("#################### Global Hizalama (Farklı Skorlar) ####################\n")
    for i in range(len(sequences)):
        for j in range(i + 1, len(sequences)):
            score, aligned_seq1, aligned_seq2 = global_alignment(sequences[i], sequences[j], match_score, mismatch_score, gap_penalty, output_file)
            f.write(f"Diziler {i+1} ve {j+1} arasındaki global hizalama skoru (farklı skorlar): {score}\n")
            f.write(f"{aligned_seq1}\n")
            f.write(f"{aligned_seq2}\n")
        f.write("\n")

def local_alignment(seq1, seq2, match_score, mismatch_score, gap_penalty, output_file=None):
    alignment = pairwise2.align.localms(seq1, seq2, match_score, mismatch_score, gap_penalty, gap_penalty)[0]
    result = f"Local hizalama skoru: {alignment.score}\n"
    result += f"{alignment.seqA}\n{alignment.seqB}\n"
    if output_file:
        with open(output_file, 'a') as f:
            f.write(result)
    return alignment.score, alignment.seqA, alignment.seqB

# Varsayılan skorlar ile local hizalama
print("\n\nVarsayılan skorlar d şıkkı için b şıkkının denemeleri. Varsayılan skorlar ile.\n\n")
match_score = 2
mismatch_score = -1
gap_penalty = -2

with open(output_file, 'a') as f:
    f.write("#################### Local Hizalama (Varsayılan Skorlar) ####################\n")
    for i in range(len(sequences)):
        for j in range(i + 1, len(sequences)):
            score, aligned_seq1, aligned_seq2 = local_alignment(sequences[i], sequences[j], match_score, mismatch_score, gap_penalty, output_file)
            f.write(f"Diziler {i+1} ve {j+1} arasındaki local hizalama skoru: {score}\n")
            f.write(f"{aligned_seq1}\n")
            f.write(f"{aligned_seq2}\n")
        f.write("\n")

# Farklı skorlar ile local hizalama
print("\n\nVarsayılan skorlar d şıkkı için b şıkkının denemeleri. Farklı skorlar ile.\n\n")
match_score = 3
mismatch_score = -2
gap_penalty = -3

with open(output_file, 'a') as f:
    f.write("#################### Local Hizalama (Farklı Skorlar) ####################\n")
    for i in range(len(sequences)):
        for j in range(i + 1, len(sequences)):
            score, aligned_seq1, aligned_seq2 = local_alignment(sequences[i], sequences[j], match_score, mismatch_score, gap_penalty, output_file)
            f.write(f"Diziler {i+1} ve {j+1} arasındaki local hizalama skoru (farklı skorlar): {score}\n")
            f.write(f"{aligned_seq1}\n")
            f.write(f"{aligned_seq2}\n")
        f.write("\n")

print("Çıktılar 'output.txt' dosyasına kaydedildi.")
