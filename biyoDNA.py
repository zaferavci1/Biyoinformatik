from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqUtils import gc_fraction
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
import numpy as np
from Bio import AlignIO
from Bio.Align.Applications import ClustalOmegaCommandline

# Dosya yolu
output_file = "/Users/zaferavci/Documents/GitHub/Biyoinformatik/outputDNA.txt"

# Fasta dosyasını oku
records = list(SeqIO.parse("/Users/zaferavci/Documents/GitHub/Biyoinformatik/DNA.fasta", "fasta"))

# ATGC sayıları
topA = 0
topT = 0
topG = 0
topC = 0

# Çıktıları dosyaya yaz
with open(output_file, 'w') as f:

    # a sorusu
    f.write("########################## a sorusu\n")
    for record in records:
        sequence = str(record.seq)
        f.write(f"Dizi ID: {record.id}\n")
        f.write(f"A sayısı: {sequence.count('A')}\n")
        f.write(f"T sayısı: {sequence.count('T')}\n")
        f.write(f"G sayısı: {sequence.count('G')}\n")
        f.write(f"C sayısı: {sequence.count('C')}\n")
        
        # ATGC toplamını hesapla
        topA += sequence.count('A')
        topT += sequence.count('T')
        topG += sequence.count('G')
        topC += sequence.count('C')

    f.write("\nToplam A sayısı: {}\n".format(topA))
    f.write("Toplam T sayısı: {}\n".format(topT))
    f.write("Toplam G sayısı: {}\n".format(topG))
    f.write("Toplam C sayısı: {}\n".format(topC))

    # b sorusu
    f.write("\n########################## b sorusu\n")
    rna_sequences = []
    for record in records:
        rna_sequence = record.seq.transcribe()
        rna_sequences.append(rna_sequence)
        f.write(f"RNA Sequence for {record.id}: {rna_sequence}\n")  # RNA dizisini yazdır

    # c sorusu
    f.write("\n########################## c sorusu\n")
    protein_sequences = []
    for record in records:
        protein_sequence = record.seq.translate()
        protein_sequences.append(protein_sequence)
        f.write(f"Protein Sequence for {record.id}: {protein_sequence}\n")  # Protein dizisini yazdır

    # d sorusu
    f.write("\n########################## d sorusu\n")
    with open("/Users/zaferavci/Documents/GitHub/Biyoinformatik/protein_sequences.fasta", "w") as protein_file:
        for i, protein_sequence in enumerate(protein_sequences):
            protein_file.write(f">Protein_{i+1}\n{protein_sequence}\n")
            f.write(f"Protein_{i+1}: {protein_sequence}\n")  # Protein dizisini output dosyasına yazdır



    # e sorusu
    f.write("\n########################## e sorusu\n")
    for record in records:
        gc_content = gc_fraction(record.seq) * 100  # gc_fraction 0 ile 1 arasında değer döndürür, %100 ile çarp
        f.write(f"Dizi ID: {record.id}, GC içeriği: {gc_content:.2f}%\n")

    # f sorusu
    f.write("\n########################## f sorusu\n")
    clustalomega_cline = ClustalOmegaCommandline(infile="/Users/zaferavci/Documents/GitHub/Biyoinformatik/DNA.fasta", 
                                                 outfile="/Users/zaferavci/Documents/GitHub/Biyoinformatik/global_alignment.fasta", 
                                                 verbose=True, auto=True)
    stdout, stderr = clustalomega_cline()
    
    # Hizalama sonucunu oku
    alignment = AlignIO.read("/Users/zaferavci/Documents/GitHub/Biyoinformatik/global_alignment.fasta", "fasta")
    
    # Hizalamayı yazdır
    f.write(str(alignment))

    # g sorusu
    f.write("\n########################## g sorusu\n")
    seq1 = records[3].seq
    seq2 = records[4].seq

    alignments = pairwise2.align.localxx(seq1, seq2)

    f.write(format_alignment(*alignments[0]))

    f.write("\n########################## h sorusu\n")

    # En uzun dizinin uzunluğunu bulun
    max_seq_length = max(len(record.seq) for record in records)

    # PFM'yi oluştur
    pfm = np.zeros((4, max_seq_length), dtype=int)

    for record in records:
        for i, nucleotide in enumerate(record.seq):
            if nucleotide == "A":
                pfm[0, i] += 1
            elif nucleotide == "C":
                pfm[1, i] += 1
            elif nucleotide == "G":
                pfm[2, i] += 1
            elif nucleotide == "T":
                pfm[3, i] += 1

    # Matrisi dosyaya yaz
    for row in pfm:
        f.write("\t".join(map(str, row)) + "\n")


    # i sorusu
    f.write("\n########################## i sorusu\n")

    ppm = pfm / len(records)

     # Matrisi dosyaya yaz
    for row in ppm:
        f.write("\t".join(map(str, row)) + "\n")

    # j sorusu
    f.write("\n########################## j sorusu\n")

    background_probs = {"A": 0.25, "C": 0.25, "G": 0.25, "T": 0.25}

    background_array = np.tile(np.array([background_probs[nucleotide] for nucleotide in "ACGT"]), (ppm.shape[1], 1)).T

    pwm = np.log2(ppm / background_array)

     # Matrisi dosyaya yaz
    for row in pwm:
        f.write("\t".join(map(str, row)) + "\n")
