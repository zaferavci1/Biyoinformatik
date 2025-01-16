from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqUtils import gc_fraction
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
import numpy as np
from Bio import AlignIO
from Bio.Align.Applications import ClustalOmegaCommandline

records = list(SeqIO.parse("dna.fasta", "fasta"))

# ATGC sayıları
topA = 0
topT = 0
topG = 0
topC = 0

# records listesindeki her bir kaydı döngüyle işlemeye başlıyoruz
for record in records:
    sequence = str(record.seq)
    print(f"Dizi ID: {record.id}")
    print(f"A sayısı: {sequence.count('A')}")
    print(f"T sayısı: {sequence.count('T')}")
    print(f"G sayısı: {sequence.count('G')}")
    print(f"C sayısı: {sequence.count('C')}")
    
    # ATGC toplamını hesaplıyoruz
    topA += sequence.count('A')
    topT += sequence.count('T')
    topG += sequence.count('G')
    topC += sequence.count('C')

# Toplam ATGC sayısını yazdırıyoruz
print(f"\nToplam A sayısı: {topA}")
print(f"Toplam T sayısı: {topT}")
print(f"Toplam G sayısı: {topG}")
print(f"Toplam C sayısı: {topC}")

# b
rna_sequences = []
for record in records:
  rna_sequence = record.seq.transcribe()
  rna_sequences.append(rna_sequence)


protein_sequences = []
for record in records:
  protein_sequence = record.seq.translate()
  protein_sequences.append(protein_sequence)

with open("protein_sequences.fasta", "w") as f:
  for i, protein_sequence in enumerate(protein_sequences):
    f.write(f">Protein_{i+1}\n{protein_sequence}\n")


for record in records:
  gc_content = gc_fraction(record.seq) * 100  # gc_fraction 0 ile 1 arasında değer döndürür, %100 ile çarpın
  print(f"Dizi ID: {record.id}, GC içeriği: {gc_content:.2f}%")



clustalomega_cline = ClustalOmegaCommandline(infile="dna.fasta", outfile="global_alignment.fasta", verbose=True, auto=True)
stdout, stderr = clustalomega_cline()


# Hizalama sonucunu oku
alignment = AlignIO.read("global_alignment.fasta", "fasta")

# Hizalamayı yazdır
print(alignment)

# İki dizi seçin (örnek olarak ilk iki dizi)
seq1 = records[0].seq
seq2 = records[1].seq

# Yerel hizalamayı hesapla
alignments = pairwise2.align.localxx(seq1, seq2)

# En iyi hizalamayı yazdır
print(format_alignment(*alignments[0]))



# En uzun dizinin uzunluğunu bulun
max_seq_length = max(len(record.seq) for record in records)

# PFM'yi başlatın
pfm = np.zeros((4, max_seq_length), dtype=int)

# Her bir pozisyondaki nükleotidleri sayın
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

# PFM'yi yazdırın
print(pfm)

# PPM'yi hesaplayın
ppm = pfm / len(records)

# PPM'yi yazdırın
print(ppm)

background_probs = {"A": 0.25, "C": 0.25, "G": 0.25, "T": 0.25}

# Arka plan olasılıkları dizisini ppm ile aynı şekle getirin
background_array = np.tile(np.array([background_probs[nucleotide] for nucleotide in "ACGT"]), (ppm.shape[1], 1)).T

# PWM'yi hesaplayın
pwm = np.log2(ppm / background_array)

# PWM'yi yazdırın
print(pwm)

