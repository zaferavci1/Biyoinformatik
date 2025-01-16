from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqUtils import gc_fraction
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
from Bio import AlignIO
from Bio.Align.Applications import ClustalOmegaCommandline


def read_fasta(fasta_file):
  """
  FASTA dosyasından protein dizilerini okur.

  Args:
    fasta_file: FASTA dosyasının yolu.

  Returns:
    Protein dizilerinin bir listesi.
  """
  sequences = []
  with open(fasta_file, "r") as f:
    for record in SeqIO.parse(f, "fasta"):
      sequences.append(str(record.seq))
  return sequences

sequences = read_fasta("protein.fasta")  # "proteinler.fasta" dosya adını kendi dosyanızla değiştirin

def find_motif(sequence, motif):
  """
  Bir protein dizisinde belirli bir motifi arar.

  Args:
    sequence: Protein dizisi.
    motif: Aranacak motif.

  Returns:
    Motifin bulunduğu konumların bir listesi.
  """
  positions = []
  for i in range(len(sequence) - len(motif) + 1):
    if sequence[i:i+len(motif)] == motif:
      positions.append(i)
  return positions

motif = "AGGTG"  # Motif örneği
for sequence in sequences:
  positions = find_motif(sequence, motif)
  if positions:
    print(f"Motif '{motif}' bulundu: {positions}")




def global_alignment(seq1, seq2, match_score, mismatch_score, gap_penalty):
  """
  İki protein dizisi arasında global hizalama yapar.

  Args:
    seq1: İlk protein dizisi.
    seq2: İkinci protein dizisi.
    match_score: Eşleşme skoru.
    mismatch_score: Uyuşmazlık skoru.
    gap_penalty: Boşluk cezası.

  Returns:
    Hizalama skorunu ve hizalanmış dizileri içeren bir tuple.
  """
  alignment = pairwise2.align.globalms(seq1, seq2, match_score, mismatch_score, gap_penalty, gap_penalty)[0]
  return alignment.score, alignment.seqA, alignment.seqB

# Varsayılan skorlar
match_score = 2
mismatch_score = -1
gap_penalty = -2

for i in range(len(sequences)):
  for j in range(i + 1, len(sequences)):
    score, aligned_seq1, aligned_seq2 = global_alignment(sequences[i], sequences[j], match_score, mismatch_score, gap_penalty)
    print(f"Diziler {i+1} ve {j+1} arasındaki global hizalama skoru: {score}")
    print(aligned_seq1)
    print(aligned_seq2)

# Farklı skorlar ile hizalama (örnek)
match_score = 3
mismatch_score = -2
gap_penalty = -3

for i in range(len(sequences)):
  for j in range(i + 1, len(sequences)):
    score, aligned_seq1, aligned_seq2 = global_alignment(sequences[i], sequences[j], match_score, mismatch_score, gap_penalty)
    print(f"Diziler {i+1} ve {j+1} arasındaki global hizalama skoru (farklı skorlar): {score}")
    print(aligned_seq1)
    print(aligned_seq2)



def local_alignment(seq1, seq2, match_score, mismatch_score, gap_penalty):
  """
  İki protein dizisi arasında local hizalama yapar.

  Args:
    seq1: İlk protein dizisi.
    seq2: İkinci protein dizisi.
    match_score: Eşleşme skoru.
    mismatch_score: Uyuşmazlık skoru.
    gap_penalty: Boşluk cezası.

  Returns:
    Hizalama skorunu ve hizalanmış dizileri içeren bir tuple.
  """
  alignment = pairwise2.align.localms(seq1, seq2, match_score, mismatch_score, gap_penalty, gap_penalty)[0]
  return alignment.score, alignment.seqA, alignment.seqB

# Varsayılan skorlar
match_score = 2
mismatch_score = -1
gap_penalty = -2

for i in range(len(sequences)):
  for j in range(i + 1, len(sequences)):
    score, aligned_seq1, aligned_seq2 = local_alignment(sequences[i], sequences[j], match_score, mismatch_score, gap_penalty)
    print(f"Diziler {i+1} ve {j+1} arasındaki local hizalama skoru: {score}")
    print(aligned_seq1)
    print(aligned_seq2)

# Farklı skorlar ile hizalama (örnek)
match_score = 3
mismatch_score = -2
gap_penalty = -3

for i in range(len(sequences)):
  for j in range(i + 1, len(sequences)):
    score, aligned_seq1, aligned_seq2 = local_alignment(sequences[i], sequences[j], match_score, mismatch_score, gap_penalty)
    print(f"Diziler {i+1} ve {j+1} arasındaki local hizalama skoru (farklı skorlar): {score}")
    print(aligned_seq1)
    print(aligned_seq2)