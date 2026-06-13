#!/usr/bin/env python3

from Bio import SeqIO

# Archivos
fasta_in = "/Users/hugo/Desktop/vp1_parachovirus/vp1_all_complete_renamed_v2.fasta"
headers_remove = "/Users/hugo/Desktop/vp1_parachovirus/id_remove.txt"
fasta_out = "/Users/hugo/Desktop/vp1_parachovirus/vp1_all_complete_renamed_filter_v2.fasta"

# Leer headers a eliminar
with open(headers_remove) as f:
    remove_set = set(
        line.strip().replace(">", "")
        for line in f
        if line.strip()
    )

# Filtrar secuencias
kept = []

for record in SeqIO.parse(fasta_in, "fasta"):
    if record.id not in remove_set:
        kept.append(record)

# Guardar nuevo fasta
SeqIO.write(kept, fasta_out, "fasta")

print(f"Secuencias originales: {sum(1 for _ in SeqIO.parse(fasta_in, 'fasta'))}")
print(f"Secuencias eliminadas: {len(remove_set)}")
print(f"Secuencias conservadas: {len(kept)}")
print(f"Archivo guardado: {fasta_out}")re
