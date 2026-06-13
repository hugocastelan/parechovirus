
import csv

fasta_file = "/Users/hugo/Desktop/vp1_parachovirus/vp1_headers_updated_remove_v3.fasta"
output_csv = "/Users/hugo/Desktop/vp1_parachovirus/dates_treetime_v3.csv"

with open(fasta_file) as fin, open(output_csv, "w", newline="") as fout:

    writer = csv.writer(fout)
    writer.writerow(["name", "date"])

    for line in fin:
        if line.startswith(">"):

            header = line.strip()[1:]  # quitar >

            fields = header.split("|")

            # última columna = fecha
            date = fields[-1]

            writer.writerow([header, date])

print(f"Archivo generado: {output_csv}")da
