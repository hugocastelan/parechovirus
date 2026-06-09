#!/usr/bin/env python3

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import os

#input 


gbk_file = "/Users/hugo/Desktop/sequence_complete_genome.gb"

outdir = "/Users/hugo/Desktop/Parechovirus_Proteins"

os.makedirs(outdir, exist_ok=True)


# File for each protein 

proteins = {
    "VP0": [],
    "VP3": [],
    "VP1": [],
    "2A": [],
    "2B": [],
    "2C": [],
    "3A": [],
    "3B": [],
    "3C": [],
    "3D": []
}


# run inside the gbk file to gen the mat_peotide 


n = 0

for record in SeqIO.parse(gbk_file, "genbank"):

    n += 1

    accession = record.id

    for feature in record.features:

        if feature.type != "mat_peptide":
            continue

        product = feature.qualifiers.get(
            "product",
            ["unknown"]
        )[0]

        seq = feature.extract(record.seq)

        product = product.upper()

        if product == "VP0":

            proteins["VP0"].append(
                SeqRecord(
                    seq,
                    id=accession,
                    description=record.description
                )
            )

        elif product == "VP3":

            proteins["VP3"].append(
                SeqRecord(
                    seq,
                    id=accession,
                    description=record.description
                )
            )

        elif product == "VP1":

            proteins["VP1"].append(
                SeqRecord(
                    seq,
                    id=accession,
                    description=record.description
                )
            )

        elif product == "2A":

            proteins["2A"].append(
                SeqRecord(
                    seq,
                    id=accession,
                    description=record.description
                )
            )

        elif product == "2B":

            proteins["2B"].append(
                SeqRecord(
                    seq,
                    id=accession,
                    description=record.description
                )
            )

        elif product == "2C":

            proteins["2C"].append(
                SeqRecord(
                    seq,
                    id=accession,
                    description=record.description
                )
            )

        elif product == "3A":

            proteins["3A"].append(
                SeqRecord(
                    seq,
                    id=accession,
                    description=record.description
                )
            )

        elif product == "3B":

            proteins["3B"].append(
                SeqRecord(
                    seq,
                    id=accession,
                    description=record.description
                )
            )

        elif product == "3C":

            proteins["3C"].append(
                SeqRecord(
                    seq,
                    id=accession,
                    description=record.description
                )
            )

        elif product == "3D":

            proteins["3D"].append(
                SeqRecord(
                    seq,
                    id=accession,
                    description=record.description
                )
            )

print(f"\nGenomas processed: {n}")

# Write fasts 

for protein, records in proteins.items():

    outfile = os.path.join(
        outdir,
        f"{protein}.fasta"
    )

    SeqIO.write(
        records,
        outfile,
        "fasta"
    )

    print(
        f"{protein}: "
        f"{len(records)} secuencias"
    )

print("\nFinished.")
