#!/usr/bin/env python3

import os
import subprocess
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


###Files 

FRAGMENTS = "/Users/hugo/Desktop/all_genes_diferenf_size.fasta"

ALIGNMENT_DIR = "/Users/hugo/Desktop/Parechovirus_Proteins"

OUTDIR = "/Users/hugo/Desktop/parachovirus_output"

GENES = [
    "VP0",
    "VP3",
    "VP1",
    "2A",
    "2B",
    "2C",
    "3A",
    "3B",
    "3C",
    "3D"
]

MIN_PID = 70
MIN_COV = 50
MIN_LEN = 90


# OUTPUT DIRS

REFDIR = os.path.join(OUTDIR, "references")
BLASTDIR = os.path.join(OUTDIR, "blast")
GENEDIR = os.path.join(OUTDIR, "genes")
FINALDIR = os.path.join(OUTDIR, "aligned")

for d in [OUTDIR, REFDIR, BLASTDIR, GENEDIR, FINALDIR]:
    os.makedirs(d, exist_ok=True)


# Load sequences 


print("\nLoading sequences")

allseqs = SeqIO.to_dict(
    SeqIO.parse(FRAGMENTS, "fasta")
)

print("Total sequences:", len(allseqs))


# Reconstruct or build a reference 


print("\nCreating references")

valid_genes = []

for gene in GENES:

    aln_file = os.path.join(
        ALIGNMENT_DIR,
        f"{gene}.fasta"
    )

    if not os.path.exists(aln_file):

        print(f"Skipping {gene}: reference not found")
        continue

    try:

        first_record = next(
            SeqIO.parse(aln_file, "fasta")
        )

    except StopIteration:

        print(f"Skipping {gene}: empty reference")
        continue

    seq = str(first_record.seq).replace("-", "")

    outfile = os.path.join(
        REFDIR,
        f"{gene}_ref.fasta"
    )

    SeqIO.write(
        [
            SeqRecord(
                Seq(seq),
                id=gene,
                description=""
            )
        ],
        outfile,
        "fasta"
    )

    valid_genes.append(gene)


# BUILD BLAST DB

print("\nBuilding BLAST database")

db_path = os.path.join(
    OUTDIR,
    "fragment_db"
)

subprocess.run(
    [
        "makeblastdb",
        "-in",
        FRAGMENTS,
        "-dbtype",
        "nucl",
        "-out",
        db_path
    ],
    check=True
)
# Recover fragments

for gene in valid_genes:

    print(f"\nProcessing {gene}")

    query = os.path.join(
        REFDIR,
        f"{gene}_ref.fasta"
    )

    blast_out = os.path.join(
        BLASTDIR,
        f"{gene}.tsv"
    )

    subprocess.run(
        [
            "blastn", #parameters of blast 
            "-query", query,
            "-db", db_path,
            "-evalue", "1e-20",
            "-max_target_seqs", "100000",
            "-outfmt",
            "6 qseqid sseqid pident length qlen sstart send bitscore",
            "-out", blast_out
        ],
        check=True
    )

    recovered = []
    seen = set()

    with open(blast_out) as fh:

        for line in fh:

            (
                qid,
                sid,
                pid,
                aln_len,
                qlen,
                sstart,
                send,
                bitscore
            ) = line.strip().split("\t")

            pid = float(pid)
            aln_len = int(aln_len)
            qlen = int(qlen)

            cov = (aln_len / qlen) * 100

            if pid < MIN_PID:
                continue

            if cov < MIN_COV:
                continue

            if sid in seen:
                continue

            seen.add(sid)

            rec = allseqs[sid]

            sstart = int(sstart)
            send = int(send)

            start = min(sstart, send) - 1
            end = max(sstart, send)

            subseq = rec.seq[start:end]

            if sstart > send:
                subseq = subseq.reverse_complement()

            extra = len(subseq) % 3

            if extra:
                subseq = subseq[:-extra]

            if len(subseq) < MIN_LEN:
                continue

            header = rec.description

            fields = header.split("|")

            accession = fields[0]

            if len(fields) >= 5:

                virus = fields[1]
                subtype = fields[2]
                country = fields[3]
                year = fields[4]

                if "-" in year:
                    year = year.split("-")[0]

                new_header = (
                    f"{accession}|"
                    f"{virus}|"
                    f"{subtype}|"
                    f"{country}|"
                    f"{year}"
                )

            else:

                new_header = accession

            recovered.append(
                SeqRecord(
                    subseq,
                    id=new_header,
                    name="",
                    description=""
                )
            )

    gene_fasta = os.path.join(
        GENEDIR,
        f"{gene}_fragments.fasta"
    )

    SeqIO.write(
        recovered,
        gene_fasta,
        "fasta"
    )

    print(
        f"{gene}: "
        f"{len(recovered)} sequences recovered"
    )


# run MAFFT


print("\nRunning MAFFT")

for gene in valid_genes:

    ref = os.path.join(
        ALIGNMENT_DIR,
        f"{gene}.fasta"
    )

    fragments = os.path.join(
        GENEDIR,
        f"{gene}_fragments.fasta"
    )

    outfile = os.path.join(
        FINALDIR,
        f"{gene}_final.fasta"
    )

    if not os.path.exists(fragments):

        print(f"Skipping {gene}: no fragments file")
        continue

    n_frag = sum(
        1 for _ in SeqIO.parse(
            fragments,
            "fasta"
        )
    )

    if n_frag == 0:

        print(f"Skipping {gene}: 0 sequences")
        continue

    print(
        f"Aligning {gene} "
        f"({n_frag} sequences)"
    )

    cmd = (
        f"mafft --addfragments "
        f"{fragments} "
        f"{ref} "
        f"> {outfile}"
    )

    try:

        subprocess.run(
            cmd,
            shell=True,
            check=True
        )

    except subprocess.CalledProcessError:

        print(
            f"MAFFT failed for {gene}"
        )

        continue

print("\nDONE")
