#!/usr/bin/env python3

import re
import sys
import tempfile
import os
import argparse

tempdir = tempfile.gettempdir()
if not os.access(tempdir, os.W_OK):
    tempdir = os.getcwd()

cigar_re = re.compile(r'[0-9]+[MIDNSHPX=]')
indel_re = re.compile(r'[+-][0-9]+')


def is_first_read(flag):
    IS_FIRST_SEGMENT = 0x40
    return (int(flag) & IS_FIRST_SEGMENT) != 0


def sam_to_pileup(handle):

    pileup = {}
    counts = {}

    for row in handle:
        if row.startswith('@'):
            continue

        fields = row.strip().split('\t')
        if len(fields) < 11:
            continue

        qname, flag, refname, rpos, mapq, cigar, rnext, pnext, tlen, seq, qual = fields[:11]

        if cigar == '*':
            continue

        pileup.setdefault(refname, {})
        counts.setdefault(refname, 0)
        counts[refname] += 1

        is_first = is_first_read(flag)

        pos = 0
        refpos = int(rpos)

        tokens = cigar_re.findall(cigar)

        if tokens and tokens[0].endswith('S'):
            pos = int(tokens[0][:-1])
            tokens.pop(0)

        if not tokens or not tokens[0].endswith('M'):
            print("ERROR: CIGAR token after soft clip must be match interval")
            sys.exit(1)

        if refpos not in pileup[refname]:
            pileup[refname][refpos] = {'s': '', 'q': ''}

        pileup[refname][refpos]['s'] += '^' + chr(int(mapq) + 33)

        for token in tokens:
            length = int(token[:-1])

            if token.endswith('M'):
                for _ in range(length):
                    pileup[refname].setdefault(refpos, {'s': '', 'q': ''})
                    base = seq[pos] if is_first else seq[pos].lower()
                    pileup[refname][refpos]['s'] += base
                    pileup[refname][refpos]['q'] += qual[pos]
                    pos += 1
                    refpos += 1

            elif token.endswith('D'):
                pileup[refname][refpos - 1]['s'] += '-' + str(length) + (
                    'N' if is_first else 'n') * length

                for i in range(refpos, refpos + length):
                    pileup[refname].setdefault(i, {'s': '', 'q': ''})
                    pileup[refname][i]['s'] += '*'

                refpos += length

            elif token.endswith('I'):
                insert = seq[pos:pos + length]
                pileup[refname][refpos - 1]['s'] += '+' + str(length) + (
                    insert if is_first else insert.lower())
                pos += length

            elif token.endswith('S'):
                break

            else:
                print("ERROR: Unknown token in CIGAR string", token)
                sys.exit(1)

        pileup[refname][refpos - 1]['s'] += '$'

    return pileup, counts


def pileup_to_conseq(pileup, qCutoff):

    conseq = ''
    to_skip = 0
    last_pos = 0

    for pos in sorted(pileup.keys()):

        astr = pileup[pos]['s']
        qstr = pileup[pos]['q']

        if to_skip > 0:
            to_skip -= 1
            continue

        if (pos - last_pos) > 1:
            conseq += 'N' * (pos - last_pos - 1)

        last_pos = pos
        alist = []
        i = 0
        j = 0

        while i < len(astr):

            if astr[i] == '^':
                q = ord(qstr[j]) - 33
                base = astr[i + 2] if q >= qCutoff else 'N'
                alist.append(base.upper())
                i += 3
                j += 1

            elif astr[i] == '*':
                alist.append('-')
                i += 1

            elif astr[i] == '$':
                i += 1

            elif i < len(astr) - 1 and astr[i + 1] in '+-':

                m = indel_re.match(astr[i + 1:])
                indel_len = int(m.group().strip('+-'))
                left = i + 1 + len(m.group())
                insertion = astr[left:left + indel_len]

                q = ord(qstr[j]) - 33
                base = astr[i].upper() if q >= qCutoff else 'N'
                token = base + m.group() + insertion.upper()

                if astr[i + 1] == '+':
                    alist.append(token)
                else:
                    alist.append(base)

                i += len(token)
                j += 1

            else:
                q = ord(qstr[j]) - 33
                base = astr[i].upper() if q >= qCutoff else 'N'
                alist.append(base)
                i += 1
                j += 1

        if alist:
            token = max(set(alist), key=alist.count)
        else:
            token = 'N'

        if '+' in token:
            m = indel_re.findall(token)[0]
            conseq += token[0]
            if int(m.strip('+')) % 3 == 0:
                conseq += token[1 + len(m):]

        elif token == '-':
            conseq += '-'

        else:
            conseq += token

        conseq = re.sub(r'([ACGT])(---)+([ACGT])', r'\1\3', conseq)

    return conseq


def convert_fasta(handle):
    result = []
    header = None
    sequence = ''

    for line in handle:
        line = line.strip()
        if not line:
            continue
        if line.startswith('>') or line.startswith('#'):
            if header:
                result.append((header, sequence))
            header = line[1:]
            sequence = ''
        else:
            sequence += line

    if header:
        result.append((header, sequence))

    return result


def main():
    parser = argparse.ArgumentParser(
        description="Generate updated reference from SAM consensus."
    )

    parser.add_argument('sam', type=argparse.FileType('r'))
    parser.add_argument('ref', type=argparse.FileType('r'))
    parser.add_argument('out', type=argparse.FileType('w'))
    parser.add_argument('-qcut', type=int, default=10)

    args = parser.parse_args()

    refseqs = dict((h.split()[0], s)
                   for h, s in convert_fasta(args.ref))

    pileup, counts = sam_to_pileup(args.sam)

    for refname, pile in pileup.items():
        conseq = pileup_to_conseq(pile, args.qcut)

        # Nota: update_reference usa mafft externamente.
        # Mantengo igual lógica.
        print(f"{refname}: consensus length {len(conseq)}")

        args.out.write(f">{refname}\n{conseq}\n")


if __name__ == "__main__":
    main()