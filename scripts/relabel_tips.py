#!/usr/bin/env python3

"""
Given a tree where some tips have missing metadata
(e.g. genotype, subtype, country), impute these values
by propagating labels from neighbouring tips.
"""

import argparse
import sys
import random
from io import StringIO
from Bio import Phylo


###############################################################
# SAFE BRANCH LENGTH
###############################################################

def get_branch_length(node):
    """
    Return branch length.
    If branch length is None, treat it as 0.
    """

    if node.branch_length is None:
        return 0.0

    return float(node.branch_length)


###############################################################
# CLIMB DOWN A CLADE
###############################################################

def climb(tips, curnode, pathlen, cutoff):

    """
    Recursive function for traversing down a tree and collecting
    terminal nodes within the specified distance cutoff.
    """

    pathlen += get_branch_length(curnode)

    if pathlen < cutoff:

        if curnode.is_terminal():

            tips.append(
                (
                    pathlen,
                    curnode.name
                )
            )

        else:

            for child in curnode.clades:

                tips = climb(
                    tips,
                    child,
                    pathlen,
                    cutoff
                )

    return tips


###############################################################
# FIND NEAREST TIPS
###############################################################

def nearest(curnode, cutoff, parents):

    """
    Find neighbouring tips within cutoff distance.
    """

    tips = []

    pathlen = get_branch_length(curnode)

    # Root has no parent
    if curnode not in parents:
        return tips

    p = parents[curnode]

    previous_node = curnode

    while True:

        if pathlen >= cutoff:
            break

        for child in p.clades:

            # Do not travel back toward the original branch
            if child == previous_node:
                continue

            if child.is_terminal():

                this_dist = (
                    pathlen
                    + get_branch_length(child)
                )

                if this_dist < cutoff:

                    tips.append(
                        (
                            this_dist,
                            child.name
                        )
                    )

            else:

                tips.extend(
                    climb(
                        [],
                        child,
                        pathlen,
                        cutoff
                    )
                )

        # Reached root
        if p not in parents:
            break

        previous_node = p

        pathlen += get_branch_length(p)

        p = parents[p]

    return tips


###############################################################
# CONSENSUS
###############################################################

def consensus(poll):

    """
    Determine consensus label among neighbours.

    If there is a tie, randomly choose ONE STRING.

    IMPORTANT:
    random.choice() returns a string.
    random.sample(..., 1) returns a list and breaks join().
    """

    if not poll:

        raise ValueError(
            "Cannot calculate consensus from empty poll."
        )

    counts = {}

    for label in poll:

        counts[label] = (
            counts.get(label, 0)
            + 1
        )

    max_count = max(
        counts.values()
    )

    winners = [
        label
        for label, count in counts.items()
        if count == max_count
    ]

    # One winner
    if len(winners) == 1:

        return winners[0]

    # Tie: return ONE STRING
    return random.choice(
        winners
    )


###############################################################
# READ TREE SAFELY
###############################################################

def read_tree_safely(nwkfile):

    """
    Read Newick while fixing known problematic geographic names.

    Apostrophes in unquoted Newick labels can destroy the tip name.
    Example:

        Cote d'Ivoire

    Replace these with a Newick-safe form.
    """

    text = nwkfile.read()

    replacements = {

        "Cote d'Ivoire":
            "Cote_d_Ivoire",

        "Côte d'Ivoire":
            "Cote_d_Ivoire",

        "Cote_d'Ivoire":
            "Cote_d_Ivoire",

        "Côte_d'Ivoire":
            "Cote_d_Ivoire"
    }

    for old, new in replacements.items():

        text = text.replace(
            old,
            new
        )

    return Phylo.read(
        StringIO(text),
        "newick"
    )


###############################################################
# RELABEL TIPS
###############################################################

def relabel_tips(
    nwkfile,
    cutoff,
    delim,
    field,
    match="",
    ignore_case=False,
    k=1,
    verbose=False,
    **kwargs
):

    """
    Use k-nearest neighbours in the tree to impute missing values.

    Example header:

    accession|virus|genotype|country|date

    field=2 -> genotype
    """

    if ignore_case:

        def normalize(x):
            return str(x).lower()

    else:

        def normalize(x):
            return str(x)

    match_normalized = normalize(
        match
    )


    ###########################################################
    # READ TREE
    ###########################################################

    phy = read_tree_safely(
        nwkfile
    )


    ###########################################################
    # STORE PARENTS
    ###########################################################

    parents = {}

    for node in phy.find_clades(
        order="level"
    ):

        for child in node.clades:

            parents[child] = node


    ###########################################################
    # LOCATE MISSING TIPS
    ###########################################################

    missing = set()

    labels = {}

    malformed = []

    tips = phy.get_terminals()


    for tip in tips:

        if tip.name is None:

            malformed.append(
                "None"
            )

            continue

        values = tip.name.split(
            delim
        )


        #######################################################
        # CHECK HEADER STRUCTURE
        #######################################################

        if (
            field >= len(values)
            or
            field < -len(values)
        ):

            sys.stderr.write(
                f"Failed to locate field {field} "
                f"in header: {values}\n"
            )

            malformed.append(
                tip.name
            )

            continue


        val = values[field]


        #######################################################
        # MISSING LABEL
        #######################################################

        if normalize(val) == match_normalized:

            missing.add(
                tip
            )

        else:

            labels[
                tip.name
            ] = val


    ###########################################################
    # SUMMARY
    ###########################################################

    if verbose:

        sys.stderr.write(
            f"Found {len(missing)} tips "
            f"with missing labels.\n"
        )

        if malformed:

            sys.stderr.write(
                f"Found {len(malformed)} malformed "
                f"tip labels.\n"
            )


    ###########################################################
    # IMPUTE EACH MISSING TIP
    ###########################################################

    n_changed = 0

    unresolved = []


    for tip in missing:

        neighbours = nearest(
            tip,
            cutoff,
            parents
        )


        #######################################################
        # KEEP ONLY NEIGHBOURS WITH KNOWN LABELS
        #######################################################

        neighbours = [
            (distance, name)
            for distance, name
            in neighbours
            if name in labels
        ]


        #######################################################
        # ORDER NEAREST -> FARTHEST
        #######################################################

        neighbours.sort(
            key=lambda x: x[0]
        )


        #######################################################
        # K NEAREST
        #######################################################

        knn = neighbours[:k]


        #######################################################
        # NO NEIGHBOURS
        #######################################################

        if len(knn) == 0:

            sys.stderr.write(
                f"\nWARNING: no neighbours for:\n"
                f"{tip.name}\n"
                f"within cutoff={cutoff}\n\n"
            )

            unresolved.append(
                tip.name
            )

            # Do NOT kill entire script
            continue


        #######################################################
        # K > 1 -> CONSENSUS
        #######################################################

        if k > 1:

            poll = [
                labels[name]
                for _, name
                in knn
            ]

            newval = consensus(
                poll
            )


        #######################################################
        # K = 1
        #######################################################

        else:

            newval = labels[
                knn[0][1]
            ]


        #######################################################
        # SAFETY: MAKE ABSOLUTELY SURE newval IS STRING
        #######################################################

        if isinstance(
            newval,
            list
        ):

            if len(newval) == 0:

                unresolved.append(
                    tip.name
                )

                continue

            newval = str(
                newval[0]
            )

        else:

            newval = str(
                newval
            )


        #######################################################
        # UPDATE TIP LABEL
        #######################################################

        values = tip.name.split(
            delim
        )

        old_name = tip.name

        values[field] = newval

        tip.name = delim.join(
            values
        )

        n_changed += 1


        #######################################################
        # ADD NEWLY IMPUTED LABEL
        #
        # This lets later missing tips use it if encountered.
        #######################################################

        labels[
            tip.name
        ] = newval


        #######################################################
        # VERBOSE
        #######################################################

        if verbose:

            sys.stderr.write(
                f"{old_name} -> "
                f"{tip.name}\n"
            )


    ###########################################################
    # FINAL SUMMARY
    ###########################################################

    if verbose:

        sys.stderr.write(
            "\n"
            + "=" * 70
            + "\n"
        )

        sys.stderr.write(
            "RELABEL SUMMARY\n"
        )

        sys.stderr.write(
            "=" * 70
            + "\n"
        )

        sys.stderr.write(
            f"Missing tips initially: "
            f"{len(missing)}\n"
        )

        sys.stderr.write(
            f"Successfully relabeled: "
            f"{n_changed}\n"
        )

        sys.stderr.write(
            f"Unresolved: "
            f"{len(unresolved)}\n"
        )

        sys.stderr.write(
            f"Malformed headers: "
            f"{len(malformed)}\n"
        )

        if unresolved:

            sys.stderr.write(
                "\nUNRESOLVED IDS:\n"
            )

            for name in unresolved:

                sys.stderr.write(
                    name
                    + "\n"
                )

        if malformed:

            sys.stderr.write(
                "\nMALFORMED HEADERS:\n"
            )

            for name in malformed:

                sys.stderr.write(
                    name
                    + "\n"
                )


    return phy


###############################################################
# MAIN
###############################################################

if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description=__doc__
    )


    parser.add_argument(
        "nwkfile",
        type=argparse.FileType("r"),
        help="Path to Newick tree."
    )


    parser.add_argument(
        "-o",
        "--outfile",
        type=argparse.FileType("w"),
        default=sys.stdout,
        help=(
            "Output Newick file. "
            "Default: stdout."
        )
    )


    parser.add_argument(
        "-d",
        "--delim",
        type=str,
        default="_",
        help=(
            "Character separating fields "
            "in tip labels."
        )
    )


    parser.add_argument(
        "-f",
        "--field",
        type=int,
        default=-1,
        help=(
            "0-based index of metadata "
            "field to impute."
        )
    )


    parser.add_argument(
        "-m",
        "--match",
        type=str,
        default="",
        help=(
            "Value representing missing "
            "metadata."
        )
    )


    parser.add_argument(
        "-i",
        "--ignore_case",
        action="store_true",
        default=False,
        help="Ignore case when matching."
    )


    parser.add_argument(
        "--cutoff",
        type=float,
        default=0.1,
        help=(
            "Maximum phylogenetic "
            "distance."
        )
    )


    parser.add_argument(
        "-k",
        type=int,
        default=1,
        help=(
            "Number of nearest neighbours "
            "used for imputation."
        )
    )


    parser.add_argument(
        "--format",
        default="newick",
        help="Tree format."
    )


    parser.add_argument(
        "--verbose",
        action="store_true",
        help="Print detailed information."
    )


    args = parser.parse_args()


    ###########################################################
    # RUN
    ###########################################################

    tree = relabel_tips(
        **vars(args)
    )


    ###########################################################
    # WRITE TREE
    ###########################################################

    Phylo.write(
        tree,
        args.outfile,
        format="newick"
    )