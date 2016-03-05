#!/usr/bin python
import sys
import os
import re
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

__author__ = 'albrecht.naumann@gmx.de'

"""
The script extract subsequences from a given fasta file by a given
coordinate txt. All extracted subsequences resulting in a output fasta file.
"""


def is_valid_file(s_file_path, s_file_type):
    """return True if file exists"""
    return (os.path.isfile(s_file_path) and
            s_file_path.lower().endswith(s_file_type))


def parse_coordinate_line(line):
    """
    parse a line of the coordinate file and return a coordinate dictionary
    """

    # return if line is empty or a comment
    l = re.sub(re.compile(r'\s+'), '', line)
    if l == '' or l.startswith("#"):
        return

    pttn = re.compile(r"(.+)\s+([0-9]+)\s+([0-9]+)")
    m = re.match(pttn, line)

    if not m:
        raise Exception("ERROR coordinate file error in line:\n\t%s"
                        % (line)
                        )
    # create coordinate dictionary
    c = {
        "id": m.group(1),
        "start": int(m.group(2)),
        "end": int(m.group(3))
    }

    if c["end"] < c["start"]:
        raise Exception(
            "coordinate file ERROR: end coordinate is larger then the" +
            " start coordinate in line:\n\t%s"
            % (line)
        )

    return c


def main(argv):
    fasta_file = False
    coordinate_file = False
    output_file = False

    args = ' '.join(argv)

    for arg in ['-s', '-c', '-o']:
        regex = re.compile(r"%s\s+(\S+)" % (arg))
        m = re.search(regex, args)

        if not m:
            sys.exit(
                "ERROR: Missing parameter\n" +
                "read README.md for more information")

        if arg == '-s' and is_valid_file(m.group(1), '.fasta'):
            fasta_file = m.group(1)
        elif arg == '-c' and is_valid_file(m.group(1), '.txt'):
            coordinate_file = m.group(1)
        elif arg == '-o':
            output_file = m.group(1)

    # is a fasta file passed
    if not fasta_file:
        sys.exit(
            "ERROR: Missing fasta file\nread README.md for more information")
    if not coordinate_file:
        sys.exit(
            "ERROR: Missing coordinate " +
            "file\nread README.md for more information")
    if not output_file:
        sys.exit(
            "ERROR: Missing output fasta " +
            "file\nread README.md for more information")

    coordinates = {}
    # read given coordinate file and save each coordinate
    with open(coordinate_file) as fh:
        for line in fh:
            # save sequence ID start and end in coordinates
            try:
                result = parse_coordinate_line(line)
                if not result:
                    continue

                if not result['id'] in coordinates:
                    coordinates[result['id']] = []

                c = {"start": result['start'], "end": result['end']}

                if c not in coordinates[result['id']]:
                    coordinates[result['id']].append(c)

            except Exception as e:
                sys.exit(e)

    out_sequences = []

    # read given fasta file and process each sequence record
    for seq_record in SeqIO.parse(fasta_file, "fasta"):

        if seq_record.id in coordinates:
            for c in coordinates[seq_record.id]:
                # extract fragment sequence
                fragment = Seq(str(seq_record.seq)[c["start"] - 1:c["end"]])
                # fragment description
                description = "fragment start: %s end: %s" % (
                    str(c["start"]), str(c["end"]))
                # create output record
                out_record = SeqRecord(
                    fragment,
                    id=seq_record.id,
                    description="<%s>" % (description)
                )
                # append output record to output sequences
                out_sequences.append(out_record)

                fh = open(output_file, "w")
                SeqIO.write(out_sequences, fh, "fasta")
                fh.close()

if __name__ == "__main__":
    main(sys.argv)
    print("done")
