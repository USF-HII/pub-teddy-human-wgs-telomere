#!/usr/bin/env python3

import argparse
import os
import sys

def main(argv):
    parser = argparse.ArgumentParser(os.path.basename(__file__))

    parser.add_argument("--decoder", required=True)
    parser.add_argument("--separator", default="\t")
    parser.add_argument("--column", default=1, type=int)
    parser.add_argument("--override-phase", action="store_true")
    parser.add_argument("fname")


    args = parser.parse_args(argv[1:])

    unmask_col = args.column - 1

    decoder_db = {}

    with open(args.decoder) as f:
        for line in f:
            k, v = line.rstrip("\n").split(":")
            decoder_db[v] = k

    entry_db = {}
    subject_ids = set()

    with open(args.fname) as f:
        header_fields = next(f).rstrip("\n").split(args.separator)

        for line in f:
            fields = line.rstrip("\n").split(args.separator)

            masked_barcode = fields[unmask_col][:-3] # strip the _P1 or _P2

            phase = fields[unmask_col][-1] # get the 1 or 2 for phase

            subject_id = decoder_db[masked_barcode]

            output = args.separator.join(fields[:unmask_col] + [subject_id] + fields[unmask_col + 1:])

            subject_ids.add(subject_id)

            entry_db.setdefault(phase, {})[subject_id] = output

    print(args.separator.join(header_fields[:unmask_col] + ["subject_id"] + header_fields[unmask_col + 1:]))

    for subject_id in sorted(subject_ids):
        for phase in sorted(entry_db, reverse=True):
            if subject_id in entry_db[phase]:
                print(entry_db[phase][subject_id])

                if args.override_phase:
                    break

if __name__ == "__main__":
    main(sys.argv)
