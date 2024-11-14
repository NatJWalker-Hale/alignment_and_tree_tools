#! /usr/bin/env python3


"""
script to blast nucleotide orthologues. Assumes one best hit locus per DB by concatenating subject
range in the case of multiple hits per query
"""


import os
import sys
import argparse
import subprocess
from operator import itemgetter

from pyfaidx import Fasta

import sequence as sq


def make_blast_db(dbf: str):
    """
    make a blast database from the FASTA-formatted dbf
    """
    sys.stderr.write("making blastdb\n")
    cmd = ["makeblastdb", "-in", dbf, "-dbtype", "nucl", "-parse_seqids"]
    print(subprocess.list2cmdline(cmd))
    subprocess.run(cmd, shell=False, check=True)


def blast_db(query: str , dbf: str, nt: int=2) -> str:
    """
    blast the queries in bait against blastdb in dbf
    """
    sys.stderr.write("searching database with blastn\n")
    blastout = f"{dbf}.blastn.outfmt6"
    cmd = ["blastn", "-query", query, "-db", dbf, "-num_threads", str(nt), "-out",
           f"{dbf}.blastn.outfmt6", "-outfmt",
           "6 qseqid sseqid bitscore evalue pident sstart send qstart qend qcovs"]
    print(subprocess.list2cmdline(cmd))
    try:
        result = subprocess.run(cmd, shell=False, check=True, capture_output=True, text=True)
        print(f"BLAST completed successfully. Return code: {result.returncode}")

        if os.path.exists(blastout) and os.path.getsize(blastout) > 0:
            print(f"Output file '{blastout}' created successfully.")
        else:
            print(f"Output file '{blastout}' is empty or does not exist.")
            return None

        return blastout
    except subprocess.CalledProcessError as e:
        print(f"Error running BLAST: {e}")
        print(f"BLAST stderr output:\n{e.stderr}")
        return None


def parse_blastn_out(outf: str) -> tuple:
    """
    parse the hits in outf. Combine the start and end ranges of multiple hits (assumes 1 locus per
    db). In the case of multiple queries, return the hit(s) with the best (average) bitscore(s)
    """
    sys.stderr.write(f"parsing hits in {outf}\n")
    hits = {}
    with open(outf, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip().split("\t")
            # fields are qseqid sseqid bitscore evalue pident sstart send qstart qend qcovs
            sstart = int(line[5])
            send = int(line[6])
            qstart = int(line[7])
            qend = int(line[8])
            bscore = float(line[2])
            evalue = float(line[3])
            pident = float(line[4])
            qcovs = float(line[9])
            # if send < sstart:
            #     send = int(line[5])
            #     sstart = int(line[6])
            try:
                hits[line[0]].append([line[1], bscore, evalue, pident, sstart, send, qstart, qend,
                                      qcovs])
            except KeyError:
                hits[line[0]] = [[line[1], bscore, evalue, pident, sstart, send, qstart, qend,
                                      qcovs]]
    filt_hits = {}  # key is qseqid, value is ['sseqid', grade, sstart, send]
    for q, hitfo in hits.items():  # hitfo is hits-info
        if len(hitfo) > 1:
            hitfo = sorted(hitfo, key=lambda x: x[6])  # sort hits by query start
            evalue = sum(x[2] for x in hitfo) / len(hitfo)
            pident = sum(x[3] for x in hitfo) / len(hitfo)
            qcovs = sum(x[8] for x in hitfo) / len(hitfo)
            i = 0
            while i < len(hitfo) - 1:
                current_qend = hitfo[i][7]
                next_qstart = hitfo[i+1][6]
                if current_qend > next_qstart:
                    overlap = current_qend - next_qstart
                    if hitfo[i][4] < hitfo[i][5]:
                        hitfo[i+1][4] += overlap + 1
                    else:
                        hitfo[i][5] += overlap + 1
                    # in order to make this work we just have to reverse the operation in reversed
                    # cases, that is, make the adjustment to the END of the first match, not the
                    # START of the second
                i += 1
            sstarts = [x[4] for x in hitfo]
            sends = [x[5] for x in hitfo]
            # hit_len = send - sstart
            grade = ((50 * qcovs/100) +
                     (25 * max(0, 1-(evalue/(10**-20)))) +
                     (25 * max(0, ((pident-50)/50))))
            # grade is like from geneious: 50, 25, 25 weighted av of query cover, evalue, and pid
            filt_hits[q] = [hitfo[0][0], grade, sstarts, sends]
        else:
            evalue = hitfo[0][2]
            pident = hitfo[0][3]
            qcovs = hitfo[0][8]
            sstart = hitfo[0][4]
            send = hitfo[0][5]
            # hit_len = send - sstart
            grade = ((50 * qcovs/100) +
                     (25 * max(0, 1-(evalue/(10**-20)))) +
                     (25 * max(0, ((pident-50)/50))))
            filt_hits[q] = [hitfo[0][0], grade, [sstart], [send]]
    # print(filt_hits)
    sort_hits = dict(sorted(filt_hits.items(), key=itemgetter(1), reverse=True))
    best = next(iter(sort_hits.items()))
    return best


def search_dbs(query: str, db_dir: str, out_dir: str):
    """
    main function
    """
    if "/" in query:
        name = query.split("/")[-1].split(".")[0]
    else:
        name = query.split(".")[0]

    dblist = []
    for dirpath, _, filenames in os.walk(db_dir):
        for f in filenames:
            if f.endswith(".fa"):
                dblist.append(os.path.abspath(os.path.join(dirpath, f)))
    print(dblist)
    outfile = f"{os.path.abspath(out_dir)}/{name}.blastn.fa"
    if os.path.isfile(outfile):
        os.remove(outfile)  # prevent appending partial file
    out_seqs = {}
    for db in dblist:
        blastdbsuf = [".ndb", ".nhr", ".nin", ".nog", ".nos", ".not",
                      ".nsq", ".ntf", ".nto"]
        for s in blastdbsuf:
            if not os.path.isfile(db + s):
                make_blast_db(db)
        blastout = blast_db(query, db)
        if blastout is None:
            sys.stderr.write(f"no hits for {db}, skipping\n")
            continue
        hit = parse_blastn_out(blastout)
        sid = hit[1][0]
        sys.stderr.write(f"model query is {hit}\n")
        sstart = hit[1][2]
        print(sstart)
        send = hit[1][3]
        print(send)
        db_seqs = Fasta(db)
        if (len(sstart) > 1) & (len(send) > 1):
            coords = sorted(zip(sstart, send), key=lambda x: x[1])
            seq = ""
            for i in coords:
                if i[0] > i[1]:
                    seq += db_seqs[sid][i[1]-1:i[0]].seq
                else:
                    seq += db_seqs[sid][i[0]-1:i[1]].seq
        else:
            if sstart[0] > send[0]:
                seq = db_seqs[sid][send[0]-1:sstart[0]].seq
            else:
                seq = db_seqs[sid][sstart[0]-1:send[0]].seq
        samp = sid.split("_")[0]
        out_seqs[samp] = seq
    sq.write_fasta(out_seqs, outfile)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="script to parse blastn hits of nucleotide \
                                     queries from a supermatrix of orthologues. Assumes one \
                                     best hit locus per DB by concatenating subject range in \
                                     the case of multiple hits per query")
    parser.add_argument("query", help="FASTA-formatted query sequences")
    parser.add_argument("db", help="Directory containing FASTAs to search")
    parser.add_argument("out", help="Output directory")
    args = parser.parse_args(sys.argv[1:] or ["--help"])

    search_dbs(args.query, args.db, args.out)
