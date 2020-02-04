#!/usr/bin/env python

import argparse
import os
import re
import sys
import os
from Bio.Blast import NCBIXML
from Bio.Seq import Seq
from Bio import SeqIO


parser = argparse.ArgumentParser()
parser.add_argument("--fasta", "-f", help="Inform the fasta file containing your sequence(s) of interest. If the"
                                            "sgRNA is going to be designed to target an ORF, the genes present in the"
                                            "fasta file must be on the non-template strand, 5' > 3'.", required=True)
parser.add_argument("--genome", "-g", help="Inform the fasta file containing the genome of the organism. It will be"
                                           "used to look for possible off-targets. If not informed, the program will"
                                           "not look for any off-targets")
parser.add_argument("--outdir", "-o", help="The output directory where your files will be stored. By default, they will"
                                           "be stored in CRISPRi directory.")
args = parser.parse_args()


def find_off_targets(genome, direc, outdir):
    cmd_blastn = 'blastn -query %s/%s/sgRNAs.fasta -task blastn-short -subject %s ' \
                 '-outfmt 5 -qcov_hsp_perc 100 -out %s/%s/blasted_sgrnas.xml' % (outdir, direc, genome, outdir, direc)
    os.system(cmd_blastn)
    if os.stat("%s/%s/blasted_sgrnas.xml" % (outdir, direc)).st_size != 0:
        off_handle = open("%s/%s/blasted_sgrnas.xml" % (outdir, direc), 'r')
        blast_records = NCBIXML.parse(off_handle)
        with open("%s/%s/off_targets.xlsx" % (outdir, direc), 'w') as off:
            results = ["sgRNA_name\tsgRNA_Seq\tTarget_Seq\tStart\tEnd\tOff-target_number\n"]
            for record in blast_records:
                for alignment in record.alignments:
                    i = 0
                    for hsp in alignment.hsps:
                        if hsp.identities >= (hsp.align_length-2):
                            i += 1
                            results.append(str(record.query)[-20:] + "\t" + hsp.query + "\t" + hsp.sbjct + "\t" +
                                           str(hsp.sbjct_start) + "\t" + str(hsp.sbjct_end) + "\t" + "%s" % i + "\n")

            off.writelines(results)
        Cmd_rm = 'rm %s/%s/blasted_sgrnas.xml' % (outdir, direc)
        os.system(Cmd_rm)


def cas9(f, outdir):
    PAMs = ['CTTCT', 'ATTCT', 'TTTCT', 'CTTCC', 'GTTCT', 'TTTCC', 'ATGCT', 'CTCCT', 'ATCCT', 'TTGCT', 'GTTCC', 'ATTCC',
            'CTGCT', 'TTCCT', 'GTCCT', 'CTCCC', 'ATCCC', 'TTCCC', 'GTGCT', 'GTCCC', 'ATGCC', 'CTGCC', 'TTGCC', 'GTGCC']
    Fold_rep = ['216.7', '216.2', '158.1', '145.2', '120.5', '110.5', '84.6', '82.2', '64.7', '53.4', '51.5', '47.3',
                '42.2', '38.5', '25.5', '24.7', '24.2', '12.3', '11.9', '7.9', '6.7', '4.0', '3.3', '2.7', '1.3']
    folder_number = 0
    for seq_record in SeqIO.parse(f, "fasta"):
        folder_number += 1
        protospacers = []
        proto_fasta = []
        # parses through each entry of a fasta file and adds to an object the complement of said entry
        rev = seq_record.seq
        # rev2 = rev.reverse_complement()
        stringseq = str(rev)
        # splits the entry by this specific PAM sequence and removes the first element
        # the removed element is not relevant as it corresponds either to the PAM itself or to a sequence 3' to the PAM
        # in the template strand
        j = 0
        for pam in PAMs:
            seq_elements = re.split(r'%s' % pam, stringseq)[1:]
            # looks for elements with length bigger than 22
            for oligo in seq_elements:
                if len(oligo) >= 22:
                    # generates a complement for sequences with length > 22
                    comples = Seq(oligo)
                    # complement = comples.complement()
                    k = 21
                    # if the last nucleotide of the sequence is not a T or a C, it will extend the sequence by 1
                    # nucleotide until it finds either a T or a C. If none is found, the sequence will not be written.
                    j += 1
                    while k < 27:
                        if len(comples) > k:
                            if comples[k] == "T" or comples[k] == "C":
                                # writes to protospacers object the entry's name and the fold repression reached with
                                # that PAM
                                i = PAMs.index(pam)
                                # adds a CCCA to 5' extremity of forward sequence and a AAAC to 3' extremity of reverse
                                # sequence
                                # writes to protospacers only the nucleotides between position 2 and last position
                                # containing T or C, thus removing the first two nucleotides, which are still part of
                                # the PAM
                                forward = comples[2:k+1]
                                ind = seq_record.seq.find(forward)
                                final = ind+len(forward)
                                range_proto = "%s:" % ind + "%s" % (final-1)
                                gene_l = len(seq_record.seq)
                                protospacers.append("> " + seq_record.id + "|sgRNA %s Fold repression: %s | Range in Gene: %s | Gene length: %s" % (j, Fold_rep[i], range_proto, gene_l))
                                protospacers.append("F: 5' GGGA" + forward.reverse_complement() + " 3'\nR: 5' AAAC" +
                                                    comples[2:k + 1] + " 3'")
                                proto_fasta.append(">" + seq_record.id + "_sgRNA_%s" % j + "\n" + forward.reverse_complement())
                                break
                            if comples[k] != "T" or comples[k] != "C":
                                if len(comples) > k:
                                    k += 1
                                else:
                                    break
                        else:
                            break
        # writes results stored in protospacers object to a .fasta file
        position = str(seq_record.id).find("gene")
        end_position = str(seq_record.id).find("xref=")
        if position != -1 and end_position != -1:
            folder_name = str(seq_record.id[position-1:end_position-5])
        elif position == -1 and end_position != -1:
            folder_name = str(seq_record.id[end_position-24:end_position-5])
        elif end_position == -1 and position != -1:
            folder_name = str(seq_record.id[position-1:position+25])
        elif position == -1 and end_position == -1:
            folder_name = str(seq_record.id[15:40])
        else:
            folder_name = "ERR"
        folder_name = folder_name.replace("(", "_")
        folder_name = folder_name.replace(")", "_")
        cmd_direc = "mkdir %s/%s_%s" % (outdir, folder_number, folder_name)
        os.system(cmd_direc)
        if os.path.exists("%s/%s_%s/" % (outdir, folder_number, folder_name)):
            with open('%s/%s_%s/primers.txt' % (outdir, folder_number, folder_name), 'w') as filehandle, open('%s/%s_%s/sgRNAs.fasta' % (outdir, folder_number, folder_name), 'w') as fa:
                filehandle.writelines("%s \n" % item for item in protospacers)
                fa.writelines("%s \n" % entry for entry in proto_fasta)
            if args.genome is not None:
                find_off_targets(args.genome, "%s_%s" % (folder_number, folder_name), outdir)


if __name__ == '__main__':
    outdir = "CRISPRi"
    if args.outdir is not None:
        Cmd_dir = 'mkdir %s' % args.outdir
        outdir = "%s" % args.outdir
        os.system(Cmd_dir)
    else:
        Cmd_dir = 'mkdir CRISPRi'
        os.system(Cmd_dir)
    cas9(args.fasta, outdir)



