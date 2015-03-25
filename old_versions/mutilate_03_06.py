#!/usr/bin/env python2.7
from __future__ import division,print_function
from Bio import SeqIO
from Bio.Alphabet import generic_dna
import argparse
import gzip
import urllib
import collections
import re
import sys
#from parser_modules import read_fasta_with_quality, read_fasta
""" 
Author : Arjun Arkal Rao
Affiliation : UCSC BME
File : /inside/depot4/users/arjun/temp_dir/parse_refflat_enspep.py 

Program info can be found in the docstring of the main function
"""

def file_type(filename):
    """
    This module is used to open an input file in the appropriate method 
    depending on whether it is a regular file, a gzipped file, or a url.
    """
    if filename.startswith("http") or filename.startswith("www") or \
            filename.startswith("ftp"):
        return urllib.urlopen(filename)
    elif filename.endswith('.gz'):
        return gzip.GzipFile(filename,'r')
    else:
        return open(filename, 'r') 

def parse_ensGene_file(ens_file):
    """
    This function accepts an ensgene file (obtained from 
    http://genome.ucsc.edu/cgi-bin/hgTables).  The columns are as follows:
    [0]  - bin
    [1]  - name
    [2]  - chrom
    [3]  - strand
    [4]  - txStart
    [5]  - txEnd
    [6]  - cdsStart
    [7]  - cdsEnd
    [8]  - exonCount
    [9]  - exonStarts
    [10] - exonEnds
    [11] - score
    [12] - name2
    [13] - cdsStartStat
    [14] - cdsEndStat
    [15] - exonFrames
    """
    parsed_object=collections.Counter()
    for line in ens_file:
        if line.startswith("#"):
            continue
        line=line.strip().split("\t")
        assert len(line)==16, "%s doesn't have 15 columns of data. " %line[1]
        ensgene=line[1]
        parsed_object[ensgene]=collections.Counter()
        parsed_object[ensgene]["chr"] = line[2]
        parsed_object[ensgene]["strand"] = line[3]
        parsed_object[ensgene]["cds_start"] = int(line[6])
        parsed_object[ensgene]["cds_end"] = int(line[7])
        parsed_object[ensgene]["exon_starts"] = [int(x) for x in 
                                                 line[9].split(",")[:-1]]
        parsed_object[ensgene]["exon_ends"] = [int(x) for x in 
                                               line[10].split(",")[:-1]]
        parsed_object[ensgene]["exonframes"] = [int(x) for x in 
                                                line[15].split(",")[:-1]]
        assert len(parsed_object[ensgene]["exon_starts"]) == \
                  len(parsed_object[ensgene]["exonframes"]), "%s had different\
                    number of entries in exon_starts and exonframes" %ensgene
        assert len(parsed_object[ensgene]["exon_ends"]) == \
                   len(parsed_object[ensgene]["exonframes"]), "%s had different\
                    number of entries in exon_ends and exonframes" %ensgene
    return parsed_object


def parse_codon_scheme(codonfile):
    """
    This function accepts a codon file with the following columns
    [0] - AA - One letter amino acid representation E.g. V
    [1] - 3LN - 3 letter amino acid name            E.g. Val
    [2] - Amino_Acid_name                           E.g. Valine
    [3] - codons - comma separated list of codons   E.g. GTA,GTC,GTG,GTT

    The columns need not be in the same order as above, however the must have
    same column names on a line starting with the "#" character.  In the absence
    of a header, it will default to the above.  The codons will be stored in a
    Counter object with the codon as key and AA as value.
    """
    AA_index,codon_index=1,3 # Default
    codon_table=collections.Counter()
    for line in codonfile:
        if line.startswith("#"):
            line=line.strip().split()
            line[0]=line[0][1:] # Strip leading "#"
            for i,x in enumerate(line):
                if x == "AA": AA_index = i
                if x == "codons": codon_index = i
            continue
        line=line.strip().split()
        for c in line[codon_index].split(","):
            codon_table[c]=line[AA_index]
    return codon_table


def parse_vcf_file(vcf_file):
    """
    This function will accept a vcf file and return the results in the form of a
    dict of dicts of the form:
                    CHROMOSOME_1
                       |----------|
                       |        POSITION_1
                       |          |    |-------(REF, ALT)
                       |          |
                       |        POSITION_2
                       |               |-------(REF, ALT)
                    CHROMOSOME_2

    """
    parsed_vcf = collections.Counter()
    for line in vcf_file:
        if line.startswith("#"):
            continue
        line=line.strip().split()
        if parsed_vcf[line[0]] == 0:
            parsed_vcf[line[0]] = collections.Counter()
        parsed_vcf[line[0]][line[1]] = (line[3], line[4])
    return vcf_file


def main(args):
    """
    This is the main function for the program.  It accepts numerous arguments
    from the user including
    in_ens  = input ensGene file to translate mutations in genomic space to 
              protein space.
    in_gen  = Input genomic fasta file to ascertain mutations from.
    in_prot = Input proteomeic fasta file for sanity checking the process.
    chrom   = A list of chromosomes to be processed by the script. This allows 
              for distributing the workload over multiple nodes on the cluster.
    codon   = A file containing the amino acid coding scheme of interest.  This
              has to be in a specific format.  Refer parse_codon_scheme().
    in_vcf  = An input vcf file containing the mutations for the sample.  The 
              chromosomes in the file must match those provided in --chroms.

    The main function parses the ensGene and proteomic fasta files into Counter
    objects to pass on to later modules. The input genomic fasta is processed
    one chromosome at a time to reduce the memory footprint of the program.
    """
    parser = argparse.ArgumentParser(description= main.__doc__ )
    parser.add_argument('--in_ens', dest='in_ens', type=file_type, 
                        help='Input ensGene file', 
                        default="/inside/depot4/users/arjun/tools/transgene/" +
                        "ensgene.txt" , required=False)
    parser.add_argument('--in_genome', dest='in_gen', type=file_type, 
                        help='Input genomic fasta file', 
                        default="/inside/depot4/users/arjun/tools/" +
                        "bwa_indexes/hg19_M_rCRS.fa", required=False)
    parser.add_argument('--in_proteome', dest='in_prot', type=file_type, 
                        help='Input proteomic fastsa file', 
                        default="/inside/depot4/users/arjun/tools/" + 
                        "transcriptome_ref/enspep.fa", required=False)
    parser.add_argument('--chrom', dest='chrom', type=str,
                        help='chromosome to process.', default="all")
    parser.add_argument('--codon', dest='codon', type=file_type,
                        help='Input codon file',
                        default="/inside/depot4/users/arjun/tools/transgene/" +
                        "codon_list" , required=False)
    parser.add_argument('--in_vcf', dest='in_vcf', type=file_type, 
                        help='Input vcf file.  Currently only SNVs.', 
                        default="None", required=True)
    params=parser.parse_args()
    all_chroms=["".join(["chr",str(x)]) for x in range(1,23)+["X","Y","M_rCRS"]]
    if params.chrom=="all":
        chrom_in=all_chroms
    else:
        chrom_in=params.chrom.strip().split(",")
    # Read in the codon scheme
    print("PROGRESS: Reading in codon scheme...", end="",
          file=sys.stderr)
    codon_scheme=parse_codon_scheme(params.codon)
    params.codon.close()
    print("\tDONE", file=sys.stderr)
    # Read in the ensGene file.
    print("PROGRESS: Reading in ensGene file to memory...", end="", 
          file=sys.stderr)
    parsed_file=parse_ensGene_file(params.in_ens)
    params.in_ens.close()
    print("\tDONE", file=sys.stderr)
    # Read in the proteome for sanity checking results.
    print("PROGRESS: Reading in proteome file to memory...", end="", 
          file=sys.stderr)
    prot_sequences = SeqIO.parse(params.in_prot,'fasta',
            alphabet=generic_protein)
    for protein in prot_sequences:
        prot[fasta.id]=list(fasta.seq.tostring())
    params.in_prot.close()
    print("\tDONE", file=sys.stderr)
    # Read in the mutations to memory
    print("PROGRESS: Reading in the mutations...", end="",
          file=sys.stderr)
    mutations=parse_vcf_file(params.in_vcf)
    params.in_vcf.close()
    assert set(mutations.keys()).issubset(chrom_in), "VCF contains mutations" +\
         " in chromosomes not provided in --chroms. The current chromosomes" + \
         " are %s" ", ".join(chrom_in)
    print("\tDONE", file=sys.stderr)
    # Read in the genomic fasta one sequence at a time and process on the fly.
    print("PROGRESS:Processing fasta file now", file=sys.stderr) 
    fasta_sequences = SeqIO.parse(params.in_gen,'fasta', 
            alphabet=generic_dna)
    for fasta in fasta_sequences:
        if name not in chrom_in:
            continue
        name,seq=fasta.id,list(fasta.seq.tostring())
        print("PROGRESS: Loaded ",name," to memory. Processing...",end="",
          file=sys.stderr)
        # Processing happens here
        print("\tDONE", file=sys.stderr)
    params.in_gen.close()


if __name__ == "__main__" :
    sys.exit(main(sys.argv))
