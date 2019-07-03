# -*- coding: utf-8 -*-
"""
Created on Fri May 24 11:01:40 2019

@author: Maria Pessoa Monteiro
"""
from CAI import CAI
from Bio import Entrez, SeqIO
from time import sleep
import csv

def choose_gene(CAI_reference_table, ipgfile="ipg.csv"):
    """
    Choose_gene accepts two arguments, the name of the file containing RSCU
    values and the name of a csv file containing the result of an identical
    protein EFetch.
    The function returns a dictionary with the protein ANs as
    keys and a tuple containing the corresponding DNA ANs, CAI, start
    positions and stop positions as values.
    The format expected for the RSCU file is a text file with a header and
    two columns, the first column with codons and the second with RSCU values,
    delimited by a single space.
    """
    with open(CAI_reference_table) as file:
        reader = csv.reader(file, delimiter=" ")
        rows = list(reader)
        rows = rows[1:]

        tri = []
        RSCU = []
        for col in rows:
            tri.append(col[0])
            RSCU.append(float(col[1]))
        RSCU_sharp = dict(zip(tri, RSCU))

    with open(ipgfile) as ipgfile:
        reader = csv.reader(ipgfile)
        ipgrows = list(reader)[1:]

    dnaACC = []
    protACC = []
    strand = []
    start = []
    stop = []
    for col in ipgrows:
        if len(col[2]) > 3:
            dnaACC.append(col[2])
            protACC.append(col[6])
            start.append(col[3])
            stop.append(col[4])
            if col[5] == "+":
                strand.append(1)
            else:
                strand.append(2)
    all_dna = []
    i = 0
    while i < len(dnaACC):
        dna_fetch = Entrez.efetch(db="nucleotide", id=dnaACC[i],
                                  rettype='gbwithparts', retmode='text',
                                  seq_start=start[i], seq_stop=stop[i],
                                  strand=strand[i])
        sleep(1)
        seq_obj = SeqIO.read(dna_fetch, "gb")
        all_dna.append(str(seq_obj.seq))
        i += 1

    CAI_list = []
    final_dnaACC = []
    final_protACC = []
    final_start = []
    final_stop = []
    i = 0
    while i < len(all_dna):
        if len(all_dna[i]) % 3 == 0:
            CAI_list.append(round(CAI(all_dna[i], RSCUs=RSCU_sharp), 3))
            final_dnaACC.append(dnaACC[i])
            final_protACC.append(protACC[i])
            final_start.append(start[i])
            final_stop.append(stop[i])
        else:
            print(f"Not a multiple of 3: {dnaACC[i]}")
        i += 1
    values = zip(final_dnaACC, CAI_list, final_start, final_stop)
    CAI_dict = dict(zip(final_protACC, values))
    return CAI_dict
