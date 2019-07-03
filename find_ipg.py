# -*- coding: utf-8 -*-
"""
Created on Fri May 24 11:01:40 2019

@author: Maria Pessoa Monteiro
"""
import csv
from Bio import Entrez

def find_ipg(protACC_list, ipgfilename="ipg.csv"):
    """
    Find_ipg accepts a protein accession number list as an argument and uses
    EFetch to retrieve identical proteins. A csv file, by default named
    ipg.csv, is saved with the results.
    The function returns a dictionary with the keys being the corresponding
    DNA ANs and the values a tuple with the start position, stop position
    and coding strand (plus = 2, minus = 1), in that order.
    """
    protACC = " ".join(protACC_list)
    ipg = Entrez.efetch(db="protein", id=protACC, rettype='ipg',
                        retmode='text')

    fetch_reader = csv.reader(ipg, delimiter="\t")
    with open(ipgfilename, "w", newline='') as out_handle:
        writer = csv.writer(out_handle)
        writer.writerows(fetch_reader)

    with open(ipgfilename) as ipgfile:
        reader = csv.reader(ipgfile)
        rows = list(reader)
        rows = rows[1:]

    dna = []
    strand = []
    start = []
    stop = []

    for col in rows:
        if len(col[2]) > 3:
            dna.append(col[2])
            start.append(col[3])
            stop.append(col[4])
            if col[5] == "+":
                strand.append(1)
            else:
                strand.append(2)
    values = zip(start, stop, strand)

    return dict(zip(dna, values))
