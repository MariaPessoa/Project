# -*- coding: utf-8 -*-
"""
Created on Wed May  8 17:19:43 2019

@author: Maria Pessoa Monteiro
"""
from Bio import SeqIO, Entrez

def find_gene(prot_acc):
    """
    Find_gene accepts a protein AN as an argument. The function
    returns the AN of the matching gene.
    """
    results = []
    prot_fetch = Entrez.efetch(db="protein", id=prot_acc, rettype="gb",
                               retmode="text")
    record = SeqIO.read(prot_fetch, "gb")

    for feat in record.features:
        if feat.type == "CDS" and 'coded_by' in feat.qualifiers.keys():
            gene_acc = str(feat.qualifiers['coded_by']).split(":")
            gene_acc[0] = gene_acc[0].strip("[").strip("'")
            results.append(f"DNA AN: {gene_acc[0]}")

    return results
