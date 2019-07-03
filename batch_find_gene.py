# -*- coding: utf-8 -*-
"""
Created on Mon May 13 18:07:11 2019

@author: Maria Pessoa Monteiro
"""
from Bio import Entrez, SeqIO

def batch_find_gene(protACC_list):
    """
    Batch_find_gene takes a protein accession number list as the only
    argument and returns the corresponding DNA accession number in a list.
    """
    protACC = " ".join(protACC_list)
    dna_acc = []
    prot_fetch = Entrez.efetch(db="protein", id=protACC, rettype="gb",
                               retmode="text")
    prot_fetch_record = list(SeqIO.parse(prot_fetch, "gb"))

    i = 0
    while i < len(prot_fetch_record):
        for feat in prot_fetch_record[i].features:
            if feat.type == "CDS" and 'coded_by' in feat.qualifiers.keys():
                gene_acc = str(feat.qualifiers['coded_by']).split(":")
                if "comp" not in gene_acc:
                    gene_acc[0] = gene_acc[0].strip("[").strip("'")
                else:
                    gene_acc[0] = gene_acc[0].strip("[").strip("]")[11:]
                if gene_acc not in dna_acc:
                    dna_acc.append(gene_acc[0])
        i += 1
    return dna_acc
