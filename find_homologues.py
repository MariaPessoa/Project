# -*- coding: utf-8 -*-
"""
Created on Mon May 13 17:08:57 2019

@author: Maria Pessoa Monteiro
"""
from Bio.Blast import NCBIWWW, NCBIXML

def find_homologues(protACC, max_number=10, filename="blast.xml"):
    """
    Find_homologues takes a protein accession number as required argument,
    and an optional max_number of results argument, default set to 10,
    and does a protein BLAST. The function returns the accession numbers
    of the BLAST proteins.
    """
    result_handle = NCBIWWW.qblast("blastp","nr", protACC,
                                   hitlist_size=max_number)
    with open(filename, "w") as out_handle:
        out_handle.write(result_handle.read())
    with open(filename) as file:
        blast_record = NCBIXML.read(file)
        protACC = []
        for rec in blast_record.alignments:
            protACC.append(rec.accession)
    return protACC
