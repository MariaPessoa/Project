# -*- coding: utf-8 -*-
"""
Created on Wed May  8 17:19:43 2019

@author: Maria Pessoa Monteiro
"""
from Bio import SeqIO, Entrez
from Bio.Seq import translate
from time import sleep

def early_find_gene(prot_acc,max_number=10):
    """
    early_find_gene accepts a protein AN and an optional argument, the maximum number
    of results (set to 10 by default), find_gene returns a list containing the
    protein description, the total number of results, a list with their
    corresponding ANs and, if found, the AN of the matching gene.
    """
    results = []
    prot_fetch = Entrez.efetch(db="protein", id=prot_acc, rettype="gb",
                               retmode="text")
    sleep(1)
    record = SeqIO.read(prot_fetch, "gb")
    results.append(f"Protein description: {record.description}")

    start_index = record.description.index("[")
    query = record.description[0:(start_index-1)] + "[All fields]"
    organism = record.description[(start_index+1):record.description.index("]")] + "[Organism]"
    search_term = query + " AND " + organism + ' AND (biomol_genomic[PROP] AND ("0"[SLEN] : "20000"[SLEN]))'

    for feat in record.features:
        if feat.type == "CDS" and 'coded_by' in feat.qualifiers.keys():
            gene_acc = str(feat.qualifiers['coded_by']).split(":")
            gene_acc[0] = gene_acc[0].strip("[").strip("'")
            results.append(f"DNA AN: {gene_acc[0]}")

    dna_search = Entrez.esearch(db="nucleotide", term=search_term,
                                idtype="acc", retmax=max_number)
    dna_search_record = Entrez.read(dna_search)
    sleep(1)

    if int(dna_search_record['Count']) != 0:
        results.append(f"Result Count: {dna_search_record['Count']}")
        results.append(f"Accession numbers: {dna_search_record['IdList']}")
    else:
        results.append("No DNA results")
        return results

    dna_acc = " ".join(dna_search_record["IdList"])
    dna_fetch = Entrez.efetch(db="nucleotide", id=dna_acc, rettype="gb",
                              retmode="text")
    dna_fetch_record = list(SeqIO.parse(dna_fetch, "gb"))
    print(len(dna_fetch_record[0].seq))

    for feat in dna_fetch_record[0].features:
        if feat.type == "CDS":
            cds_location = dna_fetch_record[0].seq[int(feat.location.start):int(feat.location.end)]
            if translate(cds_location, table=1, to_stop=True) == record.seq:
                results.append(f"Match found: {dna_fetch_record[0].id}")
    return results
