# -*- coding: utf-8 -*-
"""
Created on Sat Jun  1 14:19:15 2019

@author: Maria Pessoa Monteiro
"""
from find_gene import find_gene
from batch_find_gene import batch_find_gene
from find_homologues import find_homologues
from find_ipg import find_ipg
from choose_gene import choose_gene

def test_find_gene():
    """
    test_early_find_gene should return the matching DNA accession
    number (M22259.1) of protein AAA34866.
    """
    expect_gene = find_gene("AAA34866")
    assert expect_gene == ["DNA AN: M22259.1"]

def test_batch_find_gene():
    """
    test_batch_find_gene should return both gene's corresponding DNA accession
    numbers.
    """
    expect_batch = batch_find_gene(["AAA34866", "GAX67478"])
    assert expect_batch == ['M22259.1', 'BEGW01000001.1']

def test_find_ipg():
    """
    test_find_ipg should return the identical proteins of both proteins
    provided. GenBank may update, possibly leading to different identical
    protein results.
    """
    expect_ipg = find_ipg(["AAA34866", "GAX67478"], ipgfilename="testipg.csv")
    assert expect_ipg == {'M22259.1': ('387', '2066', 1),
                          'NC_001136.10': ('270222', '271901', 2),
                          'NM_001180165.1': ('1', '1680', 1),
                          'X05062.1': ('391', '2070', 1),
                          'X95644.1': ('8447', '10126', 2),
                          'Z74154.1': ('246', '1925', 2),
                          'BK006938.2': ('270222', '271901', 2),
                          'CM004297.1': ('288824', '290503', 2),
                          'LBMA01000004.1': ('288824', '290503', 2),
                          'HB870545.1': ('387', '2066', 1),
                          'HC927954.1': ('387', '2066', 1),
                          'BEGW01000001.1': ('1223469', '1225145', 1),
                          'CP004724.2': ('264356', '266032', 2),
                          'CP004740.2': ('253800', '255476', 2),
                          'CP004750.2': ('264165', '265841', 2),
                          'DG000040.1': ('258839', '260515', 2)}

def test_find_homologues():
    """
    test_find_homologues should return only the 3 results, as specified.
    GenBank may update, possibly leading to different BLAST results.
    """
    expect_blast1 = find_homologues("AAA34866",max_number=3,
                                    filename="testblast.xml")
    assert expect_blast1 == ['NP_010177', 'AJU66839', 'AJU64042']

def test_choose_gene():
    """
    test_choose_gene should return the DNA accession numbers, codon adaptation
    indexes, start positions and stop positions of the corresponding first two
    proteins in testipg.csv with available DNA accession numbers.
    """
    expect_CAI = choose_gene("sharp_yeast.txt",ipgfile="testipg.csv")
    assert expect_CAI == {'AAA34866.1': ('M22259.1', 0.148, '387', '2066'),
                         'NP_010177.1': ('NM_001180165.1', 0.148, '1', '1680')}
