# -*- coding: utf-8 -*-
"""
Created on Sat Jun 22 14:07:14 2019

@author: Maria Pessoa Monteiro
"""
from time import sleep
from find_homologues import find_homologues
from find_ipg import find_ipg
from choose_gene import choose_gene

protACC = find_homologues("NP_036345.2", filename="chain_blast.xml")
sleep(1)
find_ipg(protACC, ipgfilename="chain_ipg.csv")
sleep(1)
print(choose_gene("sharp_ecoli.txt", ipgfile="chain_ipg.csv"))
