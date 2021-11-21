#!/usr/bin/env python
# coding: utf-8

# # TP6 Grupo Nº 5
# Perez Ernesto Rafael rafaelperezctes@gmail.com Sofia Erdozain sofierdozain@gmail.com
#Created on Wed Nov  3 14:02:59 2021

#@author: Grupo Nº 5


import Bio.SeqIO as bsio
import matplotlib.pyplot as plt
import numpy as np
import Alineamientos
archivos_agrupados = Alineamientos.archivos_agrupados
archivos_MSA = []

############################# 3: Alineamientos MSA
"""
Alineo grupos, usando ClustalO Windows
alto riesgo E1
alto riesgo E2
alto riesgo E7
alto riesgo L1
alto riesgo L2
bajo riesgo E1
bajo riesgo E2
bajo riesgo E7
bajo riesgo L1
bajo riesgo L2
no_determinado riesgo E1
no_determinado riesgo E2
no_determinado riesgo E7
no_determinado riesgo L1
no_determinado riesgo L2

Obtengo_arboles por cada uno usando lubreria Phylo 

"""
archivo = 'BD/high_riskE1.fasta '


#http://www.clustal.org/omega/


from Bio import AlignIO
from Bio.Align.Applications import ClustalOmegaCommandline 
import subprocess
from subprocess import PIPE, Popen

def CorrerClustalOmega():
    for archivo in archivos_agrupados:
        archivo_entrada = archivo
        archivo_salida = archivo.replace(".fasta", "_MSA.phylip")
        
        # Si se tiene instalado descomentar y correr:
        #clustalomega_cline = ClustalOmegaCommandline(infile = archivo_entrada, outfile = archivo_salida, outfmt = 'phylip', verbose = True, auto = False)
        #print(clustalomega_cline)
        
        # Por practicidad lo corremos portable en windows
        clustal='"clustal-omega-1.2.2-win64/clustalo.exe"' + ' -i '+ archivo_entrada + ' -o '+ archivo_salida +' --outfmt phylip -v --force'
        result = subprocess.run(clustal, stdout=PIPE)
        archivos_MSA.append(archivo_salida)
CorrerClustalOmega()

def LeerMSA():
    for archivo in archivos_MSA:
        alignments = list(AlignIO.parse(archivo, "phylip"))
        for alignment in alignments:
            print(alignment)
            print()
    
LeerMSA()

#https://biopython.org/wiki/Phylo    
from Bio import Phylo
from Bio.Phylo.TreeConstruction import DistanceCalculator
from Bio import AlignIO
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor


def ContruirArboles():
    for archivo in archivos_MSA:  
        print("Creando arbol para:",archivo)
        # construyo arbol Filogenetico
        # paso 1 leo arcchivos
        aln = AlignIO.read(archivo, 'phylip')
        #print(aln)
        # calculo de distancia con matriz default(dna_models, protein_models, models)
        calculator = DistanceCalculator('identity')
        dm = calculator.get_distance(aln)
        #print(dm)
        # No se bien
        constructor = DistanceTreeConstructor(calculator, 'nj')
        tree = constructor.build_tree(aln)
        tree = constructor.upgma(dm)
        tree.ladderize()  # Imprime imagen
        fig=Phylo.draw(tree)
        #print(tree) # Imprime imagen texto
        #Phylo.draw_ascii(tree) # Imprime por consola arbol

ContruirArboles()

Alineamientos.Limpiar(archivos_agrupados)