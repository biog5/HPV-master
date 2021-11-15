#!/usr/bin/env python
# coding: utf-8

# # TP6 Grupo Nº 5
# Perez Ernesto Rafael rafaelperezctes@gmail.com Sofia Erdozain sofierdozain@gmail.com
#Created on Wed Nov  3 14:02:59 2021

#@author: Grupo Nº 5


import Bio.SeqIO as bsio
import matplotlib.pyplot as plt
import numpy as np

############################# 1: Descargas
# llamo al modulo de carga de BD
import TPEspecial

unspecified_risk = ["HPV30", "HPV34", "HPV39", "HPV40", "HPV53", "HPV57", "HPV59", "HPV61", "HPV62", "HPV64", "HPV66", "HPV67", "HPV68", "HPV69"]
low_risk = ["HPV6", "HPV11", "HPV16", "HPV18", "HPV31", "HPV33", "HPV35", "HPV42", "HPV43", "HPV44", "HPV45", "HPV51", "HPV52", "HPV74"]
high_risk = ["HPV16", "HPV18", "HPV6", "HPV11", "HPV31", "HPV34", "HPV33", "HPV35", "HPV39", "HPV42", "HPV44", "HPV45", "HPV51", "HPV52", "HPV56", "HPV58", "HPV66"]


# Compararemos cada protenia (e1,e2,e7,l1,l2) de cada cepa de alto riesgo
# contra cada protenia (e1,e2,e7,l1,l2) de cepas de low and unspecified riesgo
# Para esto dividiremos en 3 partes: 
# 1-Agrupación  
# 2-blast
# 3-Analisis


############################# 2: Agrupacion
# Leer
genes = ["E1","E2","E7","L1","L2"]
def Leer(proteina, salida):
    archivo= 'BD/'+proteina+'.fasta'
    proteinas = list(bsio.parse(archivo, 'fasta'))
    bsio.write(proteinas, salida, 'fasta')

def BuscarNombres(genoma, gen):
    for nombre, k in TPEspecial.proteinas.items():
        if (genoma in nombre) and (gen in nombre) :
            return nombre

def Recorrer():
    # Primero grabo las protenias ojetivo de alto riesgo
    for genoma in high_risk: # nombre del genoma
        for gen in genes:
            salida_nombre ='BD/'+"high_risk"+gen+'.fasta'
            output_file = open(salida_nombre, "a")
            archivo=BuscarNombres(str(genoma), str(gen))
            if archivo != None:                
                Leer(archivo, output_file)
                output_file.close() 
    
    # Segundo grabo las protenias ojetivo de bajo riesgo
    for genoma in low_risk: # nombre del genoma
        for gen in genes:
            salida_nombre ='BD/'+"low_risk"+gen+'.fasta'
            output_file = open(salida_nombre, "a")
            archivo=BuscarNombres(str(genoma), str(gen))
            if archivo != None:                
                Leer(archivo, output_file)
                output_file.close()  
    
    # Tercero grabo las protenias ojetivo de no especificado riesgo
    for genoma in unspecified_risk: # nombre del genoma
        for gen in genes:
            salida_nombre ='BD/'+"unspecified_risk"+gen+'.fasta'
            output_file = open(salida_nombre, "a")
            archivo=BuscarNombres(str(genoma), str(gen))
            if archivo != None:                
                Leer(archivo, output_file)
                output_file.close()  
Recorrer()


############################# 3: Blasp



"""
Ejemeplo E1 
# blastp -query BD/high_riskE1.fasta -subject BD/low_riskE1.fasta -out E1_high_vs_E1_low.blast
# blastp -query BD/high_riskE2.fasta -subject BD/low_riskE2.fasta -out E1_high_vs_E2_low.blast
# blastp -query BD/high_riskE7.fasta -subject BD/low_riskE2.fasta -out E1_high_vs_E2_low.blast
# blastp -query BD/high_riskL1.fasta -subject BD/low_riskL1.fasta -out E1_high_vs_L1_low.blast
# blastp -query BD/high_riskL2.fasta -subject BD/low_riskL2.fasta -out E1_high_vs_L2_low.blast


"""

