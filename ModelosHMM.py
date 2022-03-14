#!/usr/bin/env python
# coding: utf-8

# # TP6 Grupo Nº 5
# Perez Ernesto Rafael rafaelperezctes@gmail.com Sofia Erdozain sofierdozain@gmail.com
# Created on Wed Nov  3 14:02:59 2021

# @author: Grupo Nº 5

import numpy as np
import Bio.SearchIO as bpio
import CargaBD
import Main as principal
import Utils
import sys

# archivos_agrupados = Alineamientos.archivos_agrupados
archivos_MSA = []

############################# 3: Análisis HMM
"""
#MAFFT

#!/bin/bash
echo “----------------PIPELINE biog5--------------”

echo Obeniendo semillas para cada E1-alto riesgo
mafft high_riskE1.fasta > high_riskE1_seed_msa.fasta 
echo Obeniendo modelo para cada E1-alto riesgo
hmmbuild high_riskE1_modelo.hmm high_riskE1_seed_msa.fasta
echo Obeniendo secuencias de E1-bajo riesgo que concidan con cada modelo E1-alto riesgo
hmmsearch high_riskE1_modelo.hmm low_riskE1.fasta > high_riskE1_vs_low_riskE1_.out
echo Obeniendo secuencias de E1-no determinado riesgo que concidan con cada modelo E1-alto riesgo
hmmsearch high_riskE1_modelo.hmm unspecified_riskE1.fasta > high_riskE1_vs_unspecified_riskE1_.out


echo Obeniendo semillas para cada E2-alto riesgo
mafft high_riskE2.fasta > high_riskE2_seed_msa.fasta
echo Obeniendo modelo para cada E2-alto riesgo
hmmbuild high_riskE2_modelo.hmm high_riskE2_seed_msa.fasta
echo Obeniendo secuencias de E2-bajo riesgo que concidan con cada modelo E2-alto riesgo
hmmsearch high_riskE2_modelo.hmm low_riskE2.fasta > high_riskE2_vs_low_riskE2_.out
echo Obeniendo secuencias de E2-no determinado riesgo que concidan con cada modelo E2-alto riesgo
hmmsearch high_riskE2_modelo.hmm unspecified_riskE2.fasta > high_riskE2_vs_unspecified_riskE2_.out

echo Obeniendo semillas para cada E7-alto riesgo
mafft high_riskE7.fasta > high_riskE7_seed_msa.fasta
echo Obeniendo modelo para cada E7-alto riesgo
hmmbuild high_riskE7_modelo.hmm high_riskE7_seed_msa.fasta
echo Obeniendo secuencias de E7-bajo riesgo que concidan con cada modelo E7-alto riesgo
hmmsearch high_riskE7_modelo.hmm low_riskE7.fasta > high_riskE7_vs_low_riskE7_.out
echo Obeniendo secuencias de E7-no determinado riesgo que concidan con cada modelo E7-alto riesgo
hmmsearch high_riskE7_modelo.hmm unspecified_riskE7.fasta > high_riskE7_vs_unspecified_riskE7_.out

echo Obeniendo semillas para cada L1-alto riesgo
mafft high_riskL1.fasta > high_riskL1_seed_msa.fasta
echo Obeniendo modelo para cada L1-alto riesgo
hmmbuild high_riskL1_modelo.hmm high_riskL1_seed_msa.fasta
echo Obeniendo secuencias de L1-bajo riesgo que concidan con cada modelo L1-alto riesgo
hmmsearch high_riskL1_modelo.hmm low_riskL1.fasta > high_riskL1_vs_low_riskL1_.out
echo Obeniendo secuencias de L1-no determinado riesgo que concidan con cada modelo L1-alto riesgo
hmmsearch high_riskL1_modelo.hmm unspecified_riskL1.fasta > high_riskL1_vs_unspecified_riskL1_.out

echo Obeniendo semillas para cada L2-alto riesgo
mafft high_riskL2.fasta > high_riskL2_seed_msa.fasta
echo Obeniendo modelo para cada L2-alto riesgo
hmmbuild high_riskL2_modelo.hmm high_riskL2_seed_msa.fasta
echo Obeniendo secuencias de L2-bajo riesgo que concidan con cada modelo L2-alto riesgo
hmmsearch high_riskL2_modelo.hmm low_riskL2.fasta > high_riskL2_vs_low_riskL2_.out
echo Obeniendo secuencias de L2-no determinado riesgo que concidan con cada modelo L2-alto riesgo
hmmsearch high_riskL2_modelo.hmm unspecified_riskL2.fasta > high_riskL2_vs_unspecified_riskL2_.out


"""
archivos_resultados = ['high_riskE1_vs_low_riskE1_.out', 'high_riskE1_vs_unspecified_riskE1_.out',
                       'high_riskE2_vs_low_riskE2_.out', 'high_riskE2_vs_unspecified_riskE2_.out',
                       'high_riskE7_vs_low_riskE7_.out', 'high_riskE7_vs_unspecified_riskE7_.out',
                       'high_riskL1_vs_low_riskL1_.out', 'high_riskL1_vs_unspecified_riskL1_.out',
                       'high_riskL2_vs_low_riskL2_.out', 'high_riskL2_vs_unspecified_riskL2_.out']

archivo0 = resultados = list(bpio.read("BD/HMM/" + archivos_resultados[0], 'hmmer3-text'))


def Leer(archivo):
    resultados = list(bpio.read("BD/HMM/" + archivo, 'hmmer3-text'))
    return (resultados)


comparacionesHMM = {}


def Recorrer(resultados):
    c = 0
    lista_aux = []
    for secuencias in resultados:
        for dominios in secuencias:
            cepa2 = "HPV" + (str(dominios.hit_description)[-4:-1]).replace(" ", "")
            lista_aux.append([cepa2, int(dominios.bitscore)])
            c += 1
    return lista_aux


def Comparar():
    for archivo in archivos_resultados:
        resultados = Leer(archivo)
        lista = Recorrer(resultados)
        clave = archivo.replace("_", " ").replace(".out", "").replace("risk", "")
        leyenda = "Clasificacion: " + clave[:4] + " Proteina:" + clave[4:7] + clave[7:11] + "Clasificacion: " + clave[
                                                                                                                11:-3] + " Proteina: " + clave[
                                                                                                                                         -3:]
        print(leyenda)
        comparacionesHMM[clave] = lista
    print()


def AnalizarA():
    print("2: Relaciones evolutivas con mejores scores encontradas por cada modelo:")
    for objetivo, k in comparacionesHMM.items():
        if len(k) > 0:
            matriz_aux = np.array(k)
            eje_x = matriz_aux[:, 0]
            eje_y = list(map(int, matriz_aux[:, 1]))
            print("Modelo HMM de Alto Riesgo proteina ", objetivo[5:8])
            print("|______ Mejor score - ", "Cepa:", eje_x[0], "Riesgo:", objetivo[11:-3], "Proteina:", objetivo[-3:])
            print()

# Mutaciones

# Paso 1 obtener las secuencias alineadas en grupos por proteinas E1-alto riesgo
import pandas as pd
from Bio import AlignIO

align = AlignIO.read("BD/HMM/high_riskE1_seed_msa.fasta", 'fasta')

# Paso 2 pasarlo a tablas
align_pd = pd.DataFrame(align)
filas = align_pd.shape[0]

# Paso 3 calculo de porcentajes
porcentajes = {}
for columna in align_pd:
    # obtengo todos los elementos de una columna
    lista_elementos = align_pd[columna].value_counts().index
    # print(lista_elementos)
    cant_hist = align_pd[columna].value_counts()
    valores = [list(cant_hist.index), list((cant_hist.values * 100) / cant_hist.sum())]
    porcentajes[columna] = valores


def Main():
    print("### Recuerde que debe tener los resultados del pipeline big5 en la carpeta /BD/HMM/ ###")
    print()
    print("1: Se compararan: Proteinas de cepas de alto riesgo vs proteinas de otros riesgos, espere por favor ...")
    Comparar()
    AnalizarA()
    Menu()
# Main()

def Menu():
    print()
    print("### Modulo: Alineamiento multiples ###")
    print()
    print("Que desea realizar:")
    print("1- Listar Proteinas almacenadas")
    print("2- Listar Genomas almacenadas")
    print("3- Realizar las comparaciones de modelos HMM")
    print("4- Volver al menu principal")
    print("5- Salir")
    CargaBD.LeerInicio()
    opcion_principal = int(input("Ingrese una opción: "))
    if opcion_principal == 1:
        Utils.ListarProteinas()
        Menu()
    if opcion_principal == 2:
        Utils.ListarGenomas()
        Menu()
    if opcion_principal == 3:
        Main()
    if opcion_principal == 4:
        principal.Menu()
    if opcion_principal == 5:
        print("Gracias por utilizar BIOG5")
        sys.exit()
# Menu()
