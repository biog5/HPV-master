#!/usr/bin/env python
# coding: utf-8

# # TP6 Grupo Nº 5
# Perez Ernesto Rafael rafaelperezctes@gmail.com Sofia Erdozain sofierdozain@gmail.com
# Created on Wed Nov  3 14:02:59 2021

# @author: Grupo Nº 5


import Bio.SeqIO as bsio
import matplotlib.pyplot as plt
import numpy as np
import Utils
import Main as principal
import sys
import CargaBD
#  https://www.papilocare.com/copia-de-caracteristicas

def LeerBD():
    if not (('unspecified_risk' and "low_risk" and "high_risk") in globals()):
        global unspecified_risk
        global low_risk
        global high_risk
        CargaBD.LeerInicio()
        unspecified_risk = CargaBD.unspecified
        low_risk = CargaBD.low_risk 
        high_risk = CargaBD.high_risk


LeerBD()

############# 2: Alineamientos Blast
# Compararemos cada protenia (e1,e2,e7,l1,l2) de cada cepa de alto riesgo
# contra cada protenia (e1,e2,e7,l1,l2) de cepas de low and unspecified riesgo
# Para esto dividiremos en 3 partes: 
# 2.1-Agrupación  
# 2.2-blast
# 2.3-Analisis

archivos_temp = []
archivos_agrupados = CargaBD.archivos_agrupados
archivos_MSA = []

############################# 2.1: Agrupacion
# Leer
genes = ["E1", "E2", "E7", "L1", "L2"]

def LeerYGrabar(proteina, salida):
    archivo = 'BD/' + proteina + '.fasta'
    proteinas = list(bsio.parse(archivo, 'fasta'))
    bsio.write(proteinas, salida, 'fasta')


def BuscarNombres(genoma, gen):
    for nombre, k in CargaBD.proteinas.items():
        #tamanio1=len(nombre)
        #tamanio2=len(genoma+gen)+1
        intermedio=nombre[len(genoma):-(len(gen))]
        if (genoma in nombre) and (gen in nombre) :
            if len(intermedio)<=1 and (intermedio.isdigit()==False):
                return nombre

def Limpiar(archivos):
    from os import remove
    from os import path
    for i in archivos:
        if path.exists(i):
            # print("Eliminando: ",i)
            remove(i)

def AgruparRiesgos():
    # Primero grabo las protenias ojetivo de alto riesgo
    global archivos_agrupados
    archivos_agrupados = []  # Aunque lo tengo, lo uso por si actualizo el codigo y necesito mas clasificaciones
    for genoma in high_risk:  # nombre del genoma
        for gen in genes:
            salida_nombre = 'BD/' + "high_risk" + gen + '.fasta'
            archivo = BuscarNombres(str(genoma), str(gen))
            if archivo != None:
                archivo_salida = open(salida_nombre, "a+")
                LeerYGrabar(archivo, archivo_salida)
                if salida_nombre not in archivos_agrupados:
                    archivos_agrupados.append(salida_nombre)
                archivo_salida.close()

                # Segundo grabo las protenias ojetivo de bajo riesgo
    for genoma in low_risk:  # nombre del genoma
        for gen in genes:
            salida_nombre = 'BD/' + "low_risk" + gen + '.fasta'
            archivo = BuscarNombres(str(genoma), str(gen))
            if archivo != None:
                archivo_salida = open(salida_nombre, "a+")
                LeerYGrabar(archivo, archivo_salida)
                if salida_nombre not in archivos_agrupados:
                    archivos_agrupados.append(salida_nombre)
                archivo_salida.close()

                # Tercero grabo las protenias ojetivo de no especificado riesgo
    for genoma in unspecified_risk:  # nombre del genoma
        for gen in genes:
            salida_nombre = 'BD/' + "unspecified_risk" + gen + '.fasta'
            archivo = BuscarNombres(str(genoma), str(gen))
            if archivo != None:
                archivo_salida = open(salida_nombre, "a+")
                LeerYGrabar(archivo, archivo_salida)
                if salida_nombre not in archivos_agrupados:
                    archivos_agrupados.append(salida_nombre)
                archivo_salida.close()
            ############################# 2.2: Blastp


import subprocess
from subprocess import PIPE, Popen
from Bio.Blast import NCBIWWW, NCBIXML
from Bio import SeqIO
from io import StringIO
import Bio.SearchIO as bpio

"""
Intalar Blast Windows:
guia: https://www.ncbi.nlm.nih.gov/books/NBK52637/
release: https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.10.0/
configurar la variable de entorno:BLASTDB_LMDB_MAP_SIZE=1000000

Ejemplo parametros
-perc_identity 75 #identidad
-qcov_hsp_perc 55 #cobertura querry 
-outfmt "7  #formato
"""

comparaciones_full = [{}, {}]

def Recorrer(resultados):
    c = 0
    lista_aux = []
    for query in resultados:
        cepa_query = str(query.description)[-4:-1]
        for hit in query:
            cepa2 = "HPV" + (str(hit.description)[-4:-1]).replace(" ", "")
            for hsp in hit:
                #print(cepa2, hsp.bitscore)
                #if cepa2=="HPV68a":
                    #opcion_principal = int(input("Ingrese una opción: "))
                lista_aux.append([cepa2, int(hsp.bitscore)])
                c += 1
    return cepa_query, lista_aux


def AlinearRiesgos():
    for gen in genes:
        archivo_low = "BD/low_risk" + gen + '.fasta'
        archivo_unspecified = "BD/unspecified_risk" + gen + '.fasta'
        for genoma in high_risk:  # por cada genoma y gen de alto rieago
            archivo_entrada = BuscarNombres(str(genoma), str(gen))
            if archivo_entrada != None:
                archivo_high = 'BD/' + archivo_entrada + '.fasta'
                archivo_salida = "BD/" + genoma + gen + "_high_vs_" + gen + "_low.blast"
                archivo_salida2 = "BD/" + genoma + gen + "_high_vs_" + gen + "_unspecified.blast"
                result = subprocess.run(
                    "blastp -query " + archivo_high + " -subject " + archivo_low + " -out " + archivo_salida,
                    stdout=PIPE)
                print("Corriendo BLASTP ---> Cepa:", genoma, " Proteina:", gen, " Contra grupo: Bajo riego ")
                result2 = subprocess.run(
                    "blastp -query " + archivo_high + " -subject " + archivo_unspecified + " -out " + archivo_salida2,
                    stdout=PIPE)
                print("Corriendo BLASTP ---> Cepa:", genoma, " Proteina:", gen, " Contra grupo: No especificado riego ")
                archivos_temp.append(archivo_salida)
                archivos_temp.append(archivo_salida2)
                resultados = list(bpio.parse(archivo_salida, 'blast-text'))
                resultados2 = list(bpio.parse(archivo_salida2, 'blast-text'))
                cepa_query, lista_aux = Recorrer(resultados)
                cepa_query2, lista_aux2 = Recorrer(resultados2)
                comparaciones_full[0][archivo_entrada] = lista_aux
                comparaciones_full[1][archivo_entrada] = lista_aux2

    print("Calculos terminados se analizaran resultados: ")


############################# 2.3: Análisis

def GraficarA(objetivo, x, y, proteina, cepa):
    plt.figure(figsize=(12, 6))
    print("Graficando y guardando analisis", " Cepa: ", cepa, " Proteina: ", proteina)
    plt.bar(x, y)
    plt.xlabel("Cepas de bajo riesgo")
    plt.ylabel("BitScore")
    titulo = "Análisis del virus del Papiloma Humano - " + "Cepa: " + cepa + " Proteina: " + proteina
    titulo2 = "Análisis del virus del Papiloma Humano - " + "Cepa " + cepa + " Proteina " + proteina
    plt.title(titulo)
    archivo_grafica = "IMG/" + titulo2 + ".png"
    plt.savefig(archivo_grafica)


def AnalizarB(leyenda, E1, E2, E7, L1, L2):
    print()
    for i in ['E1', 'E2', 'E7', 'L1', 'L2']:
        print("Relaciones evolutivas con mejores scores encontradas para la proteina", i)
        print("Cepa", leyenda, "Riesgo ______ Cepa Alto Riesgo")
        for bajo, altos in eval(i).items():
            print(bajo)
            for j in altos:
                print("|______", j)
        print()


def AnalizarA():
    E1 = {}
    E2 = {}
    E7 = {}
    L1 = {}
    L2 = {}
    leyendas = ["Bajo", "No especificado"]
    for count, comparaciones in enumerate(comparaciones_full):
        for objetivo, k in comparaciones.items():
            if len(k) > 0:
                matriz_aux = np.array(k)
                eje_x = matriz_aux[:, 0]
                eje_y = list(map(int, matriz_aux[:, 1]))
                proteina = str(objetivo[-2:])
                cepa = str(objetivo[:-2])
                valor = eval(proteina).get(eje_x[0])
                if valor != None:
                    valor.append(cepa)
                else:
                    valor = [cepa]
                eval(proteina)[eje_x[0]] = valor
                # GraficarA(objetivo,eje_x,eje_y, proteina, cepa)
        AnalizarB(leyendas[count], E1, E2, E7, L1, L2)


# Limpiar(archivos_agrupados) # Solo una vez al iniciar y lo dejo xq lo necesita MSA
# Limpiar(archivos_temp)


def Main():
    print("1: Obteniendo datos y agrupandolos por riesgo, espere por favor ...")
    Limpiar(archivos_agrupados)  # Solo una vez al iniciar y lo dejo xq lo necesita MSA
    # Revisar que los archivos no se esten regrabando xq darian error
    AgruparRiesgos()
    print("2: Se alineran genes-proteinas contra grupos de riesgo, espere por favor ...")
    AlinearRiesgos()
    AnalizarA()
    Limpiar(archivos_agrupados)  # Solo una vez al iniciar y lo dejo xq lo necesita MSA
    Limpiar(archivos_temp)
    Menu()


def Menu():
    print()
    print("### Modulo: Alineamiento de a pares ###")
    print()
    print("Que desea realizar:")
    print("1- Listar Proteinas almacenadas")
    print("2- Listar Genomas almacenadas")
    print("3- Realizar alineamiento de a pares de proteinas con BLAST y analizar ")
    print("4- Volver al menu principal")
    print("5- Salir")

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
# Main()
# Menu()
