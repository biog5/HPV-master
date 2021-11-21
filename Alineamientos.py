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
import CargaBD

#  https://www.papilocare.com/copia-de-caracteristicas
unspecified_risk = ["HPV53", "HPV66", "HPV68", "HPV73", "HPV82"]
low_risk = ["HPV6", "HPV11",  "HPV42", "HPV40","HPV53", "HPV54", "HPV57"]
high_risk = ["HPV16", "HPV18", "HPV31", "HPV33", "HPV35", "HPV39", "HPV45", "HPV51", "HPV52", "HPV56", "HPV58"]


############# 2: Alineamientos Blast
# Compararemos cada protenia (e1,e2,e7,l1,l2) de cada cepa de alto riesgo
# contra cada protenia (e1,e2,e7,l1,l2) de cepas de low and unspecified riesgo
# Para esto dividiremos en 3 partes: 
# 2.1-Agrupación  
# 2.2-blast
# 2.3-Analisis

archivos_temp=[]
archivos_agrupados=[]
archivos_MSA=[]
############################# 2.1: Agrupacion
# Leer
genes = ["E1","E2","E7","L1","L2"]
def LeerYGrabar(proteina, salida):
    archivo= 'BD/'+proteina+'.fasta'
    proteinas = list(bsio.parse(archivo, 'fasta'))
    bsio.write(proteinas, salida, 'fasta')

def BuscarNombres(genoma, gen):
    for nombre, k in CargaBD.proteinas.items():
        if (genoma in nombre) and (gen in nombre) :
            return nombre

def AgruparRiesgos():
    # Primero grabo las protenias ojetivo de alto riesgo
    for genoma in high_risk: # nombre del genoma
        for gen in genes:
            salida_nombre ='BD/'+"high_risk"+gen+'.fasta'
            archivo=BuscarNombres(str(genoma), str(gen))
            if archivo != None:
                archivo_salida = open(salida_nombre, "a+")                
                LeerYGrabar(archivo, archivo_salida)
                archivos_agrupados.append(salida_nombre)
                archivo_salida.close() 
    
    # Segundo grabo las protenias ojetivo de bajo riesgo
    for genoma in low_risk: # nombre del genoma
        for gen in genes:
            salida_nombre ='BD/'+"low_risk"+gen+'.fasta'
            archivo=BuscarNombres(str(genoma), str(gen))
            if archivo != None: 
                archivo_salida = open(salida_nombre, "a+")
                LeerYGrabar(archivo, archivo_salida)
                archivos_agrupados.append(salida_nombre)
                archivo_salida.close()  
    
    # Tercero grabo las protenias ojetivo de no especificado riesgo
    for genoma in unspecified_risk: # nombre del genoma
        for gen in genes:
            salida_nombre ='BD/'+"unspecified_risk"+gen+'.fasta'
            archivo=BuscarNombres(str(genoma), str(gen))
            if archivo != None:                
                archivo_salida = open(salida_nombre, "a+")
                LeerYGrabar(archivo, archivo_salida)
                archivos_agrupados.append(salida_nombre)
                archivo_salida.close()  

#Revisar que los archivos no se esten regrabando xq darian error
AgruparRiesgos()

       

############################# 2.2: Blastp
import subprocess
from subprocess import PIPE, Popen
from Bio.Blast import NCBIWWW,NCBIXML
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

comparaciones={}
def Recorrer(resultados):
    c=0
    lista_aux=[]
    for query in resultados: 
        cepa_query = str(query.description)[-4:-1]
        for hit in query:
            cepa2 = "HPV"+(str(hit.description)[-4:-1]).replace(" ", "")
            for hsp in hit:
                print(cepa2, hsp.bitscore)
                lista_aux.append([cepa2, int(hsp.bitscore)])
                c+=1
    return cepa_query, lista_aux

def AlinearRiesgos():                
    for gen in genes:
        archivo_low = "BD/low_risk"+gen+'.fasta'
        for genoma in high_risk: # por cada genoma y gen de alto rieago
            archivo_entrada = BuscarNombres(str(genoma), str(gen))
            if archivo_entrada != None:  
                archivo_high ='BD/'+archivo_entrada+'.fasta'
                archivo_salida = "BD/"+genoma+gen+"_high_vs_"+gen+"_low.blast"
                result = subprocess.run("blastp -query "+ archivo_high+" -subject "+ archivo_low + " -out "+archivo_salida, stdout=PIPE)
                print("Corriendo BLASTP contra lista de bajo riesgo")
                print("Cepa: ",genoma," Gen: ",gen, " Salida:", archivo_salida)
                archivos_temp.append(archivo_salida)
                resultados=list(bpio.parse(archivo_salida, 'blast-text'))
                cepa_query, lista_aux= Recorrer(resultados)
                comparaciones[archivo_entrada]=lista_aux
            
AlinearRiesgos()       

def Limpiar(archivos):
    from os import remove
    from os import path
    for i in archivos:
        if path.exists(i):
            print("Eliminando: ",i)
            remove(i)           
Limpiar(archivos_temp)
#Limpiar(archivos_agrupados) # Lo dejo xq lo necesita MSA


############################# 2.3: Análisis

def GraficarA(objetivo,x,y,proteina, cepa):
    plt.figure(figsize = (12, 6))
    print("Graficando y guardando analisis", " Cepa: ",cepa," Proteina: ",proteina)
    plt.bar(x,y)
    plt.xlabel("Cepas de bajo riesgo")
    plt.ylabel("BitScore")
    titulo="Análisis del virus del Papiloma Humano - "+"Cepa: "+cepa+" Proteina: "+proteina
    titulo2="Análisis del virus del Papiloma Humano - "+"Cepa "+cepa+" Proteina "+proteina
    plt.title(titulo)
    plt.savefig("IMG/"+titulo2+".png")


def AnalizarB(E1,E2,E7,L1,L2):
    for i in ['E1','E2','E7','L1','L2']:
        print("Arbol Filogentico, nivel 1, para la proteina",i)
        print("Cepa Bajo Riesgo ______ Cepa Alto Riesgo")
        for bajo,altos in eval(i).items(): 
            print(bajo)
            for j in altos:
                print("|______",j)
        print()
    

def AnalizarA():
    E1 = {}
    E2 = {}
    E7 = {}
    L1 = {}
    L2 = {}
    for objetivo,k in comparaciones.items():       
        if len(k)>0:
            matriz_aux = np.array(k)
            eje_x = matriz_aux[:,0]
            eje_y = list(map(int, matriz_aux[:,1]))
            proteina=str(objetivo[-2:]) 
            cepa = str(objetivo[:-2])
            valor = eval(proteina).get(eje_x[0])
            if valor != None:
                valor.append(cepa)
            else:
                valor = [cepa]
            eval(proteina)[eje_x[0]] = valor
            GraficarA(objetivo,eje_x,eje_y, proteina, cepa)
    AnalizarB(E1,E2,E7,L1,L2)
    
AnalizarA()
