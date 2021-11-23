#!/usr/bin/env python
# coding: utf-8

# # TP6 Grupo Nº 5
# Perez Ernesto Rafael rafaelperezctes@gmail.com Sofia Erdozain sofierdozain@gmail.com
#Created on Wed Nov  3 14:02:59 2021

#@author: Grupo Nº 5


import Bio.SeqIO as bsio
import matplotlib.pyplot as plt
import numpy as np


print("Modulo princial")

############################# 1: Descargas
# Modulo de carga de BD
print("Carga de BD")
import CargaBD

############################# 2: Alineamientos
print("Alineamientos")
# Modulo Alineamientos usando blast
import Alineamientos

############################# 3: Alineamientos MSA
print("Alineamientos")
# Modulo Alineamientos usando clustalOmega y libreria de creacion de arboles filogeneticos
import AlineamientosMSA

############################# 4: Modelos HMM
print("Modelos HMM")
# Modulo Alineamientos usando clustalOmega y libreria de creacion de arboles filogeneticos
import ModelosHMM