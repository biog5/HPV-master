# HPV-master

# Ernesto Rafael, Perez rafaelperezctes@gmail.com 
# Sofia, Erdozain sofierdozain@gmail.com

Colaboradores: Aldana Claps aldanaclaps@gmail.com Ana Luz Garcia garcianaluz@gmail.com
# Análisis de proteinas y genomas de cepas del virus del papiloma humano

El sistema Biog5 para análisis, incia como trabajo final de la materia Bioinformática Avanzada de la UBA. Estamos en la primera etapa del desarrollo, pero se planifica a futuro incluir, múltiples módulos, integrantes de diversas ramas de las ciencias en participación comunitaria para el desarrollo de un suite de herramientas que esperamos sirvan de ayuda para la Investigación y Desarrollo de tratamientos en patologías generaras por el virus del papiloma humano.

# Desarrollo etapas:

Main 

!["E0"](https://github.com/biog5/HPV-master/blob/main/IMG/FrontEnd/A-R0.png)

1- Alineamiento: 
<pre>
  1.1 Agrupamiento
  
  1.2 Se utiliza en primera instacia Blastp para alinear las proteinas de cepas de alto riesgo contra las demas cepas
  
  1.3. Análisis de resultados 
</pre>

Menu

!["E1-1"](https://github.com/biog5/HPV-master/blob/main/IMG/FrontEnd/A-R1-1.png)

Lista de proteinas 

!["E1-2"](https://github.com/biog5/HPV-master/blob/main/IMG/FrontEnd/A-R1-2.png)

Lista de genomas

!["E1-3"](https://github.com/biog5/HPV-master/blob/main/IMG/FrontEnd/A-R1-3.png)

1.1 Salida impresa

!["E1-4"](https://github.com/biog5/HPV-master/blob/main/IMG/FrontEnd/A-R1-4.png)

1.2 Salida menu

!["E1-5"](https://github.com/biog5/HPV-master/blob/main/IMG/FrontEnd/A-R1-5.png)

1.3 Salida gráfica

!["E1-6"](https://github.com/biog5/HPV-master/blob/main/IMG/FrontEnd/A-R1-6.png)

2- Alineamiento MSA
<pre>
  2.1 Se realiza un alineamiento multiple (MSA) con Clustal Omega de todas las proteinas aprupadas por riesgos 
  
  2.2 Se generan y visualizan árboles filogenticos asociados a cada proteina
</pre>

2.1 Interfaz 

!["E2-0"](https://github.com/biog5/HPV-master/blob/main/IMG/FrontEnd/A-R2-0.png)

2.2 Árboles Filogenticos

!["E2-1"](https://github.com/biog5/HPV-master/blob/main/IMG/FrontEnd/A-R2-1.png)

3- Modelado HMM: Se utilza MAFFT para aliniamiento multiple y HMMER para modelos HMM
<pre>
  3.1 Se utiliza el pipeline biog5 y se generan modelos de altos riesgo para poteinas E1, E2, E7, L1, L2
  
  3.2 Se comparan las cepas con cada modelo

  3.3 Se visualizan los resultados de forma impresa y gráfica
</pre>

3.2 Salida impresa

!["E3-0"](https://github.com/biog5/HPV-master/blob/main/IMG/FrontEnd/A-R3-0.png)

3.2 Salida gráfica

!["E3-1"](https://github.com/biog5/HPV-master/blob/main/IMG/FrontEnd/A-R3-1.png)

4- Variantes en proteinas
<pre>
  4.1 Se realiza un alineamiento multiple por proteina y riesgo
  
  4.2 Se filtran solo los aminoacidos con un porcentaje de conservación ingresado

  4.3 Se visualizan los resultados de forma impresa y gráfica
</pre>

Sub menu
 
!["E4-0"](https://github.com/biog5/HPV-master/blob/main/IMG/FrontEnd/A-R4-0.png)

Salida impresa

!["E4-1"](https://github.com/biog5/HPV-master/blob/main/IMG/FrontEnd/A-R4-1.png)

Salida gráfica

!["E4-2"](https://github.com/biog5/HPV-master/blob/main/IMG/FrontEnd/A-R4-2.png)

5- Variantes en genomas
<pre>
  4.1 Se realiza un alineamiento multiple de genomas completos dos o más
  
  4.2 Se filtran solo los aminoacidos con un porcentaje de conservación ingresado

  4.3 Se visualizan los resultados de forma impresa y gráfica
</pre>

!["E5"](https://github.com/biog5/HPV-master/blob/main/IMG/FrontEnd/A-R5.png)

6- Estructuras NCBI
<pre>
  5.1 Se buscan y alinean secuencias de proteinas contra estructuras del NCBI (pdb)
  
  5.2 Se descargan las de interes selecionadas

  5.3 Se visualizan en un navegador
</pre>

5.1 Se busca y alinea una proteina de interes

!["E6-0"](https://github.com/biog5/HPV-master/blob/main/IMG/FrontEnd/A-R6-0.png)

!["E6-1"](https://github.com/biog5/HPV-master/blob/main/IMG/FrontEnd/A-R6-1.png)

5.2 Se descarga la estructura

!["E6-2"](https://github.com/biog5/HPV-master/blob/main/IMG/FrontEnd/A-R6-2.png)

5.3 Se visualiza

!["E6-3"](https://github.com/biog5/HPV-master/blob/main/IMG/FrontEnd/A-R6-3.png)

!["E6-4"](https://github.com/biog5/HPV-master/blob/main/IMG/FrontEnd/A-R6-4.png)

7- Configurar riesgos
<pre>
  7.1 Permite visualizar calsificaciones actuales
  
  7.2 Modificar, eliminar, agregar o restaurar nuevas cepas a las Clasificaciones

  7.3 Modificar, eliminar, agregar o restaurar nuevas cepas a las Base de Datos 
</pre>

 7.1 Listado actual

!["E7-0"](https://github.com/biog5/HPV-master/blob/main/IMG/FrontEnd/A-R7-0.png)

 7.2 Ejemplo agregar una cepa a una clasificacion
 
!["E7-1"](https://github.com/biog5/HPV-master/blob/main/IMG/FrontEnd/A-R7-1.png)

!["E7-2"](https://github.com/biog5/HPV-master/blob/main/IMG/FrontEnd/A-R7-2.png)

 7.2 Ejemplo agregar una nueva cepa
 
!["E7-3"](https://github.com/biog5/HPV-master/blob/main/IMG/FrontEnd/A-R7-3.png)

!["E7-4"](https://github.com/biog5/HPV-master/blob/main/IMG/FrontEnd/A-R7-4.png)

!["E7-5"](https://github.com/biog5/HPV-master/blob/main/IMG/FrontEnd/A-R7-5.png)


Version 1.4.0
Veanos en nuestro canal de Youtube: @biog5
