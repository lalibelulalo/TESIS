# Resultados

## Clado Calothrix

El clado calothrix contiene 6 especies  y es de interes ya que segun la filogenia estan estrechamente relacionadas y muestra un cambio en el palindromo mas abundante, pasando de **GCGATCGC** a **TGGCGCCA** (Figure \@ref(fig:FIG12)).
<div class="figure" style="text-align: center">
<img src="Clados/Callothrix_clade/figures/Calothrix_Octanuc_FrecObs_sel32_filogenia_HIG.png" alt="**Filogenia anotada del clado Calothrix.** En esta imagen se muestra un cambio abrupto en la Frecuencia observada de **GCGCATCGC** en las especies PCC\_6303, NIES-3974 y 336-3." width="100%" />
<p class="caption">(\#fig:FIG12)**Filogenia anotada del clado Calothrix.** En esta imagen se muestra un cambio abrupto en la Frecuencia observada de **GCGCATCGC** en las especies PCC\_6303, NIES-3974 y 336-3.</p>
</div>

### Red de transiciones

Para hacer mas visual la reconstrucción, construimos una red de las transiciones entre los estados ancestrales. Esto lo hicimos en r usando la función ```Create_Transition_Table()```:


```r
source("ASR_Orth_Functions/NodeAndEdges.R")
Nodes.Edges <- Create_Transition_Table(SitesTable = "Clados/Callothrix_clade/PALINDROMES/GCGATCGC/336-3/Orthologues_Palindrome_sites.txt",
                                EvolutionModel = "F81",
                                Method = "bayes",
                                Phylogeny = "Clados/Callothrix_clade/SpeciesTree_rooted.txt",
                                OrthoPath = "Clados/Callothrix_clade/PALINDROMES/GCGATCGC/336-3/Only_ORTHOLOGUES/")

```

Posteriormente creamos la red usando la función ```Create_Network()```:



y visualizamos dicha red .



Para visualizar la red usamos la paqueteria ```networkD3```. Hicimos 2 figuras, la (Figura \@ref(fig:FIG13)) muestra la red como una conexión de nodos a través de vertices con un grosor proporcional al numero de veces que ocurrió cada transición. En dicha red podemos ver algunos nodos con bordes muy gruesos como **GCAATTGC**, **GCAATCGC**, **GCAATAGC**, **GCGATTGC** (Tabla \@ref(tab:TAB6)).



<div class="figure" style="text-align: center">
<img src="Resultados_files/figure-epub3/FIG13-1.png" alt="**Red de todas las transiciones del clado Calothrix.** En esta imagen se muestra una red en la que cada nodo es un octanucleótido el cual esta unido a otro nodo por un vertice. Dicho vertice tiene un grosor proporcional al numero de veces que dicha transición ocurrió en la reconstrucción ancestral de sitios."  />
<p class="caption">(\#fig:FIG13)**Red de todas las transiciones del clado Calothrix.** En esta imagen se muestra una red en la que cada nodo es un octanucleótido el cual esta unido a otro nodo por un vertice. Dicho vertice tiene un grosor proporcional al numero de veces que dicha transición ocurrió en la reconstrucción ancestral de sitios.</p>
</div>

En la (Figura \@ref(fig:FIG14)) podemos ver las transiciones de una forma mas ordenada, con el numero de ocurrencias y la dirección en la que ocurrieron.
<div class="figure" style="text-align: center">
<img src="Resultados_files/figure-epub3/FIG14-1.png" alt="**Red de todas las transiciones del clado Calothrix.** En esta imagen se muestra la red  de la Figura \@ref(fig:FIG13) de una forma mas visual y con el numero de veces que ocurrio cada transición, asi como la dirección en la que ocurrió."  />
<p class="caption">(\#fig:FIG14)**Red de todas las transiciones del clado Calothrix.** En esta imagen se muestra la red  de la Figura \@ref(fig:FIG13) de una forma mas visual y con el numero de veces que ocurrio cada transición, asi como la dirección en la que ocurrió.</p>
</div>

<table class=" lightable-paper table table-striped table-condensed" style='font-family: "Arial Narrow", arial, helvetica, sans-serif; margin-left: auto; margin-right: auto; font-size: 11px; width: auto !important; margin-left: auto; margin-right: auto;'>
<caption style="font-size: initial !important;">(\#tab:TAB6)Transiciones entre nodos.</caption>
 <thead>
  <tr>
   <th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;"> from </th>
   <th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;"> to </th>
   <th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;"> weight </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: black !important;">GCGATTGC</span> </td>
   <td style="text-align:left;font-weight: bold;vertical-align: top !important;" rowspan="9"> <span style="     color: red !important;">GCGATCGC</span> </td>
   <td style="text-align:left;"> <span style="display: block; padding: 0 4px; border-radius: 4px; background-color: #ffa500">95</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: blue !important;">GCAATTGC</span> </td>
   
   <td style="text-align:left;"> <span style="display: block; padding: 0 4px; border-radius: 4px; background-color: #ffb01f">84</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: black !important;">GCTATCGC</span> </td>
   
   <td style="text-align:left;"> <span style="display: block; padding: 0 4px; border-radius: 4px; background-color: #ffc04e">68</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: black !important;">GCTATTGC</span> </td>
   
   <td style="text-align:left;"> <span style="display: block; padding: 0 4px; border-radius: 4px; background-color: #ffc151">67</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: black !important;">GCAATCGC</span> </td>
   
   <td style="text-align:left;"> <span style="display: block; padding: 0 4px; border-radius: 4px; background-color: #ffc65f">62</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: black !important;">--------</span> </td>
   
   <td style="text-align:left;"> <span style="display: block; padding: 0 4px; border-radius: 4px; background-color: #ffcc71">56</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: black !important;">GCGACCGC</span> </td>
   
   <td style="text-align:left;"> <span style="display: block; padding: 0 4px; border-radius: 4px; background-color: #ffcf79">53</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: black !important;">GCCATCGC</span> </td>
   
   <td style="text-align:left;"> <span style="display: block; padding: 0 4px; border-radius: 4px; background-color: #ffd588">48</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: black !important;">GCGACCGT</span> </td>
   
   <td style="text-align:left;"> <span style="display: block; padding: 0 4px; border-radius: 4px; background-color: #ffd68b">47</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: red !important;">GCGATCGC</span> </td>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: black !important;">GCGATTGC</span> </td>
   <td style="text-align:left;"> <span style="display: block; padding: 0 4px; border-radius: 4px; background-color: #ffd890">45</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: black !important;">GCGATAGC</span> </td>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: red !important;">GCGATCGC</span> </td>
   <td style="text-align:left;"> <span style="display: block; padding: 0 4px; border-radius: 4px; background-color: #ffd993">44</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: blue !important;">GCAATTGC</span> </td>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: black !important;">GCGATTGC</span> </td>
   <td style="text-align:left;"> <span style="display: block; padding: 0 4px; border-radius: 4px; background-color: #ffda96">43</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: black !important;">GCCATTGC</span> </td>
   <td style="text-align:left;font-weight: bold;vertical-align: top !important;" rowspan="2"> <span style="     color: red !important;">GCGATCGC</span> </td>
   <td style="text-align:left;"> <span style="display: block; padding: 0 4px; border-radius: 4px; background-color: #ffdd9f">40</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: black !important;">GCAATAGC</span> </td>
   
   <td style="text-align:left;"> <span style="display: block; padding: 0 4px; border-radius: 4px; background-color: #ffdfa5">38</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: red !important;">GCGATCGC</span> </td>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: black !important;">GCAATCGC</span> </td>
   <td style="text-align:left;"> <span style="display: block; padding: 0 4px; border-radius: 4px; background-color: #ffe3b0">34</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;font-weight: bold;vertical-align: top !important;" rowspan="2"> <span style="     color: black !important;">GCTATTGC</span> </td>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: black !important;">GCTATCGC</span> </td>
   <td style="text-align:left;"> <span style="display: block; padding: 0 4px; border-radius: 4px; background-color: #ffe5b6">32</span> </td>
  </tr>
  <tr>
   
   <td style="text-align:left;font-weight: bold;"> <span style="     color: black !important;">GCGATTGC</span> </td>
   <td style="text-align:left;"> <span style="display: block; padding: 0 4px; border-radius: 4px; background-color: #ffe6b9">31</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: black !important;">GTGATCGC</span> </td>
   <td style="text-align:left;font-weight: bold;vertical-align: top !important;" rowspan="3"> <span style="     color: red !important;">GCGATCGC</span> </td>
   <td style="text-align:left;"> <span style="display: block; padding: 0 4px; border-radius: 4px; background-color: #ffe8bf">29</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: black !important;">GTGACCGC</span> </td>
   
   <td style="text-align:left;"> <span style="display: block; padding: 0 4px; border-radius: 4px; background-color: #ffecca">25</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: black !important;">GTGACCGT</span> </td>
   
   <td style="text-align:left;"> <span style="display: block; padding: 0 4px; border-radius: 4px; background-color: #ffedcd">24</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: black !important;">GCGACCGT</span> </td>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: black !important;">GCGACCGC</span> </td>
   <td style="text-align:left;"> <span style="display: block; padding: 0 4px; border-radius: 4px; background-color: #ffefd3">22</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: blue !important;">GCAATTGC</span> </td>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: black !important;">GCAATCGC</span> </td>
   <td style="text-align:left;"> <span style="display: block; padding: 0 4px; border-radius: 4px; background-color: #ffefd3">22</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: black !important;">GCGATCGT</span> </td>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: red !important;">GCGATCGC</span> </td>
   <td style="text-align:left;"> <span style="display: block; padding: 0 4px; border-radius: 4px; background-color: #ffefd3">22</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: red !important;">GCGATCGC</span> </td>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: blue !important;">GCAATTGC</span> </td>
   <td style="text-align:left;"> <span style="display: block; padding: 0 4px; border-radius: 4px; background-color: #ffefd3">22</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: black !important;">GCTATAGC</span> </td>
   <td style="text-align:left;font-weight: bold;vertical-align: top !important;" rowspan="3"> <span style="     color: red !important;">GCGATCGC</span> </td>
   <td style="text-align:left;"> <span style="display: block; padding: 0 4px; border-radius: 4px; background-color: #fff0d6">21</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: black !important;">GCGTTCGC</span> </td>
   
   <td style="text-align:left;"> <span style="display: block; padding: 0 4px; border-radius: 4px; background-color: #fff0d6">21</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: black !important;">GCGATAGA</span> </td>
   
   <td style="text-align:left;"> <span style="display: block; padding: 0 4px; border-radius: 4px; background-color: #fff0d6">21</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: black !important;">GCTATTGC</span> </td>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: blue !important;">GCAATTGC</span> </td>
   <td style="text-align:left;"> <span style="display: block; padding: 0 4px; border-radius: 4px; background-color: #fff2dc">19</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: black !important;">GCCATAGC</span> </td>
   <td style="text-align:left;font-weight: bold;vertical-align: top !important;" rowspan="2"> <span style="     color: red !important;">GCGATCGC</span> </td>
   <td style="text-align:left;"> <span style="display: block; padding: 0 4px; border-radius: 4px; background-color: #fff3df">18</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: black !important;">GTGATAGA</span> </td>
   
   <td style="text-align:left;"> <span style="display: block; padding: 0 4px; border-radius: 4px; background-color: #fff3df">18</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: blue !important;">GCAATTGC</span> </td>
   <td style="text-align:left;font-weight: bold;vertical-align: top !important;" rowspan="2"> <span style="     color: black !important;">GCTATTGC</span> </td>
   <td style="text-align:left;"> <span style="display: block; padding: 0 4px; border-radius: 4px; background-color: #fff3df">18</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;font-weight: bold;vertical-align: top !important;" rowspan="4"> <span style="     color: red !important;">GCGATCGC</span> </td>
   
   <td style="text-align:left;"> <span style="display: block; padding: 0 4px; border-radius: 4px; background-color: #fff4e2">17</span> </td>
  </tr>
  <tr>
   
   <td style="text-align:left;font-weight: bold;"> <span style="     color: black !important;">GTGATCGC</span> </td>
   <td style="text-align:left;"> <span style="display: block; padding: 0 4px; border-radius: 4px; background-color: #fff4e2">17</span> </td>
  </tr>
  <tr>
   
   <td style="text-align:left;font-weight: bold;"> <span style="     color: black !important;">GCGATAGC</span> </td>
   <td style="text-align:left;"> <span style="display: block; padding: 0 4px; border-radius: 4px; background-color: #fff4e2">17</span> </td>
  </tr>
  <tr>
   
   <td style="text-align:left;font-weight: bold;"> <span style="     color: black !important;">GCGACCGC</span> </td>
   <td style="text-align:left;"> <span style="display: block; padding: 0 4px; border-radius: 4px; background-color: #fff4e2">17</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: black !important;">GCAATAGC</span> </td>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: black !important;">GCAATCGC</span> </td>
   <td style="text-align:left;"> <span style="display: block; padding: 0 4px; border-radius: 4px; background-color: #fff5e4">16</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: black !important;">GCCATTGC</span> </td>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: black !important;">GCGATTGC</span> </td>
   <td style="text-align:left;"> <span style="display: block; padding: 0 4px; border-radius: 4px; background-color: #fff6e7">15</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;font-weight: bold;vertical-align: top !important;" rowspan="2"> <span style="     color: red !important;">GCGATCGC</span> </td>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: black !important;">GCGATCGA</span> </td>
   <td style="text-align:left;"> <span style="display: block; padding: 0 4px; border-radius: 4px; background-color: #fff6e7">15</span> </td>
  </tr>
  <tr>
   
   <td style="text-align:left;font-weight: bold;"> <span style="     color: black !important;">GCGATCGT</span> </td>
   <td style="text-align:left;"> <span style="display: block; padding: 0 4px; border-radius: 4px; background-color: #fff6e7">15</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: black !important;">GTGACCGT</span> </td>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: black !important;">GCGACCGT</span> </td>
   <td style="text-align:left;"> <span style="display: block; padding: 0 4px; border-radius: 4px; background-color: #fff7ea">14</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: black !important;">GCCATTGC</span> </td>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: black !important;">GCCATCGC</span> </td>
   <td style="text-align:left;"> <span style="display: block; padding: 0 4px; border-radius: 4px; background-color: #fff7ea">14</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: red !important;">GCGATCGC</span> </td>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: black !important;">GCGAACGC</span> </td>
   <td style="text-align:left;"> <span style="display: block; padding: 0 4px; border-radius: 4px; background-color: #fff7ea">14</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: black !important;">GCGGTCGC</span> </td>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: red !important;">GCGATCGC</span> </td>
   <td style="text-align:left;"> <span style="display: block; padding: 0 4px; border-radius: 4px; background-color: #fff7ea">14</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: black !important;">GCAATAGC</span> </td>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: blue !important;">GCAATTGC</span> </td>
   <td style="text-align:left;"> <span style="display: block; padding: 0 4px; border-radius: 4px; background-color: #fff8ed">13</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: red !important;">GCGATCGC</span> </td>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: black !important;">GCTATCGC</span> </td>
   <td style="text-align:left;"> <span style="display: block; padding: 0 4px; border-radius: 4px; background-color: #fff8ed">13</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: black !important;">GCGATCGA</span> </td>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: red !important;">GCGATCGC</span> </td>
   <td style="text-align:left;"> <span style="display: block; padding: 0 4px; border-radius: 4px; background-color: #fff8ed">13</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: black !important;">GTGACCGC</span> </td>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: black !important;">GCGACCGC</span> </td>
   <td style="text-align:left;"> <span style="display: block; padding: 0 4px; border-radius: 4px; background-color: #fff9f0">12</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;font-weight: bold;vertical-align: top !important;" rowspan="2"> <span style="     color: red !important;">GCGATCGC</span> </td>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: black !important;">GCGTTCGC</span> </td>
   <td style="text-align:left;"> <span style="display: block; padding: 0 4px; border-radius: 4px; background-color: #fff9f0">12</span> </td>
  </tr>
  <tr>
   
   <td style="text-align:left;font-weight: bold;"> <span style="     color: black !important;">GCGGTCGC</span> </td>
   <td style="text-align:left;"> <span style="display: block; padding: 0 4px; border-radius: 4px; background-color: #fff9f0">12</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: black !important;">GCAATAGC</span> </td>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: black !important;">GCGATAGC</span> </td>
   <td style="text-align:left;"> <span style="display: block; padding: 0 4px; border-radius: 4px; background-color: #fff9f0">12</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: black !important;">GCGGTTGC</span> </td>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: red !important;">GCGATCGC</span> </td>
   <td style="text-align:left;"> <span style="display: block; padding: 0 4px; border-radius: 4px; background-color: #fffaf3">11</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: red !important;">GCGATCGC</span> </td>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: black !important;">GCGAGCGC</span> </td>
   <td style="text-align:left;"> <span style="display: block; padding: 0 4px; border-radius: 4px; background-color: #fffaf3">11</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: black !important;">ACGTTCGC</span> </td>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: red !important;">GCGATCGC</span> </td>
   <td style="text-align:left;"> <span style="display: block; padding: 0 4px; border-radius: 4px; background-color: #fffaf3">11</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: black !important;">GTGACCGT</span> </td>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: black !important;">GTGACCGC</span> </td>
   <td style="text-align:left;"> <span style="display: block; padding: 0 4px; border-radius: 4px; background-color: #fffaf3">11</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: black !important;">GCTATAGC</span> </td>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: black !important;">GCGATAGC</span> </td>
   <td style="text-align:left;"> <span style="display: block; padding: 0 4px; border-radius: 4px; background-color: #fffaf3">11</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: red !important;">GCGATCGC</span> </td>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: black !important;">GCGCTCGC</span> </td>
   <td style="text-align:left;"> <span style="display: block; padding: 0 4px; border-radius: 4px; background-color: #fffaf3">11</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: black !important;">GCCATTGC</span> </td>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: black !important;">GCTATTGC</span> </td>
   <td style="text-align:left;"> <span style="display: block; padding: 0 4px; border-radius: 4px; background-color: #fffbf6">10</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: black !important;">GCGAGCGC</span> </td>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: red !important;">GCGATCGC</span> </td>
   <td style="text-align:left;"> <span style="display: block; padding: 0 4px; border-radius: 4px; background-color: #fffbf6">10</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;font-weight: bold;vertical-align: top !important;" rowspan="3"> <span style="     color: red !important;">GCGATCGC</span> </td>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: black !important;">GCGATCAC</span> </td>
   <td style="text-align:left;"> <span style="display: block; padding: 0 4px; border-radius: 4px; background-color: #fffbf6">10</span> </td>
  </tr>
  <tr>
   
   <td style="text-align:left;font-weight: bold;"> <span style="     color: black !important;">GCGTTAGC</span> </td>
   <td style="text-align:left;"> <span style="display: block; padding: 0 4px; border-radius: 4px; background-color: #fffbf6">10</span> </td>
  </tr>
  <tr>
   
   <td style="text-align:left;font-weight: bold;"> <span style="     color: black !important;">GTGATTGC</span> </td>
   <td style="text-align:left;"> <span style="display: block; padding: 0 4px; border-radius: 4px; background-color: #fffbf6">10</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: black !important;">GCAGTTGC</span> </td>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: black !important;">GCGGTTGC</span> </td>
   <td style="text-align:left;"> <span style="display: block; padding: 0 4px; border-radius: 4px; background-color: #fffbf6">10</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: black !important;">GAGACCGC</span> </td>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: red !important;">GCGATCGC</span> </td>
   <td style="text-align:left;"> <span style="display: block; padding: 0 4px; border-radius: 4px; background-color: #fffcf9">9</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: red !important;">GCGATCGC</span> </td>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: black !important;">GCGATAAC</span> </td>
   <td style="text-align:left;"> <span style="display: block; padding: 0 4px; border-radius: 4px; background-color: #fffcf9">9</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: black !important;">GAGATAGA</span> </td>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: red !important;">GCGATCGC</span> </td>
   <td style="text-align:left;"> <span style="display: block; padding: 0 4px; border-radius: 4px; background-color: #fffcf9">9</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: red !important;">GCGATCGC</span> </td>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: black !important;">GCCATCGC</span> </td>
   <td style="text-align:left;"> <span style="display: block; padding: 0 4px; border-radius: 4px; background-color: #fffcf9">9</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: black !important;">GCGACAGA</span> </td>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: black !important;">GCGATAGA</span> </td>
   <td style="text-align:left;"> <span style="display: block; padding: 0 4px; border-radius: 4px; background-color: #fffcf9">9</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: black !important;">GCTATAGC</span> </td>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: black !important;">GCTATCGC</span> </td>
   <td style="text-align:left;"> <span style="display: block; padding: 0 4px; border-radius: 4px; background-color: #fffcf9">9</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: red !important;">GCGATCGC</span> </td>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: black !important;">ACGATCGC</span> </td>
   <td style="text-align:left;"> <span style="display: block; padding: 0 4px; border-radius: 4px; background-color: #fffcf9">9</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: black !important;">GTGATAGA</span> </td>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: black !important;">GCGATAGA</span> </td>
   <td style="text-align:left;"> <span style="display: block; padding: 0 4px; border-radius: 4px; background-color: #fffcf9">9</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: black !important;">GCGATAAA</span> </td>
   <td style="text-align:left;font-weight: bold;vertical-align: top !important;" rowspan="2"> <span style="     color: red !important;">GCGATCGC</span> </td>
   <td style="text-align:left;"> <span style="display: block; padding: 0 4px; border-radius: 4px; background-color: #fffcf9">9</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: black !important;">GCGCTCGC</span> </td>
   
   <td style="text-align:left;"> <span style="display: block; padding: 0 4px; border-radius: 4px; background-color: #fffcf9">9</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: red !important;">GCGATCGC</span> </td>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: black !important;">GCGATTGT</span> </td>
   <td style="text-align:left;"> <span style="display: block; padding: 0 4px; border-radius: 4px; background-color: #fffcf9">9</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: black !important;">GCTGTTGC</span> </td>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: red !important;">GCGATCGC</span> </td>
   <td style="text-align:left;"> <span style="display: block; padding: 0 4px; border-radius: 4px; background-color: #fffdfc">8</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: black !important;">GCCATCGC</span> </td>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: black !important;">GCTATCGC</span> </td>
   <td style="text-align:left;"> <span style="display: block; padding: 0 4px; border-radius: 4px; background-color: #fffdfc">8</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: black !important;">GCTATCGC</span> </td>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: black !important;">GCAATCGC</span> </td>
   <td style="text-align:left;"> <span style="display: block; padding: 0 4px; border-radius: 4px; background-color: #fffdfc">8</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: black !important;">GAGACCGT</span> </td>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: red !important;">GCGATCGC</span> </td>
   <td style="text-align:left;"> <span style="display: block; padding: 0 4px; border-radius: 4px; background-color: #fffdfc">8</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: black !important;">GCTATAGC</span> </td>
   <td style="text-align:left;font-weight: bold;vertical-align: top !important;" rowspan="2"> <span style="     color: black !important;">GCTATTGC</span> </td>
   <td style="text-align:left;"> <span style="display: block; padding: 0 4px; border-radius: 4px; background-color: #fffdfc">8</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: black !important;">GCTGTTGC</span> </td>
   
   <td style="text-align:left;"> <span style="display: block; padding: 0 4px; border-radius: 4px; background-color: #fffdfc">8</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: black !important;">GCGAACGC</span> </td>
   <td style="text-align:left;font-weight: bold;vertical-align: top !important;" rowspan="3"> <span style="     color: red !important;">GCGATCGC</span> </td>
   <td style="text-align:left;"> <span style="display: block; padding: 0 4px; border-radius: 4px; background-color: #fffdfc">8</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: black !important;">ACGATCGC</span> </td>
   
   <td style="text-align:left;"> <span style="display: block; padding: 0 4px; border-radius: 4px; background-color: #fffdfc">8</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: black !important;">ACAATTGC</span> </td>
   
   <td style="text-align:left;"> <span style="display: block; padding: 0 4px; border-radius: 4px; background-color: #fffdfc">8</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: black !important;">GCCATAGC</span> </td>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: black !important;">GCGATAGC</span> </td>
   <td style="text-align:left;"> <span style="display: block; padding: 0 4px; border-radius: 4px; background-color: #fffdfc">8</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: black !important;">GTGATTGC</span> </td>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: red !important;">GCGATCGC</span> </td>
   <td style="text-align:left;"> <span style="display: block; padding: 0 4px; border-radius: 4px; background-color: #ffffff">7</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: black !important;">GCGGTTGC</span> </td>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: black !important;">GCGATTGC</span> </td>
   <td style="text-align:left;"> <span style="display: block; padding: 0 4px; border-radius: 4px; background-color: #ffffff">7</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: black !important;">GCTTTAGC</span> </td>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: red !important;">GCGATCGC</span> </td>
   <td style="text-align:left;"> <span style="display: block; padding: 0 4px; border-radius: 4px; background-color: #ffffff">7</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: red !important;">GCGATCGC</span> </td>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: black !important;">GCAGTTGC</span> </td>
   <td style="text-align:left;"> <span style="display: block; padding: 0 4px; border-radius: 4px; background-color: #ffffff">7</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: black !important;">GCGATAGC</span> </td>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: black !important;">GCGATTGC</span> </td>
   <td style="text-align:left;"> <span style="display: block; padding: 0 4px; border-radius: 4px; background-color: #ffffff">7</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: black !important;">GTAACCGC</span> </td>
   <td style="text-align:left;font-weight: bold;vertical-align: top !important;" rowspan="2"> <span style="     color: red !important;">GCGATCGC</span> </td>
   <td style="text-align:left;"> <span style="display: block; padding: 0 4px; border-radius: 4px; background-color: #ffffff">7</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: black !important;">GTGATAGC</span> </td>
   
   <td style="text-align:left;"> <span style="display: block; padding: 0 4px; border-radius: 4px; background-color: #ffffff">7</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: black !important;">GCAGTTGC</span> </td>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: blue !important;">GCAATTGC</span> </td>
   <td style="text-align:left;"> <span style="display: block; padding: 0 4px; border-radius: 4px; background-color: #ffffff">7</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: black !important;">GCGACCGA</span> </td>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: red !important;">GCGATCGC</span> </td>
   <td style="text-align:left;"> <span style="display: block; padding: 0 4px; border-radius: 4px; background-color: #ffffff">7</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;font-weight: bold;vertical-align: top !important;" rowspan="3"> <span style="     color: red !important;">GCGATCGC</span> </td>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: black !important;">GCCATTGC</span> </td>
   <td style="text-align:left;"> <span style="display: block; padding: 0 4px; border-radius: 4px; background-color: #ffffff">7</span> </td>
  </tr>
  <tr>
   
   <td style="text-align:left;font-weight: bold;"> <span style="     color: black !important;">CCGATCGC</span> </td>
   <td style="text-align:left;"> <span style="display: block; padding: 0 4px; border-radius: 4px; background-color: #ffffff">7</span> </td>
  </tr>
  <tr>
   
   <td style="text-align:left;font-weight: bold;"> <span style="     color: black !important;">GCGATTGG</span> </td>
   <td style="text-align:left;"> <span style="display: block; padding: 0 4px; border-radius: 4px; background-color: #ffffff">7</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: black !important;">GCAATTAC</span> </td>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: blue !important;">GCAATTGC</span> </td>
   <td style="text-align:left;"> <span style="display: block; padding: 0 4px; border-radius: 4px; background-color: #ffffff">7</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: black !important;">GCGAACGT</span> </td>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: red !important;">GCGATCGC</span> </td>
   <td style="text-align:left;"> <span style="display: block; padding: 0 4px; border-radius: 4px; background-color: #ffffff">7</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: black !important;">GCTATCGC</span> </td>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: black !important;">GCCATCGC</span> </td>
   <td style="text-align:left;"> <span style="display: block; padding: 0 4px; border-radius: 4px; background-color: #ffffff">7</span> </td>
  </tr>
</tbody>
</table>



#### Transiciones entre Nodo 9 y Nodo 10
Para entender más como es que se gana o se pierden los sitios palindrómicos revisamos la transición en tre los nodos 9 y 10. Esto es porque es esta transicion de nodos la que separa a los dos subclados entre los que hay una repentino cambio de abundancia de sitios palindrómicos (Figura \@ref(fig:FIG15)).

<div class="figure">
<img src="Resultados_files/figure-epub3/FIG15-1.png" alt="**Filogenia del clado Calotrix**. En rojo y azul se muestran los subclados unidos (en verde) por la transición entre los nodos 9 y 10. "  />
<p class="caption">(\#fig:FIG15)**Filogenia del clado Calotrix**. En rojo y azul se muestran los subclados unidos (en verde) por la transición entre los nodos 9 y 10. </p>
</div>

Para hacer esto filtramos los datos de la red para mostrar unicamente las transiciones que se dieron entre los nodos 9 y 10 e hicimos las mismas figuras.
En la (Figura \@ref(fig:FIG16)) se muestra la red como una conexión de nodos a través de vertices con un grosor proporcional al numero de veces que ocurrió cada transición. En la (Figura \@ref(fig:FIG17)) podemos ver las transiciones de una forma mas ordenada, con el numero de ocurrencias y la dirección en la que ocurrieron.



<div class="figure" style="text-align: center">
<img src="Resultados_files/figure-epub3/FIG16-1.png" alt="**Red de las transiciones entre los Nodos 9 y 10 del clado Calothrix.** En esta imagen se muestra una red en la que cada nodo es un octanucleótido el cual esta unido a otro nodo por un vertice. Dicho vertice tiene un grosor proporcional al numero de veces que dicha transición ocurrió en la reconstrucción ancestral de sitios."  />
<p class="caption">(\#fig:FIG16)**Red de las transiciones entre los Nodos 9 y 10 del clado Calothrix.** En esta imagen se muestra una red en la que cada nodo es un octanucleótido el cual esta unido a otro nodo por un vertice. Dicho vertice tiene un grosor proporcional al numero de veces que dicha transición ocurrió en la reconstrucción ancestral de sitios.</p>
</div>

<div class="figure" style="text-align: center">
<img src="Resultados_files/figure-epub3/FIG17-1.png" alt="**Red de las transiciones entre los nodos 9 y 10 del clado Calothrix.** En esta imagen se muestra la red  de la Figura \@ref(fig:FIG16) de una forma mas visual y con el numero de veces que ocurrio cada transición, asi como la dirección en la que ocurrió."  />
<p class="caption">(\#fig:FIG17)**Red de las transiciones entre los nodos 9 y 10 del clado Calothrix.** En esta imagen se muestra la red  de la Figura \@ref(fig:FIG16) de una forma mas visual y con el numero de veces que ocurrio cada transición, asi como la dirección en la que ocurrió.</p>
</div>

### Mutaciones en los codones

Para entender como es que se van ganando o perdiendo los sitios palindrómicos hicimos un análisis del tipo mutaciones de los sitios. Esto lo hicimos viendo en que marco de lectura se encontraba cada nodo y revisando la secuencia de aminoacidos que codificaban. En la (Figura \@ref(fig:FIG18)) mostramos 3 gráficos que indican la abundancia de los peptidos codificados por los sitios palindrómicos de acuerdo al marco de lectura en el que se encuentran.

<div class="figure" style="text-align: center">
<img src="Resultados_files/figure-epub3/FIG18-1.png" alt="**Abundancia de peptidos por cada nodo segun el marco de lectura.**." width="50%" /><img src="Resultados_files/figure-epub3/FIG18-2.png" alt="**Abundancia de peptidos por cada nodo segun el marco de lectura.**." width="50%" /><img src="Resultados_files/figure-epub3/FIG18-3.png" alt="**Abundancia de peptidos por cada nodo segun el marco de lectura.**." width="50%" />
<p class="caption">(\#fig:FIG18)**Abundancia de peptidos por cada nodo segun el marco de lectura.**.</p>
</div>

En la (Figura \@ref(fig:FIG19)) mostramos 3 gráficos que indican la abundancia del tipo de mutaciones que hay en cada nodo de acuerdo al marco de lectura. Lo sitios de mutaciones mostrados pueden ser de los siguientes tipos:

* Conservative (la secuencia de AA cambió pero tiene similitud de acuerdo al score de BLOSUM62)
* ConservativeNoSiteMut (la secuencia de AA cambió pero tiene similitud de acuerdo al score de BLOSUM62. Sin embargo, el sitio no sufrió mutaciones)
* Deletion (La secuencia de AA tiene sufrio 1 o mas deleciones)
* NoMutation (La secuencia de AA no sufrio mutaciones)
* NoSynonym (La secuencia de AA cambió)
* NoSynonymNoSiteMut (La secuencia de AA cambió. Sin embargo, el sitio no sufrió mutaciones.)
* Synonym (El sitio sufrió mutaciones. Sin embargo, la secuencia de AA no cambió.)

<div class="figure" style="text-align: center">
<img src="Resultados_files/figure-epub3/FIG19-1.png" alt="**Abundancia del tipo de sustitución por cada nodo segun el marco de lectura.**." width="50%" /><img src="Resultados_files/figure-epub3/FIG19-2.png" alt="**Abundancia del tipo de sustitución por cada nodo segun el marco de lectura.**." width="50%" /><img src="Resultados_files/figure-epub3/FIG19-3.png" alt="**Abundancia del tipo de sustitución por cada nodo segun el marco de lectura.**." width="50%" />
<p class="caption">(\#fig:FIG19)**Abundancia del tipo de sustitución por cada nodo segun el marco de lectura.**.</p>
</div>

#### Análisis de sitios en los cuales su ancestro era HIP1

Para tratar de entender como es que los sitios HIP1 se pierden hicimos un análisis unicamente en en las transiciones en las que el nodo ancestral tenia un sitio HIP1.

En la (Figura \@ref(fig:FIG20)) mostramos 3 gráficos que indican la frecuencia del tipo de sustituciones que hubo para estos casos para cada nodo en cada uno de los marcos de lectura.

En la (Figura \@ref(fig:FIG21)) mostramos 3 gráficos (uno por cada marco de lectura) que indican la frecuencia de las mutaciones en cada uno de los 8 nucleótidos del sitio HIP.

En la (Figura \@ref(fig:FIG22)) mostramos 3 gráficos (uno por cada marco de lectura) que indican la frecuencia del tipo sustitucion de bases.

<div class="figure" style="text-align: center">
<img src="Resultados_files/figure-epub3/FIG20-1.png" alt="**Abundancia del tipo de sustitución por cada nodo segun el marco de lectura. Unicamente para transiciones en los que el nodo ancestral era un sitio HIP1.**" width="50%" /><img src="Resultados_files/figure-epub3/FIG20-2.png" alt="**Abundancia del tipo de sustitución por cada nodo segun el marco de lectura. Unicamente para transiciones en los que el nodo ancestral era un sitio HIP1.**" width="50%" /><img src="Resultados_files/figure-epub3/FIG20-3.png" alt="**Abundancia del tipo de sustitución por cada nodo segun el marco de lectura. Unicamente para transiciones en los que el nodo ancestral era un sitio HIP1.**" width="50%" />
<p class="caption">(\#fig:FIG20)**Abundancia del tipo de sustitución por cada nodo segun el marco de lectura. Unicamente para transiciones en los que el nodo ancestral era un sitio HIP1.**</p>
</div>

<div class="figure" style="text-align: center">
<img src="Resultados_files/figure-epub3/FIG21-1.png" alt="**Frecuencia de las mutaciones de cada nucleótido del sitio HIP para cada nodo segun el marco de lectura.**." width="50%" /><img src="Resultados_files/figure-epub3/FIG21-2.png" alt="**Frecuencia de las mutaciones de cada nucleótido del sitio HIP para cada nodo segun el marco de lectura.**." width="50%" /><img src="Resultados_files/figure-epub3/FIG21-3.png" alt="**Frecuencia de las mutaciones de cada nucleótido del sitio HIP para cada nodo segun el marco de lectura.**." width="50%" />
<p class="caption">(\#fig:FIG21)**Frecuencia de las mutaciones de cada nucleótido del sitio HIP para cada nodo segun el marco de lectura.**.</p>
</div>

<div class="figure" style="text-align: center">
<img src="Resultados_files/figure-epub3/FIG22-1.png" alt="**Frecuencia del tipo de sustituciónes de base en los sitios HIP para cada marco de lectura**." width="50%" /><img src="Resultados_files/figure-epub3/FIG22-2.png" alt="**Frecuencia del tipo de sustituciónes de base en los sitios HIP para cada marco de lectura**." width="50%" /><img src="Resultados_files/figure-epub3/FIG22-3.png" alt="**Frecuencia del tipo de sustituciónes de base en los sitios HIP para cada marco de lectura**." width="50%" />
<p class="caption">(\#fig:FIG22)**Frecuencia del tipo de sustituciónes de base en los sitios HIP para cada marco de lectura**.</p>
</div>

#### Análisis de sitios en los cuales solo el nodo actual tiene HIP1

Para tratar de entender como es que los sitios HIP se ganan, hicimos un analisis unicamente en las transiciones en las que el nodo actual tenia un sitio HIP1. 

En la Figura \@ref(fig:FIG23) mostramos 3 gráficos (uno por cada marco de lectura) que indican la frecuencia del tipo de sustituciones que hubo para estos casos para cada nodo en cada uno de los marcos de lectura.


#### Análisis de sitios en los cuales solo el nodo actual tiene HIP1

Para tratar de entender como es que los sitios HIP se ganan, hicimos un analisis unicamente en las transiciones en las que el nodo actual tenia un sitio HIP1. 

En la Figura \@ref(fig:FIG23) mostramos 3 gráficos (uno por cada marco de lectura) que indican la frecuencia del tipo de sustituciones que hubo para estos casos para cada nodo en cada uno de los marcos de lectura.

<div class="figure" style="text-align: center">
<img src="Resultados_files/figure-epub3/FIG23-1.png" alt="**Abundancia del tipo de sustitución por cada nodo segun el marco de lectura. Unicamente para transiciones en los que el nodo actual era un sitio HIP1.**." width="50%" /><img src="Resultados_files/figure-epub3/FIG23-2.png" alt="**Abundancia del tipo de sustitución por cada nodo segun el marco de lectura. Unicamente para transiciones en los que el nodo actual era un sitio HIP1.**." width="50%" /><img src="Resultados_files/figure-epub3/FIG23-3.png" alt="**Abundancia del tipo de sustitución por cada nodo segun el marco de lectura. Unicamente para transiciones en los que el nodo actual era un sitio HIP1.**." width="50%" />
<p class="caption">(\#fig:FIG23)**Abundancia del tipo de sustitución por cada nodo segun el marco de lectura. Unicamente para transiciones en los que el nodo actual era un sitio HIP1.**.</p>
</div>








































































## Clado A18-40

El clado A18-40 fue de particular interés ya que entre las especies **Synechococcus_sp_A18−40** y **Synechococcus sp RS9915** hay un cambio abrupto del palíndromo con mayor OE. Más aun, los genomas de ambas especies son muy parecidos y ambas especias son cercanas segun la filogenia (Figura \@ref(fig:FIG7)).
<div class="figure" style="text-align: center">
<img src="figures/A18-40_Octanuc_OE_sel32_filogenia_HIG.png" alt="**Filogenia anotada del clado A18-40.** En esta imagen se muestra un cambio abrupto en la tasa OE de la especie Synechococcus sp A18-40." width="100%" />
<p class="caption">(\#fig:FIG7)**Filogenia anotada del clado A18-40.** En esta imagen se muestra un cambio abrupto en la tasa OE de la especie Synechococcus sp A18-40.</p>
</div>

Para saber que tan parecidos eran los genomas hicimos dos análisis de sintenia, uno de enfocado en los ortólogos (Figura \@ref(fig:FIG8)) y otro enfocado en el genoma (Figura \@ref(fig:FIG9)).

<div class="figure" style="text-align: center">
<img src="figures/circos.png" alt="**Sintenía de ortólogos entre las especies Synechococcus sp A18-40 y Synechococcus sp RS99150.** En esta imagen se hizo un análisis de sintenia de ortólogos para ver que tan parecidos eran los genomas. En azul se muestra la especie Synechococcus sp A18-40 y en verde Synechococcus sp RS99150." width="75%" />
<p class="caption">(\#fig:FIG8)**Sintenía de ortólogos entre las especies Synechococcus sp A18-40 y Synechococcus sp RS99150.** En esta imagen se hizo un análisis de sintenia de ortólogos para ver que tan parecidos eran los genomas. En azul se muestra la especie Synechococcus sp A18-40 y en verde Synechococcus sp RS99150.</p>
</div>

<div class="figure" style="text-align: center">
<img src="figures/out.png" alt="**Sintenia de DNA entre especies Synechococcus sp A18−40 y Synechococcus sp RS99150.** En esta imagen se hizo un análisis de sintenia de ortólogos para ver que tan parecidos eran los genomas. En el eje X se muestra la especie Synechococcus sp A18−40 y en el eje Y Synechococcus sp RS99150." width="800" />
<p class="caption">(\#fig:FIG9)**Sintenia de DNA entre especies Synechococcus sp A18−40 y Synechococcus sp RS99150.** En esta imagen se hizo un análisis de sintenia de ortólogos para ver que tan parecidos eran los genomas. En el eje X se muestra la especie Synechococcus sp A18−40 y en el eje Y Synechococcus sp RS99150.</p>
</div>

### CGTTAACG

El palíndromo con la tasa OE mas alta es **CGTTAACG** y lo tiene la especie **Synechococcus sp A18−40** (Tabla \@ref(tab:TAB1)) con un conteo de 112 sitios para dicho palindromo mientras que en las demás especies oscila entre 3 y 15 sitios. 


  
<table>
<caption>(\#tab:TAB1)Conteo del palíndromo **CGTTAACG** en el clado **A18-40**</caption>
 <thead>
  <tr>
   <th style="text-align:left;"> spp </th>
   <th style="text-align:left;"> palindrome </th>
   <th style="text-align:right;"> obs </th>
   <th style="text-align:right;"> markov3 </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> Synechococcus_sp_A18-46_1 </td>
   <td style="text-align:left;"> CGTTAACG </td>
   <td style="text-align:right;"> 10 </td>
   <td style="text-align:right;"> 9.19 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Synechococcus_sp_A18-40 </td>
   <td style="text-align:left;"> CGTTAACG </td>
   <td style="text-align:right;"> 112 </td>
   <td style="text-align:right;"> 9.71 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Synechococcus_sp_WH_8103 </td>
   <td style="text-align:left;"> CGTTAACG </td>
   <td style="text-align:right;"> 9 </td>
   <td style="text-align:right;"> 9.14 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Synechococcus_sp_RS9915 </td>
   <td style="text-align:left;"> CGTTAACG </td>
   <td style="text-align:right;"> 10 </td>
   <td style="text-align:right;"> 8.58 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Synechococcus_sp_A15-28 </td>
   <td style="text-align:left;"> CGTTAACG </td>
   <td style="text-align:right;"> 3 </td>
   <td style="text-align:right;"> 7.58 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Parasynechococcus_marenigrum_WH_8102 </td>
   <td style="text-align:left;"> CGTTAACG </td>
   <td style="text-align:right;"> 15 </td>
   <td style="text-align:right;"> 9.35 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Synechococcus_sp_A15-24 </td>
   <td style="text-align:left;"> CGTTAACG </td>
   <td style="text-align:right;"> 11 </td>
   <td style="text-align:right;"> 7.99 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Synechococcus_sp_BOUM118 </td>
   <td style="text-align:left;"> CGTTAACG </td>
   <td style="text-align:right;"> 11 </td>
   <td style="text-align:right;"> 8.24 </td>
  </tr>
</tbody>
</table>



#### Ubicación de sitios CGTTAACG en los genomas del clado A18-40
Para tratar de entender la distribución del palindromo en el genoma buscamos la ubicación de cada sitio y analizamos la distancia entre cada uno de ellos (Tabla \@ref(tab:TAB2)).
En dicho análisis pudimos observar que había 101 sitios que se encontraban entre repeticiones de 340 nuclótidos (columna 3 de la Tabla \@ref(tab:TAB2)).



<table>
<caption>(\#tab:TAB2)**Ubicación de los sitios CGTTACG.** La tabla muestra las priemras 15 lineas de la tabla. La primera columna se muestra el numero de sitio. La segunda columna muestra el intervalo en el que se encuentra el palíndromo. La tercera columna muestra a cuantos nucléotidos se encuentra el ultimo sitio. La cuarta columna indica la diferencia entre la distancia del ultimo palindromo y la distancia del siguiente.</caption>
 <thead>
  <tr>
   <th style="text-align:right;"> SiteNumber </th>
   <th style="text-align:right;"> Interval </th>
   <th style="text-align:right;"> Dist2NextPal </th>
   <th style="text-align:right;"> DifBetDist </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 315323:315331 </td>
   <td style="text-align:right;"> 315323 </td>
   <td style="text-align:right;"> -315323 </td>
  </tr>
  <tr>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 522737:522745 </td>
   <td style="text-align:right;"> 207406 </td>
   <td style="text-align:right;"> 107917 </td>
  </tr>
  <tr>
   <td style="text-align:right;"> 3 </td>
   <td style="text-align:right;"> 594946:594954 </td>
   <td style="text-align:right;"> 72201 </td>
   <td style="text-align:right;"> 135205 </td>
  </tr>
  <tr>
   <td style="text-align:right;"> 4 </td>
   <td style="text-align:right;"> 683860:683868 </td>
   <td style="text-align:right;"> 88906 </td>
   <td style="text-align:right;"> -16705 </td>
  </tr>
  <tr>
   <td style="text-align:right;"> 5 </td>
   <td style="text-align:right;"> 893774:893782 </td>
   <td style="text-align:right;"> 209906 </td>
   <td style="text-align:right;"> -121000 </td>
  </tr>
  <tr>
   <td style="text-align:right;"> 6 </td>
   <td style="text-align:right;"> 894122:894130 </td>
   <td style="text-align:right;"> 340 </td>
   <td style="text-align:right;"> 209566 </td>
  </tr>
  <tr>
   <td style="text-align:right;"> 7 </td>
   <td style="text-align:right;"> 894470:894478 </td>
   <td style="text-align:right;"> 340 </td>
   <td style="text-align:right;"> 0 </td>
  </tr>
  <tr>
   <td style="text-align:right;"> 8 </td>
   <td style="text-align:right;"> 894818:894826 </td>
   <td style="text-align:right;"> 340 </td>
   <td style="text-align:right;"> 0 </td>
  </tr>
  <tr>
   <td style="text-align:right;"> 9 </td>
   <td style="text-align:right;"> 895166:895174 </td>
   <td style="text-align:right;"> 340 </td>
   <td style="text-align:right;"> 0 </td>
  </tr>
  <tr>
   <td style="text-align:right;"> 10 </td>
   <td style="text-align:right;"> 895514:895522 </td>
   <td style="text-align:right;"> 340 </td>
   <td style="text-align:right;"> 0 </td>
  </tr>
  <tr>
   <td style="text-align:right;"> 11 </td>
   <td style="text-align:right;"> 895862:895870 </td>
   <td style="text-align:right;"> 340 </td>
   <td style="text-align:right;"> 0 </td>
  </tr>
  <tr>
   <td style="text-align:right;"> 12 </td>
   <td style="text-align:right;"> 896210:896218 </td>
   <td style="text-align:right;"> 340 </td>
   <td style="text-align:right;"> 0 </td>
  </tr>
  <tr>
   <td style="text-align:right;"> 13 </td>
   <td style="text-align:right;"> 896558:896566 </td>
   <td style="text-align:right;"> 340 </td>
   <td style="text-align:right;"> 0 </td>
  </tr>
  <tr>
   <td style="text-align:right;"> 14 </td>
   <td style="text-align:right;"> 896906:896914 </td>
   <td style="text-align:right;"> 340 </td>
   <td style="text-align:right;"> 0 </td>
  </tr>
  <tr>
   <td style="text-align:right;"> 15 </td>
   <td style="text-align:right;"> 897254:897262 </td>
   <td style="text-align:right;"> 340 </td>
   <td style="text-align:right;"> 0 </td>
  </tr>
</tbody>
</table>



Para entender un poco más esta secuencia hicimos un blast el cual arrojó que dicha repetcion de 340 nucleotidos era un motivo **SWM_repeat** el cual se encuentra aaltamente repetido en una proteína de la superficie celular requerida para la movilidad (**QNJ16559.1**). Además, segun la secuencia de aminoácidos, el palindromo esta segmentado en 2 partes. La primera mitad corresponde a **TTA ACG** en el primer y segundo codon  y **CG** en el ultimo codón del motivo **SWM_repeat** (Figure \@ref(fig:FIG10)).

<div class="figure" style="text-align: center">
<img src="Clados/Clado_A18-40/SWM_repeat_AA.png" alt="**Traducción del motivo SWM\_repeat.** En esta imagen se muestra la secuencia traducida del motivo SWM\_repeat la cual se encuentra en la especie Synechococcus sp A18-40." width="100%" />
<p class="caption">(\#fig:FIG10)**Traducción del motivo SWM\_repeat.** En esta imagen se muestra la secuencia traducida del motivo SWM\_repeat la cual se encuentra en la especie Synechococcus sp A18-40.</p>
</div>
Debido a que estos 101 sitios solo estaban presentes en la especie  **Synechococcus sp A18−40** se concluyo que la tasa elevada de OE solo se debia a que tenia presente dicha proteína, ya que si se omitía del conteo, las tasas OE eran homogeneas en todo el clado.

### ATGCGCAT

El segundo palíndromo con la tasa OE mas alta es **ATGCGCAT** en la misma especie **Synechococcus sp A18−40** (Tabla \@ref(tab:TAB3)) con un conteo de 135 sitios para dicho palindromo mientras que en las demás especies oscila entre 22 y 31 sitios. 


  
<table>
<caption>(\#tab:TAB3)Conteo del palíndromo ATGCGCAT en el clado A18-40</caption>
 <thead>
  <tr>
   <th style="text-align:left;"> spp </th>
   <th style="text-align:left;"> palindrome </th>
   <th style="text-align:right;"> obs </th>
   <th style="text-align:right;"> markov3 </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> Synechococcus_sp_A18-46_1 </td>
   <td style="text-align:left;"> ATGCGCAT </td>
   <td style="text-align:right;"> 31 </td>
   <td style="text-align:right;"> 21.14 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Synechococcus_sp_A18-40 </td>
   <td style="text-align:left;"> ATGCGCAT </td>
   <td style="text-align:right;"> 135 </td>
   <td style="text-align:right;"> 21.67 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Synechococcus_sp_WH_8103 </td>
   <td style="text-align:left;"> ATGCGCAT </td>
   <td style="text-align:right;"> 28 </td>
   <td style="text-align:right;"> 20.73 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Synechococcus_sp_RS9915 </td>
   <td style="text-align:left;"> ATGCGCAT </td>
   <td style="text-align:right;"> 31 </td>
   <td style="text-align:right;"> 20.50 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Synechococcus_sp_A15-28 </td>
   <td style="text-align:left;"> ATGCGCAT </td>
   <td style="text-align:right;"> 22 </td>
   <td style="text-align:right;"> 19.56 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Parasynechococcus_marenigrum_WH_8102 </td>
   <td style="text-align:left;"> ATGCGCAT </td>
   <td style="text-align:right;"> 28 </td>
   <td style="text-align:right;"> 20.60 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Synechococcus_sp_A15-24 </td>
   <td style="text-align:left;"> ATGCGCAT </td>
   <td style="text-align:right;"> 30 </td>
   <td style="text-align:right;"> 19.93 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Synechococcus_sp_BOUM118 </td>
   <td style="text-align:left;"> ATGCGCAT </td>
   <td style="text-align:right;"> 29 </td>
   <td style="text-align:right;"> 19.87 </td>
  </tr>
</tbody>
</table>



#### Ubicación de sitios ATGCGCAT en los genomas del clado A18-40
Al analizar la distribución del palindromo en el genoma pudimos observar que, al igual que CGTTAACG, había 101 sitios que se encontraban entre repeticiones de 340 nuclótidos (columna 5 de la Tabla \@ref(tab:TAB4))



<table>
<caption>(\#tab:TAB4)**Ubicación de los sitios ATGCGCAT.** La primera columna se muestra el numero de sitio. La segunda columna muestra el intervalo en el que se encuentra el palíndromo. La tercera columna muestra a cuantos nucléotidos se encuentra el ultimo sitio. La cuarta columna indica la diferencia entre la distancia del ultimo palindromo y la distancia del siguiente.</caption>
 <thead>
  <tr>
   <th style="text-align:right;"> SiteNumber </th>
   <th style="text-align:right;"> Interval </th>
   <th style="text-align:right;"> Dist2NextPal </th>
   <th style="text-align:right;"> DifBetDist </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 24777:24785 </td>
   <td style="text-align:right;"> 24777 </td>
   <td style="text-align:right;"> -24777 </td>
  </tr>
  <tr>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 91754:91762 </td>
   <td style="text-align:right;"> 66969 </td>
   <td style="text-align:right;"> -42192 </td>
  </tr>
  <tr>
   <td style="text-align:right;"> 3 </td>
   <td style="text-align:right;"> 106816:106824 </td>
   <td style="text-align:right;"> 15054 </td>
   <td style="text-align:right;"> 51915 </td>
  </tr>
  <tr>
   <td style="text-align:right;"> 4 </td>
   <td style="text-align:right;"> 303729:303737 </td>
   <td style="text-align:right;"> 196905 </td>
   <td style="text-align:right;"> -181851 </td>
  </tr>
  <tr>
   <td style="text-align:right;"> 5 </td>
   <td style="text-align:right;"> 415357:415365 </td>
   <td style="text-align:right;"> 111620 </td>
   <td style="text-align:right;"> 85285 </td>
  </tr>
  <tr>
   <td style="text-align:right;"> 6 </td>
   <td style="text-align:right;"> 420623:420631 </td>
   <td style="text-align:right;"> 5258 </td>
   <td style="text-align:right;"> 106362 </td>
  </tr>
  <tr>
   <td style="text-align:right;"> 7 </td>
   <td style="text-align:right;"> 540975:540983 </td>
   <td style="text-align:right;"> 120344 </td>
   <td style="text-align:right;"> -115086 </td>
  </tr>
  <tr>
   <td style="text-align:right;"> 8 </td>
   <td style="text-align:right;"> 540987:540995 </td>
   <td style="text-align:right;"> 4 </td>
   <td style="text-align:right;"> 120340 </td>
  </tr>
  <tr>
   <td style="text-align:right;"> 9 </td>
   <td style="text-align:right;"> 571743:571751 </td>
   <td style="text-align:right;"> 30748 </td>
   <td style="text-align:right;"> -30744 </td>
  </tr>
  <tr>
   <td style="text-align:right;"> 10 </td>
   <td style="text-align:right;"> 694673:694681 </td>
   <td style="text-align:right;"> 122922 </td>
   <td style="text-align:right;"> -92174 </td>
  </tr>
  <tr>
   <td style="text-align:right;"> 11 </td>
   <td style="text-align:right;"> 708637:708645 </td>
   <td style="text-align:right;"> 13956 </td>
   <td style="text-align:right;"> 108966 </td>
  </tr>
  <tr>
   <td style="text-align:right;"> 12 </td>
   <td style="text-align:right;"> 760389:760397 </td>
   <td style="text-align:right;"> 51744 </td>
   <td style="text-align:right;"> -37788 </td>
  </tr>
  <tr>
   <td style="text-align:right;"> 13 </td>
   <td style="text-align:right;"> 867192:867200 </td>
   <td style="text-align:right;"> 106795 </td>
   <td style="text-align:right;"> -55051 </td>
  </tr>
  <tr>
   <td style="text-align:right;"> 14 </td>
   <td style="text-align:right;"> 893783:893791 </td>
   <td style="text-align:right;"> 26583 </td>
   <td style="text-align:right;"> 80212 </td>
  </tr>
  <tr>
   <td style="text-align:right;"> 15 </td>
   <td style="text-align:right;"> 894131:894139 </td>
   <td style="text-align:right;"> 340 </td>
   <td style="text-align:right;"> 26243 </td>
  </tr>
</tbody>
</table>



Al hacer un blast el cual arrojó que dicha repetcion de 340 nucleotidos era el mismo motivo **SWM_repeat** de la misma proteína de la superficie celular requerida para la movilidad (**QNJ16559.1**). Al revisar la ubicacion del palíndromo pudimos notar que las mismas cantidades de los palíndromos **ATGCGCAT** y **CGTTAAGC** se debe a que estan a un nucleótido de distancia el uno del otro (Figure \@ref(fig:FIG11)). 

<div class="figure" style="text-align: center">
<img src="Clados/Clado_A18-40/SWM_3.png" alt="**Traducción del motivo SWM\_repeat.** En esta imagen se muestra la secuencia traducida del motivo SWM\_repeat. Subrayado en amarillo y rojo se señalan los palíndromos **CGTTAAGC** y **ATGCGCAT** respectivamente, los cuales se encuentran a un nucleotido de distancia." width="1253" />
<p class="caption">(\#fig:FIG11)**Traducción del motivo SWM\_repeat.** En esta imagen se muestra la secuencia traducida del motivo SWM\_repeat. Subrayado en amarillo y rojo se señalan los palíndromos **CGTTAAGC** y **ATGCGCAT** respectivamente, los cuales se encuentran a un nucleotido de distancia.</p>
</div>

