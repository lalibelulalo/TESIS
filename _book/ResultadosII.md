# Sección II {-}

## Clado Calothrix

El clado calothrix contiene 6 especies  y es de interes ya que segun la filogenia estan estrechamente relacionadas y muestra un cambio en el palindromo mas abundante, pasando de **GCGATCGC** a **TGGCGCCA** (Figure \@ref(fig:FIG12)).
<div class="figure" style="text-align: center">
<img src="Clados/Callothrix_clade/figures/Calothrix_Octanuc_FrecObs_sel32_filogenia_HIG.png" alt="**Filogenia anotada del clado Calothrix.** En esta imagen se muestra un cambio abrupto en la Frecuencia observada de **GCGCATCGC** en las especies NIES-4105, NIES-4071 y PCC\_7716." width="100%" />
<p class="caption">(\#fig:FIG12)**Filogenia anotada del clado Calothrix.** En esta imagen se muestra un cambio abrupto en la Frecuencia observada de **GCGCATCGC** en las especies NIES-4105, NIES-4071 y PCC\_7716.</p>
</div>

### Conjunto de sitios HIP1 usando la especie 336-3 como referencia

#### Red de transiciones 

Para hacer mas visual la reconstrucción, construimos una red de las transiciones entre los estados ancestrales. Esto lo hicimos en r usando la función ```Create_Transition_Table()```:


```r
source("ASR_Orth_Functions/NodeAndEdges.R")
Nodes.Edges <- Create_Transition_Table_No_Fit(SitesTable = "Clados/Callothrix_clade/PALINDROMES/GCGATCGC/336-3/Orthologues_Palindrome_sites.txt",
                                EvolutionModel = "F81",
                                Method = "bayes",
                                Phylogeny = "Clados/Callothrix_clade/SpeciesTree_rooted.txt",
                                OrthoPath = "Clados/Callothrix_clade/PALINDROMES/GCGATCGC/336-3/Only_ORTHOLOGUES/")

```

Posteriormente creamos la red usando la función ```Create_Network()```:



y visualizamos dicha red .



Para visualizar la red usamos la paqueteria ```networkD3```. Hicimos 2 figuras, la (Figura \@ref(fig:FIG13)) muestra la red como una conexión de nodos a través de vertices con un grosor proporcional al numero de veces que ocurrió cada transición. En dicha red podemos ver algunos nodos con bordes muy gruesos como **GCGATTGC**, **GCAATTGC**, **GCTATCGC**, **GCTATTGC** (Tabla \@ref(tab:TAB6)).



<div class="figure" style="text-align: center">
<img src="ResultadosII_files/figure-epub3/FIG13-1.png" alt="**Red de todas las transiciones del clado Calothrix.** En esta imagen se muestra una red en la que cada nodo es un octanucleótido el cual esta unido a otro nodo por un vertice. Dicho vertice tiene un grosor proporcional al numero de veces que dicha transición ocurrió en la reconstrucción ancestral de sitios."  />
<p class="caption">(\#fig:FIG13)**Red de todas las transiciones del clado Calothrix.** En esta imagen se muestra una red en la que cada nodo es un octanucleótido el cual esta unido a otro nodo por un vertice. Dicho vertice tiene un grosor proporcional al numero de veces que dicha transición ocurrió en la reconstrucción ancestral de sitios.</p>
</div>

En la (Figura \@ref(fig:FIG14)) podemos ver las transiciones de una forma mas ordenada, con el numero de ocurrencias y la dirección en la que ocurrieron.
<div class="figure" style="text-align: center">
<img src="ResultadosII_files/figure-epub3/FIG14-1.png" alt="**Red de todas las transiciones del clado Calothrix.** En esta imagen se muestra la red  de la Figura \@ref(fig:FIG13) de una forma mas visual y con el numero de veces que ocurrio cada transición, asi como la dirección en la que ocurrió."  />
<p class="caption">(\#fig:FIG14)**Red de todas las transiciones del clado Calothrix.** En esta imagen se muestra la red  de la Figura \@ref(fig:FIG13) de una forma mas visual y con el numero de veces que ocurrio cada transición, asi como la dirección en la que ocurrió.</p>
</div>

<div style="border: 1px solid #ddd; padding: 0px; overflow-y: scroll; height:200px; overflow-x: scroll; width:500px; "><table class=" lightable-paper table table-striped table-condensed" style='font-family: "Arial Narrow", arial, helvetica, sans-serif; margin-left: auto; margin-right: auto; font-size: 11px; width: auto !important; margin-left: auto; margin-right: auto;'>
<caption style="font-size: initial !important;">(\#tab:TAB6)Transiciones entre nodos.</caption>
 <thead>
  <tr>
   <th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;position: sticky; top:0; background-color: #FFFFFF;"> from </th>
   <th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;position: sticky; top:0; background-color: #FFFFFF;"> to </th>
   <th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;position: sticky; top:0; background-color: #FFFFFF;"> weight </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: black !important;">GCGATTGC</span> </td>
   <td style="text-align:left;font-weight: bold;vertical-align: top !important;" rowspan="11"> <span style="     color: red !important;">GCGATCGC</span> </td>
   <td style="text-align:left;"> <span style="display: block; padding: 0 4px; border-radius: 4px; background-color: #ffa500">139</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: black !important;">GCTATTGC</span> </td>
   
   <td style="text-align:left;"> <span style="display: block; padding: 0 4px; border-radius: 4px; background-color: #ffc55c">91</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: blue !important;">GCAATTGC</span> </td>
   
   <td style="text-align:left;"> <span style="display: block; padding: 0 4px; border-radius: 4px; background-color: #ffc55c">91</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: black !important;">GCAATCGC</span> </td>
   
   <td style="text-align:left;"> <span style="display: block; padding: 0 4px; border-radius: 4px; background-color: #ffc760">89</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: black !important;">GCTATCGC</span> </td>
   
   <td style="text-align:left;"> <span style="display: block; padding: 0 4px; border-radius: 4px; background-color: #ffcd73">79</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: black !important;">GCGACCGC</span> </td>
   
   <td style="text-align:left;"> <span style="display: block; padding: 0 4px; border-radius: 4px; background-color: #ffd68d">66</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: black !important;">GCGATAGC</span> </td>
   
   <td style="text-align:left;"> <span style="display: block; padding: 0 4px; border-radius: 4px; background-color: #ffd68d">66</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: black !important;">--------</span> </td>
   
   <td style="text-align:left;"> <span style="display: block; padding: 0 4px; border-radius: 4px; background-color: #ffdb9a">59</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: black !important;">GCCATCGC</span> </td>
   
   <td style="text-align:left;"> <span style="display: block; padding: 0 4px; border-radius: 4px; background-color: #ffdc9e">57</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: black !important;">GTGATCGC</span> </td>
   
   <td style="text-align:left;"> <span style="display: block; padding: 0 4px; border-radius: 4px; background-color: #ffe3b1">47</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: black !important;">GCGACCGT</span> </td>
   
   <td style="text-align:left;"> <span style="display: block; padding: 0 4px; border-radius: 4px; background-color: #ffe4b3">46</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: blue !important;">GCAATTGC</span> </td>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: black !important;">GCGATTGC</span> </td>
   <td style="text-align:left;"> <span style="display: block; padding: 0 4px; border-radius: 4px; background-color: #ffe5b7">44</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: black !important;">GCCATTGC</span> </td>
   <td style="text-align:left;font-weight: bold;vertical-align: top !important;" rowspan="3"> <span style="     color: red !important;">GCGATCGC</span> </td>
   <td style="text-align:left;"> <span style="display: block; padding: 0 4px; border-radius: 4px; background-color: #ffe7bb">42</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: black !important;">GCAATAGC</span> </td>
   
   <td style="text-align:left;"> <span style="display: block; padding: 0 4px; border-radius: 4px; background-color: #ffe9c1">39</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: black !important;">GCGATCGT</span> </td>
   
   <td style="text-align:left;"> <span style="display: block; padding: 0 4px; border-radius: 4px; background-color: #ffebc8">35</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;font-weight: bold;vertical-align: top !important;" rowspan="2"> <span style="     color: black !important;">GCTATTGC</span> </td>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: black !important;">GCGATTGC</span> </td>
   <td style="text-align:left;"> <span style="display: block; padding: 0 4px; border-radius: 4px; background-color: #ffedcc">33</span> </td>
  </tr>
  <tr>
   
   <td style="text-align:left;font-weight: bold;"> <span style="     color: black !important;">GCTATCGC</span> </td>
   <td style="text-align:left;"> <span style="display: block; padding: 0 4px; border-radius: 4px; background-color: #ffedce">32</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: black !important;">GCGATAGA</span> </td>
   <td style="text-align:left;font-weight: bold;vertical-align: top !important;" rowspan="5"> <span style="     color: red !important;">GCGATCGC</span> </td>
   <td style="text-align:left;"> <span style="display: block; padding: 0 4px; border-radius: 4px; background-color: #ffedce">32</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: black !important;">GCGTTCGC</span> </td>
   
   <td style="text-align:left;"> <span style="display: block; padding: 0 4px; border-radius: 4px; background-color: #fff1d8">27</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: black !important;">GTGACCGC</span> </td>
   
   <td style="text-align:left;"> <span style="display: block; padding: 0 4px; border-radius: 4px; background-color: #fff1d8">27</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: black !important;">GCGATCGA</span> </td>
   
   <td style="text-align:left;"> <span style="display: block; padding: 0 4px; border-radius: 4px; background-color: #fff2da">26</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: black !important;">GCGGTCGC</span> </td>
   
   <td style="text-align:left;"> <span style="display: block; padding: 0 4px; border-radius: 4px; background-color: #fff2dc">25</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: black !important;">GCGACCGT</span> </td>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: black !important;">GCGACCGC</span> </td>
   <td style="text-align:left;"> <span style="display: block; padding: 0 4px; border-radius: 4px; background-color: #fff4e0">23</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: red !important;">GCGATCGC</span> </td>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: black !important;">GCGATTGC</span> </td>
   <td style="text-align:left;"> <span style="display: block; padding: 0 4px; border-radius: 4px; background-color: #fff4e0">23</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: blue !important;">GCAATTGC</span> </td>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: black !important;">GCAATCGC</span> </td>
   <td style="text-align:left;"> <span style="display: block; padding: 0 4px; border-radius: 4px; background-color: #fff4e0">23</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: black !important;">GTGACCGT</span> </td>
   <td style="text-align:left;font-weight: bold;vertical-align: top !important;" rowspan="4"> <span style="     color: red !important;">GCGATCGC</span> </td>
   <td style="text-align:left;"> <span style="display: block; padding: 0 4px; border-radius: 4px; background-color: #fff4e2">22</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: black !important;">GCGAACGC</span> </td>
   
   <td style="text-align:left;"> <span style="display: block; padding: 0 4px; border-radius: 4px; background-color: #fff5e3">21</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: black !important;">GCCATAGC</span> </td>
   
   <td style="text-align:left;"> <span style="display: block; padding: 0 4px; border-radius: 4px; background-color: #fff6e5">20</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: black !important;">ACGATCGC</span> </td>
   
   <td style="text-align:left;"> <span style="display: block; padding: 0 4px; border-radius: 4px; background-color: #fff6e5">20</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: black !important;">GCCATTGC</span> </td>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: black !important;">GCCATCGC</span> </td>
   <td style="text-align:left;"> <span style="display: block; padding: 0 4px; border-radius: 4px; background-color: #fff7e9">18</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: black !important;">GTGATAGA</span> </td>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: red !important;">GCGATCGC</span> </td>
   <td style="text-align:left;"> <span style="display: block; padding: 0 4px; border-radius: 4px; background-color: #fff7e9">18</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: blue !important;">GCAATTGC</span> </td>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: black !important;">GCTATTGC</span> </td>
   <td style="text-align:left;"> <span style="display: block; padding: 0 4px; border-radius: 4px; background-color: #fff7e9">18</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: black !important;">GCTATTGC</span> </td>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: blue !important;">GCAATTGC</span> </td>
   <td style="text-align:left;"> <span style="display: block; padding: 0 4px; border-radius: 4px; background-color: #fff7e9">18</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: black !important;">GCCATTGC</span> </td>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: black !important;">GCGATTGC</span> </td>
   <td style="text-align:left;"> <span style="display: block; padding: 0 4px; border-radius: 4px; background-color: #fff8eb">17</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: black !important;">GCAATAGC</span> </td>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: black !important;">GCAATCGC</span> </td>
   <td style="text-align:left;"> <span style="display: block; padding: 0 4px; border-radius: 4px; background-color: #fff8ed">16</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: black !important;">GCGATCAC</span> </td>
   <td style="text-align:left;font-weight: bold;vertical-align: top !important;" rowspan="2"> <span style="     color: red !important;">GCGATCGC</span> </td>
   <td style="text-align:left;"> <span style="display: block; padding: 0 4px; border-radius: 4px; background-color: #fff8ed">16</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: black !important;">GCTATAGC</span> </td>
   
   <td style="text-align:left;"> <span style="display: block; padding: 0 4px; border-radius: 4px; background-color: #fff8ed">16</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: black !important;">GCAATAGC</span> </td>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: black !important;">GCGATAGC</span> </td>
   <td style="text-align:left;"> <span style="display: block; padding: 0 4px; border-radius: 4px; background-color: #fff9ef">15</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: black !important;">GCGCTCGC</span> </td>
   <td style="text-align:left;font-weight: bold;vertical-align: top !important;" rowspan="3"> <span style="     color: red !important;">GCGATCGC</span> </td>
   <td style="text-align:left;"> <span style="display: block; padding: 0 4px; border-radius: 4px; background-color: #fff9ef">15</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: black !important;">GCGAGCGC</span> </td>
   
   <td style="text-align:left;"> <span style="display: block; padding: 0 4px; border-radius: 4px; background-color: #fffaf1">14</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: black !important;">GTGATTGC</span> </td>
   
   <td style="text-align:left;"> <span style="display: block; padding: 0 4px; border-radius: 4px; background-color: #fffaf1">14</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: black !important;">GCGACAGA</span> </td>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: black !important;">GCGATAGA</span> </td>
   <td style="text-align:left;"> <span style="display: block; padding: 0 4px; border-radius: 4px; background-color: #fffaf1">14</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: black !important;">GAAATCGC</span> </td>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: red !important;">GCGATCGC</span> </td>
   <td style="text-align:left;"> <span style="display: block; padding: 0 4px; border-radius: 4px; background-color: #fffaf1">14</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: black !important;">GCCATTGC</span> </td>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: black !important;">GCTATTGC</span> </td>
   <td style="text-align:left;"> <span style="display: block; padding: 0 4px; border-radius: 4px; background-color: #fffaf3">13</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: black !important;">GTGACCGT</span> </td>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: black !important;">GCGACCGT</span> </td>
   <td style="text-align:left;"> <span style="display: block; padding: 0 4px; border-radius: 4px; background-color: #fffaf3">13</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: black !important;">GCGGTTGC</span> </td>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: red !important;">GCGATCGC</span> </td>
   <td style="text-align:left;"> <span style="display: block; padding: 0 4px; border-radius: 4px; background-color: #fffaf3">13</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: red !important;">GCGATCGC</span> </td>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: black !important;">GTGATCGC</span> </td>
   <td style="text-align:left;"> <span style="display: block; padding: 0 4px; border-radius: 4px; background-color: #fffaf3">13</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: black !important;">GAGATCGC</span> </td>
   <td style="text-align:left;font-weight: bold;vertical-align: top !important;" rowspan="3"> <span style="     color: red !important;">GCGATCGC</span> </td>
   <td style="text-align:left;"> <span style="display: block; padding: 0 4px; border-radius: 4px; background-color: #fffaf3">13</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: black !important;">ACGTTCGC</span> </td>
   
   <td style="text-align:left;"> <span style="display: block; padding: 0 4px; border-radius: 4px; background-color: #fffaf3">13</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: black !important;">GCGATAAA</span> </td>
   
   <td style="text-align:left;"> <span style="display: block; padding: 0 4px; border-radius: 4px; background-color: #fffaf3">13</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: black !important;">GCAATAGC</span> </td>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: blue !important;">GCAATTGC</span> </td>
   <td style="text-align:left;"> <span style="display: block; padding: 0 4px; border-radius: 4px; background-color: #fffbf5">12</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: black !important;">GTGACCGC</span> </td>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: black !important;">GCGACCGC</span> </td>
   <td style="text-align:left;"> <span style="display: block; padding: 0 4px; border-radius: 4px; background-color: #fffbf5">12</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: black !important;">ACAATTGC</span> </td>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: red !important;">GCGATCGC</span> </td>
   <td style="text-align:left;"> <span style="display: block; padding: 0 4px; border-radius: 4px; background-color: #fffbf5">12</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: black !important;">GCTATAGC</span> </td>
   <td style="text-align:left;font-weight: bold;vertical-align: top !important;" rowspan="2"> <span style="     color: black !important;">GCGATAGC</span> </td>
   <td style="text-align:left;"> <span style="display: block; padding: 0 4px; border-radius: 4px; background-color: #fffbf5">12</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: black !important;">GCCATAGC</span> </td>
   
   <td style="text-align:left;"> <span style="display: block; padding: 0 4px; border-radius: 4px; background-color: #fffbf5">12</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: black !important;">GCAGTTGC</span> </td>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: red !important;">GCGATCGC</span> </td>
   <td style="text-align:left;"> <span style="display: block; padding: 0 4px; border-radius: 4px; background-color: #fffcf7">11</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: black !important;">GTGATAGA</span> </td>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: black !important;">GCGATAGA</span> </td>
   <td style="text-align:left;"> <span style="display: block; padding: 0 4px; border-radius: 4px; background-color: #fffcf7">11</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;font-weight: bold;vertical-align: top !important;" rowspan="2"> <span style="     color: red !important;">GCGATCGC</span> </td>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: black !important;">ACGATCGC</span> </td>
   <td style="text-align:left;"> <span style="display: block; padding: 0 4px; border-radius: 4px; background-color: #fffcf7">11</span> </td>
  </tr>
  <tr>
   
   <td style="text-align:left;font-weight: bold;"> <span style="     color: black !important;">GCGACCGC</span> </td>
   <td style="text-align:left;"> <span style="display: block; padding: 0 4px; border-radius: 4px; background-color: #fffcf7">11</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: black !important;">GTGACCGT</span> </td>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: black !important;">GTGACCGC</span> </td>
   <td style="text-align:left;"> <span style="display: block; padding: 0 4px; border-radius: 4px; background-color: #fffcf7">11</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: black !important;">GAGACCGC</span> </td>
   <td style="text-align:left;font-weight: bold;vertical-align: top !important;" rowspan="2"> <span style="     color: red !important;">GCGATCGC</span> </td>
   <td style="text-align:left;"> <span style="display: block; padding: 0 4px; border-radius: 4px; background-color: #fffcf9">10</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: black !important;">GCGACCGA</span> </td>
   
   <td style="text-align:left;"> <span style="display: block; padding: 0 4px; border-radius: 4px; background-color: #fffcf9">10</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: black !important;">GCTATCGC</span> </td>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: black !important;">GCCATCGC</span> </td>
   <td style="text-align:left;"> <span style="display: block; padding: 0 4px; border-radius: 4px; background-color: #fffcf9">10</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: black !important;">GCGATAAC</span> </td>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: red !important;">GCGATCGC</span> </td>
   <td style="text-align:left;"> <span style="display: block; padding: 0 4px; border-radius: 4px; background-color: #fffcf9">10</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: black !important;">GCTATAGC</span> </td>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: black !important;">GCTATTGC</span> </td>
   <td style="text-align:left;"> <span style="display: block; padding: 0 4px; border-radius: 4px; background-color: #fffcf9">10</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: red !important;">GCGATCGC</span> </td>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: black !important;">GCAATCGC</span> </td>
   <td style="text-align:left;"> <span style="display: block; padding: 0 4px; border-radius: 4px; background-color: #fffcf9">10</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: black !important;">GCGATTGT</span> </td>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: red !important;">GCGATCGC</span> </td>
   <td style="text-align:left;"> <span style="display: block; padding: 0 4px; border-radius: 4px; background-color: #fffcf9">10</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: black !important;">GCTATAGC</span> </td>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: black !important;">GCTATCGC</span> </td>
   <td style="text-align:left;"> <span style="display: block; padding: 0 4px; border-radius: 4px; background-color: #fffdfb">9</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: black !important;">GCAATCGT</span> </td>
   <td style="text-align:left;font-weight: bold;vertical-align: top !important;" rowspan="2"> <span style="     color: red !important;">GCGATCGC</span> </td>
   <td style="text-align:left;"> <span style="display: block; padding: 0 4px; border-radius: 4px; background-color: #fffdfb">9</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: black !important;">GCGATTGA</span> </td>
   
   <td style="text-align:left;"> <span style="display: block; padding: 0 4px; border-radius: 4px; background-color: #fffdfb">9</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: red !important;">GCGATCGC</span> </td>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: black !important;">GCGATAGC</span> </td>
   <td style="text-align:left;"> <span style="display: block; padding: 0 4px; border-radius: 4px; background-color: #fffdfb">9</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: black !important;">GCGATGGC</span> </td>
   <td style="text-align:left;font-weight: bold;vertical-align: top !important;" rowspan="6"> <span style="     color: red !important;">GCGATCGC</span> </td>
   <td style="text-align:left;"> <span style="display: block; padding: 0 4px; border-radius: 4px; background-color: #fffdfb">9</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: black !important;">GCAATATC</span> </td>
   
   <td style="text-align:left;"> <span style="display: block; padding: 0 4px; border-radius: 4px; background-color: #fffdfb">9</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: black !important;">GCGTTCTC</span> </td>
   
   <td style="text-align:left;"> <span style="display: block; padding: 0 4px; border-radius: 4px; background-color: #fffdfb">9</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: black !important;">GTGATAGC</span> </td>
   
   <td style="text-align:left;"> <span style="display: block; padding: 0 4px; border-radius: 4px; background-color: #fffdfb">9</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: black !important;">GCAATTTC</span> </td>
   
   <td style="text-align:left;"> <span style="display: block; padding: 0 4px; border-radius: 4px; background-color: #fffdfb">9</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: black !important;">GCGAACGT</span> </td>
   
   <td style="text-align:left;"> <span style="display: block; padding: 0 4px; border-radius: 4px; background-color: #fffdfb">9</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: black !important;">GCAGTTGC</span> </td>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: black !important;">GCGGTTGC</span> </td>
   <td style="text-align:left;"> <span style="display: block; padding: 0 4px; border-radius: 4px; background-color: #fffdfb">9</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: black !important;">GCTGTTGC</span> </td>
   <td style="text-align:left;font-weight: bold;vertical-align: top !important;" rowspan="4"> <span style="     color: red !important;">GCGATCGC</span> </td>
   <td style="text-align:left;"> <span style="display: block; padding: 0 4px; border-radius: 4px; background-color: #fffefd">8</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: black !important;">GAGACCGT</span> </td>
   
   <td style="text-align:left;"> <span style="display: block; padding: 0 4px; border-radius: 4px; background-color: #fffefd">8</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: black !important;">GTAACCGC</span> </td>
   
   <td style="text-align:left;"> <span style="display: block; padding: 0 4px; border-radius: 4px; background-color: #fffefd">8</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: black !important;">GAGATAGA</span> </td>
   
   <td style="text-align:left;"> <span style="display: block; padding: 0 4px; border-radius: 4px; background-color: #fffefd">8</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: black !important;">GCAGTTGC</span> </td>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: blue !important;">GCAATTGC</span> </td>
   <td style="text-align:left;"> <span style="display: block; padding: 0 4px; border-radius: 4px; background-color: #fffefd">8</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: black !important;">GAAATTGC</span> </td>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: black !important;">GAAATCGC</span> </td>
   <td style="text-align:left;"> <span style="display: block; padding: 0 4px; border-radius: 4px; background-color: #fffefd">8</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: black !important;">GCTGTTGC</span> </td>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: black !important;">GCTATTGC</span> </td>
   <td style="text-align:left;"> <span style="display: block; padding: 0 4px; border-radius: 4px; background-color: #fffefd">8</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: black !important;">GCGATCAA</span> </td>
   <td style="text-align:left;font-weight: bold;vertical-align: top !important;" rowspan="2"> <span style="     color: red !important;">GCGATCGC</span> </td>
   <td style="text-align:left;"> <span style="display: block; padding: 0 4px; border-radius: 4px; background-color: #fffefd">8</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: black !important;">GCTTTAGC</span> </td>
   
   <td style="text-align:left;"> <span style="display: block; padding: 0 4px; border-radius: 4px; background-color: #fffefd">8</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: blue !important;">GCAATTGC</span> </td>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: black !important;">GCAATAGC</span> </td>
   <td style="text-align:left;"> <span style="display: block; padding: 0 4px; border-radius: 4px; background-color: #fffefd">8</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: black !important;">GTGATCGT</span> </td>
   <td style="text-align:left;font-weight: bold;vertical-align: top !important;" rowspan="4"> <span style="     color: red !important;">GCGATCGC</span> </td>
   <td style="text-align:left;"> <span style="display: block; padding: 0 4px; border-radius: 4px; background-color: #ffffff">7</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: black !important;">GTGACCGA</span> </td>
   
   <td style="text-align:left;"> <span style="display: block; padding: 0 4px; border-radius: 4px; background-color: #ffffff">7</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: black !important;">GTAATCGC</span> </td>
   
   <td style="text-align:left;"> <span style="display: block; padding: 0 4px; border-radius: 4px; background-color: #ffffff">7</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: black !important;">GCGGTAGC</span> </td>
   
   <td style="text-align:left;"> <span style="display: block; padding: 0 4px; border-radius: 4px; background-color: #ffffff">7</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: black !important;">GCGGTTGC</span> </td>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: black !important;">GCGATTGC</span> </td>
   <td style="text-align:left;"> <span style="display: block; padding: 0 4px; border-radius: 4px; background-color: #ffffff">7</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: red !important;">GCGATCGC</span> </td>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: black !important;">GCGAACGC</span> </td>
   <td style="text-align:left;"> <span style="display: block; padding: 0 4px; border-radius: 4px; background-color: #ffffff">7</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: black !important;">GCAATTGA</span> </td>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: black !important;">GCGATTGA</span> </td>
   <td style="text-align:left;"> <span style="display: block; padding: 0 4px; border-radius: 4px; background-color: #ffffff">7</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: black !important;">GCAATTAC</span> </td>
   <td style="text-align:left;font-weight: bold;vertical-align: top !important;" rowspan="3"> <span style="     color: red !important;">GCGATCGC</span> </td>
   <td style="text-align:left;"> <span style="display: block; padding: 0 4px; border-radius: 4px; background-color: #ffffff">7</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: black !important;">GCGATCGG</span> </td>
   
   <td style="text-align:left;"> <span style="display: block; padding: 0 4px; border-radius: 4px; background-color: #ffffff">7</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: black !important;">ACTATCGC</span> </td>
   
   <td style="text-align:left;"> <span style="display: block; padding: 0 4px; border-radius: 4px; background-color: #ffffff">7</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: black !important;">GCAATGGC</span> </td>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: black !important;">GCGATGGC</span> </td>
   <td style="text-align:left;"> <span style="display: block; padding: 0 4px; border-radius: 4px; background-color: #ffffff">7</span> </td>
  </tr>
</tbody>
</table></div>



#### Transiciones entre Nodo 9 y Nodo 10

Para entender más como es que se gana o se pierden los sitios palindrómicos revisamos la transición en tre los nodos 9 y 10. Esto es porque es esta transicion de nodos la que separa a los dos subclados entre los que hay una repentino cambio de abundancia de sitios palindrómicos (Figura \@ref(fig:FIG15)).

<div class="figure">
<img src="ResultadosII_files/figure-epub3/FIG15-1.png" alt="**Filogenia del clado Calotrix**. En rojo y azul se muestran los subclados unidos (en verde) por la transición entre los nodos 9 y 10. "  />
<p class="caption">(\#fig:FIG15)**Filogenia del clado Calotrix**. En rojo y azul se muestran los subclados unidos (en verde) por la transición entre los nodos 9 y 10. </p>
</div>

Para hacer esto filtramos los datos de la red para mostrar unicamente las transiciones que se dieron entre los nodos 9 y 10 e hicimos las mismas figuras.
En la (Figura \@ref(fig:FIG16)) se muestra la red como una conexión de nodos a través de vertices con un grosor proporcional al numero de veces que ocurrió cada transición. En la (Figura \@ref(fig:FIG17)) podemos ver las transiciones de una forma mas ordenada, con el numero de ocurrencias y la dirección en la que ocurrieron.



<div class="figure" style="text-align: center">
<img src="ResultadosII_files/figure-epub3/FIG16-1.png" alt="**Red de las transiciones entre los Nodos 9 y 10 del clado Calothrix.** En esta imagen se muestra una red en la que cada nodo es un octanucleótido el cual esta unido a otro nodo por un vertice. Dicho vertice tiene un grosor proporcional al numero de veces que dicha transición ocurrió en la reconstrucción ancestral de sitios."  />
<p class="caption">(\#fig:FIG16)**Red de las transiciones entre los Nodos 9 y 10 del clado Calothrix.** En esta imagen se muestra una red en la que cada nodo es un octanucleótido el cual esta unido a otro nodo por un vertice. Dicho vertice tiene un grosor proporcional al numero de veces que dicha transición ocurrió en la reconstrucción ancestral de sitios.</p>
</div>

<div class="figure" style="text-align: center">
<img src="ResultadosII_files/figure-epub3/FIG17-1.png" alt="**Red de las transiciones entre los nodos 9 y 10 del clado Calothrix.** En esta imagen se muestra la red  de la Figura \@ref(fig:FIG16) de una forma mas visual y con el numero de veces que ocurrio cada transición, asi como la dirección en la que ocurrió."  />
<p class="caption">(\#fig:FIG17)**Red de las transiciones entre los nodos 9 y 10 del clado Calothrix.** En esta imagen se muestra la red  de la Figura \@ref(fig:FIG16) de una forma mas visual y con el numero de veces que ocurrio cada transición, asi como la dirección en la que ocurrió.</p>
</div>

#### Mutaciones en los codones

Para entender como es que se van ganando o perdiendo los sitios palindrómicos hicimos un análisis del tipo mutaciones de los sitios. Esto lo hicimos viendo en que marco de lectura se encontraba cada nodo y revisando la secuencia de aminoacidos que codificaban. En la (Figura \@ref(fig:FIG18)) mostramos 3 gráficos que indican la abundancia de los peptidos codificados por los sitios palindrómicos de acuerdo al marco de lectura en el que se encuentran. En esta figura podemos observar que el marco de lectura es el que contiene la mayoria de los sitios

<div class="figure" style="text-align: center">
<img src="ResultadosII_files/figure-epub3/FIG18-1.png" alt="**Abundancia de peptidos por cada nodo segun el marco de lectura.**." width="110%" /><img src="ResultadosII_files/figure-epub3/FIG18-2.png" alt="**Abundancia de peptidos por cada nodo segun el marco de lectura.**." width="110%" /><img src="ResultadosII_files/figure-epub3/FIG18-3.png" alt="**Abundancia de peptidos por cada nodo segun el marco de lectura.**." width="110%" />
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
<img src="ResultadosII_files/figure-epub3/FIG19-1.png" alt="**Abundancia del tipo de sustitución por cada nodo segun el marco de lectura.**." width="110%" /><img src="ResultadosII_files/figure-epub3/FIG19-2.png" alt="**Abundancia del tipo de sustitución por cada nodo segun el marco de lectura.**." width="110%" /><img src="ResultadosII_files/figure-epub3/FIG19-3.png" alt="**Abundancia del tipo de sustitución por cada nodo segun el marco de lectura.**." width="110%" />
<p class="caption">(\#fig:FIG19)**Abundancia del tipo de sustitución por cada nodo segun el marco de lectura.**.</p>
</div>

#### Análisis de sitios en los cuales su ancestro era HIP1

Para tratar de entender como es que los sitios HIP1 se pierden hicimos un análisis unicamente en en las transiciones en las que el nodo ancestral tenia un sitio HIP1.

En la (Figura \@ref(fig:FIG20)) mostramos 3 gráficos que indican la frecuencia del tipo de sustituciones que hubo para estos casos para cada nodo en cada uno de los marcos de lectura.

<div class="figure" style="text-align: center">
<img src="ResultadosII_files/figure-epub3/FIG20-1.png" alt="**Abundancia del tipo de sustitución por cada nodo segun el marco de lectura. Unicamente para transiciones en los que el nodo ancestral era un sitio HIP1.**" width="110%" /><img src="ResultadosII_files/figure-epub3/FIG20-2.png" alt="**Abundancia del tipo de sustitución por cada nodo segun el marco de lectura. Unicamente para transiciones en los que el nodo ancestral era un sitio HIP1.**" width="110%" /><img src="ResultadosII_files/figure-epub3/FIG20-3.png" alt="**Abundancia del tipo de sustitución por cada nodo segun el marco de lectura. Unicamente para transiciones en los que el nodo ancestral era un sitio HIP1.**" width="110%" />
<p class="caption">(\#fig:FIG20)**Abundancia del tipo de sustitución por cada nodo segun el marco de lectura. Unicamente para transiciones en los que el nodo ancestral era un sitio HIP1.**</p>
</div>

En la (Figura \@ref(fig:FIG21)) mostramos 3 gráficos (uno por cada marco de lectura) que indican la frecuencia de las mutaciones en cada uno de los 8 nucleótidos del sitio HIP.

<div class="figure" style="text-align: center">
<img src="ResultadosII_files/figure-epub3/FIG21-1.png" alt="**Frecuencia de las mutaciones de cada nucleótido del sitio HIP para cada nodo segun el marco de lectura.**." width="110%" /><img src="ResultadosII_files/figure-epub3/FIG21-2.png" alt="**Frecuencia de las mutaciones de cada nucleótido del sitio HIP para cada nodo segun el marco de lectura.**." width="110%" /><img src="ResultadosII_files/figure-epub3/FIG21-3.png" alt="**Frecuencia de las mutaciones de cada nucleótido del sitio HIP para cada nodo segun el marco de lectura.**." width="110%" />
<p class="caption">(\#fig:FIG21)**Frecuencia de las mutaciones de cada nucleótido del sitio HIP para cada nodo segun el marco de lectura.**.</p>
</div>

En la (Figura \@ref(fig:FIG22)) mostramos 3 gráficos (uno por cada marco de lectura) que indican la frecuencia del tipo sustitucion de bases.

<div class="figure" style="text-align: center">
<img src="ResultadosII_files/figure-epub3/FIG22-1.png" alt="**Frecuencia del tipo de sustituciónes de base en los sitios HIP para cada marco de lectura**." width="110%" /><img src="ResultadosII_files/figure-epub3/FIG22-2.png" alt="**Frecuencia del tipo de sustituciónes de base en los sitios HIP para cada marco de lectura**." width="110%" /><img src="ResultadosII_files/figure-epub3/FIG22-3.png" alt="**Frecuencia del tipo de sustituciónes de base en los sitios HIP para cada marco de lectura**." width="110%" />
<p class="caption">(\#fig:FIG22)**Frecuencia del tipo de sustituciónes de base en los sitios HIP para cada marco de lectura**.</p>
</div>

#### Análisis de sitios en los cuales solo el nodo actual tiene HIP1

Para tratar de entender como es que los sitios HIP se ganan, hicimos un analisis unicamente en las transiciones en las que el nodo actual tenia un sitio HIP1. 

En la Figura \@ref(fig:FIG23) mostramos 3 gráficos (uno por cada marco de lectura) que indican la frecuencia del tipo de sustituciones que hubo para estos casos para cada nodo en cada uno de los marcos de lectura.

<div class="figure" style="text-align: center">
<img src="ResultadosII_files/figure-epub3/FIG23-1.png" alt="**Abundancia del tipo de sustitución por cada nodo segun el marco de lectura. Unicamente para transiciones en los que el nodo actual era un sitio HIP1.**." width="110%" /><img src="ResultadosII_files/figure-epub3/FIG23-2.png" alt="**Abundancia del tipo de sustitución por cada nodo segun el marco de lectura. Unicamente para transiciones en los que el nodo actual era un sitio HIP1.**." width="110%" /><img src="ResultadosII_files/figure-epub3/FIG23-3.png" alt="**Abundancia del tipo de sustitución por cada nodo segun el marco de lectura. Unicamente para transiciones en los que el nodo actual era un sitio HIP1.**." width="110%" />
<p class="caption">(\#fig:FIG23)**Abundancia del tipo de sustitución por cada nodo segun el marco de lectura. Unicamente para transiciones en los que el nodo actual era un sitio HIP1.**.</p>
</div>



### Conjunto de sitios HIP1 usando la especie NIES-3974 como referencia

#### Red de transiciones



visualizamos dicha red .



Para visualizar la red usamos la paqueteria ```networkD3```. Hicimos 2 figuras, la (Figura \@ref(fig:FIG44)) muestra la red como una conexión de nodos a través de vertices con un grosor proporcional al numero de veces que ocurrió cada transición.



<div class="figure" style="text-align: center">
<img src="ResultadosII_files/figure-epub3/FIG44-1.png" alt="**Red de todas las transiciones del clado Calothrix.** En esta imagen se muestra una red en la que cada nodo es un octanucleótido el cual esta unido a otro nodo por un vertice. Dicho vertice tiene un grosor proporcional al numero de veces que dicha transición ocurrió en la reconstrucción ancestral de sitios. Esta red es una represntación de aquellos sitios que se encuentran en la especie NIES-3974."  />
<p class="caption">(\#fig:FIG44)**Red de todas las transiciones del clado Calothrix.** En esta imagen se muestra una red en la que cada nodo es un octanucleótido el cual esta unido a otro nodo por un vertice. Dicho vertice tiene un grosor proporcional al numero de veces que dicha transición ocurrió en la reconstrucción ancestral de sitios. Esta red es una represntación de aquellos sitios que se encuentran en la especie NIES-3974.</p>
</div>

En la (Figura \@ref(fig:FIG45)) podemos ver las transiciones de una forma mas ordenada, con el numero de ocurrencias y la dirección en la que ocurrieron.
<div class="figure" style="text-align: center">
<img src="ResultadosII_files/figure-epub3/FIG45-1.png" alt="**Red de todas las transiciones del clado Calothrix.** En esta imagen se muestra la red  de la Figura \@ref(fig:FIG44) de una forma mas visual y con el numero de veces que ocurrio cada transición, asi como la dirección en la que ocurrió."  />
<p class="caption">(\#fig:FIG45)**Red de todas las transiciones del clado Calothrix.** En esta imagen se muestra la red  de la Figura \@ref(fig:FIG44) de una forma mas visual y con el numero de veces que ocurrio cada transición, asi como la dirección en la que ocurrió.</p>
</div>

#### Transiciones entre Nodo 9 y Nodo 10

Para entender más como es que se gana o se pierden los sitios palindrómicos revisamos la transición en tre los nodos 9 y 10. Esto es porque es esta transicion de nodos la que separa a los dos subclados entre los que hay una repentino cambio de abundancia de sitios palindrómicos (Figura \@ref(fig:FIG15)).

Para hacer esto filtramos los datos de la red para mostrar unicamente las transiciones que se dieron entre los nodos 9 y 10 e hicimos las mismas figuras.
En la (Figura \@ref(fig:FIG46)) se muestra la red como una conexión de nodos a través de vertices con un grosor proporcional al numero de veces que ocurrió cada transición. En la (Figura \@ref(fig:FIG47)) podemos ver las transiciones de una forma mas ordenada, con el numero de ocurrencias y la dirección en la que ocurrieron.



<div class="figure" style="text-align: center">
<img src="ResultadosII_files/figure-epub3/FIG46-1.png" alt="**Red de las transiciones entre los Nodos 9 y 10 del clado Calothrix.** En esta imagen se muestra una red en la que cada nodo es un octanucleótido el cual esta unido a otro nodo por un vertice. Dicho vertice tiene un grosor proporcional al numero de veces que dicha transición ocurrió en la reconstrucción ancestral de sitios."  />
<p class="caption">(\#fig:FIG46)**Red de las transiciones entre los Nodos 9 y 10 del clado Calothrix.** En esta imagen se muestra una red en la que cada nodo es un octanucleótido el cual esta unido a otro nodo por un vertice. Dicho vertice tiene un grosor proporcional al numero de veces que dicha transición ocurrió en la reconstrucción ancestral de sitios.</p>
</div>

<div class="figure" style="text-align: center">
<img src="ResultadosII_files/figure-epub3/FIG47-1.png" alt="**Red de las transiciones entre los nodos 9 y 10 del clado Calothrix.** En esta imagen se muestra la red  de la Figura \@ref(fig:FIG16) de una forma mas visual y con el numero de veces que ocurrio cada transición, asi como la dirección en la que ocurrió."  />
<p class="caption">(\#fig:FIG47)**Red de las transiciones entre los nodos 9 y 10 del clado Calothrix.** En esta imagen se muestra la red  de la Figura \@ref(fig:FIG16) de una forma mas visual y con el numero de veces que ocurrio cada transición, asi como la dirección en la que ocurrió.</p>
</div>

#### Mutaciones en los codones

Para entender como es que se van ganando o perdiendo los sitios palindrómicos hicimos un análisis del tipo mutaciones de los sitios. Esto lo hicimos viendo en que marco de lectura se encontraba cada nodo y revisando la secuencia de aminoacidos que codificaban. En la (Figura \@ref(fig:FIG48)) mostramos 3 gráficos que indican la abundancia de los peptidos codificados por los sitios palindrómicos de acuerdo al marco de lectura en el que se encuentran. En esta figura podemos observar que el marco de lectura es el que contiene la mayoria de los sitios

<div class="figure" style="text-align: center">
<img src="ResultadosII_files/figure-epub3/FIG48-1.png" alt="**Abundancia de peptidos por cada nodo segun el marco de lectura.**." width="110%" /><img src="ResultadosII_files/figure-epub3/FIG48-2.png" alt="**Abundancia de peptidos por cada nodo segun el marco de lectura.**." width="110%" /><img src="ResultadosII_files/figure-epub3/FIG48-3.png" alt="**Abundancia de peptidos por cada nodo segun el marco de lectura.**." width="110%" />
<p class="caption">(\#fig:FIG48)**Abundancia de peptidos por cada nodo segun el marco de lectura.**.</p>
</div>

En la (Figura \@ref(fig:FIG49)) mostramos 3 gráficos que indican la abundancia del tipo de mutaciones que hay en cada nodo de acuerdo al marco de lectura. Lo sitios de mutaciones mostrados pueden ser de los siguientes tipos:

* Conservative (la secuencia de AA cambió pero tiene similitud de acuerdo al score de BLOSUM62)
* ConservativeNoSiteMut (la secuencia de AA cambió pero tiene similitud de acuerdo al score de BLOSUM62. Sin embargo, el sitio no sufrió mutaciones)
* Deletion (La secuencia de AA tiene sufrio 1 o mas deleciones)
* NoMutation (La secuencia de AA no sufrio mutaciones)
* NoSynonym (La secuencia de AA cambió)
* NoSynonymNoSiteMut (La secuencia de AA cambió. Sin embargo, el sitio no sufrió mutaciones.)
* Synonym (El sitio sufrió mutaciones. Sin embargo, la secuencia de AA no cambió.)

<div class="figure" style="text-align: center">
<img src="ResultadosII_files/figure-epub3/FIG49-1.png" alt="**Abundancia del tipo de sustitución por cada nodo segun el marco de lectura.**." width="110%" /><img src="ResultadosII_files/figure-epub3/FIG49-2.png" alt="**Abundancia del tipo de sustitución por cada nodo segun el marco de lectura.**." width="110%" /><img src="ResultadosII_files/figure-epub3/FIG49-3.png" alt="**Abundancia del tipo de sustitución por cada nodo segun el marco de lectura.**." width="110%" />
<p class="caption">(\#fig:FIG49)**Abundancia del tipo de sustitución por cada nodo segun el marco de lectura.**.</p>
</div>

#### Análisis de sitios en los cuales su ancestro era HIP1

Para tratar de entender como es que los sitios HIP1 se pierden hicimos un análisis unicamente en en las transiciones en las que el nodo ancestral tenia un sitio HIP1.

En la (Figura \@ref(fig:FIG50)) mostramos 3 gráficos que indican la frecuencia del tipo de sustituciones que hubo para estos casos para cada nodo en cada uno de los marcos de lectura.

<div class="figure" style="text-align: center">
<img src="ResultadosII_files/figure-epub3/FIG50-1.png" alt="**Abundancia del tipo de sustitución por cada nodo segun el marco de lectura. Unicamente para transiciones en los que el nodo ancestral era un sitio HIP1.**" width="110%" /><img src="ResultadosII_files/figure-epub3/FIG50-2.png" alt="**Abundancia del tipo de sustitución por cada nodo segun el marco de lectura. Unicamente para transiciones en los que el nodo ancestral era un sitio HIP1.**" width="110%" /><img src="ResultadosII_files/figure-epub3/FIG50-3.png" alt="**Abundancia del tipo de sustitución por cada nodo segun el marco de lectura. Unicamente para transiciones en los que el nodo ancestral era un sitio HIP1.**" width="110%" />
<p class="caption">(\#fig:FIG50)**Abundancia del tipo de sustitución por cada nodo segun el marco de lectura. Unicamente para transiciones en los que el nodo ancestral era un sitio HIP1.**</p>
</div>

En la (Figura \@ref(fig:FIG51)) mostramos 3 gráficos (uno por cada marco de lectura) que indican la frecuencia de las mutaciones en cada uno de los 8 nucleótidos del sitio HIP.

<div class="figure" style="text-align: center">
<img src="ResultadosII_files/figure-epub3/FIG51-1.png" alt="**Frecuencia de las mutaciones de cada nucleótido del sitio HIP para cada nodo segun el marco de lectura.**." width="110%" /><img src="ResultadosII_files/figure-epub3/FIG51-2.png" alt="**Frecuencia de las mutaciones de cada nucleótido del sitio HIP para cada nodo segun el marco de lectura.**." width="110%" /><img src="ResultadosII_files/figure-epub3/FIG51-3.png" alt="**Frecuencia de las mutaciones de cada nucleótido del sitio HIP para cada nodo segun el marco de lectura.**." width="110%" />
<p class="caption">(\#fig:FIG51)**Frecuencia de las mutaciones de cada nucleótido del sitio HIP para cada nodo segun el marco de lectura.**.</p>
</div>

En la (Figura \@ref(fig:FIG52)) mostramos 3 gráficos (uno por cada marco de lectura) que indican la frecuencia del tipo sustitucion de bases.

<div class="figure" style="text-align: center">
<img src="ResultadosII_files/figure-epub3/FIG52-1.png" alt="**Frecuencia del tipo de sustituciónes de base en los sitios HIP para cada marco de lectura**." width="110%" /><img src="ResultadosII_files/figure-epub3/FIG52-2.png" alt="**Frecuencia del tipo de sustituciónes de base en los sitios HIP para cada marco de lectura**." width="110%" /><img src="ResultadosII_files/figure-epub3/FIG52-3.png" alt="**Frecuencia del tipo de sustituciónes de base en los sitios HIP para cada marco de lectura**." width="110%" />
<p class="caption">(\#fig:FIG52)**Frecuencia del tipo de sustituciónes de base en los sitios HIP para cada marco de lectura**.</p>
</div>

#### Análisis de sitios en los cuales solo el nodo actual tiene HIP1

Para tratar de entender como es que los sitios HIP se ganan, hicimos un analisis unicamente en las transiciones en las que el nodo actual tenia un sitio HIP1. 

En la Figura \@ref(fig:FIG53) mostramos 3 gráficos (uno por cada marco de lectura) que indican la frecuencia del tipo de sustituciones que hubo para estos casos para cada nodo en cada uno de los marcos de lectura.

<div class="figure" style="text-align: center">
<img src="ResultadosII_files/figure-epub3/FIG53-1.png" alt="**Abundancia del tipo de sustitución por cada nodo segun el marco de lectura. Unicamente para transiciones en los que el nodo actual era un sitio HIP1.**." width="110%" /><img src="ResultadosII_files/figure-epub3/FIG53-2.png" alt="**Abundancia del tipo de sustitución por cada nodo segun el marco de lectura. Unicamente para transiciones en los que el nodo actual era un sitio HIP1.**." width="110%" /><img src="ResultadosII_files/figure-epub3/FIG53-3.png" alt="**Abundancia del tipo de sustitución por cada nodo segun el marco de lectura. Unicamente para transiciones en los que el nodo actual era un sitio HIP1.**." width="110%" />
<p class="caption">(\#fig:FIG53)**Abundancia del tipo de sustitución por cada nodo segun el marco de lectura. Unicamente para transiciones en los que el nodo actual era un sitio HIP1.**.</p>
</div>


### Conjunto de sitios HIP1 usando la especie PCC_6303 como referencia

#### Red de transiciones



visualizamos dicha red .



Para visualizar la red usamos la paqueteria ```networkD3```. Hicimos 2 figuras, la (Figura \@ref(fig:FIG54)) muestra la red como una conexión de nodos a través de vertices con un grosor proporcional al numero de veces que ocurrió cada transición.



<div class="figure" style="text-align: center">
<img src="ResultadosII_files/figure-epub3/FIG54-1.png" alt="**Red de todas las transiciones del clado Calothrix.** En esta imagen se muestra una red en la que cada nodo es un octanucleótido el cual esta unido a otro nodo por un vertice. Dicho vertice tiene un grosor proporcional al numero de veces que dicha transición ocurrió en la reconstrucción ancestral de sitios. Esta red es una represntación de aquellos sitios que se encuentran en la especie NIES-3974."  />
<p class="caption">(\#fig:FIG54)**Red de todas las transiciones del clado Calothrix.** En esta imagen se muestra una red en la que cada nodo es un octanucleótido el cual esta unido a otro nodo por un vertice. Dicho vertice tiene un grosor proporcional al numero de veces que dicha transición ocurrió en la reconstrucción ancestral de sitios. Esta red es una represntación de aquellos sitios que se encuentran en la especie NIES-3974.</p>
</div>

En la (Figura \@ref(fig:FIG55)) podemos ver las transiciones de una forma mas ordenada, con el numero de ocurrencias y la dirección en la que ocurrieron.
<div class="figure" style="text-align: center">
<img src="ResultadosII_files/figure-epub3/FIG55-1.png" alt="**Red de todas las transiciones del clado Calothrix.** En esta imagen se muestra la red  de la Figura \@ref(fig:FIG54) de una forma mas visual y con el numero de veces que ocurrio cada transición, asi como la dirección en la que ocurrió."  />
<p class="caption">(\#fig:FIG55)**Red de todas las transiciones del clado Calothrix.** En esta imagen se muestra la red  de la Figura \@ref(fig:FIG54) de una forma mas visual y con el numero de veces que ocurrio cada transición, asi como la dirección en la que ocurrió.</p>
</div>

#### Transiciones entre Nodo 9 y Nodo 10

Para entender más como es que se gana o se pierden los sitios palindrómicos revisamos la transición en tre los nodos 9 y 10. Esto es porque es esta transicion de nodos la que separa a los dos subclados entre los que hay una repentino cambio de abundancia de sitios palindrómicos (Figura \@ref(fig:FIG15)).

Para hacer esto filtramos los datos de la red para mostrar unicamente las transiciones que se dieron entre los nodos 9 y 10 e hicimos las mismas figuras.
En la (Figura \@ref(fig:FIG56)) se muestra la red como una conexión de nodos a través de vertices con un grosor proporcional al numero de veces que ocurrió cada transición. En la (Figura \@ref(fig:FIG57)) podemos ver las transiciones de una forma mas ordenada, con el numero de ocurrencias y la dirección en la que ocurrieron.



<div class="figure" style="text-align: center">
<img src="ResultadosII_files/figure-epub3/FIG56-1.png" alt="**Red de las transiciones entre los Nodos 9 y 10 del clado Calothrix.** En esta imagen se muestra una red en la que cada nodo es un octanucleótido el cual esta unido a otro nodo por un vertice. Dicho vertice tiene un grosor proporcional al numero de veces que dicha transición ocurrió en la reconstrucción ancestral de sitios."  />
<p class="caption">(\#fig:FIG56)**Red de las transiciones entre los Nodos 9 y 10 del clado Calothrix.** En esta imagen se muestra una red en la que cada nodo es un octanucleótido el cual esta unido a otro nodo por un vertice. Dicho vertice tiene un grosor proporcional al numero de veces que dicha transición ocurrió en la reconstrucción ancestral de sitios.</p>
</div>

<div class="figure" style="text-align: center">
<img src="ResultadosII_files/figure-epub3/FIG57-1.png" alt="**Red de las transiciones entre los nodos 9 y 10 del clado Calothrix.** En esta imagen se muestra la red  de la Figura \@ref(fig:FIG56) de una forma mas visual y con el numero de veces que ocurrio cada transición, asi como la dirección en la que ocurrió."  />
<p class="caption">(\#fig:FIG57)**Red de las transiciones entre los nodos 9 y 10 del clado Calothrix.** En esta imagen se muestra la red  de la Figura \@ref(fig:FIG56) de una forma mas visual y con el numero de veces que ocurrio cada transición, asi como la dirección en la que ocurrió.</p>
</div>

#### Mutaciones en los codones

Para entender como es que se van ganando o perdiendo los sitios palindrómicos hicimos un análisis del tipo mutaciones de los sitios. Esto lo hicimos viendo en que marco de lectura se encontraba cada nodo y revisando la secuencia de aminoacidos que codificaban. En la (Figura \@ref(fig:FIG58)) mostramos 3 gráficos que indican la abundancia de los peptidos codificados por los sitios palindrómicos de acuerdo al marco de lectura en el que se encuentran. En esta figura podemos observar que el marco de lectura es el que contiene la mayoria de los sitios

<div class="figure" style="text-align: center">
<img src="ResultadosII_files/figure-epub3/FIG58-1.png" alt="**Abundancia de peptidos por cada nodo segun el marco de lectura.**." width="110%" /><img src="ResultadosII_files/figure-epub3/FIG58-2.png" alt="**Abundancia de peptidos por cada nodo segun el marco de lectura.**." width="110%" /><img src="ResultadosII_files/figure-epub3/FIG58-3.png" alt="**Abundancia de peptidos por cada nodo segun el marco de lectura.**." width="110%" />
<p class="caption">(\#fig:FIG58)**Abundancia de peptidos por cada nodo segun el marco de lectura.**.</p>
</div>

En la (Figura \@ref(fig:FIG59)) mostramos 3 gráficos que indican la abundancia del tipo de mutaciones que hay en cada nodo de acuerdo al marco de lectura. Lo sitios de mutaciones mostrados pueden ser de los siguientes tipos:

* Conservative (la secuencia de AA cambió pero tiene similitud de acuerdo al score de BLOSUM62)
* ConservativeNoSiteMut (la secuencia de AA cambió pero tiene similitud de acuerdo al score de BLOSUM62. Sin embargo, el sitio no sufrió mutaciones)
* Deletion (La secuencia de AA tiene sufrio 1 o mas deleciones)
* NoMutation (La secuencia de AA no sufrio mutaciones)
* NoSynonym (La secuencia de AA cambió)
* NoSynonymNoSiteMut (La secuencia de AA cambió. Sin embargo, el sitio no sufrió mutaciones.)
* Synonym (El sitio sufrió mutaciones. Sin embargo, la secuencia de AA no cambió.)

<div class="figure" style="text-align: center">
<img src="ResultadosII_files/figure-epub3/FIG59-1.png" alt="**Abundancia del tipo de sustitución por cada nodo segun el marco de lectura.**." width="110%" /><img src="ResultadosII_files/figure-epub3/FIG59-2.png" alt="**Abundancia del tipo de sustitución por cada nodo segun el marco de lectura.**." width="110%" /><img src="ResultadosII_files/figure-epub3/FIG59-3.png" alt="**Abundancia del tipo de sustitución por cada nodo segun el marco de lectura.**." width="110%" />
<p class="caption">(\#fig:FIG59)**Abundancia del tipo de sustitución por cada nodo segun el marco de lectura.**.</p>
</div>

#### Análisis de sitios en los cuales su ancestro era HIP1

Para tratar de entender como es que los sitios HIP1 se pierden hicimos un análisis unicamente en en las transiciones en las que el nodo ancestral tenia un sitio HIP1.

En la (Figura \@ref(fig:FIG60)) mostramos 3 gráficos que indican la frecuencia del tipo de sustituciones que hubo para estos casos para cada nodo en cada uno de los marcos de lectura.

<div class="figure" style="text-align: center">
<img src="ResultadosII_files/figure-epub3/FIG60-1.png" alt="**Abundancia del tipo de sustitución por cada nodo segun el marco de lectura. Unicamente para transiciones en los que el nodo ancestral era un sitio HIP1.**" width="110%" /><img src="ResultadosII_files/figure-epub3/FIG60-2.png" alt="**Abundancia del tipo de sustitución por cada nodo segun el marco de lectura. Unicamente para transiciones en los que el nodo ancestral era un sitio HIP1.**" width="110%" /><img src="ResultadosII_files/figure-epub3/FIG60-3.png" alt="**Abundancia del tipo de sustitución por cada nodo segun el marco de lectura. Unicamente para transiciones en los que el nodo ancestral era un sitio HIP1.**" width="110%" />
<p class="caption">(\#fig:FIG60)**Abundancia del tipo de sustitución por cada nodo segun el marco de lectura. Unicamente para transiciones en los que el nodo ancestral era un sitio HIP1.**</p>
</div>

En la (Figura \@ref(fig:FIG61)) mostramos 3 gráficos (uno por cada marco de lectura) que indican la frecuencia de las mutaciones en cada uno de los 8 nucleótidos del sitio HIP.

<div class="figure" style="text-align: center">
<img src="ResultadosII_files/figure-epub3/FIG61-1.png" alt="**Frecuencia de las mutaciones de cada nucleótido del sitio HIP para cada nodo segun el marco de lectura.**." width="110%" /><img src="ResultadosII_files/figure-epub3/FIG61-2.png" alt="**Frecuencia de las mutaciones de cada nucleótido del sitio HIP para cada nodo segun el marco de lectura.**." width="110%" /><img src="ResultadosII_files/figure-epub3/FIG61-3.png" alt="**Frecuencia de las mutaciones de cada nucleótido del sitio HIP para cada nodo segun el marco de lectura.**." width="110%" />
<p class="caption">(\#fig:FIG61)**Frecuencia de las mutaciones de cada nucleótido del sitio HIP para cada nodo segun el marco de lectura.**.</p>
</div>

En la (Figura \@ref(fig:FIG62)) mostramos 3 gráficos (uno por cada marco de lectura) que indican la frecuencia del tipo sustitucion de bases.

<div class="figure" style="text-align: center">
<img src="ResultadosII_files/figure-epub3/FIG62-1.png" alt="**Frecuencia del tipo de sustituciónes de base en los sitios HIP para cada marco de lectura**." width="110%" /><img src="ResultadosII_files/figure-epub3/FIG62-2.png" alt="**Frecuencia del tipo de sustituciónes de base en los sitios HIP para cada marco de lectura**." width="110%" /><img src="ResultadosII_files/figure-epub3/FIG62-3.png" alt="**Frecuencia del tipo de sustituciónes de base en los sitios HIP para cada marco de lectura**." width="110%" />
<p class="caption">(\#fig:FIG62)**Frecuencia del tipo de sustituciónes de base en los sitios HIP para cada marco de lectura**.</p>
</div>

#### Análisis de sitios en los cuales solo el nodo actual tiene HIP1

Para tratar de entender como es que los sitios HIP se ganan, hicimos un analisis unicamente en las transiciones en las que el nodo actual tenia un sitio HIP1. 

En la Figura \@ref(fig:FIG63) mostramos 3 gráficos (uno por cada marco de lectura) que indican la frecuencia del tipo de sustituciones que hubo para estos casos para cada nodo en cada uno de los marcos de lectura.

<div class="figure" style="text-align: center">
<img src="ResultadosII_files/figure-epub3/FIG63-1.png" alt="**Abundancia del tipo de sustitución por cada nodo segun el marco de lectura. Unicamente para transiciones en los que el nodo actual era un sitio HIP1.**." width="110%" /><img src="ResultadosII_files/figure-epub3/FIG63-2.png" alt="**Abundancia del tipo de sustitución por cada nodo segun el marco de lectura. Unicamente para transiciones en los que el nodo actual era un sitio HIP1.**." width="110%" /><img src="ResultadosII_files/figure-epub3/FIG63-3.png" alt="**Abundancia del tipo de sustitución por cada nodo segun el marco de lectura. Unicamente para transiciones en los que el nodo actual era un sitio HIP1.**." width="110%" />
<p class="caption">(\#fig:FIG63)**Abundancia del tipo de sustitución por cada nodo segun el marco de lectura. Unicamente para transiciones en los que el nodo actual era un sitio HIP1.**.</p>
</div>


### Conjunto de sitios HIP1 unicos en las especies 336-3, NIES-3974 y PCC_6303

#### Red de transiciones

Para aumentar el numero de "experimentos", buscamos todos los sitios HIP1 UNICOS existentes en el subclado que contiene a las especies 336-3, NIES-3974 y PCC_6303.

Posteriormente creamos la red usando la función ```Create_Network()```:



y visualizamos dicha red .



Para visualizar la red usamos la paqueteria ```networkD3```. Hicimos 2 figuras, la (Figura \@ref(fig:FIG34)) muestra la red como una conexión de nodos a través de vertices con un grosor proporcional al numero de veces que ocurrió cada transición.



<div class="figure" style="text-align: center">
<img src="ResultadosII_files/figure-epub3/FIG34-1.png" alt="**Red de todas las transiciones del clado Calothrix.** En esta imagen se muestra una red en la que cada nodo es un octanucleótido el cual esta unido a otro nodo por un vertice. Dicho vertice tiene un grosor proporcional al numero de veces que dicha transición ocurrió en la reconstrucción ancestral de sitios. Esta red es una represntación de aquellos sitios que son unicos en las especies 336-3, NIES-3974 y PCC\_6303"  />
<p class="caption">(\#fig:FIG34)**Red de todas las transiciones del clado Calothrix.** En esta imagen se muestra una red en la que cada nodo es un octanucleótido el cual esta unido a otro nodo por un vertice. Dicho vertice tiene un grosor proporcional al numero de veces que dicha transición ocurrió en la reconstrucción ancestral de sitios. Esta red es una represntación de aquellos sitios que son unicos en las especies 336-3, NIES-3974 y PCC\_6303</p>
</div>

En la (Figura \@ref(fig:FIG35)) podemos ver las transiciones de una forma mas ordenada, con el numero de ocurrencias y la dirección en la que ocurrieron.
<div class="figure" style="text-align: center">
<img src="ResultadosII_files/figure-epub3/FIG35-1.png" alt="**Red de todas las transiciones del clado Calothrix.** En esta imagen se muestra la red  de la Figura \@ref(fig:FIG34) de una forma mas visual y con el numero de veces que ocurrio cada transición, asi como la dirección en la que ocurrió."  />
<p class="caption">(\#fig:FIG35)**Red de todas las transiciones del clado Calothrix.** En esta imagen se muestra la red  de la Figura \@ref(fig:FIG34) de una forma mas visual y con el numero de veces que ocurrio cada transición, asi como la dirección en la que ocurrió.</p>
</div>

#### Transiciones entre Nodo 9 y Nodo 10

En la (Figura \@ref(fig:FIG36)) se muestra la red como una conexión de nodos a través de vertices con un grosor proporcional al numero de veces que ocurrió cada transición. En la (Figura \@ref(fig:FIG37)) podemos ver las transiciones de una forma mas ordenada, con el numero de ocurrencias y la dirección en la que ocurrieron.



<div class="figure" style="text-align: center">
<img src="ResultadosII_files/figure-epub3/FIG36-1.png" alt="**Red de las transiciones entre los Nodos 9 y 10 del clado Calothrix.** En esta imagen se muestra una red en la que cada nodo es un octanucleótido el cual esta unido a otro nodo por un vertice. Dicho vertice tiene un grosor proporcional al numero de veces que dicha transición ocurrió en la reconstrucción ancestral de sitios."  />
<p class="caption">(\#fig:FIG36)**Red de las transiciones entre los Nodos 9 y 10 del clado Calothrix.** En esta imagen se muestra una red en la que cada nodo es un octanucleótido el cual esta unido a otro nodo por un vertice. Dicho vertice tiene un grosor proporcional al numero de veces que dicha transición ocurrió en la reconstrucción ancestral de sitios.</p>
</div>

<div class="figure" style="text-align: center">
<img src="ResultadosII_files/figure-epub3/FIG37-1.png" alt="**Red de las transiciones entre los nodos 9 y 10 del clado Calothrix.** En esta imagen se muestra la red  de la Figura \@ref(fig:FIG16) de una forma mas visual y con el numero de veces que ocurrio cada transición, asi como la dirección en la que ocurrió."  />
<p class="caption">(\#fig:FIG37)**Red de las transiciones entre los nodos 9 y 10 del clado Calothrix.** En esta imagen se muestra la red  de la Figura \@ref(fig:FIG16) de una forma mas visual y con el numero de veces que ocurrio cada transición, asi como la dirección en la que ocurrió.</p>
</div>

#### Mutaciones en los codones

Para entender como es que se van ganando o perdiendo los sitios palindrómicos hicimos un análisis del tipo mutaciones de los sitios. Esto lo hicimos viendo en que marco de lectura se encontraba cada nodo y revisando la secuencia de aminoacidos que codificaban. En la (Figura \@ref(fig:FIG38)) mostramos 3 gráficos que indican la abundancia de los peptidos codificados por los sitios palindrómicos de acuerdo al marco de lectura en el que se encuentran. En esta figura podemos observar que el marco de lectura es el que contiene la mayoria de los sitios

<div class="figure" style="text-align: center">
<img src="ResultadosII_files/figure-epub3/FIG38-1.png" alt="**Abundancia de peptidos por cada nodo segun el marco de lectura.**." width="110%" /><img src="ResultadosII_files/figure-epub3/FIG38-2.png" alt="**Abundancia de peptidos por cada nodo segun el marco de lectura.**." width="110%" /><img src="ResultadosII_files/figure-epub3/FIG38-3.png" alt="**Abundancia de peptidos por cada nodo segun el marco de lectura.**." width="110%" />
<p class="caption">(\#fig:FIG38)**Abundancia de peptidos por cada nodo segun el marco de lectura.**.</p>
</div>

En la (Figura \@ref(fig:FIG39)) mostramos 3 gráficos que indican la abundancia del tipo de mutaciones que hay en cada nodo de acuerdo al marco de lectura. Lo sitios de mutaciones mostrados pueden ser de los siguientes tipos:

* Conservative (la secuencia de AA cambió pero tiene similitud de acuerdo al score de BLOSUM62)
* ConservativeNoSiteMut (la secuencia de AA cambió pero tiene similitud de acuerdo al score de BLOSUM62. Sin embargo, el sitio no sufrió mutaciones)
* Deletion (La secuencia de AA tiene sufrio 1 o mas deleciones)
* NoMutation (La secuencia de AA no sufrio mutaciones)
* NoSynonym (La secuencia de AA cambió)
* NoSynonymNoSiteMut (La secuencia de AA cambió. Sin embargo, el sitio no sufrió mutaciones.)
* Synonym (El sitio sufrió mutaciones. Sin embargo, la secuencia de AA no cambió.)

<div class="figure" style="text-align: center">
<img src="ResultadosII_files/figure-epub3/FIG39-1.png" alt="**Abundancia del tipo de sustitución por cada nodo segun el marco de lectura.**." width="110%" /><img src="ResultadosII_files/figure-epub3/FIG39-2.png" alt="**Abundancia del tipo de sustitución por cada nodo segun el marco de lectura.**." width="110%" /><img src="ResultadosII_files/figure-epub3/FIG39-3.png" alt="**Abundancia del tipo de sustitución por cada nodo segun el marco de lectura.**." width="110%" />
<p class="caption">(\#fig:FIG39)**Abundancia del tipo de sustitución por cada nodo segun el marco de lectura.**.</p>
</div>

#### Análisis de sitios en los cuales su ancestro era HIP1

Para tratar de entender como es que los sitios HIP1 se pierden hicimos un análisis unicamente en en las transiciones en las que el nodo ancestral tenia un sitio HIP1.

En la (Figura \@ref(fig:FIG40)) mostramos 3 gráficos que indican la frecuencia del tipo de sustituciones que hubo para estos casos para cada nodo en cada uno de los marcos de lectura.

<div class="figure" style="text-align: center">
<img src="ResultadosII_files/figure-epub3/FIG40-1.png" alt="**Abundancia del tipo de sustitución por cada nodo segun el marco de lectura. Unicamente para transiciones en los que el nodo ancestral era un sitio HIP1.**" width="110%" /><img src="ResultadosII_files/figure-epub3/FIG40-2.png" alt="**Abundancia del tipo de sustitución por cada nodo segun el marco de lectura. Unicamente para transiciones en los que el nodo ancestral era un sitio HIP1.**" width="110%" /><img src="ResultadosII_files/figure-epub3/FIG40-3.png" alt="**Abundancia del tipo de sustitución por cada nodo segun el marco de lectura. Unicamente para transiciones en los que el nodo ancestral era un sitio HIP1.**" width="110%" />
<p class="caption">(\#fig:FIG40)**Abundancia del tipo de sustitución por cada nodo segun el marco de lectura. Unicamente para transiciones en los que el nodo ancestral era un sitio HIP1.**</p>
</div>

En la (Figura \@ref(fig:FIG41)) mostramos 3 gráficos (uno por cada marco de lectura) que indican la frecuencia de las mutaciones en cada uno de los 8 nucleótidos del sitio HIP.

<div class="figure" style="text-align: center">
<img src="ResultadosII_files/figure-epub3/FIG41-1.png" alt="**Frecuencia de las mutaciones de cada nucleótido del sitio HIP para cada nodo segun el marco de lectura.**." width="110%" /><img src="ResultadosII_files/figure-epub3/FIG41-2.png" alt="**Frecuencia de las mutaciones de cada nucleótido del sitio HIP para cada nodo segun el marco de lectura.**." width="110%" /><img src="ResultadosII_files/figure-epub3/FIG41-3.png" alt="**Frecuencia de las mutaciones de cada nucleótido del sitio HIP para cada nodo segun el marco de lectura.**." width="110%" />
<p class="caption">(\#fig:FIG41)**Frecuencia de las mutaciones de cada nucleótido del sitio HIP para cada nodo segun el marco de lectura.**.</p>
</div>

En la (Figura \@ref(fig:FIG42)) mostramos 3 gráficos (uno por cada marco de lectura) que indican la frecuencia del tipo sustitucion de bases.

<div class="figure" style="text-align: center">
<img src="ResultadosII_files/figure-epub3/FIG42-1.png" alt="**Frecuencia del tipo de sustituciónes de base en los sitios HIP para cada marco de lectura**." width="110%" /><img src="ResultadosII_files/figure-epub3/FIG42-2.png" alt="**Frecuencia del tipo de sustituciónes de base en los sitios HIP para cada marco de lectura**." width="110%" /><img src="ResultadosII_files/figure-epub3/FIG42-3.png" alt="**Frecuencia del tipo de sustituciónes de base en los sitios HIP para cada marco de lectura**." width="110%" />
<p class="caption">(\#fig:FIG42)**Frecuencia del tipo de sustituciónes de base en los sitios HIP para cada marco de lectura**.</p>
</div>

#### Análisis de sitios en los cuales solo el nodo actual tiene HIP1

Para tratar de entender como es que los sitios HIP se ganan, hicimos un analisis unicamente en las transiciones en las que el nodo actual tenia un sitio HIP1. 

En la Figura \@ref(fig:FIG43) mostramos 3 gráficos (uno por cada marco de lectura) que indican la frecuencia del tipo de sustituciones que hubo para estos casos para cada nodo en cada uno de los marcos de lectura.

<div class="figure" style="text-align: center">
<img src="ResultadosII_files/figure-epub3/FIG43-1.png" alt="**Abundancia del tipo de sustitución por cada nodo segun el marco de lectura. Unicamente para transiciones en los que el nodo actual era un sitio HIP1.**." width="110%" /><img src="ResultadosII_files/figure-epub3/FIG43-2.png" alt="**Abundancia del tipo de sustitución por cada nodo segun el marco de lectura. Unicamente para transiciones en los que el nodo actual era un sitio HIP1.**." width="110%" /><img src="ResultadosII_files/figure-epub3/FIG43-3.png" alt="**Abundancia del tipo de sustitución por cada nodo segun el marco de lectura. Unicamente para transiciones en los que el nodo actual era un sitio HIP1.**." width="110%" />
<p class="caption">(\#fig:FIG43)**Abundancia del tipo de sustitución por cada nodo segun el marco de lectura. Unicamente para transiciones en los que el nodo actual era un sitio HIP1.**.</p>
</div>


### TGGCGCCA

#### Red de transiciones

Para hacer mas visual la reconstrucción, construimos una red de las transiciones entre los estados ancestrales. Esto lo hicimos en r usando la función ```Create_Transition_Table()```:


```r
source("ASR_Orth_Functions/NodeAndEdges.R")
Nodes.Edges <- Create_Transition_Table(SitesTable = "Clados/Callothrix_clade/PALINDROMES/TGGCGCCA/PCC_7716/Orthologues_Palindrome_sites.txt",
                                EvolutionModel = "F81",
                                Method = "bayes",
                                Phylogeny = "Clados/Callothrix_clade/SpeciesTree_rooted.txt",
                                OrthoPath = "Clados/Callothrix_clade/PALINDROMES/TGGCGCCA/PCC_7716/Only_ORTHOLOGUES/")

```

Posteriormente creamos la red usando la función ```Create_Network()```:



y visualizamos dicha red .



Para visualizar la red usamos la paqueteria ```networkD3```. Hicimos 2 figuras, la (Figura \@ref(fig:FIG24)) muestra la red como una conexión de nodos a través de vertices con un grosor proporcional al numero de veces que ocurrió cada transición. En dicha red podemos ver algunos nodos con bordes muy gruesos como **GCAATTGC**, **GCAATCGC**, **GCAATAGC**, **GCGATTGC** (Tabla \@ref(tab:TAB62)).



<div class="figure" style="text-align: center">
<img src="ResultadosII_files/figure-epub3/FIG24-1.png" alt="**Red de todas las transiciones del clado Calothrix.** En esta imagen se muestra una red en la que cada nodo es un octanucleótido el cual esta unido a otro nodo por un vertice. Dicho vertice tiene un grosor proporcional al numero de veces que dicha transición ocurrió en la reconstrucción ancestral de sitios."  />
<p class="caption">(\#fig:FIG24)**Red de todas las transiciones del clado Calothrix.** En esta imagen se muestra una red en la que cada nodo es un octanucleótido el cual esta unido a otro nodo por un vertice. Dicho vertice tiene un grosor proporcional al numero de veces que dicha transición ocurrió en la reconstrucción ancestral de sitios.</p>
</div>

En la (Figura \@ref(fig:FIG25)) podemos ver las transiciones de una forma mas ordenada, con el numero de ocurrencias y la dirección en la que ocurrieron.
<div class="figure" style="text-align: center">
<img src="ResultadosII_files/figure-epub3/FIG25-1.png" alt="**Red de todas las transiciones del clado Calothrix.** En esta imagen se muestra la red  de la Figura \@ref(fig:FIG24) de una forma mas visual y con el numero de veces que ocurrio cada transición, asi como la dirección en la que ocurrió."  />
<p class="caption">(\#fig:FIG25)**Red de todas las transiciones del clado Calothrix.** En esta imagen se muestra la red  de la Figura \@ref(fig:FIG24) de una forma mas visual y con el numero de veces que ocurrio cada transición, asi como la dirección en la que ocurrió.</p>
</div>

<table class=" lightable-paper table table-striped table-condensed" style='font-family: "Arial Narrow", arial, helvetica, sans-serif; margin-left: auto; margin-right: auto; font-size: 11px; width: auto !important; margin-left: auto; margin-right: auto;'>
<caption style="font-size: initial !important;">(\#tab:TAB62)Transiciones entre nodos.</caption>
 <thead>
  <tr>
   <th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;"> from </th>
   <th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;"> to </th>
   <th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;"> weight </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;font-weight: bold;vertical-align: top !important;" rowspan="2"> <span style="     color: red !important;">TGGCGCCA</span> </td>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: black !important;">TGGTGCCA</span> </td>
   <td style="text-align:left;"> <span style="display: block; padding: 0 4px; border-radius: 4px; background-color: #ffa500">54</span> </td>
  </tr>
  <tr>
   
   <td style="text-align:left;font-weight: bold;"> <span style="     color: black !important;">TGGCACCA</span> </td>
   <td style="text-align:left;"> <span style="display: block; padding: 0 4px; border-radius: 4px; background-color: #ffcb6c">34</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: black !important;">CGGCGCCA</span> </td>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: red !important;">TGGCGCCA</span> </td>
   <td style="text-align:left;"> <span style="display: block; padding: 0 4px; border-radius: 4px; background-color: #ffd487">29</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;font-weight: bold;vertical-align: top !important;" rowspan="4"> <span style="     color: red !important;">TGGCGCCA</span> </td>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: black !important;">TGGAGCCA</span> </td>
   <td style="text-align:left;"> <span style="display: block; padding: 0 4px; border-radius: 4px; background-color: #ffd487">29</span> </td>
  </tr>
  <tr>
   
   <td style="text-align:left;font-weight: bold;"> <span style="     color: black !important;">TGGCGCAA</span> </td>
   <td style="text-align:left;"> <span style="display: block; padding: 0 4px; border-radius: 4px; background-color: #ffd892">27</span> </td>
  </tr>
  <tr>
   
   <td style="text-align:left;font-weight: bold;"> <span style="     color: black !important;">TGGCGACA</span> </td>
   <td style="text-align:left;"> <span style="display: block; padding: 0 4px; border-radius: 4px; background-color: #ffda97">26</span> </td>
  </tr>
  <tr>
   
   <td style="text-align:left;font-weight: bold;"> <span style="     color: black !important;">TGGTGCAA</span> </td>
   <td style="text-align:left;"> <span style="display: block; padding: 0 4px; border-radius: 4px; background-color: #ffe4b3">21</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: black !important;">TGGCGCCG</span> </td>
   <td style="text-align:left;font-weight: bold;vertical-align: top !important;" rowspan="2"> <span style="     color: red !important;">TGGCGCCA</span> </td>
   <td style="text-align:left;"> <span style="display: block; padding: 0 4px; border-radius: 4px; background-color: #ffe6b8">20</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: black !important;">--------</span> </td>
   
   <td style="text-align:left;"> <span style="display: block; padding: 0 4px; border-radius: 4px; background-color: #ffedce">16</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;font-weight: bold;vertical-align: top !important;" rowspan="2"> <span style="     color: red !important;">TGGCGCCA</span> </td>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: black !important;">TAGCACCA</span> </td>
   <td style="text-align:left;"> <span style="display: block; padding: 0 4px; border-radius: 4px; background-color: #ffedce">16</span> </td>
  </tr>
  <tr>
   
   <td style="text-align:left;font-weight: bold;"> <span style="     color: black !important;">AGGCGCCA</span> </td>
   <td style="text-align:left;"> <span style="display: block; padding: 0 4px; border-radius: 4px; background-color: #ffefd3">15</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: black !important;">TGGCGCCT</span> </td>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: red !important;">TGGCGCCA</span> </td>
   <td style="text-align:left;"> <span style="display: block; padding: 0 4px; border-radius: 4px; background-color: #fff1d9">14</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;font-weight: bold;vertical-align: top !important;" rowspan="5"> <span style="     color: red !important;">TGGCGCCA</span> </td>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: black !important;">TGGCGCCG</span> </td>
   <td style="text-align:left;"> <span style="display: block; padding: 0 4px; border-radius: 4px; background-color: #fff3de">13</span> </td>
  </tr>
  <tr>
   
   <td style="text-align:left;font-weight: bold;"> <span style="     color: black !important;">TGGCGCTA</span> </td>
   <td style="text-align:left;"> <span style="display: block; padding: 0 4px; border-radius: 4px; background-color: #fff5e3">12</span> </td>
  </tr>
  <tr>
   
   <td style="text-align:left;font-weight: bold;"> <span style="     color: black !important;">TTGCGCCA</span> </td>
   <td style="text-align:left;"> <span style="display: block; padding: 0 4px; border-radius: 4px; background-color: #fff5e3">12</span> </td>
  </tr>
  <tr>
   
   <td style="text-align:left;font-weight: bold;"> <span style="     color: black !important;">CGGCGCCA</span> </td>
   <td style="text-align:left;"> <span style="display: block; padding: 0 4px; border-radius: 4px; background-color: #fff5e3">12</span> </td>
  </tr>
  <tr>
   
   <td style="text-align:left;font-weight: bold;"> <span style="     color: black !important;">TGGTGCTA</span> </td>
   <td style="text-align:left;"> <span style="display: block; padding: 0 4px; border-radius: 4px; background-color: #fff5e3">12</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: black !important;">TGGTGCTA</span> </td>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: black !important;">TGGTGCCA</span> </td>
   <td style="text-align:left;"> <span style="display: block; padding: 0 4px; border-radius: 4px; background-color: #fff7e9">11</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: red !important;">TGGCGCCA</span> </td>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: black !important;">TAGCGCCA</span> </td>
   <td style="text-align:left;"> <span style="display: block; padding: 0 4px; border-radius: 4px; background-color: #fff7e9">11</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: black !important;">AGGCGCCA</span> </td>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: red !important;">TGGCGCCA</span> </td>
   <td style="text-align:left;"> <span style="display: block; padding: 0 4px; border-radius: 4px; background-color: #fff7e9">11</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;font-weight: bold;vertical-align: top !important;" rowspan="4"> <span style="     color: red !important;">TGGCGCCA</span> </td>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: black !important;">TGGCTCCA</span> </td>
   <td style="text-align:left;"> <span style="display: block; padding: 0 4px; border-radius: 4px; background-color: #fff7e9">11</span> </td>
  </tr>
  <tr>
   
   <td style="text-align:left;font-weight: bold;"> <span style="     color: black !important;">TGGCGCGA</span> </td>
   <td style="text-align:left;"> <span style="display: block; padding: 0 4px; border-radius: 4px; background-color: #fff7e9">11</span> </td>
  </tr>
  <tr>
   
   <td style="text-align:left;font-weight: bold;"> <span style="     color: black !important;">TGGGGCAA</span> </td>
   <td style="text-align:left;"> <span style="display: block; padding: 0 4px; border-radius: 4px; background-color: #fffbf4">9</span> </td>
  </tr>
  <tr>
   
   <td style="text-align:left;font-weight: bold;"> <span style="     color: black !important;">TGGAGCAA</span> </td>
   <td style="text-align:left;"> <span style="display: block; padding: 0 4px; border-radius: 4px; background-color: #fffdf9">8</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: black !important;">TGGTGCCA</span> </td>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: black !important;">TGGTGCTA</span> </td>
   <td style="text-align:left;"> <span style="display: block; padding: 0 4px; border-radius: 4px; background-color: #ffffff">7</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;font-weight: bold;vertical-align: top !important;" rowspan="2"> <span style="     color: red !important;">TGGCGCCA</span> </td>
   <td style="text-align:left;font-weight: bold;"> <span style="     color: black !important;">TGGCGCCT</span> </td>
   <td style="text-align:left;"> <span style="display: block; padding: 0 4px; border-radius: 4px; background-color: #ffffff">7</span> </td>
  </tr>
  <tr>
   
   <td style="text-align:left;font-weight: bold;"> <span style="     color: black !important;">TTGCACCA</span> </td>
   <td style="text-align:left;"> <span style="display: block; padding: 0 4px; border-radius: 4px; background-color: #ffffff">7</span> </td>
  </tr>
</tbody>
</table>



#### Transiciones entre Nodo 9 y Nodo 10
Para entender más como es que se gana o se pierden los sitios palindrómicos revisamos la transición en tre los nodos 9 y 10. Esto es porque es esta transicion de nodos la que separa a los dos subclados entre los que hay una repentino cambio de abundancia de sitios palindrómicos (Figura \@ref(fig:FIG15)).


Para hacer esto filtramos los datos de la red para mostrar unicamente las transiciones que se dieron entre los nodos 9 y 10 e hicimos las mismas figuras.
En la (Figura \@ref(fig:FIG26)) se muestra la red como una conexión de nodos a través de vertices con un grosor proporcional al numero de veces que ocurrió cada transición. En la (Figura \@ref(fig:FIG27)) podemos ver las transiciones de una forma mas ordenada, con el numero de ocurrencias y la dirección en la que ocurrieron.



<div class="figure" style="text-align: center">
<img src="ResultadosII_files/figure-epub3/FIG26-1.png" alt="**Red de las transiciones entre los Nodos 9 y 10 del clado Calothrix.** En esta imagen se muestra una red en la que cada nodo es un octanucleótido el cual esta unido a otro nodo por un vertice. Dicho vertice tiene un grosor proporcional al numero de veces que dicha transición ocurrió en la reconstrucción ancestral de sitios."  />
<p class="caption">(\#fig:FIG26)**Red de las transiciones entre los Nodos 9 y 10 del clado Calothrix.** En esta imagen se muestra una red en la que cada nodo es un octanucleótido el cual esta unido a otro nodo por un vertice. Dicho vertice tiene un grosor proporcional al numero de veces que dicha transición ocurrió en la reconstrucción ancestral de sitios.</p>
</div>

<div class="figure" style="text-align: center">
<img src="ResultadosII_files/figure-epub3/FIG27-1.png" alt="**Red de las transiciones entre los nodos 9 y 10 del clado Calothrix.** En esta imagen se muestra la red  de la Figura \@ref(fig:FIG26) de una forma mas visual y con el numero de veces que ocurrio cada transición, asi como la dirección en la que ocurrió."  />
<p class="caption">(\#fig:FIG27)**Red de las transiciones entre los nodos 9 y 10 del clado Calothrix.** En esta imagen se muestra la red  de la Figura \@ref(fig:FIG26) de una forma mas visual y con el numero de veces que ocurrio cada transición, asi como la dirección en la que ocurrió.</p>
</div>

#### Mutaciones en los codones

Para entender como es que se van ganando o perdiendo los sitios palindrómicos hicimos un análisis del tipo mutaciones de los sitios. Esto lo hicimos viendo en que marco de lectura se encontraba cada nodo y revisando la secuencia de aminoacidos que codificaban. En la (Figura \@ref(fig:FIG28)) mostramos 3 gráficos que indican la abundancia de los peptidos codificados por los sitios palindrómicos de acuerdo al marco de lectura en el que se encuentran.

<div class="figure" style="text-align: center">
<img src="ResultadosII_files/figure-epub3/FIG28-1.png" alt="**Abundancia de peptidos por cada nodo segun el marco de lectura.**." width="110%" /><img src="ResultadosII_files/figure-epub3/FIG28-2.png" alt="**Abundancia de peptidos por cada nodo segun el marco de lectura.**." width="110%" /><img src="ResultadosII_files/figure-epub3/FIG28-3.png" alt="**Abundancia de peptidos por cada nodo segun el marco de lectura.**." width="110%" />
<p class="caption">(\#fig:FIG28)**Abundancia de peptidos por cada nodo segun el marco de lectura.**.</p>
</div>

En la (Figura \@ref(fig:FIG29)) mostramos 3 gráficos que indican la abundancia del tipo de mutaciones que hay en cada nodo de acuerdo al marco de lectura. Lo sitios de mutaciones mostrados pueden ser de los siguientes tipos:

* Conservative (la secuencia de AA cambió pero tiene similitud de acuerdo al score de BLOSUM62)
* ConservativeNoSiteMut (la secuencia de AA cambió pero tiene similitud de acuerdo al score de BLOSUM62. Sin embargo, el sitio no sufrió mutaciones)
* Deletion (La secuencia de AA tiene sufrio 1 o mas deleciones)
* NoMutation (La secuencia de AA no sufrio mutaciones)
* NoSynonym (La secuencia de AA cambió)
* NoSynonymNoSiteMut (La secuencia de AA cambió. Sin embargo, el sitio no sufrió mutaciones.)
* Synonym (El sitio sufrió mutaciones. Sin embargo, la secuencia de AA no cambió.)

<div class="figure" style="text-align: center">
<img src="ResultadosII_files/figure-epub3/FIG29-1.png" alt="**Abundancia del tipo de sustitución por cada nodo segun el marco de lectura.**." width="110%" /><img src="ResultadosII_files/figure-epub3/FIG29-2.png" alt="**Abundancia del tipo de sustitución por cada nodo segun el marco de lectura.**." width="110%" /><img src="ResultadosII_files/figure-epub3/FIG29-3.png" alt="**Abundancia del tipo de sustitución por cada nodo segun el marco de lectura.**." width="110%" />
<p class="caption">(\#fig:FIG29)**Abundancia del tipo de sustitución por cada nodo segun el marco de lectura.**.</p>
</div>

#### Análisis de sitios en los cuales su ancestro era TGGCGCCA

Para tratar de entender como es que los sitios HIP1 se pierden hicimos un análisis unicamente en en las transiciones en las que el nodo ancestral tenia un sitio TGGCGCCA.

En la (Figura \@ref(fig:FIG30)) mostramos 3 gráficos que indican la frecuencia del tipo de sustituciones que hubo para estos casos para cada nodo en cada uno de los marcos de lectura.

En la (Figura \@ref(fig:FIG31)) mostramos 3 gráficos (uno por cada marco de lectura) que indican la frecuencia de las mutaciones en cada uno de los 8 nucleótidos del sitio HIP.

En la (Figura \@ref(fig:FIG32)) mostramos 3 gráficos (uno por cada marco de lectura) que indican la frecuencia del tipo sustitucion de bases.

<div class="figure" style="text-align: center">
<img src="ResultadosII_files/figure-epub3/FIG30-1.png" alt="**Abundancia del tipo de sustitución por cada nodo segun el marco de lectura. Unicamente para transiciones en los que el nodo ancestral era un sitio TGGCGCCA.**" width="110%" /><img src="ResultadosII_files/figure-epub3/FIG30-2.png" alt="**Abundancia del tipo de sustitución por cada nodo segun el marco de lectura. Unicamente para transiciones en los que el nodo ancestral era un sitio TGGCGCCA.**" width="110%" /><img src="ResultadosII_files/figure-epub3/FIG30-3.png" alt="**Abundancia del tipo de sustitución por cada nodo segun el marco de lectura. Unicamente para transiciones en los que el nodo ancestral era un sitio TGGCGCCA.**" width="110%" />
<p class="caption">(\#fig:FIG30)**Abundancia del tipo de sustitución por cada nodo segun el marco de lectura. Unicamente para transiciones en los que el nodo ancestral era un sitio TGGCGCCA.**</p>
</div>

<div class="figure" style="text-align: center">
<img src="ResultadosII_files/figure-epub3/FIG31-1.png" alt="**Frecuencia de las mutaciones de cada nucleótido del sitio TGGCGCCA para cada nodo segun el marco de lectura.**." width="110%" /><img src="ResultadosII_files/figure-epub3/FIG31-2.png" alt="**Frecuencia de las mutaciones de cada nucleótido del sitio TGGCGCCA para cada nodo segun el marco de lectura.**." width="110%" /><img src="ResultadosII_files/figure-epub3/FIG31-3.png" alt="**Frecuencia de las mutaciones de cada nucleótido del sitio TGGCGCCA para cada nodo segun el marco de lectura.**." width="110%" />
<p class="caption">(\#fig:FIG31)**Frecuencia de las mutaciones de cada nucleótido del sitio TGGCGCCA para cada nodo segun el marco de lectura.**.</p>
</div>

<div class="figure" style="text-align: center">
<img src="ResultadosII_files/figure-epub3/FIG32-1.png" alt="**Frecuencia del tipo de sustituciónes de base en los sitios TGGCGCCA para cada marco de lectura**." width="110%" /><img src="ResultadosII_files/figure-epub3/FIG32-2.png" alt="**Frecuencia del tipo de sustituciónes de base en los sitios TGGCGCCA para cada marco de lectura**." width="110%" /><img src="ResultadosII_files/figure-epub3/FIG32-3.png" alt="**Frecuencia del tipo de sustituciónes de base en los sitios TGGCGCCA para cada marco de lectura**." width="110%" />
<p class="caption">(\#fig:FIG32)**Frecuencia del tipo de sustituciónes de base en los sitios TGGCGCCA para cada marco de lectura**.</p>
</div>

#### Análisis de sitios en los cuales solo el nodo actual tiene TGGCGCCA

Para tratar de entender como es que los sitios TGGCGCCA se ganan, hicimos un analisis unicamente en las transiciones en las que el nodo actual tenia un sitio TGGCGCCA. 

En la Figura \@ref(fig:FIG33) mostramos 3 gráficos (uno por cada marco de lectura) que indican la frecuencia del tipo de sustituciones que hubo para estos casos para cada nodo en cada uno de los marcos de lectura.

<div class="figure" style="text-align: center">
<img src="ResultadosII_files/figure-epub3/FIG33-1.png" alt="**Abundancia del tipo de sustitución por cada nodo segun el marco de lectura. Unicamente para transiciones en los que el nodo actual era un sitio TGGCGCCA.**." width="110%" /><img src="ResultadosII_files/figure-epub3/FIG33-2.png" alt="**Abundancia del tipo de sustitución por cada nodo segun el marco de lectura. Unicamente para transiciones en los que el nodo actual era un sitio TGGCGCCA.**." width="110%" /><img src="ResultadosII_files/figure-epub3/FIG33-3.png" alt="**Abundancia del tipo de sustitución por cada nodo segun el marco de lectura. Unicamente para transiciones en los que el nodo actual era un sitio TGGCGCCA.**." width="110%" />
<p class="caption">(\#fig:FIG33)**Abundancia del tipo de sustitución por cada nodo segun el marco de lectura. Unicamente para transiciones en los que el nodo actual era un sitio TGGCGCCA.**.</p>
</div>


