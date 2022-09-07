# SelectiveSweeps
Recopilació del codi de les pràctiques curriculars per descriure els sweeps selectius en la història evolutiva de l'espècie humana

# Mapa del repositori
Tots els codis estan dins de la carpeta Scripts. Cal destacar que molts dels codis s'han d'executar dins del servidor d'andromeda, i que contenen paths a arxius que estan a la carpeta del meu usuari del servidor. Tots els codis estan comentats quan hi ha un path a un arxiu del meu usuari, i caldria executar el codi des d'allà o bé modificar-lo lleugerament per obrir els arxius necessaris des d'alguna altra ubicació. També s'ha de poder accedir a la base de dades d'SQL del servidor, on hi ha les taules necessàries amb informació de les variants. 

## globalFunctions.R
Un arxiu de codi que es necessita executar per poder funcionar la majoria de codis restants. Es crida a la majoria de paquets d'R que s'utilitzen en tots els anàlisis, i es creen alguns arxius amb informació de les poblacions i les metapoblacions del 1000 Genomes Project. També es defineixen funcions per poder accedir a les taules d'SQL del servidor amb informació d'iSAFE, GEVA, etc. 
![Geographical-distribution-of-26-populations-in-the-1000-Genomes-Project-The-major](https://user-images.githubusercontent.com/112875889/188856397-63e4e934-79da-4d51-b7ed-18cd22b0027c.png)

## AgeOfVariantsPlot
Codis per la representació de les edats de les variants de regions del genoma especificades (iSAFE no significatiu vs. significatiu) en boxplots.
### - PlotGEVA.R
Representa les variants d'UNA regió amb les edats de GEVA. Es fa l'anàlisi de la variància entre les poblacions d'una metapoblació amb Kruskal-Wallis, si és possible.
Exemple amb la regió del gen EDAR:
![image](https://user-images.githubusercontent.com/112875889/188856645-cf53c0bd-d6ee-45ba-a078-c400641803bc.png)

### - PlotRelate.R
Representa les variants d'UNA regió amb les edats de Relate. Es fa l'anàlisi de la variància entre les poblacions d'una metapoblació amb Kruskal-Wallis, si és possible. 
Exemple amb la regió del gen EDAR:
![image](https://user-images.githubusercontent.com/112875889/188856759-93ade5e1-2f95-406f-bae0-059b7a623074.png)

### - PlotsGEVA.R
Representa el mateix gràfic amb les edats de GEVA però per tantes regions com s'especifiqui (únicament regions preestablertes de PopHumanScan). Si només es diu un o diversos cromosomes, es representaran totes les regions que es puguin d'aquests cromsomes. 
Exemple amb el cromosoma 19:
![image](https://user-images.githubusercontent.com/112875889/188857403-b3b945a9-926f-4f42-9f0e-6913fa1a688c.png)

### - PlotsRelate.R
Representa el mateix gràfic amb les edats de Relate però per tantes regions com s'especifiqui (únicament regions preestablertes de PopHumanScan). Si només es diu un o diversos cromosomes, es representaran totes les regions que es puguin d'aquests cromsomes. 
Exemple amb el cromosoma 19:
![image](https://user-images.githubusercontent.com/112875889/188857565-feb574aa-0867-4e68-89b7-d6fb8d6e160d.png)

## PHSRegionsAnalysis
### - Taula_Grafics_GeneOntology.R

Hi ha el codi per tot l'anàlisi efectuat a partir de les regions de PopHumanScan. De totes les variants de les taules, és la que compara les poblacions d'una metapoblació amb Kruskal-Wallis (després de fer test de Levene per demostrar que les variàncies són homogènies), per tal de comprovar que hi ha una distribució uniforme de les edats de les variants amb iSAFE significatiu entre poblacions, i així calcular l'edat estimada del sweep selectiu per cada regió si és possible. 
També hi ha el codi per les diverses representacions, així com per l'anàlisi de Gene Ontology, tot i que l'anàlisi de Gene Ontology més complert és el que es fa per les regions de l'article d'iHS (Johnson and Voight, 2018). 

### - DiscardedRegions.R
Aquest codi serveix per analitzar amb més profunditat les regions de PopHumanScan que no es poden analitzar a la taula per diversos motius. També es comparen en termes de iHS i nSL les regions descartades i no descartades de PHS per observar possibles diferències (tot i que al final no es va trobar cap resultat important en aquest sentit). 
![image](https://user-images.githubusercontent.com/112875889/188866179-e6f3f126-fd37-4dfd-b63a-a220dc083203.png)


## ihsRegionsAnalysis
### - Taula_Grafics_GeneOntology.R

S'efectua un anàlisi molt similar al de les regions de PopHumanScan, però partint de les regions descrites en l'article de Johnson i Voight (2018). Alguns processos són lleugerament diferents, ja que en aquest cas, les regions candidates de selecció natural s'han identificat per a una població determinada, i per tant no té sentit fer comparacions entre les diferents poblacions. La distribució i la generació de les taules és relativament diferent, però l'anàlisi en sí és anàleg. 
Per cada regió, s'utilitzen les 5 variants amb iSAFE més elevat per estimar l'edat de l'esdeveniment de selecció, i es classifiquen les regions en intervals de 1000 generacions per fer l'anàlisi de Gene Ontology. Es filtren els resultats de Gene Ontology causats per clusters de gens amb funcions relacionades en una mateixa regió, que són la majoria, per obtenir resultats causats per gens independents.

### - ResultsExample.md
Exemple dels resultats obtinguts al portar a terme aquest anàlisis, amb les imatges obtingudes per clarificar l'explicació. 
![image](https://user-images.githubusercontent.com/112875889/188861981-5a742286-7d18-403d-aebf-5a0ab9bf57c5.png)
