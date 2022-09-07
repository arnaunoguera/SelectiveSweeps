# SelectiveSweeps
Recopilació del codi de les pràctiques curriculars per descriure els sweeps selectius en la història evolutiva de l'espècie humana

# Mapa del repositori

## Scripts
### - globalFunctions.R
Un arxiu de codi que es necessita executar per poder funcionar la majoria de codis restants. Es crida a la majoria de paquets d'R que s'utilitzen en tots els anàlisis, i es creen alguns arxius amb informació de les poblacions i les metapoblacions del 1000 Genomes Project. També es defineixen funcions per poder accedir a les taules d'SQL del servidor amb informació d'iSAFE, GEVA, etc. 

### - AgeOfVariantsPlot
Codis per la representació de les edats de les variants de regions del genoma especificades (iSAFE no significatiu vs. significatiu) en boxplots.
#### PlotGEVA.R
Representa les variants d'UNA regió amb les edats de GEVA. Es fa l'anàlisi de la variància entre les poblacions d'una metapoblació amb Kruskal-Wallis, si és possible. 
#### PlotRelate.R
Representa les variants d'UNA regió amb les edats de Relate. Es fa l'anàlisi de la variància entre les poblacions d'una metapoblació amb Kruskal-Wallis, si és possible. 
