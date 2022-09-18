#Codi per analitzar les regions de PopHumanScan on no es pot estimar l'edat d'un esdeveniment de selecció natural per qualsevol de les causes
#i per comparar aquestes regions amb les analitzades en termes d'iHS i nSL

#És necessària la taula de resultats, resultats filtrats i de PopHumanScan,
#així com les funcions dels diagrames de Venn de l'arxiu Taula_Grafics_GeneOntology.R d'aquesta mateixa carpeta

##############################################################################
###################### Taula de regions descartades ##########################
############################################################################## 


#Aquesta funció serveix per filtrar les regions que NO se seleccionen amb l'anàlisi 
variantage_table_discarded <- function(resultats, resultats_filtrats, PHS = phs) {
    #Em quedo amb les regions que no es trobin a la taula filtrada
    PHS <- PHS %>% select(c("GeneID", "sweepGroup", "chr", "start", "end", "filterTime", "rawMetapopulation")) %>% filter(!sweepGroup %in% resultats_filtrats$ID) %>%
    #Afegeixo la columna causa
        mutate(Cause = resultats %>% filter(!ID %in% resultats_filtrats$ID) %>% pull(Cause))
    #Els valors NA de Cause (són els que sí que tenen variants significatives però o no es pot fer ANOVA o mostra diferències significatives) els canvio
    for (n in 1:nrow(PHS)) {
        if (is.na(PHS$Cause[n])) {
            id <- PHS$sweepGroup[n]
            taula <- resultats[ID == id, c('AFR_test','EUR_test','EAS_test','SAS_test')]
            #Miro si les 4 columnes tenen valors NA o None, que no s'ha fet cap ANOVA
            if (all(is.na(taula) | taula == 'None')) {
                PHS$Cause[n] <- 'No test'
            } else if (all(is.na(taula) | taula == 'None' | taula == "Non-homogeneous")) {
                #En aquest cas, significa que no s'ha fet cap test, i les metapoblacions que en podrien haver fet, no podien per l'heterogeneitat de les variancies
                PHS$Cause[n] <- 'Non-homogeneous'
            } else {
                #Si no és el segon cas, singnifica que s'ha pogut fer algun ANOVA pero que ha donat significatiu
                PHS$Cause[n] <- 'Kruskal signif'
            }
        }
    }
    return(PHS)
}
discarded_phs <- variantage_table_discarded(resultats, resultats_filtrats)

###############################################################################
############################ Diagrama de Venn #################################
###############################################################################

#Diagrama de Venn per les regions descartades amb l'anàlisi fet segons la/les metapoblació que s'ha detectat amb popHumanScan
venn_phs_pops(discarded_phs)



###############################################################################
##################### Causa de les regions descartades ########################
###############################################################################

#Representació de les causes per les quals es descarten les altres regions de PHS en l'anàlisi
options(repr.plot.width = 25, repr.plot.height = 10, warn = 1)
discarded_phs %>% ggplot(mapping = aes(y = Cause, fill = Cause)) + geom_bar(width = 0.7, show.legend = FALSE) +
    labs(x= '# of regions', y = 'Cause', title = 'Cause of discarded regions') + 
    theme(text = element_text(size = 25)) + scale_y_discrete(labels= c("No variants" = "No variants in the region", "No test" = "No Kruskal-Wallis\n test performed",
                                                                       "N.S." = "No variants with\nsignificant iSAFE", "Kruskal signif" = "Significant differences\namong populations",
                                                                       "Non-homogeneous" = "Non-homogeneous variances"),
                                                            limits = c("N.S.","No test","No variants", "Kruskal signif", "Non-homogeneous")) +
    stat_count(geom = "text", colour = "black", size = 6, aes(label = ..count..),position=position_stack(vjust=0.5)) + scale_fill_brewer(palette="Set2")


###############################################################################
################# Anàlisi d'iHS de les regions descartades ####################
###############################################################################

#Afegiré a les taules la informació de ihs de la regió: iHS màxim i mínim i el nombre de variants amb la regió amb iHS >2 o <2 per comparar entre regions seleccionades o no

#Aquesta funció farà un histograma de la distribució dels valors d'iHS per comparar entre regions seleccionades o no amb l'anàlisi 
perfil_ihs <- function(taula, threshold_ihs = 2,
                       AFR = AFRpops, EAS = EASpops, EUR = EURpops, SAS = SASpops) {
    #Afegeixo a la taula columnes per guardar tota la informació necessària: 4 columnes per població
    taula <- taula %>% select(GeneID, sweepGroup, chr, start, end, Cause) %>% 
        mutate('ACB__max_ihs' = as.numeric(NA), 'ASW__max_ihs' = as.numeric(NA), 'ESN__max_ihs' = as.numeric(NA), 'GWD__max_ihs' = as.numeric(NA),
               'LWK__max_ihs' = as.numeric(NA), 'MSL__max_ihs' = as.numeric(NA), 'YRI__max_ihs' = as.numeric(NA), 'CDX__max_ihs' = as.numeric(NA),
               'CHB__max_ihs' = as.numeric(NA), 'CHS__max_ihs' = as.numeric(NA), 'JPT__max_ihs' = as.numeric(NA), 'KHV__max_ihs' = as.numeric(NA), 
               'CEU__max_ihs' = as.numeric(NA), 'FIN__max_ihs' = as.numeric(NA), 'GBR__max_ihs' = as.numeric(NA), 'IBS__max_ihs' = as.numeric(NA),
               'TSI__max_ihs' = as.numeric(NA), 'BEB__max_ihs' = as.numeric(NA), 'GIH__max_ihs' = as.numeric(NA), 'ITU__max_ihs' = as.numeric(NA),
               'PJL__max_ihs' = as.numeric(NA), 'STU__max_ihs' = as.numeric(NA), 'ACB__min_ihs' = as.numeric(NA), 'ASW__min_ihs' = as.numeric(NA),
               'ESN__min_ihs' = as.numeric(NA), 'GWD__min_ihs' = as.numeric(NA), 'LWK__min_ihs' = as.numeric(NA), 'MSL__min_ihs' = as.numeric(NA),
               'YRI__min_ihs' = as.numeric(NA), 'CDX__min_ihs' = as.numeric(NA), 'CHB__min_ihs' = as.numeric(NA), 'CHS__min_ihs' = as.numeric(NA),
               'JPT__min_ihs' = as.numeric(NA), 'KHV__min_ihs' = as.numeric(NA), 'CEU__min_ihs' = as.numeric(NA), 'FIN__min_ihs' = as.numeric(NA),
               'GBR__min_ihs' = as.numeric(NA), 'IBS__min_ihs' = as.numeric(NA), 'TSI__min_ihs' = as.numeric(NA), 'BEB__min_ihs' = as.numeric(NA),
               'GIH__min_ihs' = as.numeric(NA), 'ITU__min_ihs' = as.numeric(NA), 'PJL__min_ihs' = as.numeric(NA), 'STU__min_ihs' = as.numeric(NA),
               'ACB__npos_ihs' = as.numeric(NA), 'ASW__npos_ihs' = as.numeric(NA), 'ESN__npos_ihs' = as.numeric(NA), 'GWD__npos_ihs' = as.numeric(NA),
               'LWK__npos_ihs' = as.numeric(NA), 'MSL__npos_ihs' = as.numeric(NA), 'YRI__npos_ihs' = as.numeric(NA), 'CDX__npos_ihs' = as.numeric(NA),
               'CHB__npos_ihs' = as.numeric(NA), 'CHS__npos_ihs' = as.numeric(NA), 'JPT__npos_ihs' = as.numeric(NA), 'KHV__npos_ihs' = as.numeric(NA), 
               'CEU__npos_ihs' = as.numeric(NA), 'FIN__npos_ihs' = as.numeric(NA), 'GBR__npos_ihs' = as.numeric(NA), 'IBS__npos_ihs' = as.numeric(NA),
               'TSI__npos_ihs' = as.numeric(NA), 'BEB__npos_ihs' = as.numeric(NA), 'GIH__npos_ihs' = as.numeric(NA), 'ITU__npos_ihs' = as.numeric(NA),
               'PJL__npos_ihs' = as.numeric(NA), 'STU__npos_ihs' = as.numeric(NA), 'ACB__nneg_ihs' = as.numeric(NA), 'ASW__nneg_ihs' = as.numeric(NA),
               'ESN__nneg_ihs' = as.numeric(NA), 'GWD__nneg_ihs' = as.numeric(NA), 'LWK__nneg_ihs' = as.numeric(NA), 'MSL__nneg_ihs' = as.numeric(NA),
               'YRI__nneg_ihs' = as.numeric(NA), 'CDX__nneg_ihs' = as.numeric(NA), 'CHB__nneg_ihs' = as.numeric(NA), 'CHS__nneg_ihs' = as.numeric(NA),
               'JPT__nneg_ihs' = as.numeric(NA), 'KHV__nneg_ihs' = as.numeric(NA), 'CEU__nneg_ihs' = as.numeric(NA), 'FIN__nneg_ihs' = as.numeric(NA),
               'GBR__nneg_ihs' = as.numeric(NA), 'IBS__nneg_ihs' = as.numeric(NA), 'TSI__nneg_ihs' = as.numeric(NA), 'BEB__nneg_ihs' = as.numeric(NA),
               'GIH__nneg_ihs' = as.numeric(NA), 'ITU__nneg_ihs' = as.numeric(NA), 'PJL__nneg_ihs' = as.numeric(NA), 'STU__nneg_ihs' = as.numeric(NA))
    #Per cada regió, calcularem els valors mitjans d'iHS de les poblacions que es pugui
    for (n in 1:nrow(taula)) {
        chr <- gsub('chr', '', taula$chr[n])
        inici <- taula$start[n]
        final <- taula$end[n]
        #Obtinc la taula d'iHS per les variants de la regió
        iHS <- sqlToDf(chr, inici, final, 'ihs')
        #Si no hi hagués cap variant, ho deixo tot en NA i passo a la següent iteració (potser hauria de posar 0 a les npos i nneg però no sé)
        if (nrow(iHS) == 0) {
            next
        }
        #Si hi ha variants, faig un pivot_longer per tenir cada variant de cada població en un registre i elimino els valors NA
        iHS <- iHS %>% pivot_longer(cols = 'ACB':'YRI', names_to = 'pop', values_to = 'ihs') %>% filter(!is.na(ihs)) %>%
            #llavors, elimino les poblacions americanes, agrupo per població i en calculo el mínim, màxim i el nombre de variants que superin els valors límits 
            filter(!pop %in% c('CLM', 'MXL', 'PEL', 'PUR')) %>% group_by(pop) %>% summarise(max_ihs = max(ihs), min_ihs = min(ihs),
                                                                                            npos_ihs = sum(ihs >= threshold_ihs), nneg_ihs = sum(ihs <= - threshold_ihs), .groups = 'drop')
        #Ara he de guardar cadascun dels 4 valors a la taula principal a la població corresponent
        for (p in 1:nrow(iHS)) {
            pop_ihs <- iHS$pop[p]
            #Per cada fila, he de guardar els 4 valors:
            for (valor in c('max_ihs', 'min_ihs', 'npos_ihs', 'nneg_ihs')) {
                columna <- paste(pop_ihs, valor, sep='__')
                taula[n, columna] <- iHS[p, valor]
            }
        }
        #Ara ja tenim les dades d'iHS per cada població que es pot da la taula principal
    }
    #Faig un pivot longer per tenir cada població i regió (i dada) en un registre i filtro les files on no hi hagi cap dada d'iHS en aquella població i regió
    #Pero realment el que m'interessa és tenir 4 columnes per les dades i cada poblacioó en una fila. Separaré en 2 columnes pop.data
    taula <- taula %>% pivot_longer(cols = 'ACB__max_ihs':'STU__nneg_ihs', names_to = 'pop_data', values_to = 'valor') %>% filter(!is.na(valor)) %>%
        separate(pop_data, sep = "__", into = c("pop", "data")) %>% 
        #Llavors faig un pivot_wider per tenir cada dada en una columna i un registre per regió i població; i afegeixo la columna metapop per saber a quina metapoblació es correspon
        # i afegeixo una columna que miri si se seleccionen més els al·lels derivats (npos_ihs>nneg:ihs) o els ancestrals (npos_ihs < nneg_ihs)
        pivot_wider(names_from = data, values_from = valor) %>% mutate(metapop = case_when(pop %in% AFR ~ 'AFR', pop %in% EUR ~ 'EUR', pop %in% EAS ~ 'EAS', pop %in% SAS ~ 'SAS', TRUE ~ 'AMR'),
                                                                       selected_alleles = as.character(NA))
    for (n in 1:nrow(taula)) {
        if (taula[n, 'npos_ihs'] > taula[n, 'nneg_ihs']) {
            taula[n, 'selected_alleles'] <- 'Derived'
        } else if (taula[n, 'npos_ihs'] < taula[n, 'nneg_ihs']) {
            taula[n, 'selected_alleles'] <- 'Ancestral'
        } else if (taula[n, 'npos_ihs'] == 0) {
            taula[n, 'selected_alleles'] <- 'None'
        }
    }
    return(taula)
}

#S0aplica tant a la taula de les regions filtrades com les descartades per comparar
ihs_filtrats <- perfil_ihs(resultats_filtrats_phs %>% mutate(Cause = 0))
ihs_discarded <- perfil_ihs(discarded_phs)

#Impressió dels resultats:
for (causa in unique(ihs_discarded$Cause)) {
    ihs_temp <- ihs_discarded %>% filter(Cause == causa)
    cat('\nCAUSA: ', causa, '\n')
    cat("- Valor màxim d'iHS\n")
    print(summary(ihs_temp$max_ihs))
    cat("- Valor mínim d'iHS\n")
    print(summary(ihs_temp$min_ihs))
    cat("- Nombre de variants amb iHS > 2\n")
    print(summary(ihs_temp$npos_ihs))
    cat("- Nombre de variants amb iHS < 2\n")
    print(summary(ihs_temp$nneg_ihs))
    cat("- Al·lels seleccionats\n")
    print(ihs_temp %>% group_by(selected_alleles) %>% summarise(percent= n()/nrow(ihs_temp)*100, .groups = 'drop'))
    rm(ihs_temp)
}

#Ara ja es pot representar l'histograma, amb aquestes dades
#Ajunto les dos taules, afegint abans una columna que diferencii de quina venien
taula <- rbind(ihs_filtrats, ihs_discarded) 
options(repr.plot.width = 25, repr.plot.height = 7, warn = 1)
grafic <- taula %>% ggplot(mapping = aes(x = Cause, fill = selected_alleles)) + geom_bar() + 
    labs(x= 'Selected regions', y = '# of regions', title = 'Selected alleles according to significant (>2) iHS scores') +
    facet_wrap(~ Cause, ncol = 6,  scales = 'free', labeller = labeller(Cause = c('0' = 'Selected', 'Kruskal signif' = 'Significant Kruskal-Wallis', 'N.S.' = 'No variants with significant iSAFE', 
                                                                                 'No test' = 'No test performed', 'No variants' = 'No variants', 'Non-homogeneous' = 'Non-homogeneous'))) +
    theme(text = element_text(size = 20))
    #scale_fill_manual(values = metapopsData$colors, name = 'Population')  +
    #Faig servir facet_grid2 de ggh4x perquè deixa 'alliberar' l'eix y, mentre que facet_grid normal no deixa que canviï en una mateixa fila, que és el que necessito (amb l'argument independent)
grafic

###############################################################################
################# Anàlisi d'nSL de les regions descartades ####################
###############################################################################

#Afegiré a les taules la informació de nsl de la regió: nSL màxim i mínim i el nombre de variants amb la regió amb nSL >2 o <2 per comparar entre regions seleccionades o no

#Aquesta funció farà un histograma de la distribució dels valors d'nSL per comparar entre regions seleccionades o no amb l'anàlisi 
perfil_nsl <- function(taula, threshold_nsl = 2,
                       AFR = AFRpops, EAS = EASpops, EUR = EURpops, SAS = SASpops) {
    #Afegeixo a la taula columnes per guardar tota la informació necessària: 4 columnes per població
    taula <- taula %>% select(GeneID, sweepGroup, chr, start, end, Cause) %>% 
        mutate('ACB__max_nsl' = as.numeric(NA), 'ASW__max_nsl' = as.numeric(NA), 'ESN__max_nsl' = as.numeric(NA), 'GWD__max_nsl' = as.numeric(NA),
               'LWK__max_nsl' = as.numeric(NA), 'MSL__max_nsl' = as.numeric(NA), 'YRI__max_nsl' = as.numeric(NA), 'CDX__max_nsl' = as.numeric(NA),
               'CHB__max_nsl' = as.numeric(NA), 'CHS__max_nsl' = as.numeric(NA), 'JPT__max_nsl' = as.numeric(NA), 'KHV__max_nsl' = as.numeric(NA), 
               'CEU__max_nsl' = as.numeric(NA), 'FIN__max_nsl' = as.numeric(NA), 'GBR__max_nsl' = as.numeric(NA), 'IBS__max_nsl' = as.numeric(NA),
               'TSI__max_nsl' = as.numeric(NA), 'BEB__max_nsl' = as.numeric(NA), 'GIH__max_nsl' = as.numeric(NA), 'ITU__max_nsl' = as.numeric(NA),
               'PJL__max_nsl' = as.numeric(NA), 'STU__max_nsl' = as.numeric(NA), 'ACB__min_nsl' = as.numeric(NA), 'ASW__min_nsl' = as.numeric(NA),
               'ESN__min_nsl' = as.numeric(NA), 'GWD__min_nsl' = as.numeric(NA), 'LWK__min_nsl' = as.numeric(NA), 'MSL__min_nsl' = as.numeric(NA),
               'YRI__min_nsl' = as.numeric(NA), 'CDX__min_nsl' = as.numeric(NA), 'CHB__min_nsl' = as.numeric(NA), 'CHS__min_nsl' = as.numeric(NA),
               'JPT__min_nsl' = as.numeric(NA), 'KHV__min_nsl' = as.numeric(NA), 'CEU__min_nsl' = as.numeric(NA), 'FIN__min_nsl' = as.numeric(NA),
               'GBR__min_nsl' = as.numeric(NA), 'IBS__min_nsl' = as.numeric(NA), 'TSI__min_nsl' = as.numeric(NA), 'BEB__min_nsl' = as.numeric(NA),
               'GIH__min_nsl' = as.numeric(NA), 'ITU__min_nsl' = as.numeric(NA), 'PJL__min_nsl' = as.numeric(NA), 'STU__min_nsl' = as.numeric(NA),
               'ACB__npos_nsl' = as.numeric(NA), 'ASW__npos_nsl' = as.numeric(NA), 'ESN__npos_nsl' = as.numeric(NA), 'GWD__npos_nsl' = as.numeric(NA),
               'LWK__npos_nsl' = as.numeric(NA), 'MSL__npos_nsl' = as.numeric(NA), 'YRI__npos_nsl' = as.numeric(NA), 'CDX__npos_nsl' = as.numeric(NA),
               'CHB__npos_nsl' = as.numeric(NA), 'CHS__npos_nsl' = as.numeric(NA), 'JPT__npos_nsl' = as.numeric(NA), 'KHV__npos_nsl' = as.numeric(NA), 
               'CEU__npos_nsl' = as.numeric(NA), 'FIN__npos_nsl' = as.numeric(NA), 'GBR__npos_nsl' = as.numeric(NA), 'IBS__npos_nsl' = as.numeric(NA),
               'TSI__npos_nsl' = as.numeric(NA), 'BEB__npos_nsl' = as.numeric(NA), 'GIH__npos_nsl' = as.numeric(NA), 'ITU__npos_nsl' = as.numeric(NA),
               'PJL__npos_nsl' = as.numeric(NA), 'STU__npos_nsl' = as.numeric(NA), 'ACB__nneg_nsl' = as.numeric(NA), 'ASW__nneg_nsl' = as.numeric(NA),
               'ESN__nneg_nsl' = as.numeric(NA), 'GWD__nneg_nsl' = as.numeric(NA), 'LWK__nneg_nsl' = as.numeric(NA), 'MSL__nneg_nsl' = as.numeric(NA),
               'YRI__nneg_nsl' = as.numeric(NA), 'CDX__nneg_nsl' = as.numeric(NA), 'CHB__nneg_nsl' = as.numeric(NA), 'CHS__nneg_nsl' = as.numeric(NA),
               'JPT__nneg_nsl' = as.numeric(NA), 'KHV__nneg_nsl' = as.numeric(NA), 'CEU__nneg_nsl' = as.numeric(NA), 'FIN__nneg_nsl' = as.numeric(NA),
               'GBR__nneg_nsl' = as.numeric(NA), 'IBS__nneg_nsl' = as.numeric(NA), 'TSI__nneg_nsl' = as.numeric(NA), 'BEB__nneg_nsl' = as.numeric(NA),
               'GIH__nneg_nsl' = as.numeric(NA), 'ITU__nneg_nsl' = as.numeric(NA), 'PJL__nneg_nsl' = as.numeric(NA), 'STU__nneg_nsl' = as.numeric(NA))
    #Per cada regió, calcularem els valors mitjans d'sS de les poblacions que es pugui
    for (n in 1:nrow(taula)) {
        chr <- gsub('chr', '', taula$chr[n])
        inici <- taula$start[n]
        final <- taula$end[n]
        #Obtinc la taula d'nSL per les variants de la regió
        nSL <- sqlToDf(chr, inici, final, 'nsl')
        #Si no hi hagués cap variant, ho deixo tot en NA i passo a la següent iteració (potser hauria de posar 0 a les npos i nneg però no sé)
        if (nrow(nSL) == 0) {
            next
        }
        #Si hi ha varianrs, faig un pivot_longer per tenir cada variant de cada població en un registre i elimino els valors NA
        nSL <- nSL %>% pivot_longer(cols = 'ACB':'YRI', names_to = 'pop', values_to = 'nsl') %>% filter(!is.na(nsl)) %>%
            #llavors, elimino les poblacions americanes, agrupo per població i en calculo el mínim, màxim i el nombre de variants que superin els valors límits 
            filter(!pop %in% c('CLM', 'MXL', 'PEL', 'PUR')) %>% group_by(pop) %>% summarise(max_nsl = max(nsl), min_nsl = min(nsl),
                                                                                            npos_nsl = sum(nsl >= threshold_nsl), nneg_nsl = sum(nsl <= - threshold_nsl), .groups = 'drop')
        #Ara he de guardar cadascun dels 4 valors a la taula principal a la població corresponent
        for (p in 1:nrow(nSL)) {
            pop_nsl <- nSL$pop[p]
            #Per cada fila, he de guardar els 4 valors:
            for (valor in c('max_nsl', 'min_nsl', 'npos_nsl', 'nneg_nsl')) {
                columna <- paste(pop_nsl, valor, sep='__')
                taula[n, columna] <- nSL[p, valor]
            }
        }
        #Ara ja tenim les dades d'nSL per cada població que es pot da la taula principal
    }
    #Faig un pivot longer per tenir cada població i regió (i dada) en un registre i filtro les files on no hi hagi cap dada d'nSL en aquella població i regió
    #Pero realment el que m'interessa és tenir 4 columnes per les dades i cada poblacioó en una fila. Separaré en 2 columnes pop.data
    taula <- taula %>% pivot_longer(cols = 'ACB__max_nsl':'STU__nneg_nsl', names_to = 'pop_data', values_to = 'valor') %>% filter(!is.na(valor)) %>%
        separate(pop_data, sep = "__", into = c("pop", "data")) %>% 
        #Llavors faig un pivot_wider per tenir cada dada en una columna i un registre per regió i població; i afegeixo la columna metapop per saber a quina metapoblació es correspon
        pivot_wider(names_from = data, values_from = valor) %>% mutate(metapop = case_when(pop %in% AFR ~ 'AFR', pop %in% EUR ~ 'EUR', pop %in% EAS ~ 'EAS', pop %in% SAS ~ 'SAS', TRUE ~ 'AMR'),
                                                                       selected_alleles = as.character(NA))
    for (n in 1:nrow(taula)) {
        if (taula[n, 'npos_nsl'] > taula[n, 'nneg_nsl']) {
            taula[n, 'selected_alleles'] <- 'Derived'
        } else if (taula[n, 'npos_nsl'] < taula[n, 'nneg_nsl']) {
            taula[n, 'selected_alleles'] <- 'Ancestral'
        } else if (taula[n, 'npos_nsl'] == 0) {
            taula[n, 'selected_alleles'] <- 'None'
        }
    }
    return(taula)
}

nsl_filtrats <- perfil_nsl(resultats_filtrats_phs %>% mutate(Cause = 0))
nsl_discarded <- perfil_nsl(discarded_phs)

#Impressió dels resultats
for (causa in unique(nsl_discarded$Cause)) {
    nsl_temp <- nsl_discarded %>% filter(Cause == causa)
    cat('\nCAUSA: ', causa, '\n')
    cat("- Valor màxim d'nSL\n")
    print(summary(nsl_temp$max_nsl))
    cat("- Valor mínim d'nSL\n")
    print(summary(nsl_temp$min_nsl))
    cat("- Nombre de variants amb nSL > 2\n")
    print(summary(nsl_temp$npos_nsl))
    cat("- Nombre de variants amb nSL < 2\n")
    print(summary(nsl_temp$nneg_nsl))
    cat("- Al·lels seleccionats\n")
    print(nsl_temp %>% group_by(selected_alleles) %>% summarise(percent= n()/nrow(nsl_temp)*100, .groups = 'drop'))
    rm(nsl_temp)
}

#Ara ja es pot representar l'histograma, amb aquestes dades. Ajunto les dos taules
taula <- rbind(nsl_filtrats, nsl_discarded)
options(repr.plot.width = 25, repr.plot.height = 7, warn = 1)
grafic <- taula %>% ggplot(mapping = aes(x = Cause, fill = selected_alleles)) + geom_bar() + 
    labs(x= 'Selected regions', y = '# of regions · population', title = 'Selected alleles according to significant (>2) nSL scores') +
    facet_wrap(~ Cause, ncol = 6,  scales = 'free', labeller = labeller(Cause = c('0' = 'Selected', 'Kruskal signif' = 'Significant Kruskal-Wallis', 'N.S.' = 'No variants with significant iSAFE', 
                                                                                 'No test' = 'No test performed', 'No variants' = 'No variants', 'Non-homogeneous' = 'Non-homogeneous'))) +
    theme(text = element_text(size = 20))
    #scale_fill_manual(values = metapopsData$colors, name = 'Population')  +
    #Faig servir facet_grid2 de ggh4x perquè deixa 'alliberar' l'eix y, mentre que facet_grid normal no deixa que canviï en una mateixa fila, que és el que necessito (amb l'argument independent)
rm(taula)
grafic



