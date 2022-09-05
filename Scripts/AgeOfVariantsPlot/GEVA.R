#S'ha d'executar a andromeda després s'haver executat globalFunctions.R per poder accedir a les dades necessàries
#S'obtenen les edats estimades de GEVA
###S'especifica una regió i s'obté un gràfic on es representa, per totes les poblacions que tinguin al menys una variant amb iSAFE significatiu a la regió, un boxplot amb l'edat estimada de les variants amb iSAFE no significatiu vs. les que tenen iSAFE significatiu
#Dins de les variants amb iSAFE significatiu, s'escullen les 15 amb iSAFE superior i es descarten la resta. Aquest número es pot canviar (p. ex.: 5). 
#Compara la distribució d'edats de les variants significatives de les poblacions d'una metapoblació amb un test de Kruskal-Wallis i en representa el p-valor.


###Gràfic per plotejar una regió concreta on es mira si l'edat estimada de les variants amb iSAFE significatiu és coherent entre les poblacions d'una metapoblació AMB KRUSKAL-WILLIS i LEVENE

#S'han d'introduir també les taules atlas i isafe
diferencia_edats <- function(chrom, inici, final, pop_color = popPal, 
                             AFR = AFRpops, EAS = EASpops, EUR = EURpops, SAS = SASpops, AMR = TRUE) {
    library(ggsignif)
    #Guardo la regió que es vol representar
    region <- paste0(chrom, ':', inici, '-', final)
    #Em quedo amb la regió corresponent de les taules (i filtro per qualitat de GEVA > 0.5)
    merged_temp <- getInfoRegion(region) %>% filter(QualScore_Jnt >= 0.5) 
    #Si AMR == FALSE vol dir que s'ha d'eliminar la informació de les poblacions americanes
    if (AMR == FALSE) {
        merged_temp <- merged_temp[, -c("CLM", "MXL", "PEL", "PUR")]
    }
    #Faig el pivot longer, filtro els valors NA i classifico en iSAFE significant o no
    merged_temp <- merged_temp %>%
        pivot_longer(cols = c('ACB':'YRI'), names_to = 'pop', values_to='isafe') %>%
        filter(!is.na(isafe)) %>% mutate(grup = case_when(isafe >= 0.1 ~ 'Significant', TRUE ~ 'Nonsignificant'))
    #Això serveix per conseguir un vector amb les poblacions que tenen alguna variant significatica en iSAFE per aquesta regió
    pob_significatiu <- merged_temp %>% group_by(pop) %>% filter(grup == 'Significant') %>%
        summarise(significatiu = length(grup), .groups = 'keep') %>% collect %>% .[[1]]
    #Si no hi ha cap població significativa, s'acaba ja i no es pot representar res
    if (length(pob_significatiu) == 0) {
        print('This region cannot be represented, as no populations have variants with a significant iSAFE value')
    }
    #Filtro la taula per quedar-me només les columnes importants i les poblacions amb iSAFEs significatius
    #Una vegada filtrada, creo una columna nova que distingeixi la metapoblació segons la població i que identifiqui cada combinació de població i rsid individualment
    merged_temp <- merged_temp %>% select(rsid, pop, isafe, AgeMode_Jnt, grup) %>%
        filter(pop %in% pob_significatiu) %>% mutate(metapop = case_when(pop %in% AFR ~ 'AFR',
                                                                         pop %in% EUR ~ 'EUR',
                                                                         pop %in% EAS ~ 'EAS',
                                                                         pop %in% SAS ~ 'SAS',
                                                                         TRUE ~ 'AMR'),
                                                    rsid_pop = paste(rsid, pop, sep=':'))
    #Creo una taula on hi hagi les files corresponents als 15 valors significatius d'iSAFE més alts en cada població. SI n'hi ha menys de 15, es guarden els que hi hagi. Si hi ha empats, se'n queda > 15
    isafe15 <- merged_temp %>% filter(grup == 'Significant') %>% group_by(pop) %>% slice_max(isafe, n=15)
    #Ara de la taula principal, elimino les files que no estiguin a isafe15
    merged_temp <- merged_temp %>% filter(grup == 'Nonsignificant' | rsid_pop %in% isafe15$rsid_pop)
    #Faig una taula on posar (després) el resultat del test Kruskal-Wallis per metapoblació i la primera i última població de la metapop per representar-ho
    kruskal <- merged_temp %>% filter(isafe >= 0.1) %>% group_by(metapop) %>% arrange(pop) %>%
        summarise(primer = unique(pop)[1], ultim = unique(pop)[n_distinct(pop)], .groups = 'keep')
    #Aqui guardaré els p_valors de l'anàlisi kruskal-wallis de cada població que es pugui. No puc filtrar els iSAFE < 0.1 de la taula principal perquè els necessito pel gràfic
    p_vals <- c('AFR' = NA, 'EUR' = NA, 'EAS' = NA, 'SAS' = NA, 'AMR' = NA)
    for (i in 1:5) {
        merged_temp_temp <- merged_temp  %>% filter(grup == 'Significant' & metapop == names(p_vals)[i])
        #Si no hi ha valors per aquestes posicions, es queda el NA
        if (nrow(merged_temp_temp) == 0) {
            next
        }
        #Si no es pot fer cap anàlisi, es posa un None
        if (n_distinct(merged_temp_temp$pop) == 1 | n_distinct(merged_temp_temp$rsid) == 1) {
            p_vals[i] <- 'None'
            next
        }
        #Per poder fer el test de Kruskal-Wallis, fa falta fer un test d'homogeneitat i que no hi hagin diferències significatives entre les variàncies de les diferents poblacions
        #Faig un test de Levene amb la mediana perquè és no paramètric i no passa res perquè no segueixi una distribució normal i accedeixo al seu p-valor
        p_valor_levene <- levene.test(y = merged_temp_temp$AgeMode_Jnt, group = merged_temp_temp$pop, location = "median")[["p.value"]]
        if (is.na(p_valor_levene)) {
            p_vals[i] <- 'None'
        } else if (p_valor_levene < 0.05) {
            p_vals[i] <- 'Non-homogeneous'
            next
        }
        #Si es fa l'anàlisi, es guarda el p_valor al vector
        p_vals[i] <- round(kruskal.test(AgeMode_Jnt ~ pop, data = merged_temp_temp)[['p.value']], digits = 5)
    }
    #Afegeixo els p_valors a la taula kruskal
    kruskal <- kruskal %>% mutate(p_valor = p_vals[metapop])
    #Faig la llista de comparacions que s'ha de donar al geom_signiff
    comparisons_list <- list()
    for (n in 1:nrow(kruskal)) {
        comparisons_list[[n]] <- c(kruskal$primer[n], kruskal$ultim[n])
    }
    annotations_vector <- kruskal$p_valor
    #Faig això per tenir a pop_color només les poblacions que es representaran, i ordenades per metapoblació
    pops <- unique(merged_temp %>% arrange(metapop) %>% pull(pop))
    pop_color <- pop_color[pops]
    #Faig el gràfic. Convertir la x en factor després d'ordenar per pop i per metapop serveix per tenir les poblacions ordenades per metapoblació
    #i per la representació de geom_signif, perquè tenim la primer i última població de cada metapop alfabèticament
    options(repr.plot.width = 20, repr.plot.height = 10, warn = 1)
    grafic <- merged_temp %>% arrange(pop) %>% arrange(metapop) %>%
        ggplot(mapping = aes(x=factor(pop, levels=unique(pop)), y=AgeMode_Jnt/1000, fill=pop, color = grup, alpha = grup, by = metapop)) +
        geom_boxplot() + scale_alpha_manual(values = c(0.35, 1), name= 'iSAFE') + scale_color_manual(values = c('grey','#6b0000'), name = 'iSAFE') +
        scale_fill_manual(values = pop_color, breaks = names(pop_color), name = 'Population') +
        labs(x='', y='Estimated variant age (thousands of generations)', title = paste0('Chr ', chrom, ': ', inici, '-', final)) +
        geom_signif(comparisons = comparisons_list, annotations = annotations_vector, textsize = 7, color = 'black', alpha = 1) + theme(text = element_text(size = 20))
    return(grafic)
}
