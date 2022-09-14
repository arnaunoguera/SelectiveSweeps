#S'ha d'executar a andromeda després s'haver executat globalFunctions.R per poder accedir a les dades necessàries
#També es necessiten els fitxers amb les edats de Relate. El path del codi està al meu usuari del servidor
###S'especifica una regió i s'obté un gràfic on es representa, per totes les poblacions que tinguin al menys una variant amb iSAFE significatiu a la regió, un boxplot amb l'edat estimada de les variants amb iSAFE no significatiu vs. les que tenen iSAFE significatiu
#Dins de les variants amb iSAFE significatiu, s'escullen les 15 amb iSAFE superior i es descarten la resta. Aquest número es pot canviar (p. ex.: 5) a nvar. 
#Compara la distribució d'edats de les variants significatives de les poblacions d'una metapoblació amb un test de Kruskal-Wallis i en representa el p-valor.


diferencia_edats <- function(chr, inici, final, pop_color = popPal, 
                             AFR = AFRpops, EAS = EASpops, EUR = EURpops, SAS = SASpops, AMR = TRUE, nvar = 15) {
    library(ggsignif)
    #Obtinc la regió necessària de la taula isafe
    isafe <- sqlToDf(chr, inici, final, 'isafe')
    #Si AMR == FALSE vol dir que s'ha d'eliminar la informació de les poblacions americanes
    if (AMR == FALSE) {
        isafe <- isafe[, -c("CLM", "MXL", "PEL", "PUR")]
    }    
    #Faig el pivot longer a isafe per tenir un registre per cada població i variant, filtro els NA i classifico per si la variant és significant o no. 
    isafe <- isafe %>%
        pivot_longer(cols = c('ACB':'YRI'), names_to = 'pop', values_to='isafe') %>%
        filter(!is.na(isafe)) %>% mutate(grup = case_when(isafe >= 0.1 ~ 'Significant', TRUE ~ 'Nonsignificant'))
    #Això serveix per conseguir un vector amb les poblacions que tenen alguna variant significatica en iSAFE per aquesta regió
    pop_significatiu <- unique(isafe %>% filter(grup == 'Significant') %>% pull(pop))
    #Si no hi ha cap població significativa, s'acaba ja i no es pot representar res
    if (length(pop_significatiu) == 0) {
        return('This region cannot be represented, as no populations have variants with a significant iSAFE value')
    }
    #Filtro la taula per quedar-me només amb la informació de poblacions amb variants significatives i afegeixo a quina metapoblació pertanyen i la columna rsid_pop
    isafe <- isafe %>% filter(pop %in% pop_significatiu) %>% mutate(metapop = case_when(pop %in% AFR ~ 'AFR',
                                                                                        pop %in% EUR ~ 'EUR',
                                                                                        pop %in% EAS ~ 'EAS',
                                                                                        pop %in% SAS ~ 'SAS',
                                                                                        TRUE ~ 'AMR'),
                                                                    rsid_pop = paste(rsid, pop, sep=':'))
    #Obro la taula de les edats del cromosoma en qüestió i filtro per la regió i perquè la població sigui a pop_significatiu
    allele_ages <- fread(paste0("/home/anoguera/Data/relate_ages/allele_ages_chr", chr, '.csv')) %>% filter(BP >= inici & BP <= final) %>% select(BP, pop, lower_age, upper_age, pvalue)
    #Ara ajunto les dues taules en el merged
    merged <- merge(isafe, allele_ages, by.x = c('physicalPos', 'pop'), by.y = c('BP', 'pop'))
    #Una vegada fet el merge, em quedo amb les 15 variants més significatives per iSAFE (no ho faig abans perquè potser perdríem les variants significatives per no estar a la taula d'edats)
    #Creo una taula on hi hagi les files corresponents als 15 valors significatius d'iSAFE més alts en cada població. SI n'hi ha menys de 15, es guarden els que hi hagi. Si hi ha empats, se'n queda > 15
    isafe15 <- merged %>% select(grup, pop, isafe, rsid_pop) %>% filter(grup == 'Significant') %>% group_by(pop) %>% slice_max(isafe, n=nvar)
    #Ara de la taula principal, elimino les files que siguin significants  no estiguin a isafe15 
    merged <- merged %>% filter(grup == 'Nonsignificant' | rsid_pop %in% isafe15$rsid_pop)
    #En cas que no hi hagi cap variant significativa en tota la població després de tot, diem que no es pot representar la població
    pop_significatiu <- unique(merged %>% filter(grup == 'Significant') %>% pull(pop))
    if (length(pop_significatiu) == 0) {
        return('This region cannot be represented, as no populations have an estimated allele age for variants with a significant iSAFE value')
    }
    #Em quedo només amb les poblacions que tinguin variants significatives al final i calculo l'edat estimada de la variant, que és la mitjana de les edats superior i inferior (per mutacions neutres??)
    merged <- merged %>% filter(pop %in% pop_significatiu) %>% mutate(age = 0.5*(lower_age + upper_age))
    #Faig una taula on posar (després) el resultat del test Kruskal-Wallis per metapoblació i la primera i última població de la metapop per representar-ho
    kruskal <- merged %>% filter(grup == 'Significant') %>% group_by(metapop) %>% arrange(pop) %>%
        summarise(primer = unique(pop)[1], ultim = unique(pop)[n_distinct(pop)], .groups = 'keep')
    #Aqui guardaré els p_valors de l'anàlisi kruskal-wallis de cada població que es pugui. No puc filtrar els iSAFE < 0.1 de la taula principal perquè els necessito pel gràfic
    p_vals <- c('AFR' = NA, 'EUR' = NA, 'EAS' = NA, 'SAS' = NA, 'AMR' = NA)
    for (i in 1:5) {
        merged_temp_temp <- merged  %>% filter(grup == 'Significant' & metapop == names(p_vals)[i])
        #Si no hi ha valors per aquesta metapoblació, es queda el NA
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
        p_valor_levene <- levene.test(y = merged_temp_temp$age, group = merged_temp_temp$pop, location = "median")[["p.value"]]
        if (is.na(p_valor_levene)) {
            p_vals[i] <- 'None'
            next
        } else if (p_valor_levene < 0.05) {
            p_vals[i] <- 'Non-homogeneous'
            next
        }
        #Si es fa l'anàlisi, es guarda el p_valor al vector
        p_vals[i] <- round(kruskal.test(age ~ pop, data = merged_temp_temp)[['p.value']], digits = 5)
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
    pops <- unique(merged %>% arrange(metapop) %>% pull(pop))
    pop_color <- pop_color[pops]
    #Faig el gràfic. Convertir la x en factor després d'ordenar per pop i per metapop serveix per tenir les poblacions ordenades per metapoblació
    #i per la representació de geom_signif, perquè tenim la primer i última població de cada metapop alfabèticament
    options(repr.plot.width = 20, repr.plot.height = 10, warn = 1)
    grafic <- merged %>% arrange(pop) %>% arrange(metapop) %>%
        ggplot(mapping = aes(x=factor(pop, levels=unique(pop)), y=age/1000, fill=pop, color = grup, alpha = grup, by = metapop)) +
        geom_boxplot() + scale_color_manual(values = c('gray','#6b0000'), name= 'iSAFE') + scale_alpha_manual(values = c(0.35, 1), name= 'iSAFE') +
        scale_fill_manual(values = pop_color, breaks = names(pop_color), name = 'Population') +
        labs(x='', y='Estimated variant age (thousands of generations)', title = paste0('Chr ', chr, ': ', inici, '-', final)) +
        geom_signif(comparisons = comparisons_list, annotations = annotations_vector, textsize = 7, color = 'black', alpha = 1) + theme(text = element_text(size = 20))
    return(grafic)
}

##Un exemple d'execució d'aquesta funció és el següent:
diferencia_edats(2, 72494272, 72676317, AMR = FALSE) #EXOC6B
