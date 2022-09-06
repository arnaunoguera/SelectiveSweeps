###Gràfic de les edats estimades de variants amb iSAFE significatiu vs no significatiu, representant diversos gràfics alhora de les regions de PopHumanScan
#Idealment s'executa dins d'andromeda. S'had'haver executat globalFunctions.R abans
#S'introdueixen els inicis i finals com a dos vectors, de forma que el primer valor d'inicis es correspon al primer de finals;
#i length(inicis)=length(finals). El vector de cromosomes no ha de ser de la mateixa longitud: un cop per chr és suficient. 
#S'ha d'accedir també ales taules d'iHS i iSAFE i s'ha d'introduir la taula de PopHumanScan
#Es considera significant un valor d'iSAFE > 0.1 i es filtren les edats d'atlas amb qualitat inferior a 0.5
agemode_isafe_plot <- function(chrom = NA, inicis = NA, finals = NA, PHS, pop_color = popPal,
                               AFR = AFRpops, EAS = EASpops, EUR = EURpops, SAS = SASpops, AMR = TRUE) {
    #Miro si s'han introduit inicis/finals concrets o si serà per tots
    if (all(!is.na(chrom))) {
        PHS <- PHS %>% filter(gsub('chr','',chr) %in% chrom)
        if (all(!is.na(inicis)) & all(!is.na(finals))) {
            #Això no portarà problemes perquè totes les posicions d'inici i de final són úniques, independentment de cromosoma:
            #No s'ha d'especificar a quin cromosoma pertany una posició d'inici 
            PHS <- PHS %>% filter(start %in% inicis & end %in% finals)
        }
    }
    #Creo una data table buida
    taula <- data.table()
    #Es repeteix per cada registre de PHS que quedi després de filtrar els desitjats
    for (n in 1:nrow(PHS)) {
        pos_inici <- PHS$start[n]
        pos_final <- PHS$end[n]
        chr <- gsub('chr', '', PHS$chr[n])
        region <- paste0(chr, ':', pos_inici, '-', pos_final)
        #Creo la taula que ja son les d'iSAFE i atlas juntes i filtro per la QualScore
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
        #Filtro la taula per quedar-me només les poblacions amb iSAFEs significatius
        #Una vegada filtrada, creo una columna nova que distingeixi la regió (segons els gens que conté)
        merged_temp <- merged_temp %>% filter(pop %in% pob_significatiu) %>%
            mutate(region = paste0(PHS$GeneID[n],'\n(', region, ')'), 
                   metapop = case_when(pop %in% AFR ~ 'AFR',
                                       pop %in% EUR ~ 'EUR',
                                       pop %in% EAS ~ 'EAS',
                                       pop %in% SAS ~ 'SAS',
                                       TRUE ~ 'AMR'))
        #I en cada iteració, afegeixo la taula temporal a una taula gran
        taula <- bind_rows(taula, merged_temp)
    }
    #Faig una llista de les poblacions que es representaran per filtrar pop_colors (deixant-ho ordenat per metapoblacions)
    pops <- unique(taula %>% arrange(metapop) %>% pull(pop))
    pop_color <- pop_color[pops]
    #Faig el gràfic, separant per regió amb el facet_wrap(). El factor de la x em serveix per ordenar les poblacions com estan a la taula
    #(per metapoblació). Ordeno les poblacions tal i com estan a pop_color
    options(repr.plot.width = 25, repr.plot.height = 30, warn = 1)
    grafic <- taula %>% arrange(metapop) %>%
        ggplot(mapping = aes(x= factor(pop, levels = unique(pop)), y=AgeMode_Jnt/1000, fill=pop, color = grup)) +
        geom_boxplot() + scale_color_manual(values = c('black','#6b0000'), name= 'iSAFE') +
        scale_fill_manual(values = pop_color, breaks = names(pop_color), name = 'Population') + theme(text = element_text(size = 15)) +
        facet_wrap(~region, scales = 'free', ncol = 4) + labs(x='', y='Estimated variant age (thousands of generations)')
    #Si enlloc de facet_wrap es fa facet_grid, posant l'argument space='free', la mida de cada gràfic s'ajusta automàticament. S'ha de treure l'argument nrow
    return(grafic)
}

#Un exemple d'utilitzar la funció:
agemode_isafe_plot(chrom = 5, PHS = pophumanscan, AMR = FALSE) #Es representen totes les regions de PHS del chr 5
