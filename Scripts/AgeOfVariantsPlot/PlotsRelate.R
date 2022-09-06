####Gràfic de les edats estimades de variants amb iSAFE significatiu vs no significatiu amb edats de Relate
#Fan falta els fitxers de les edats de Relate (el path està al meu usuari del servidor) i accedir a les taules del servidor
#S'introdueixen els inicis i finals com a dos vectors, de forma que el primer valor d'inicis es correspon al primer de finals;
#i length(inicis)=length(finals). El vector de cromosomes no ha de ser de la mateixa longitud: un cop per chr és suficient. 
#Es considera significant un valor d'iSAFE > 0.1
age_isafe_plot <- function(chrom = NA, inicis = NA, finals = NA, PHS = phs, pop_color = popPal,
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
    #Ho faig per cada cromosoma per accedir un sol cop a cada taula
    for (chr in unique(chrom)) {
        chr_phs <- paste0('chr', chr)
        #Escric el path al fitxer de les edats
        path <- paste0('/home/anoguera/Data/relate_ages/allele_ages_chr', chr, '.csv')
        #Filtro la part de PHS que sigui del cromosoma pertinent
        PHS_temp <- phs %>% filter(chr == chr_phs)
        allele_ages <- fread(path) %>% select(BP, pop, lower_age, upper_age)
        #Ara repeteixo això per cada registre de PHS que hi hagi per la regió
        for (p in 1:nrow(PHS_temp)) {
            pos_inici <- PHS_temp$start[p]
            pos_final <- PHS_temp$end[p]
            region <- paste0(chr, ':', pos_inici, '-', pos_final)
            #Obtinc la taula d'iSAFE de la regió 
            isafe <- sqlToDf(chr, pos_inici, pos_final, 'isafe')
            #Si AMR == FALSE vol dir que s'ha d'eliminar la informació de les poblacions americanes
            if (AMR == FALSE) {
                isafe <- isafe[, -c("CLM", "MXL", "PEL", "PUR")]
            }
            #Faig el pivot longer, filtro els valors NA i classifico en iSAFE significant o no
            isafe <- isafe %>% pivot_longer(cols = c('ACB':'YRI'), names_to = 'pop', values_to='isafe') %>%
                filter(!is.na(isafe)) %>% mutate(grup = case_when(isafe >= 0.1 ~ 'Significant', TRUE ~ 'Nonsignificant'))
            #Consegueixo un vector amb les poblacions amb iSAFE significatiu
            pop_significatiu <- unique(isafe %>% filter(grup == 'Significant') %>% pull(pop))
            #Filtro la taula per quedar-me només les poblacions amb iSAFEs significatius
            #Una vegada filtrada, creo una columna nova que distingeixi la regió (segons els gens que conté) i la metapoblació
            isafe <- isafe %>% filter(pop %in% pop_significatiu) %>% mutate(region = paste0(PHS_temp$GeneID[p],'\n(', region, ')'), 
                                                                            metapop = case_when(pop %in% AFR ~ 'AFR',
                                                                                                pop %in% EUR ~ 'EUR',
                                                                                                pop %in% EAS ~ 'EAS',
                                                                                                pop %in% SAS ~ 'SAS',
                                                                                                TRUE ~ 'AMR'))
            #Filtro la part de relate de la regió
            relate_temp <- allele_ages %>% filter(BP >= pos_inici & BP <= pos_final)
            #Ajunto les taules de relate i isafe
            merged <- merge(isafe, relate_temp, by.x = c('physicalPos', 'pop'), by.y = c('BP', 'pop')) %>% mutate(age = 0.5*(lower_age + upper_age), rsid_pop = paste(rsid, pop, sep = ':'))
            rm(isafe)
            rm(relate_temp)
#Una vegada fet el merge, em quedo amb les 15 variants més significatives per iSAFE (no ho faig abans perquè potser perdríem les variants significatives per no estar a la taula d'edats)
#Creo una taula on hi hagi les files corresponents als 15 valors significatius d'iSAFE més alts en cada població. SI n'hi ha menys de 15, es guarden els que hi hagi. Si hi ha empats, se'n queda > 15
            isafe15 <- merged %>% select(grup, pop, isafe, rsid_pop) %>% filter(grup == 'Significant') %>% group_by(pop) %>% slice_max(isafe, n=15)
            #Ara de la taula principal, elimino les files que siguin significants  no estiguin a isafe15 
            merged <- merged %>% filter(grup == 'Nonsignificant' | rsid_pop %in% isafe15$rsid_pop)
            if (nrow(merged %>% filter(grup == 'Significant')) == 0) {
                next
            }
            pop_significatiu <- unique(merged %>% filter(grup == 'Significant') %>% pull(pop))
            merged<- merged %>% filter(pop %in% pop_significatiu)
            #Afegeixo la taula temporal a la taula gran
            taula <- bind_rows(taula, merged)
        }
        rm(allele_ages)
    }
    #Faig una llista de les poblacions que es representaran per filtrar pop_colors (deixant-ho ordenat per metapoblacions)
    pops <- unique(taula %>% arrange(metapop) %>% pull(pop))
    pop_color <- pop_color[pops]
    #Faig el gràfic, separant per regió amb el facet_wrap(). El factor de la x em serveix per ordenar les poblacions com estan a la taula
    #(per metapoblació). Ordeno les poblacions tal i com estan a pop_color
    height <- round(n_distinct(taula$region)/4 + 0.5, digits = 0) * 5 #La llargada del plot
    options(repr.plot.width = 25, repr.plot.height = height, warn = 1)
    grafic <- taula %>% arrange(metapop) %>%
        ggplot(mapping = aes(x= factor(pop, levels = unique(pop)), y=age/1000, fill=pop, color = grup)) +
        geom_boxplot() + scale_color_manual(values = c('black','#6b0000'), name= 'iSAFE') +
        scale_fill_manual(values = pop_color, breaks = names(pop_color), name = 'Population') + theme(text = element_text(size = 15)) +
        facet_wrap(~region, scales = 'free', ncol = 4) + labs(x='', y='Estimated variant age (thousands of generations)')
    #Si enlloc de facet_wrap es fa facet_grid, posant l'argument space='free', la mida de cada gràfic s'ajusta automàticament. S'ha de treure l'argument nrow
    return(grafic)
}
