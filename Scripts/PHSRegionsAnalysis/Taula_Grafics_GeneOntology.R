####################################################################################
####################################################################################
########################### Creació de la taula #################################### 
####################################################################################
####################################################################################

#Anàlisis de les regions de PHS per estimar l'edat dels esdeveniments de selecció en les diferents metapoblacions
#Es compara la distribució de les variants significatives per iSAFE de cada regió entre les poblacions de cada metapoblació
#Per poder estimar l'edat de la selecció, les distribucions:
    # - Han de ser homogènies (test de Levene). Si no, no es pot analitzar la regió
    # - No han de mostrar diferències significatives en un test de Kruskal-Wallis (ANOVA no paramètric)
#Les estimes de l'edat són de Relate, però el codi es pot adaptar per GEVA

#S'ha de tenir la taula de PopHumanScan, s'ha d'haver executat globalFunctions.R, s'ha de tenir accés a les taules d'iSAFE del servidor
#i s'ha de tenir les taules de Relate (els paths son al meu usuari del servidor)
#Important descartar les regions de PHS del cromosoma X (no hi ha dades en les altres taules)

#Aquesta funció et retorna la taula desitjada on dona les edats d'una variant significativa en una metapoblació, FILTRANT LES 15 VARIANTS AMB MAJOR ISAFE A CADA POBLACIÓ
#si no són significativament diferents les unes respecte a les altres en les diferents poblacions
#Es consideren significatives puntuacions d'iSAFE >0.1 i edats de variants amb qualitat >0.5
#No hi ha cap nombre mínim de variants significatives d'iSAFE per filtrar

variantage_table_kw <- function(pophumanscan = phs, AFR = AFRpops, EAS = EASpops, EUR = EURpops, SAS = SASpops) {
    #Valors mínims dels estadístics per filtrar (per canviar si fa falta)
    min_iSAFE <- 0.1
    #La taula haurà de tenir tantes files com regions a pophumanscan
    l <- nrow(pophumanscan)
    #Creo la taula buida on guardar tota la informació que busquem però amb totes les files ja i tot amb NA
    #A GeneID, ja puc guardar directament la columna sencera dels gens, igual que ID i metapop_PHS, que diu les poblacions on s'ha detectat
    taula <- data.table(Region = as.character(rep(NA,l)), ID = pophumanscan$sweepGroup, GeneID = pophumanscan$GeneID, metapop_PHS = pophumanscan$rawMetapopulation,
                        AFR = as.numeric(rep(NA,l)), AFR_n = as.numeric(rep(NA,l)), AFR_test = as.character(rep(NA,l)), AFR_aoh = as.numeric(rep(NA,l)),
                        EUR = as.numeric(rep(NA,l)), EUR_n = as.numeric(rep(NA,l)), EUR_test = as.character(rep(NA,l)), EUR_aoh = as.numeric(rep(NA,l)),
                        EAS = as.numeric(rep(NA,l)), EAS_n = as.numeric(rep(NA,l)), EAS_test = as.character(rep(NA,l)), EAS_aoh = as.numeric(rep(NA,l)),
                        SAS = as.numeric(rep(NA,l)), SAS_n = as.numeric(rep(NA,l)), SAS_test = as.character(rep(NA,l)), SAS_aoh = as.numeric(rep(NA,l)),
                        Cause = as.character(rep(NA,l)),
                        YRI_isafe_max = as.numeric(rep(NA, l)), YRI_age_median = as.numeric(rep(NA,l)), YRI_age_sd = as.numeric(rep(NA,l)), YRI_num = as.numeric(rep(NA,l)),
                        LWK_isafe_max = as.numeric(rep(NA, l)), LWK_age_median = as.numeric(rep(NA,l)), LWK_age_sd = as.numeric(rep(NA,l)), LWK_num = as.numeric(rep(NA,l)),
                        GWD_isafe_max = as.numeric(rep(NA, l)), GWD_age_median = as.numeric(rep(NA,l)), GWD_age_sd = as.numeric(rep(NA,l)), GWD_num = as.numeric(rep(NA,l)),
                        MSL_isafe_max = as.numeric(rep(NA, l)), MSL_age_median = as.numeric(rep(NA,l)), MSL_age_sd = as.numeric(rep(NA,l)), MSL_num = as.numeric(rep(NA,l)),
                        ESN_isafe_max = as.numeric(rep(NA, l)), ESN_age_median = as.numeric(rep(NA,l)), ESN_age_sd = as.numeric(rep(NA,l)), ESN_num = as.numeric(rep(NA,l)),
                        ACB_isafe_max = as.numeric(rep(NA, l)), ACB_age_median = as.numeric(rep(NA,l)), ACB_age_sd = as.numeric(rep(NA,l)), ACB_num = as.numeric(rep(NA,l)), 
                        ASW_isafe_max = as.numeric(rep(NA, l)), ASW_age_median = as.numeric(rep(NA,l)), ASW_age_sd = as.numeric(rep(NA,l)), ASW_num = as.numeric(rep(NA,l)), 
                        CEU_isafe_max = as.numeric(rep(NA, l)), CEU_age_median = as.numeric(rep(NA,l)), CEU_age_sd = as.numeric(rep(NA,l)), CEU_num = as.numeric(rep(NA,l)), 
                        TSI_isafe_max = as.numeric(rep(NA, l)), TSI_age_median = as.numeric(rep(NA,l)), TSI_age_sd = as.numeric(rep(NA,l)), TSI_num = as.numeric(rep(NA,l)),
                        FIN_isafe_max = as.numeric(rep(NA, l)), FIN_age_median = as.numeric(rep(NA,l)), FIN_age_sd = as.numeric(rep(NA,l)), FIN_num = as.numeric(rep(NA,l)),
                        GBR_isafe_max = as.numeric(rep(NA, l)), GBR_age_median = as.numeric(rep(NA,l)), GBR_age_sd = as.numeric(rep(NA,l)), GBR_num = as.numeric(rep(NA,l)),
                        IBS_isafe_max = as.numeric(rep(NA, l)), IBS_age_median = as.numeric(rep(NA,l)), IBS_age_sd = as.numeric(rep(NA,l)), IBS_num = as.numeric(rep(NA,l)),
                        CHB_isafe_max = as.numeric(rep(NA, l)), CHB_age_median = as.numeric(rep(NA,l)), CHB_age_sd = as.numeric(rep(NA,l)), CHB_num = as.numeric(rep(NA,l)),
                        JPT_isafe_max = as.numeric(rep(NA, l)), JPT_age_median = as.numeric(rep(NA,l)), JPT_age_sd = as.numeric(rep(NA,l)), JPT_num = as.numeric(rep(NA,l)),
                        CHS_isafe_max = as.numeric(rep(NA, l)), CHS_age_median = as.numeric(rep(NA,l)), CHS_age_sd = as.numeric(rep(NA,l)), CHS_num = as.numeric(rep(NA,l)),
                        CDX_isafe_max = as.numeric(rep(NA, l)), CDX_age_median = as.numeric(rep(NA,l)), CDX_age_sd = as.numeric(rep(NA,l)), CDX_num = as.numeric(rep(NA,l)), 
                        KHV_isafe_max = as.numeric(rep(NA, l)), KHV_age_median = as.numeric(rep(NA,l)), KHV_age_sd = as.numeric(rep(NA,l)), KHV_num = as.numeric(rep(NA,l)),
                        GIH_isafe_max = as.numeric(rep(NA, l)), GIH_age_median = as.numeric(rep(NA,l)), GIH_age_sd = as.numeric(rep(NA,l)), GIH_num = as.numeric(rep(NA,l)), 
                        PJL_isafe_max = as.numeric(rep(NA, l)), PJL_age_median = as.numeric(rep(NA,l)), PJL_age_sd = as.numeric(rep(NA,l)), PJL_num = as.numeric(rep(NA,l)),
                        BEB_isafe_max = as.numeric(rep(NA, l)), BEB_age_median = as.numeric(rep(NA,l)), BEB_age_sd = as.numeric(rep(NA,l)), BEB_num = as.numeric(rep(NA,l)),
                        STU_isafe_max = as.numeric(rep(NA, l)), STU_age_median = as.numeric(rep(NA,l)), STU_age_sd = as.numeric(rep(NA,l)), STU_num = as.numeric(rep(NA,l)), 
                        ITU_isafe_max = as.numeric(rep(NA, l)), ITU_age_median = as.numeric(rep(NA,l)), ITU_age_sd = as.numeric(rep(NA,l)), ITU_num = as.numeric(rep(NA,l)))
    #Ara faig aquesta iteració per a cada regió de pophumanscan
    last_chr <- 0 #Faig això per comprovar si el chr és el mateix que l'anterior; per només obrir un cop cada taula de Relate (va per chr)
    for (n in 1:l) {
        pos_inici <- pophumanscan$start[n]
        pos_final <- pophumanscan$end[n]
        #Així s'obté el número de cromosoma de pophumanscan, perquè està com a chr2
        chr <- gsub('chr', '', pophumanscan$chr[n])
        #Afegeixo la regió a la taula
        region <- paste0(chr, ':', pos_inici, '-', pos_final)
        taula[n,"Region"] <- region
        #Agafo la regió d'iSAFE de la taula i faig pivot longer i selecciono la informació que no sigui NA 
        isafe <- sqlToDf(chr, pos_inici, pos_final, 'isafe') %>% pivot_longer(cols = c('ACB':'YRI'), names_to = 'pop', values_to='isafe') %>%
            filter(!is.na(isafe) & (pop %in% AFR | pop %in% EUR | pop %in% SAS | pop %in% EAS))
        #Obro la taula de les edats, però si el chr és el mateix que l'anterior, ja està oberta
        if (chr != last_chr) {
            print(chr)
            path <- paste0('/home/anoguera/Data/relate_ages/allele_ages_chr', chr, '.csv')
            allele_ages <- fread(path) %>% select(BP, pop, lower_age, upper_age)
            last_chr <- chr
        }
        #Ajunto la regió corresponent de la taula allele_ages amb la taula d'isafe
        merged <- merge(isafe, allele_ages %>% filter(BP >= pos_inici & BP <= pos_final), by.x = c('physicalPos', 'pop'), by.y = c('BP', 'pop'))
        #Si no hi ha cap variant d'aquella regió a isafe o a relate, ja es pot acabar la iteració per la regió, i l'afegeixo a la taula
        if (nrow(merged) == 0) {
            taula[n,"Cause"] <- 'No variants'
            next
        }
        #Afegeixo rsid_pop i el grup segons si és significant o no el valor i la mitjana de l'edat superior i inferior de relate
        merged <- merged %>% mutate(rsid_pop = paste(rsid, pop, sep=':'),
                                    grup = case_when(isafe >= min_iSAFE ~ 1, TRUE ~ 0), 
                                    age = 0.5*(lower_age + upper_age))
        #Si alguna regió no té cap valor significatiu d'iSAFE, ho poso a la taula i passo a la següent regió
        if (merged %>% summarise(quant = sum(grup == 1)) == 0) {
            taula[n,"Cause"] <- 'N.S.'
            next
        }
        #Si hi ha valors, primer de tot filtro per quedar-me només amb els 15 valors d'iSAFE superiors en cada població
        isafe15 <- merged %>% filter(grup == 1) %>% group_by(pop) %>% slice_max(isafe, n=15)
        #Ara de la taula principal, elimino les files que no estiguin a isafe15
        merged <- merged %>% filter(grup == 0 | rsid_pop %in% isafe15$rsid_pop)
        #Si hi ha valors, calculo les coses desitjades per a cada població 
        taula_temp <- merged %>% filter(grup == 1) %>% group_by(pop) %>% summarise(isafe_max = max(isafe), age_median = median(age),
                                                                                   age_sd = sd(age), num = n(), .groups = 'keep')
        #Em quedo només amb la informació de merged de poblacions que tinguin variants significatives
        merged <- merged %>% filter(pop %in% unique(taula_temp$pop))
        #No hi ha mínim de variants amb iSAFe significatiu, així que algunes regions/poblacions no podran calcular l'sd (1 sola variant)
        #Ara a taula_temp tinc 1 fila/població i 4 columnes amb la seva informació, i ho afegeixo a la taula principal
        for (fila in 1:nrow(taula_temp)) {
            pob_actual <- taula_temp[fila,'pop']
            #Cada fila té aquests 4 valors a afegir:
            for (columna in c('isafe_max', 'age_median', 'age_sd', 'num')) {
                #Accedeixo a la posició de la taula que cal modificar
                taula[n,paste(pob_actual, columna, sep='_')] <- taula_temp[fila, columna]
            }
        }
        #Ara falta calcular si hi ha diferències significatives entre les poblacions d'una metapoblació. Itero sobre les 4 metapoblacions
        for (metapop in c('AFR', 'EUR', 'EAS', 'SAS')) {
            #Filtro una taula on només hi hagi les variants de la regió per les poblacions de la metapoblació
            #Per accedir al vector AFR, per exemple, faig get('AFR'), amb la metapop de la iteració
            merged_temp <- merged %>% filter(pop %in% get(metapop))
            #Si aquella metapoblació no té cap variant significativa, deixo el NA que hi havia a l'edat i poso que la seva n és 0
            if (nrow(merged_temp) == 0) {
                taula[n,paste0(metapop, '_n')] <- 0
                next
            }
            #Si hi ha alguna variant significativa, calculo la mitjana
            #Faig la mitjana fent la mitjana de les medianes ja calculades per població a taula_temp
            taula[n,metapop] <- round(mean(as.numeric(taula_temp %>% filter(pop %in% get(metapop)) %>% pull(age_median))), digits = 0)
            #Ara s'ha de fer un test de Kruskal-Wallis per les poblacions de la metapoblació. El que passa és que si hi ha 1 sola població, no té sentit
            #Per tant, miro quantes poblacions hi ha diferents
            cont <- n_distinct(merged_temp$pop)
            taula[n,paste0(metapop, '_n')] <- cont
            #També faré aquí un test d'homogeneitat per mirar si hi ha diferències significatives entre la variància de variannts amb iSAFE significatiu o no significatiu
            #Ho faré per totes les variants significatives i no significatives de la regió alhora perquè en principi no canvia molt. Poso el p-valor a la taula però no faig res en funció d'ell
            if (nrow(merged_temp %>% filter(grup == 0)) == 0) { #Però miro que hi hagi variants no significatives per evitar un eror que dona una regió en concret sense variants no significatives
                taula[n,paste0(metapop, '_aoh')] <- NA
            } else {
                taula[n,paste0(metapop, '_aoh')] <- levene.test(y = merged_temp$age, group = merged_temp$grup, location = "median")[["p.value"]]
            }
            merged_temp <- merged_temp %>% filter(grup == 1) #Ja no fa falta tenir les variants no significants per res
            #Si només hi ha una població, l'edat és la mitjana calculada, però s'ha de considerar que només hi ha dades per una població i no s'ha fet cap test estadístic
            #També hi ha un problema si només hi ha una sola variant significativa per població; en aquest cas, el nº files de merged_temp_temp = cont
            #Per tant, en aquest cas tampoc es pot fer el test
            if (cont == 1 | nrow(merged_temp) == cont) {
                taula[n,paste0(metapop, '_test')] <- 'None'
                next
            }
            #Per poder fer el test de Kruskal-Wallis, fa falta fer un test d'homogeneitat i que no hi hagin diferències significatives entre les variàncies de les diferents poblacions
            #Faig un test de Levene amb la mediana perquè és no paramètric i no passa res perquè no segueixi una distribució normal i accedeixo al seu p-valor
            #En alguns casos, si hi ha massa pocs valors, el p_valor de Levene és NA, així que miro aquests casos i dic que no es pot fer cap test.
            p_valor_levene <- levene.test(y = merged_temp$age, group = merged_temp$pop, location = "median")[["p.value"]]
            if (is.na(p_valor_levene)) {
                taula[n,paste0(metapop, '_test')] <- 'None'
            } else if (p_valor_levene < 0.05) {
                taula[n,paste0(metapop, '_test')] <- 'Non-homogeneous'
                next
            }
            #Finalment, si hi ha 2 o més poblacions diferents, es fa el test de Kruskal-Wallis per mirar si hi ha diferències significatives en l'edat
            #Em quedo amb el p-valor, que es troba amb [['p.value']] del test de Kruskal-Wallis
            p_valor <- kruskal.test(age ~ pop, data = merged_temp)[['p.value']]
            #Poso aquest p-valor en la columna del test
            taula[n,paste0(metapop, '_test')] <- as.character(round(as.numeric(p_valor), digits = 5))
        }
    }
    return(taula)
}


############################################################
############### Filtre dels resultats ######################
############################################################

#A aquesta funció se li dona el resultat de variantage_table i filtra aquelles regions que són d'interès
#(en el sentit que es detecti un sweep selectiu al menys a una metapoblació que concordi entre les subpoblacions)
#Es pot triar quin és el p_valor a partir del qual es consideren diferències significatives al cridar la funció
variantage_table_filter <- function(variantage_table_result, sign_level = 0.05) {
    #Filtro les files on cap de les metapoblacions tingui un resultat interessant:
    #Que s'hagi fet un test ANOVA i no hagi sortit significatiu per al menys una de les metapoblacions
    variantage_table_result <- variantage_table_result %>% filter( is.na(Cause) &
                                                                  (as.numeric(AFR_test) > sign_level | as.numeric(EUR_test) > sign_level |
                                                                   as.numeric(EAS_test) > sign_level | as.numeric(SAS_test) > sign_level))
    #Aquests as.numeric() de sobre causen warnings pero justament m'interessa que si son strings es converteixin en NA per filtrar
    return(variantage_table_result)
}

resultats_filtrats <- variantage_table_filter(resultats)


############################################################
########## Representació de les edats detectades ###########
############################################################

###HISTOGRAMA per representar les mitjanes estimades de les edats dels sweeps selectius per cada metapoblació
histogram_SweepAge_metapop <- function(resultats_filtrats) {
    #Faig que cada edat en una metapoblació passi a ser un registre, però llavors farà falta filtrar aquelles edats que siguin NA o bé que no s'hagi fet el test estadístic a la seva metapoblació
    resultats_filtrats <- resultats_filtrats %>% pivot_longer(cols = c('AFR', 'EUR', 'EAS', 'SAS'), names_to = 'metapop_kruskal', values_to = 'sweepAge') %>% 
        #He de mirar si a la població que toca, s'ha fet un test estadístic i no és significatiu
        filter(!is.na(sweepAge) & ((metapop_kruskal == 'AFR' & AFR_test != 'None' & AFR_test != 'Non-homogeneous' & AFR_test >= 0.05) |
                                   (metapop_kruskal == 'EUR' & EUR_test != 'None' & EUR_test != 'Non-homogeneous' & EUR_test >= 0.05) |
                                   (metapop_kruskal == 'EAS' & EAS_test != 'None' & EAS_test != 'Non-homogeneous' & EAS_test >= 0.05) |
                                   (metapop_kruskal == 'SAS' & SAS_test != 'None' & SAS_test != 'Non-homogeneous' & SAS_test >= 0.05)))
    options(repr.plot.width = 20, repr.plot.height = 15, warn = 1)
    histogram <- resultats_filtrats %>% ggplot(mapping = aes(x = sweepAge/1000, fill = metapop_kruskal)) + geom_histogram(binwidth = 1, color = 'black') +
        labs(x= 'Sweep Age (thousands of generations)', y = '# of regions', title = 'Distribution of sweep ages among metapopulations') +
        scale_fill_manual(values = c("#F7F14A", "#33B033", "#5691C4", "#A965BA"), name = 'Metapopulation') + theme(text = element_text(size = 20)) +
        facet_wrap(~metapop_kruskal, ncol = 1 )
    return(histogram)
}

histogram_SweepAge_metapop(resultats_filtrats)


###HISTOGRAMA per representar les mitjanes estimades de les edats dels sweeps selectius per cada població dins de les metapoblacions
histogram_SweepAge_pop <- function(resultats_filtrats, pop_color = popPal[c(1:7, 12:26)],
                                   AFR = AFRpops, EAS = EASpops, EUR = EURpops, SAS = SASpops) {
    #Selecciono les columnes de les edats en poblacions, faig que cada edat passi a ser un registre  i filtro les NA
    resultats_filtrats <- resultats_filtrats %>% select(Region, ID, YRI_age_median,LWK_age_median,GWD_age_median,MSL_age_median,ESN_age_median,ACB_age_median,ASW_age_median,
                                                        CEU_age_median,TSI_age_median,FIN_age_median,GBR_age_median,IBS_age_median,CHB_age_median,JPT_age_median,CHS_age_median,
                                                        CDX_age_median,KHV_age_median,GIH_age_median,PJL_age_median,BEB_age_median,STU_age_median,ITU_age_median) %>%
        #Trec el suffix _age_median de les columnes que ho tenen perquedar-me amb la població i faig que cada regio/població sigui un registre i filtro les NA
        rename_with(~str_remove(., '_age_median')) %>% pivot_longer(cols = 'YRI':'ITU', names_to = 'pop', values_to = 'age_median') %>% filter(!is.na(age_median)) %>%
        #afegeixo la columna metapop amb les metapoblacions
        mutate(metapop = case_when(pop %in% AFR ~ 'AFR', pop %in% EUR ~ 'EUR', pop %in% EAS ~ 'EAS', pop %in% SAS ~ 'SAS'))
    #Ara represento l'histograma
    options(repr.plot.width = 20, repr.plot.height = 15, warn = 1)
    histogram <- resultats_filtrats %>% ggplot(mapping = aes(x = age_median/1000, fill = pop)) + geom_histogram(binwidth = 5, color = 'black', position = position_dodge2(preserve = "single")) +
        labs(x= 'Sweep Age (thousands of generations)', y = '# of regions', title = 'Distribution of sweep ages among populations') + 
        scale_fill_manual(values = pop_color, name = 'Population') + theme(text = element_text(size = 20)) +
        facet_wrap(~metapop, ncol = 1 ) + scale_x_continuous(breaks = c(0:33)*10)
    return(histogram)
}

histogram_SweepAge_pop(resultats_filtrats) + coord_cartesian(xlim = c(0, 100))



######################################################################
###### Anàlisi de la variància entre variants significatives/no ######
######################################################################

#Representació del nombre de regions on s'ha detectat un sweep selectiu que tenen diferències en la variància entre variants significatives i no per iSAFE
#A mode informatiu, per comprovar que hi ha una diferència entre la variància de l'edat de les variants que tenen iSAFE significatiu i les que no
#Els resultats no són molt posotoius: ~50% de les regions mostren diferències significatives només

aohv_plot <- function(resultats_filtrats) {
    resultats_filtrats <- resultats_filtrats %>% select(Region, AFR_aoh, EUR_aoh, EAS_aoh, SAS_aoh, AFR_test, EUR_test, EAS_test, SAS_test) %>%
        pivot_longer(cols = 'AFR_aoh':'SAS_aoh', names_to = 'metapop_aoh', values_to = 'p_valor', values_drop_na = TRUE) %>%
        filter((metapop_aoh == 'AFR_aoh' & AFR_test != 'None' & AFR_test != 'Non-homogeneous' & AFR_test > 0.05) |
               (metapop_aoh == 'EUR_aoh' & EUR_test != 'None' & EUR_test != 'Non-homogeneous' & EUR_test > 0.05) |
               (metapop_aoh == 'EAS_aoh' & EAS_test != 'None' & EAS_test != 'Non-homogeneous' & EAS_test > 0.05) |
               (metapop_aoh == 'SAS_aoh' & SAS_test != 'None' & SAS_test != 'Non-homogeneous' & SAS_test > 0.05)) %>%
        mutate(significant = case_when(p_valor < 0.05 ~ '< 0.05', TRUE ~ '> 0.05'))
    options(repr.plot.width = 20, repr.plot.height = 7, warn = 1)
    grafic <- resultats_filtrats %>% ggplot(mapping = aes(x = significant, fill = metapop_aoh, color = significant)) + geom_bar() +
        labs(x= 'p-value of Levene test', y = '# of regions', title = 'Significance of Levene test for homogeneity of variance among metapopulations in filtered regions') + 
        scale_fill_manual(values = c('AFR_aoh' = "#F7F14A", 'EAS_aoh' = "#33B033", 'EUR_aoh' = "#5691C4", 'SAS_aoh' = "#A965BA"), name = 'Metapopulation', guide = 'none') +
        theme(text = element_text(size = 20)) + scale_color_manual(values = c('red', 'black')) +
        facet_grid(cols = vars(metapop_aoh), labeller = labeller(metapop_aoh = c('AFR_aoh' = 'AFR', 'EUR_aoh' = 'EUR', 'EAS_aoh' = 'EAS', 'SAS_aoh' = 'SAS')))
    return(grafic)
}
aohv_plot(resultats_filtrats)



############################################################
############## Resultats de PHS filtrats####################
############################################################

#Per quedar-nos amb les regions de PopHumanScan que hem filtrat en l'anàlisi i comparar en quines poblacions es mostra senyal, etc. 

variantage_table_phs <- function(resultats_filtrats, PHS = phs, sign_level = 0.05) {
    #Em quedo només amb les regions detectades en anàlisis anteriors
    #Afegeixo la columna per les poblacions detectades
    #Defineixo el nivell se significància per l'ANOVA, a partir del qual es considera que les edats estimades entre les diferents poblacions no són significativament diferents
    PHS <- PHS %>% select(c("GeneID", "sweepGroup", "chr", "start", "end", "filterTime", "rawMetapopulation")) %>% filter(sweepGroup %in% resultats_filtrats$ID) %>%
        mutate(filtered_metapop = as.character(NA))
    #Itero sobre cada fila dels resultats
    for (n in 1:nrow(resultats_filtrats)) {
        region <- unlist(strsplit(resultats_filtrats$Region[n], ':'))
        chr <- as.numeric(region[1])
        start <- as.numeric(unlist(strsplit(region[2], '-'))[1])
        end <- as.numeric(unlist(strsplit(region[2], '-'))[2])
        #A pops, quines poblacions mostren un sweep selectiu en moments no significativament diferents per cada regió
        pops <- c()
        for (pop in c('AFR', 'EUR', 'EAS', 'SAS')) {
            columna <- paste0(pop, '_test')
            #Mirem les columnes dels resultats on es guarda un p_valor o bé un None o un NA
            valor <- resultats_filtrats[n, eval(as.symbol(columna))]
            #Si és un p-valor que és > 0.05, afegim aquesta població a pops
            if (is.na(valor) | valor == 'None' | valor == 'Non-homogeneous') {
                next
            } else if (as.numeric(valor) >= sign_level) {
                #Afegeixo la metapoblació al vector
                pops <- append(pops, pop)
            }
        }
        #Faig un str per afegir a la taula i afegeixo les metapoblacions corresponents
        PHS[n, "filtered_metapop"] <- paste(pops, collapse = ',')
    }
    return(PHS)
}
resultats_filtrats_phs <- variantage_table_phs(resultats_filtrats)


##########################################################################
######################### Diagrames de Venn ##############################
##########################################################################

library(VennDiagram)
#Funció per dibuixar el diagrama venn  de signatures a partirt de phs
#Es dibuixa el diagrama de Venn del tipus de signatura (LD, SFS o Protein Changes) però per a les regions que es poden analitzar
venn_phs_signature <- function(PHS) {
    #Guardaré en un vector cada SweepGroup que mostri una determinada signatura, buscant una substring de filterTime amb grepl
    #LD correspon a RecentSweep, SFS correspon a MidSweep i OldSweep i MKT és Ancient
    LD <- PHS %>% filter(grepl('RecentSweep', filterTime, fixed = TRUE)) %>% pull(sweepGroup)
    SFS <- PHS %>% filter(grepl('OldSweep', filterTime, fixed = TRUE) | grepl('MidSweep', filterTime, fixed = TRUE)) %>% pull(sweepGroup)
    MKT <- PHS %>% filter(grepl('Ancient', filterTime, fixed = TRUE)) %>% pull(sweepGroup)
    #Així obtenim la llista de les regions per representar
    llista <- list(LD = LD, SFS = SFS, MKT = MKT)
    #Obtinc una taula amb les mides de les àrees a representar
    taula <- get.venn.partitions(llista, keep.elements = FALSE)
    #Construeixo el diagrama
    options(repr.plot.width = 10, repr.plot.height = 10, warn = 1)
    diagram <- draw.triple.venn(area1 = sum(taula %>% filter(LD == TRUE) %>% pull(..count..)), 
                                area2 = sum(taula %>% filter(SFS == TRUE) %>% pull(..count..)),
                                area3 = sum(taula %>% filter(MKT == TRUE) %>% pull(..count..)),
                                n12 = sum(taula %>% filter(LD == TRUE & SFS == TRUE) %>% pull(..count..)),
                                n23 = sum(taula %>% filter(MKT == TRUE & SFS == TRUE) %>% pull(..count..)),
                                n13 = sum(taula %>% filter(LD == TRUE & MKT == TRUE) %>% pull(..count..)),
                                n123 = sum(taula %>% filter(LD == TRUE & SFS == TRUE & MKT == TRUE) %>% pull(..count..)),
                                category = c('LD', 'SFS', 'Protein Changes'), euler.d = FALSE, scaled = FALSE,
                                fill = c("#33B033", "#5691C4", "#F7F14A"), lty = 'blank', alpha = 0.35, 
                                cex = 3, cat.cex = 3.5, cat.col = 'grey1', fontfamily = "sans", cat.fontfamily = "sans")
    return(diagram)
}

venn_phs_signature(resultats_filtrats_phs)


#Aquest següent representa el diagrama de Venn de les poblacions on es detecta selecció en cada regió de PHS (segons PHS, no segons l'anàlisi de la taula!)
#(en les regions que es poden analitzar)
library(VennDiagram)
#Funció per dibuixar el diagrama venn de poblacions a partirt de phs
venn_phs_pops <- function(PHS, printing = 'raw') { #Posar c('raw', 'percent') a printing si volem obtenir també els percentatges
    #Guardaré en un vector cada SweepGroup que mostri una senyal per a cada metapoblació, buscant una substring de filterTime amb grepl
    llista <- list(AFR = c(), EUR = c(), EAS = c(), SAS = c())
    for (metapop in c('AFR', 'EUR', 'EAS', 'SAS')) {
        llista[[metapop]] <- PHS %>% filter(grepl(metapop, rawMetapopulation, fixed = TRUE)) %>% pull(sweepGroup)
    }
    #Faig la llista de les interseccions
    taula <- get.venn.partitions(llista, keep.elements = FALSE)
    #Ara que tinc la llista de les regions detectades per cada població, construeixo el diagrama
    options(repr.plot.width = 15, repr.plot.height = 15, warn = 1)
    draw.quad.venn(area1 = sum(taula %>% filter(AFR == TRUE) %>% pull(..count..)),
                   area2 = sum(taula %>% filter(EUR == TRUE) %>% pull(..count..)),
                   area3 = sum(taula %>% filter(EAS == TRUE) %>% pull(..count..)),
                   area4 = sum(taula %>% filter(SAS == TRUE) %>% pull(..count..)),
                   n12 = sum(taula %>% filter(AFR == TRUE & EUR == TRUE) %>% pull(..count..)),
                   n13 = sum(taula %>% filter(AFR == TRUE & EAS == TRUE) %>% pull(..count..)),
                   n14 = sum(taula %>% filter(AFR == TRUE & SAS == TRUE) %>% pull(..count..)),
                   n23 = sum(taula %>% filter(EUR == TRUE & EAS == TRUE) %>% pull(..count..)),
                   n24 = sum(taula %>% filter(EUR == TRUE & SAS == TRUE) %>% pull(..count..)),
                   n34 = sum(taula %>% filter(EAS == TRUE & SAS == TRUE) %>% pull(..count..)),
                   n123 = sum(taula %>% filter(AFR == TRUE & EUR == TRUE & EAS == TRUE) %>% pull(..count..)),
                   n124 = sum(taula %>% filter(AFR == TRUE & EUR == TRUE & SAS == TRUE) %>% pull(..count..)),
                   n134 = sum(taula %>% filter(AFR == TRUE & EAS == TRUE & SAS == TRUE) %>% pull(..count..)),
                   n234 = sum(taula %>% filter(EUR == TRUE & EAS == TRUE & SAS == TRUE) %>% pull(..count..)),
                   n1234 = sum(taula %>% filter(AFR == TRUE & EUR == TRUE & EAS == TRUE & SAS == TRUE) %>% pull(..count..)),
                   category = c('AFR', 'EUR', 'EAS', 'SAS'), print.mode = printing,
                   fill = c("#F7F14A","#5691C4", "#33B033", "#A965BA"), lty = 'blank', alpha = 0.4, 
                   cex = 3, cat.cex = 3.5, cat.col = 'grey1', fontfamily = "sans", cat.fontfamily = "sans")
}
#Diagrama de Venn per les regions seleccionades amb l'anàlisi fet segons la/les metapoblació que s'ha detectat amb popHumanScan
venn_phs_pops(resultats_filtrats_phs)


#########################################################################
############# Comparació de les metapoblacions detectades ###############
#########################################################################

#Per les regions seleccionades, representació de les metapoblacions on PHS detecta selecció vs. les metapoblacions on hi ha un sweep selectiu sense diferències en una metapoblació
options(repr.plot.width = 20, repr.plot.height = 10, warn = 1)
taula_phs <- resultats_filtrats_phs %>% summarise(AFR = nrow(resultats_filtrats_phs %>% filter(grepl('AFR', rawMetapopulation, fixed = TRUE))),
                                              EUR = nrow(resultats_filtrats_phs %>% filter(grepl('EUR', rawMetapopulation, fixed = TRUE))),
                                              EAS = nrow(resultats_filtrats_phs %>% filter(grepl('EAS', rawMetapopulation, fixed = TRUE))),
                                              SAS = nrow(resultats_filtrats_phs %>% filter(grepl('SAS', rawMetapopulation, fixed = TRUE)))) %>% 
    pivot_longer(cols = 1:4, names_to = 'metapop', values_to = 'count') %>% mutate(font = 'phs')
taula_filter <- resultats_filtrats_phs %>% summarise(AFR = nrow(resultats_filtrats_phs %>% filter(grepl('AFR', filtered_metapop, fixed = TRUE))),
                                                     EUR = nrow(resultats_filtrats_phs %>% filter(grepl('EUR', filtered_metapop, fixed = TRUE))),
                                                     EAS = nrow(resultats_filtrats_phs %>% filter(grepl('EAS', filtered_metapop, fixed = TRUE))),
                                                     SAS = nrow(resultats_filtrats_phs %>% filter(grepl('SAS', filtered_metapop, fixed = TRUE)))) %>%
    pivot_longer(cols = 1:4, names_to = 'metapop', values_to = 'count') %>% mutate(font = 'filter')
taula <- rbind(taula_phs, taula_filter)
taula %>% ggplot(mapping = aes(x = font, y = count, fill = metapop, stat = 'identity')) + geom_col(width = 0.7, show.legend = TRUE) +
    labs(x= 'Source', y = '# of regions', title = 'Metapopulations where a selective sweep was detected in the selected regions') +
    scale_x_discrete(limits = c('phs', 'filter'), labels = c('detected by PopHumanScan', 'detected with Kruskal-Wallis test')) +
    theme(text = element_text(size = 25)) + scale_fill_manual(values = c("#F7F14A","#33B033", "#5691C4", "#A965BA"), name = 'Metapopulation')


####################################################################################
####################################################################################
#################### Gene Ontology Enrichment Analysis #############################
####################################################################################
####################################################################################

#Fet en intervals de 5000 generacions, perquè no hi ha suficients gens en les regions analitzades per poder fer intervals de 1000 generacions o menors
#L'anàlisi de Gene Ontology és molt més complert i correcte en les regions de l'article d'iHS que en les de PHS

#Funció per obtenir un vector amb els gens corresponents
genelist <- function(genevector) {
    #Obtinc tots els noms de gens que hi ha al vector inicial separats amb comes
    genevector <- unique(unlist(strsplit(genevector, ',')))
    #Elimino els elements on hi diu noncoding, trobant els seus indexs dins del vector mitjançant grep
    indexs <- grep('nonCoding', genevector)
    if (length(indexs) > 0) {
        genevector <- genevector[- indexs]
    }
    return(genevector)
}

#Aquesta funció obté, a partir de la taula de resultats filtrats, una llista amb els gens que es troben a cada regió que mostra sweeps selectius per
#cada metapoblació per cada interval determinat en generacions
metapop_genelist <- function(resultats_filtrats) {
    resultats_filtrats <- resultats_filtrats %>% #select('Region', 'GeneID', 'AFR', 'EUR', 'EAS', 'SAS', 'AFR_test', 'EUR_test', 'EAS_test', 'SAS_test') %>%
        pivot_longer(cols = c('AFR', 'EUR', 'EAS', 'SAS'),names_to = 'metapop_kruskal', values_to = 'sweepAge') %>% 
        #He de mirar si a la població que toca, s'ha fet un test estadístic i no és significatiu
        filter(!is.na(sweepAge) & ((metapop_kruskal == 'AFR' & AFR_test != 'None' & AFR_test != 'Non-homogeneous' & AFR_test >= 0.05) |
                               (metapop_kruskal == 'EUR' & EUR_test != 'None' & EUR_test != 'Non-homogeneous' & EUR_test >= 0.05) |
                               (metapop_kruskal == 'EAS' & EAS_test != 'None' & EAS_test != 'Non-homogeneous' & EAS_test >= 0.05) |
                               (metapop_kruskal == 'SAS' & SAS_test != 'None' & SAS_test != 'Non-homogeneous' & SAS_test >= 0.05))) %>%
        #Afegeixo una columna que agrupi l'edat dels sweeps
        mutate(sweepAgeGroup = as.character(NA))
    #Creo un vector amb els diferents límits entre grups que faré
    vect = c(0:14) * 5000 #Arribo fins a 70000 en intervals de 5000
    for (n in 2:15) {
        for (fila in 1:nrow(resultats_filtrats)) {
            valor <- resultats_filtrats[fila, 'sweepAge']
            #Per cada categoria de Sweep Ages que faig, miro si el valor en qüestió hi cau dins, i si és el cas ho poso a la columna sweepAgeGroup
            if (valor > vect[n-1] & valor <= vect[n]) {
                resultats_filtrats[fila, 'sweepAgeGroup'] <- paste(vect[n-1], vect[n], sep ='-')
            }
        }
    }
    #Creo una llista buida on guardar les subllistes de gens
    llista <- list()
    for (metapop in c('AFR', 'EUR', 'EAS', 'SAS')) {
        taula_temp <- resultats_filtrats %>% filter(metapop_kruskal == metapop)
        for (grup in unique(taula_temp$sweepAgeGroup)) {
            llista[[paste(metapop, grup, sep=':')]] <- c(genelist(taula_temp %>% filter(sweepAgeGroup == grup) %>% pull(GeneID)))
        }
    }
    return(llista)
}


#Ara es pot fer l'anàlisi amb gprofiler2 
library(gprofiler2)
resultats <- gost(genes_list, organism='hsapiens', evcodes= TRUE)[['result']]

#Es possible filtrar els resultats per eliminar els causats per clusters de gens (com en l'anàlisi fet en les regions d'iHS)
#Però no vaig crear el codi perquè vam deixar aquest anàlisi abans de resoldre el problema dels clusters, que en aquest cas tampoc afecta tant els resultats com en les regions d'iHS
#S'hauria de fer un procès anàleg al de l'altre anàlisi

######################################################################
################ Representació dels resultats ########################
######################################################################

#Heat Map amb els resultats de Gene Ontology

#Un llistat dels intervals possibles per ordenar l'eix
intervals <- c(0:100) * 1000
llistat <- c()
for (n in 2:101) {
    llistat[n-1] <- paste(intervals[n-1], intervals[n], sep = '-')
}
GOresult1_metapop <- GOresult1 %>% separate(query, sep = ":", into = c("metapop", "generations")) %>% select(metapop, generations, p_value, intersection_size, term_name)
#Em quedo amb els intervals (ordenats) on s'hagi detectat algo
llistat <- llistat[llistat %in% GOresult1_metapop$generations]
#Represento el plot
options(repr.plot.width = 30, repr.plot.height = 50, warn = 1)
GOresult1_metapop %>% ggplot(aes(x = generations, y = term_name)) + geom_tile(aes(fill = log10(p_value))) + 
    geom_text(aes(label = intersection_size)) +
    scale_fill_gradient(low = "red", high = "white") +
    facet_grid(rows  = vars(metapop), scale = 'free_y', space = 'free' ) +
    theme(text = element_text(size = 23), axis.text.x = element_text(angle = 90) ) + 
    labs(x='Selective sweep (generations ago)', y='', title = 'Functional enrichment analysis results for each metapopulation', fill = 'log10(p-value)') +
    scale_x_discrete(limits = llistat)
