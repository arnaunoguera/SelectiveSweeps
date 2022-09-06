#Amb les regions de l'article d'iHS (Johnson and Voight, 2018)
#Aquestes regions s'han detectat concretament per una població, així que no es pot fer comparacions entre les poblacions d'una metapoblació
#Per tant, no es fa cap test estadístic, només s'estima l'edat del sweep selectiu de cada regió si es possible
#Es necessita la taula de les regions i les taules de Relate (les dues estan amb el path al meu usuari del servidor) i accès a les dades d'iSAFE el servidor

#################################################################################
##################### Creació de la taula #######################################
#################################################################################

#Aquesta funció et retorna la taula desitjada on es dona l'edat estimada del sweep selectiu de cada regió, FILTRANT LES 5 VARIANTS AMB MAJOR ISAFE A CADA POBLACIÓ
#si no són significativament diferents les unes respecte a les altres en les diferents poblacions
#Es consideren significatives puntuacions d'iSAFE >0.1 i edats de variants amb qualitat >0.5
#No hi ha cap nombre mínim de variants significatives d'iSAFE per filtrar
sweepage_table <- function(taula_regions = fread("/home/anoguera/Data/iHS_selected_regions.csv") %>% filter(CHR != 23), AFR = AFRpops, EAS = EASpops, EUR = EURpops, SAS = SASpops) {
    #Valors mínims dels estadístics per filtrar (per canviar si fa falta)
    min_iSAFE <- 0.1
    #Afegeixo a la taula les columnes amb la informació que m'interessa: isafe max, age median, age sd i num de variants (amb iSAFE significant), i metapop i region, que ja es pot omplir
    taula_regions <- taula_regions %>% select(POP, CHR, START_POS, STOP_POS, TAG_RS) %>% arrange(CHR) %>% #Ordeno per chr per obrir cada arxiu de chr un sol cop
        mutate(isafe_max = as.numeric(NA), age_median = as.numeric(NA), age_sd = as.numeric(NA), num_var = as.numeric(NA), 
               METAPOP = case_when(POP %in% AFR ~ 'AFR', POP %in% EUR ~ 'EUR',  POP %in% EAS ~ 'EAS', POP %in% SAS ~ 'SAS', TRUE ~ 'AMR'), 
               REGION = paste(CHR, START_POS, STOP_POS, sep = ':'))
    #Ara faig aquesta iteració per a cada regió de la taula
    last_chr <- 0 #Faig això per comprovar si el chr és el mateix que l'anterior
    allele_ages <- data.table() #Per poder borrar-ho després de cada chr
    for (n in 1:nrow(taula_regions)) {
        pos_inici <- taula_regions$START_POS[n]
        pos_final <- taula_regions$STOP_POS[n]
        chr <- taula_regions$CHR[n]
        print(n)
        #En aquest cas també ens hem de quedar amb la pop actual
        pop_actual <- taula_regions$POP[n]
        #Agafo la regió d'iSAFE de la taula i em quedo només amb la població d'interés i filtro que no sigui NA
        isafe <- sqlToDf(chr, pos_inici, pos_final, 'isafe') %>% select('physicalPos', 'rsid', pop_actual) %>% rename('isafe' = pop_actual) %>% filter(!is.na(isafe))
        if (nrow(isafe) == 0) {
            taula_regions[n, 'CAUSE'] <- 'NVi' #No hi ha variants amb informació d'iSAFE per la població en la regió
            next
        }
        #Obro la taula de les edats, però si el chr és el mateix que l'anterior, ja està oberta
        if (chr != last_chr) {
            rm(allele_ages)
            print(chr)

            path <- paste0('/home/anoguera/Data/relate_ages/allele_ages_chr', chr, '.csv')
            allele_ages <- fread(path) %>% select(BP, pop, lower_age, upper_age)
            last_chr <- chr
        }
        #Ajunto la regió corresponent de la taula allele_ages per la població actual amb la taula d'isafe
        merged <- merge(isafe, allele_ages %>% filter(pop == pop_actual & BP >= pos_inici & BP <= pos_final), by.x = c('physicalPos'), by.y = c('BP'))
        #Si no hi ha cap variant d'aquella raegió a relate que solapi amb les d'iSAFE, ja es pot acabar la iteració per la regió, i l'afegeixo a la taula
        if (nrow(merged) == 0) {
            taula_regions[n,"CAUSE"] <- 'NVr' #No variants a relate
            next
        }
        #Em quedo només amb les variants significatives d'iSAFE (i em quedo amb les 5 superiors) i calculo l'edat de les variants fent la mitjana entre l'estima superior i inferior
        merged <- merged %>% filter(isafe >= min_iSAFE) %>% slice_max(isafe, n=5) %>% mutate(age = 0.5*(lower_age + upper_age)) 
        #Si alguna regió no té cap valor significatiu d'iSAFE, ho poso a la taula i passo a la següent regió
        if (nrow(merged) == 0) {
            taula_regions[n,"CAUSE"] <- 'N.S.' #No hi ha valors significatius d'iSAFE a la regió/població
            next
        }
        #Ara guardem les dades que necessitem per la població i regió en una taula
        dades <- merged %>% summarise(isafe_max = max(isafe), age_median = median(age), age_sd = sd(age), num_var = n())
        #No hi ha mínim de variants amb iSAFe significatiu, així que algunes regions/poblacions no podran calcular l'sd (1 sola variant)
        #Ara a taula_temp tinc 4 columnes amb la seva informació, i ho afegeixo a la taula principal. Cada fila té aquests 4 valors a afegir:
        for (columna in c('isafe_max', 'age_median', 'age_sd', 'num_var')) {
            #Accedeixo a la posició de la taula que cal modificar
            taula_regions[n,columna] <- dades[1, columna]
        }
    }
    rm(allele_ages)
    return(taula_regions)
}

#Executo la funció 
resultat <- sweepage_table()


#######################################################################
#################### Filtre dels resultats ############################
#######################################################################

#A aquesta funció se li dona el resultat de sweepage_table i filtra aquelles regions que són d'interès
#(en el sentit que es pugui estimar una edat per al sweep selectiu)
#Es pot triar quin és el p_valor a partir del qual es consideren diferències significatives al cridar la funció
sweepage_table_filter <- function(variantage_table_result, sign_level = 0.05) {
    #Filtro les files on no hi hagi estima
    sweepage_table_result <- variantage_table_result %>% filter(!is.na(age_median))
    return(sweepage_table_result)
}

#Executo la funció 
resultats_filtrats <- variantage_table_filte    #Valors mínims dels estadístics per filtrar (per r(resultat)

#######################################################################
############## Representació de les desviacions típiques ##############
#######################################################################

#Un boxplot de les desviacions típiques en l'estima de l'edat per veure quina precisió hi ha
options(repr.plot.width = 23, repr.plot.height = 13, warn = 1)
resultats_filtrats %>% arrange(METAPOP) %>% replace_na(list(age_sd = 0)) %>%
    ggplot(mapping = aes(x = factor(POP, levels=unique(POP)), y = age_sd, fill = POP)) + 
    geom_boxplot()+ scale_fill_manual(values = popPal) + 
    theme(text = element_text(size = 20)) + labs(x= '', y = 'sd(age)') #+ coord_cartesian(ylim = c(0, 35000))
    #Activar coord_cartesian per fer zoom a la part interessant i ignorar els outliers


#######################################################################
############## Representació de les edats detectades ##################
#######################################################################

###HISTOGRAMA per representar les mitjanes estimades de les edats dels sweeps selectius per cada metapoblació
#S'ajunten els resultats de les poblacions d'una metapoblació
options(repr.plot.width = 20, repr.plot.height = 15, warn = 1)
resultats_filtrats %>% ggplot(mapping = aes(x = age_median/1000, fill = METAPOP)) + geom_histogram(binwidth = 1, color = 'black') +
    labs(x= 'Sweep Age (thousands of generations)', y = '# of regions', title = 'Distribution of sweep ages among metapopulations') +
    scale_fill_manual(values = c("#F7F14A", "#BA5852","#33B033", "#5691C4", "#A965BA"), name = 'Metapopulation') + theme(text = element_text(size = 20)) +
    facet_wrap(~METAPOP, ncol = 1, scales = 'free_y') #+ coord_cartesian(xlim = c(0, 60))
    #Activar coord_cartesian per fer zoom a la part important de la distribució

###HISTOGRAMA per representar les mitjanes estimades de les edats dels sweeps selectius per cada població dins de les metapoblacions
options(repr.plot.width = 20, repr.plot.height = 15, warn = 1)
resultats_filtrats %>% ggplot(mapping = aes(x = age_median/1000, fill = POP)) + geom_histogram(binwidth = 2, color = 'black', position = position_dodge2(preserve = "single")) +
    labs(x= 'Sweep Age (thousands of generations)', y = '# of regions', title = 'Distribution of sweep ages among populations') + 
    scale_fill_manual(values = popPal, name = 'Population') + theme(text = element_text(size = 20)) +
    facet_wrap(~METAPOP, ncol = 1, scales = 'free_y') + scale_x_continuous(breaks = c(0:33)*10) + coord_cartesian(xlim = c(0, 60))


#######################################################################
############### Diagrama dels resultats descartats ####################
#######################################################################

#Representació de les causes per les quals es descarten les altres regions de la taula inicial en l'anàlisi
options(repr.plot.width = 25, repr.plot.height = 10, warn = 1)
resultat %>% filter(!is.na(CAUSE)) %>% ggplot(mapping = aes(y = CAUSE, fill = CAUSE)) + geom_bar(width = 0.7, show.legend = FALSE) +
    labs(x= '# of regions', y = 'Cause', title = 'Cause of discarded regions') + 
    theme(text = element_text(size = 25)) + scale_y_discrete(labels= c("NVi" = "No variants with iSAFE\ninformation in the region", "NVr" = "No variants with Relate\ninformation in the region",
                                                                       "N.S." = "No variants with\nsignificant iSAFE value"),
                                                            limits = c("N.S.","NVi", "NVr")) +
    stat_count(geom = "text", colour = "black", size = 6, aes(label = ..count..),position=position_stack(vjust=0.5)) + scale_fill_brewer(palette="Set2")


##########################################################################################################
##########################################################################################################
################################ Gene Ontology Enrichment Analysis #######################################
##########################################################################################################
##########################################################################################################

#########################################################################
############### Obtenció de les regions de cada rang ####################
#########################################################################

#Aquesta funció obté, a partir de la taula de resultats filtrats, una taula amb les regions seleccionades en cada interval determinat 
#No fa falta assegurar que les regions no solapin entre poblacions perquè gprofiler ja agafa cada gen un sol cop
#Es fa amb intervals de 1000 fins a 700000 perquè a partir d'allà, ja no hi ha regions prou condensades en cap interval com per trobar resultats interessants
metapop_genelist <- function(resultats_filtrats) {
    vect = c(0:70) * 1000 #Arribo fins a 70000 en intervals de 1000
    vector_ages <- vector()
    for (n in 2:length(vect)) {
        for (metapop in c('AFR', 'EUR', 'EAS', 'SAS', 'AMR')) {
            vector_ages <- append(vector_ages, paste0(metapop, ':', vect[n-1], '-', vect[n]))
        } #A vector_ages hi ha tots els intervals possibles
    }
    taula <- data.table(Interval = vector_ages, Regions = as.character(NA), Genes = as.character(NA)) %>% #Creo la taula on s'afegirà tot
        separate(Interval, sep = ':', into = c('Metapop', 'Interval'))
    for (n in 1:nrow(taula)) {
        metapop <- taula$Metapop[n]
        age_min <- as.numeric(unlist(strsplit(taula$Interval[n], '-'))[1])
        age_max <- as.numeric(unlist(strsplit(taula$Interval[n], '-'))[2])
        resultats_temp <- resultats_filtrats %>% filter(METAPOP == metapop & age_median >= age_min & age_median < age_max)
        if (nrow(resultats_temp) == 0) {
            next
        }
        taula[n,'Regions'] <- paste(resultats_temp$REGION, collapse =',')
    }
    return(taula %>% filter(!is.na(Regions)) %>% separate_rows(Regions, sep=',')) #Serveix per tenir una entrada per regió detectada 
}
genes_table <- metapop_genelist(resultats_filtrats) #Creo la taula a partir dels resultats



#########################################################################
######################## GO amb gprofiler2 ##############################
#########################################################################

library(gprofiler2) #Aquests paquets només em funcionaven des de la terminal, no des de JupyterNotebook
library(biomaRt)

#He d'aconseguir la llista de gens, i no de regions, perquè gprofiler treballa amb GRCh38 i les regions són amb el 37 (hg19)
#Executar a la terminal, perque es on tinc biomaRt
ensembl <- useEnsembl(biomart = 'genes', dataset = 'hsapiens_gene_ensembl', version = 'GRCh37') #Així accedeixo a la versió de hg19
print('Connected to Ensembl')
for (n in 1:nrow(genes_table)) {
    vals <- unlist(strsplit(genes_table$Regions[n], ':'))
    chr <- vals[1]
    inici <- vals[2]
    final <- vals[3]
    gens <- getBM(attributes = 'ensembl_gene_id', filters =  c('chromosome_name', 'start', 'end'), values = list(chr, inici, final), mart = ensembl) %>% pull(ensembl_gene_id)
    genes_table[n,"Genes"] <- paste(gens, collapse=',') #Amb cada iteració, obtinc els gens d'una de les regions i es col·loquen a la columna Genes de la taula
    #Els gens d'una regió s'obtenen com a str on cada gen es separa per una coma
    print(n)
}
#Elimino les files que no tenen cap gen a la regió. Aquesta taula s'ha de guardar per poder situar després els gens en les regions
genes_table <- genes_table %>% filter(Genes != '')
#Ara ajunto tots els IDs de gens d'un mateix interval/població per poder fer el Gene Ontology
#Guardo genes_table per poder mirar després si tots els gens d'un resultat sortien de la mateixa regió
GO_table <- genes_table %>% group_by(Metapop, Interval) %>% summarise(Gene_ids = paste(Genes, collapse = ','))

#Faig l'anàlisi de Gene Ontology població per població (si no, es sobrecarrega gprofiler2)  
for (metapop in c('AFR', 'EAS', 'EUR', 'SAS', 'AMR')) {
    #Obtinc una llista on cada nom és un interval i l'element és un vector amb tots els Gene IDs d'aquell interval
    llista <- as.list(GO_table %>% filter(Metapop == metapop) %>% pull(Gene_ids) %>% strsplit(',')) 
    names(llista) <- GO_table %>% filter(Metapop == metapop) %>% pull(Interval)
    #Guardo els resultats de l'anàlisi (el path és al meu usuari)
    fwrite(x = gost(llista, organism='hsapiens', evcodes= TRUE)[['result']],
           file = paste0("/home/anoguera/Data/resultatsGO_regions_ihs/", metapop, "1000top5.csv"))
    print(metapop)
    rm(llista)
}

