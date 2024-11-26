
#haps_file: input table of haplotype frequencies 
#PT_file: input table of phenotypic values 
#regimes:  subsets haplotype data by specified selective regimes (nw, n, w, c)
#generation:  subsets haplotype data by specified generation
#term_name: specify model term of interest from lm
#output_name: output file name
#dest: output directory

#outputs table showing significance of relationship between haplotype frequency and model term of interest for all loci in haps_file

get_loci <- function(haps_file, PT_file, regimes, generation, term_name, output_name, dest) {

    allhaps = read.table(haps_file, header=TRUE)
    PT_data <- read.csv(PT_file)
   
#set condition to subset haplotype table based on input regimes

    if ("nw" %in% regimes & "n" %in% regimes) {
	    condition <- "nic == 'y'"
    } 

    if ("defended" %in% regimes & "nw_undefended" %in% regimes) {
	    condition <- "(nic == 'y' & wasp == 'y') | pop == 'R14'"
    }

    if ("nw_undefended" %in% regimes & "n" %in% regimes) {
	    condition <- "nic == 'y' & defended == 'no'"
    }

    if ("n" %in% regimes & "c" %in% regimes) {
	    condition <- "wasp == 'n' & defended == 'no'"
    } 

    if ("nw_undefended" %in% regimes & "c" %in% regimes) {
	    condition <- "(nic == 'y' & wasp == 'y' & defended == 'no') | (nic == 'n' & wasp == 'n')"
	  }

    if ("nw" %in% regimes & "c" %in% regimes) {
	    condition <- "(nic == 'n' & wasp == 'n') | (nic == 'y' & wasp == 'y')"
        } 

    if ("w" %in% regimes & "c" %in% regimes) {
      condition <- "nic == 'n'"
        } 

    if ("w" %in% regimes & "nw" %in% regimes) {
	    condition <- "wasp == 'y'"
            } 
    if ("defended" %in% regimes & "c" %in% regimes) {
	    condition <- "(defended == 'yes' | (nic == 'n' & wasp == 'n'))"
	}

    
    hap_data <- data.table(allhaps)

    #separate information in sample name into separate factors 
    hap_data <- hap_data %>% separate(pool,c("nic", "wasp", "pop","gen")) 
    hap_data$chr <- as.factor(hap_data$chr)
    hap_data$nic <- as.factor(hap_data$nic)
    hap_data$wasp <- as.factor(hap_data$wasp)
    hap_data$pop <- as.factor(hap_data$pop)
    hap_data$gen <- as.factor(hap_data$gen)
    hap_data$founder <- as.factor(hap_data$founder)

    
    #average phenotypic data by population
    population_summary <- PT_data %>% group_by(regime, population, generation, parasitism_group, nicotine_dose) %>%
      summarise(survival = mean(adults), n = n(), stdev_abs = sd(adults), 
                se_abs = stdev_abs/sqrt(n))
    
    #code toxin survival as a ratio of survival on toxin/surival on control
    survival_ratios <- population_summary %>%
      group_by(regime, population, generation, parasitism_group) %>%
      summarise(high_ratio = (survival[nicotine_dose == "1.25"]/survival[nicotine_dose == "0"]),
                low_ratio = (survival[nicotine_dose == "0.75"]/survival[nicotine_dose == "0"])) 
    
    survival_ratios$population <- as.character(survival_ratios$population)
    survival_ratios$generation <- as.character(survival_ratios$generation)
    
    #create table of survival coded as ratio for each population
    ratio <- survival_ratios %>%
      ungroup() %>%
      pivot_wider(names_from = parasitism_group, values_from = c(high_ratio, low_ratio)) %>%
      mutate(pop = paste("R", population, sep = "")) %>%
      mutate(gen = paste("G", generation, sep = "")) %>%
      unite("sample", pop:gen, sep = "") %>%
      dplyr::select(sample, regime, high_ratio_NW,, high_ratio_W, low_ratio_NW, low_ratio_W) 
    
    #import survival ratio data to haplotype frequency table
    hap_data$low_ratio_NW <- ratio$low_ratio_NW[match(hap_data$pop, ratio$sample)]
    hap_data$high_ratio_NW <- ratio$high_ratio_NW[match(hap_data$pop, ratio$sample)]
    hap_data$low_ratio_W <- ratio$low_ratio_W[match(hap_data$pop, ratio$sample)]
    hap_data$high_ratio_W <- ratio$high_ratio_W[match(hap_data$pop, ratio$sample)]
    

  #label populations based on presence/absence phenotype of interest 
	hap_data <- hap_data %>%
    		mutate(defended = case_when((pop == "R1A" & gen == "16") |
						(pop == "R14" & gen == "16") |
						(pop == "R14" & gen == "20") | 
					        (pop == "R24A" & gen == "20") ~ "yes",
					         TRUE ~ "no")) 

  #where founder haplotype could not be uniquely identified, merge haplotype ID and sum total frequency of indistinguishable haplotypes
	df2 <- hap_data %>%
     		group_by(chr, pos, pop, gen, cuttree) %>%
		    summarize(nic = nic, wasp = wasp, defended = defended, low_ratio_W = low_ratio_W, 
			  high_ratio_W = high_ratio_W, low_ratio_NW = low_ratio_NW,  high_ratio_NW = high_ratio_NW, 
		    mergeFounder = paste0(founder,collapse="_"), sumfreq = sum(freq)) %>%
	   	  ungroup() %>%
	        unique() %>%
	        dplyr::select(-cuttree)

    #find sites where haplotype frequencies are <2%
	  TooLow <- df2 %>%
	    filter(sumfreq > 0.02)
		  # mergeFounder, msumfreq 

    #remove sites where founder haplotypes are below <2%
    SmallerTable <- df2 %>%
      filter(mergeFounder %in% TooLow$mergeFounder) %>%
      mutate(afreq = asin(sqrt(sumfreq)))
		
    #subset haplotype table by specified regimes, at each genomic site model haplotype frequency as a function of phenotype
    #based on input formula, perform anova, output dataframe with significance of model terms at each site
    model_sites <- SmallerTable %>% 
      filter_(condition) %>%
      filter(!is.na(afreq)) %>%
      filter(gen == "16" | gen == "20") %>%
      group_by(chr, pos) %>%
      nest() %>% 
      mutate(fit = map(data, ~ lm(formula, data = .)),
             signif = map(fit, anova)) %>%
      unnest(signif) %>%
      select(chr, pos, 'Pr(>F)') 

    #identify sites with missing data
     missing_sites <- gen16 %>% 
     	group_by(chr, pos) %>% 
     	summarise(n=n()) %>%
     	filter(n < 8)  

    #remove sites with missing data
    gen16 <- gen16[!gen16$pos %in% missing_sites$pos, , drop = FALSE]

    #run lm at one site to determine number of model terms and set dimension of dataframe to store results 
    one_site <- SmallerTable %>%
   		filter_(condition) %>%
    		filter(gen == "16" | gen == "20") %>%
		    filter(chr == "chrX" & pos == 611074) %>%
		    lm(formula, data = .)

    aov <- (tidy(anova(one_site)))
    factors <- aov[1]
    
    multiplier <- nrow(model_sites)/nrow(factors)
    factors_long <- do.call("rbind", replicate(multiplier, factors, simplify = FALSE))

    #create table of significance for each model term at each genomic locus
    model_sites <- cbind(model_sites, factors_long)  %>% na.omit()
    
    #extract only model term of interest, calculate logp
    cd_sites <- subset(model_sites, term == term_name) %>%
      rename(p = 'Pr(>F)') %>%
      mutate(logp = -log(p)/log(10)) %>%
      ungroup() %>%
      mutate(row_id = row_number())
    
    #recode chromosome names, ensure that loci fall within correct chromosomal regions
    chi2 <- cd_sites %>% filter(!is.na(logp)) %>%
      mutate(Ichr=recode(chr,'chrX'=1,'chr2L'=2,'chr2R'=3,'chr3L'=4,'chr3R'=5)) %>%
      unite("ID", chr:pos, sep = ".", remove = FALSE) %>%
      filter( (Ichr==1 & pos>331074 & pos<22971074) | 
                (Ichr==2 & pos>212275 & pos<23172275) | 
                (Ichr==3 & pos>1181678 & pos<25041678) | 
                (Ichr==4 & pos>328250 & pos<27808250) | 
                (Ichr==5 & pos>3213137 & pos<31833137)  ) 
    
    
    #extract "QTL peaks", local maxima based on logp
    hits = chi2 %>%
      mutate(rpos = 1e6*round(pos/1e6,0)) %>%
      group_by(chr,rpos) %>%
      filter(logp == max(logp)) %>%
      filter(logp>1.3)
    
    for_hits <- cd_sites %>% unite("chr_pos", chr:pos, sep = ".")
    
    #write output for log2_drop
    hits_filename = paste(dest, "/", output_name, "_hits.txt", sep = '')
    sites_filename = paste(dest, "/", output_name, "_sites.txt", sep = '')
    write.table(for_hits, sites_filename, row.names=FALSE, sep="\t", quote=FALSE)
    write.table(hits, hits_filename, row.names=FALSE, sep="\t", quote=FALSE)

}




