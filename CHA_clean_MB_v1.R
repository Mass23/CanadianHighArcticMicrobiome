# authors: Ianina Althsuler, David Touchette, Massimo Bourquin
# date: October 2024
# description: analysis of the CHA dataset (environmental and microbial data)
library(dplyr)
library(tidyverse)
library(mgcv)
library(ggplot2)
library(gratia)
library(vegan)
library(performance)
library(ape)
library(picante)
library(phytools)
library(purrr)
library(phyloseq)
library(ggrepel)
library(ggpubr)
library(microViz)

setwd('~/Desktop/CHA')

##########################################################
####### 1. data preprocessing ############################ 
##########################################################
load_microbial_data <- function(){
    # load files
    asv_table_taxonomy = read.csv('raw_data/CHA_table_silva.csv')
    rownames(asv_table_taxonomy) = asv_table_taxonomy$ASV
    asv_table_taxonomy$ASV = NULL

    phylo_tree = read.tree('raw_data/CHA_rooted-tree.nwk')
    asv_table <- asv_table_taxonomy %>% select(-c(G7.0,taxonomy))

    # create taxonomy table
    tax_table <- asv_table_taxonomy %>% select(taxonomy) %>% 
      separate(taxonomy, c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), "; ")  %>% 
      as.matrix()
    tax_table <- gsub("d__","",
                         gsub("p__","",
                              gsub("o__","",
                                   gsub("c__","",
                                        gsub("g__","",
                                             gsub("s__","",
                                                  gsub("f__","", tax_table)))))))
    
    return(list(tax=tax_table, asv=asv_table, tree=phylo_tree))
    }

log_min_half <- function(values){
    min_half = min(values[values > 0])/2
    return(log(values+min_half))
    }

load_environmental_data <- function(){
    env_data = read.csv('raw_data/metadataCHA.csv')
    env_data$Layer <- factor(env_data$Layer)
    env_data$Sample = env_data$Well
    env_data$Well = NULL
    rownames(env_data) = env_data$Sample
    
    # vegetation as 0, 1, 2 and normalise
    env_data$Vegetation[env_data$Vegetation == 'None'] = 0
    env_data$Vegetation[env_data$Vegetation == 'Sparse'] = 1
    env_data$Vegetation[env_data$Vegetation == 'Dense'] = 2
    env_data$Vegetation = as.numeric(env_data$Vegetation)
    
    #lets try with log transformation of variables
    env_data$logTemp = log_min_half(env_data$Temperature)
    env_data$logSal = log_min_half(env_data$Salinity)
    env_data$logN = log_min_half(env_data$N)
    env_data$logMoist = log_min_half(env_data$Moisture)
    env_data$logC = log_min_half(env_data$C)
    return(env_data)
}

ggrare <- function(physeq_object, step = 10, label = NULL, color = NULL, plot = TRUE, parallel = FALSE, se = TRUE) {
  x <- methods::as(phyloseq::otu_table(physeq_object), "matrix")
  if (phyloseq::taxa_are_rows(physeq_object)) { x <- t(x) }
  ## This script is adapted from vegan `rarecurve` function
  tot <- rowSums(x)
  S <- rowSums(x > 0)
  nr <- nrow(x)
  rarefun <- function(i) {
    cat(paste("rarefying sample", rownames(x)[i]), sep = "\n")
    n <- seq(1, tot[i], by = step)
    if (n[length(n)] != tot[i]) {
      n <- c(n, tot[i])
    }
    y <- vegan::rarefy(x[i, ,drop = FALSE], n, se = se)
    if (nrow(y) != 1) {
      rownames(y) <- c(".S", ".se")
      return(data.frame(t(y), Size = n, Sample = rownames(x)[i]))
    } else {
      return(data.frame(.S = y[1, ], Size = n, Sample = rownames(x)[i]))
    }
  }
  if (parallel) {
    out <- parallel::mclapply(seq_len(nr), rarefun, mc.preschedule = FALSE)
  } else {
    out <- lapply(seq_len(nr), rarefun)
  }
  df <- do.call(rbind, out)
  # Get sample data
  if (!is.null(phyloseq::sample_data(physeq_object, FALSE))) {
    sdf <- methods::as(phyloseq::sample_data(physeq_object), "data.frame")
    sdf$Sample <- rownames(sdf)
    data <- merge(df, sdf, by = "Sample")
    labels <- data.frame(x = tot, y = S, Sample = rownames(x))
    labels <- merge(labels, sdf, by = "Sample")
  }
  # Add, any custom-supplied plot-mapped variables
  if ( length(color) > 1 ) {
    data$color <- color
    names(data)[names(data) == "color"] <- deparse(substitute(color))
    color <- deparse(substitute(color))
  }
  if ( length(label) > 1 ) {
    labels$label <- label
    names(labels)[names(labels) == "label"] <- deparse(substitute(label))
    label <- deparse(substitute(label))
  }
  p <- ggplot2::ggplot(data = data,
                       ggplot2::aes_string(x = "Size",
                                           y = ".S",
                                           group = "Sample",
                                           color = color))
  p <- p + ggplot2::labs(x = "Sequence Sample Size", y = "Species Richness")
  if (!is.null(label)) {
    p <- p + ggplot2::geom_text(data = labels,
                                ggplot2::aes_string(x = "x",
                                                    y = "y",
                                                    label = label,
                                                    color = color),
                                size = 4, hjust = 0)
  }
  p <- p + ggplot2::geom_line()
  if (se) { ## add standard error if available
    p <- p +
      ggplot2::geom_ribbon(ggplot2::aes_string(ymin = ".S - .se",
                                               ymax = ".S + .se",
                                               color = NULL,
                                               fill = color),
                           alpha = 0.2)
  }
  if (plot) {
    plot(p)
  }
  invisible(p)
}

rare_curves <- function(all_data){
    #Rarefaction curve before filtering
    rare_curves_data <- ggrare(all_data, step = 1000, plot = TRUE, parallel = FALSE, se = FALSE)
    rare_curves_plot <- rare_curves_data +
      theme_classic() +
      facet_wrap(~factor(Layer, c("Upper", "Lower"))) +
      theme(legend.title = element_text(size = 14, face="bold"), 
            legend.text = element_text(size = 14),
            axis.title.x = element_text(size = 16, face="bold"),
            axis.title.y = element_text(size = 16, face="bold"),
            axis.text.x =  element_text(size = 12), 
            axis.text.y =  element_text(size = 12, angle=90, hjust = 0.5)) +
      scale_y_continuous(breaks = c(0, 3000, 6000, 9000), limits = c(0, 9000))
    
    ggsave("figures/figure_s2.pdf", rare_curves_plot, width = 18.2, height = 18.2, units = "cm")
}

datasets_to_phyloseq <- function(asv_table, tax_table, phy_tree, env_data){
    print(head(asv_table))
    print(head(tax_table))
    print(head(env_data))
    asv_table_ps <- otu_table(asv_table, taxa_are_rows = TRUE)
    env_data_ps <- sample_data(env_data)
    tax_table_ps <- tax_table(tax_table)
    phy_tree_ps <- phy_tree(phy_tree)
    all_data <- merge_phyloseq(asv_table_ps, env_data_ps, tax_table_ps, phy_tree_ps)
    
    #Filter: remove chloroplast, mitochondria, eukaryote, remove singletons
    all_data <- subset_taxa(all_data, (Order!="Chloroplast") | is.na(Order))
    all_data <- subset_taxa(all_data, (Family!="Mitochondria") | is.na(Family))
    all_data <- subset_taxa(all_data, !is.na(Kingdom) & !Kingdom %in% c("Eukaryota","", "Unassigned") & !Genus %in% c("Bacteria", "Bacteria Kingdom"))
    all_data <- filter_taxa(all_data, function(x) sum(x) > 1, TRUE)

    # Validation and taxonomy fixes
    all_data_validated <- phyloseq_validate(all_data, remove_undetected = TRUE)
    all_data_validated  <- tax_fix(all_data_validated, # fix taxonomy
                        min_length = 4,
                        unknowns = c("uncultured","Incertae_Sedis", "Unknown_Family"),
                        sep = " ", anon_unique = TRUE,
                        suffix_rank = "classified")
    
    rare_curves(all_data_validated)
    
    return(all_data)
    }

add_alpha_to_env <- function(phyloseq_data){
  alpha_data <- estimate_richness(phyloseq_data, measures = c("Observed", "Shannon"))
  alpha_data <- cbind(sample_data(phyloseq_data), alpha_data) %>% mutate(Evenness = Shannon/log(Observed))

  pd_tab = pd(samp = t(otu_table(phyloseq_data)), tree = phy_tree(phyloseq_data))
  pd_tab$Sample = rownames(t(otu_table(phyloseq_data)))
  alpha_data$PD <- map_dbl(alpha_data$Sample, function(x) pd_tab$PD[pd_tab$Sample == x])

  sample_data(phyloseq_data) <- alpha_data
  env_alpha_data <- data.frame(sample_data(phyloseq_data))
  return(env_alpha_data)
}

##########################################################
####### 2. Alpha diversity GAMs ##########################
##########################################################
compare_layers <- function(env_data_with_alpha){
summary(lm(logASV_number ~ Layer + Vegetation + Layer:Vegetation, data = env_data_with_alpha))
summary(lm(Shannon ~ Layer + Vegetation + Layer:Vegetation, data = env_data_with_alpha))
summary(lm(logPD ~ Layer + Vegetation + Layer:Vegetation, data = env_data_with_alpha))

p1 = ggplot(env_data_with_alpha, aes(x=Layer, y=logASV_number)) + geom_boxplot(outlier.shape = NA) + geom_jitter() + facet_grid(~Vegetation)
p2 = ggplot(env_data_with_alpha, aes(x=Layer, y=Shannon)) + geom_boxplot(outlier.shape = NA) + geom_jitter() + facet_grid(~Vegetation)
p3 = ggplot(env_data_with_alpha, aes(x=Layer, y=logPD)) + geom_boxplot(outlier.shape = NA) + geom_jitter() + facet_grid(~Vegetation)
}

alpha_gams <- function(env_data_with_alpha){
~~# Full model:::
  mod_nasv_full = gam(logASV_number ~ Vegetation +
                        s(pH, bs='tp', k=6) + s(logSal, bs='tp', k=6) + s(logN, bs='tp', k=6) + 
                        s(logC, bs='tp', k=6) + s(logMoist, bs='tp', k=6) + s(logTemp, bs='tp', k=6), 
                      data = CHA_metadata, method = 'REML')
gam.check(mod_nasv_full)
summary(mod_nasv_full)
draw(mod_nasv_full, residuals = T, select = c('s(pH)', 's(logC)'))

mod_shannon_full = gam(Shannon ~ Vegetation +
                         s(pH, bs='tp', k=4) + s(logSal, bs='tp', k=4) + s(logN, bs='tp', k=4) + 
                         s(logC, bs='tp', k=4) + s(logMoist, bs='tp', k=4) + s(logTemp, bs='tp', k=), 
                       data = CHA_metadata, method = 'REML')
gam.check(mod_shannon_full)
summary(mod_shannon_full)
draw(mod_shannon_full, residuals = T, select = c('s(pH)',  's(logC)', 's(logSal)', 's(logN)'))

mod_pd_full = gam(logPD ~ s(Latitude, Longitude, bs='tp', k=20, by = as.factor(Vegetation)) +
                    s(pH, bs='ts', k=4) + s(logSal, bs='ts', k=4) + s(logN, bs='ts', k=4) + 
                    s(logC, bs='ts', k=4) + s(logMoist, bs='ts', k=4) + s(logTemp, bs='ts', k=4), 
                  data = CHA_metadata, method = 'REML')
gam.check(mod_pd_full)
summary(mod_pd_full)
draw(mod_pd_full, residuals = T, select = 's(logC)')
}

##########################################################
####### 3. Tax & Tree         ############################ 
##########################################################

##########################################################
####### 4. WGNCA              ############################ 
##########################################################

##########################################################
####### 5. Beta diversity     ############################ 
##########################################################


##########################################################
####### 6. Main function      ############################ 
##########################################################
main <- function(){
  dir.create('stats')
  dir.create('figures')
  dir.create('data')
  
  # Data loading and preprocessing
  print('loading microbial data...')
  microbial_data = load_microbial_data()
  asv_table = microbial_data$asv
  tax_table = microbial_data$tax
  phy_tree = microbial_data$tree
  
  print('loading environmental data...')
  env_data = load_environmental_data()
  
  print('merging datasets into phyloseq object...')
  phyloseq_data = datasets_to_phyloseq(asv_table, tax_table, phy_tree, env_data)
  env_data_with_alpha = add_alpha_to_env(phyloseq_data)
  write.csv(env_data_with_alpha, file='data/env_with_alpha.csv', quote = F, row.names = F, sep = ',')
    
  # Alpha diversity analyses --> GAMs (4x3, Observed, Shannon, PD)
    
  # Taxonomy & Tree --> 1 figure
  
  # WGCNA --> SI plots + 1 figure
    
  # Beta-diversity
  
}

main()


