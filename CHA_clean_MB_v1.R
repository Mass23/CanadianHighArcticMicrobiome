# authors: Ianina Althsuler, David Touchette, Massimo Bourquin
# date: March 2025
# description: analysis of the CHA dataset (environmental and microbial data)
library(dplyr)
library(mgcv)
library(ggplot2); theme_set(theme_bw())
library(ggeffects)
library(gratia)
library(vegan)
library(performance)
library(ape)
library(picante)
library(phytools)
library(purrr)
library(phyloseq)
library(microViz)
library(ggrepel)
library(WGCNA)
library(tidyverse)
library(ggtree)
library(ggtreeExtra)
library(btools)
library(ggpubr)
library(ggstance)
library(geosphere)
library(ggtext)
library(reshape2)

setwd("~/CHA")

##########################################################
####### 1. data preprocessing ############################ 
##########################################################
load_microbial_data <- function(){
    # load files
    asv_table_taxonomy = read.csv('raw_data/CHA_table_silva.csv') %>% column_to_rownames(var="ASV")
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
    
    ggsave("plots/figure_s2.pdf", rare_curves_plot, width = 18.2, height = 18.2, units = "cm")
}

datasets_to_phyloseq <- function(asv_table, tax_table, phy_tree, env_data){
    rownames(env_data) = env_data$Well
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
  alpha_data$logASV_number = log(alpha_data$Observed + 1)
  alpha_data <- cbind(sample_data(phyloseq_data), alpha_data) %>% mutate(Evenness = Shannon/log(Observed))
  sample_data(phyloseq_data) <- alpha_data
  
  #Calculate Faith's phylogenetic diversity
  pd_data <- estimate_pd(phyloseq_data)
  pd_data$logPD = log(pd_data$PD + 1)
  pd_data <- pd_data %>% select(-SR)
  all_data <- cbind(as.data.frame(sample_data(phyloseq_data)), pd_data)
  sample_data(phyloseq_data) <- all_data
  
  env_alpha_data <- data.frame(sample_data(phyloseq_data))
  return(env_alpha_data)
}

##########################################################
####### 2. Alpha diversity GAMs ##########################
##########################################################
alpha_div_gams <- function(CHA_metadata){
  
  # boxplots
  sink('stats/linear_models_alpha_div.txt')
  print(summary(lm(logASV_number ~ Layer + Vegetation + Layer:Vegetation, data = CHA_metadata)))
  print(summary(lm(Shannon ~ Layer + Vegetation + Layer:Vegetation, data = CHA_metadata)))
  print(summary(lm(logPD ~ Layer + Vegetation + Layer:Vegetation, data = CHA_metadata)))
  sink()
  
  p1 = ggplot(CHA_metadata, aes(x=Layer, y=Observed)) + geom_boxplot(outlier.shape = NA) + geom_jitter() + facet_grid(~Vegetation) + stat_compare_means()
  p2 = ggplot(CHA_metadata, aes(x=Layer, y=Shannon)) + geom_boxplot(outlier.shape = NA) + geom_jitter() + facet_grid(~Vegetation) + stat_compare_means()
  p3 = ggplot(CHA_metadata, aes(x=Layer, y=PD)) + geom_boxplot(outlier.shape = NA) + geom_jitter() + facet_grid(~Vegetation) + stat_compare_means()
  p = ggarrange(p1, p2, p3, nrow=3)
  ggsave('plots/boxplots_alpha.pdf', p, width=8, height=16)
  
  # Full models
  mod_nasv_full = gam(Observed ~ Vegetation +
                          s(pH, bs='tp', k=6) + s(logSal, bs='tp', k=6) + s(logN, bs='tp', k=6) + 
                          s(logC, bs='tp', k=6) + s(logMoist, bs='tp', k=6) + s(logTemp, bs='tp', k=6), 
                        data = CHA_metadata, method = 'REML')
  gam.check(mod_nasv_full)
  sink('stats/gam_observed_summary.txt')
  print(summary(mod_nasv_full))
  sink()
  p1 = draw(mod_nasv_full, residuals = T, select = c('s(pH)', 's(logC)', 's(logN)', 's(logSal)'), nrow = 2) + theme_minimal()
  ggsave('plots/gams_obs.pdf', p1, width = 8, height = 8)
  
  mod_shannon_full = gam(Shannon ~ Vegetation +
                           s(pH, bs='tp', k=4) + s(logSal, bs='tp', k=4) + s(logN, bs='tp', k=4) + 
                           s(logC, bs='tp', k=4) + s(logMoist, bs='tp', k=4) + s(logTemp, bs='tp', k=), 
                         data = CHA_metadata, method = 'REML')
  gam.check(mod_shannon_full)
  summary(mod_shannon_full)
  p2 = draw(mod_shannon_full, residuals = T, select = c('s(pH)', 's(logC)', 's(logN)', 's(logSal)'), nrow = 2, continuous_colour = ggplot2::scale_color_manual('darkred'))
  ggsave('plots/gams_shannon.pdf', p2, width = 8, height = 8)
  
  mod_pd_full = gam(PD ~ s(Latitude, Longitude, bs='tp', k=20, by = as.factor(Vegetation)) +
                      s(pH, bs='ts', k=4) + s(logSal, bs='ts', k=4) + s(logN, bs='ts', k=4) + 
                      s(logC, bs='ts', k=4) + s(logMoist, bs='ts', k=4) + s(logTemp, bs='ts', k=4), 
                    data = CHA_metadata, method = 'REML')
  gam.check(mod_pd_full)
  summary(mod_pd_full)
  p3 = draw(mod_pd_full, residuals = T, select = c('s(pH)', 's(logC)', 's(logN)', 's(logSal)'), nrow = 2, discrete_colour = ggplot2::scale_color_manual('darkblue'))
  ggsave('plots/gams_pd.pdf', p3, width = 8, height = 8)
  
  # main figure
  ph1 = p2[[1]] + scale_color_manual(values='tomato') + ylab('Partial effect (Shannon index)')
  ph1[[11]]$title = ''
  ph1[[11]]$x = ''

  ph2 = p3[[1]] + scale_color_manual(values='tomato') + ylab("Partial effect (Faith's PD)")
  ph2[[11]]$title = ''

  ca1 = p2[[4]] + scale_color_manual(values='aquamarine2') + ylab("")
  ca1[[11]]$title = ''
  ca1[[11]]$x = ''
  ca2 = p3[[4]] + scale_color_manual(values='aquamarine2') + ylab("")
  ca2[[11]]$title = ''

  ni1 = p2[[3]] + scale_color_manual(values='darkblue') + ylab('Partial effect (Shannon index)')
  ni1[[11]]$title = ''
  ni1[[11]]$x = ''
  ni2 = p3[[3]] + scale_color_manual(values='darkblue') + ylab("Partial effect (Faith's PD)")
  ni2[[11]]$title = ''

  sa1 = p2[[2]] + scale_color_manual(values='orange') + ylab("")
  sa1[[11]]$title = ''
  sa1[[11]]$x = ''
  sa2 = p3[[2]] + scale_color_manual(values='orange') + ylab("")
  sa2[[11]]$title = ''


  ph = ggarrange(ph1, ph2, nrow = 2, ncol = 1, labels = c('A. pH', ''), align='v')
  sal = ggarrange(sa1, sa2, nrow = 2, ncol = 1, labels = c('B. Salinity', ''), align='v')
  plot.list1 <- lapply(list(ph, sal), 
                      function(p) p + theme(plot.background = element_rect(color = "black")))
  both1 = ggarrange(plotlist =  plot.list1, nrow = 1, ncol = 2)
  
  ni = ggarrange(ni1, ni2, nrow = 2, ncol = 1, labels = c('C. Nitrogen', ''), align='v')
  ca = ggarrange(ca1, ca2, nrow = 2, ncol = 1, labels = c('D. Carbon', ''), align='v')
  plot.list2 <- lapply(list(ni, ca), 
                      function(p) p + theme(plot.background = element_rect(color = "black")))
  both2 = ggarrange(plotlist =  plot.list2, nrow = 1, ncol = 2)
  
  p = ggarrange(both1, both2, ncol = 1, nrow = 2)
  ggsave('plots/fig3.pdf', p, width = 8, height = 12)

  # fun
  library(visibly)
  
  mod_shannon_contour = gam(Shannon ~ s(pH, bs='tp', k=4) + s(logSal, bs='tp', k=4), data = CHA_metadata, method = 'REML')
  pdf('plots/gam_vis_shannon_ph_sal.pdf')
  vis.gam(mod_shannon_contour, type='response', plot.type = 'contour')
  dev.off()
  
  mod_pd_contour = gam(PD ~ s(pH, bs='tp', k=4) + s(logSal, bs='tp', k=4), data = CHA_metadata, method = 'REML')
  pdf('plots/gam_vis_shannon_ph_sal.pdf')
  vis.gam(mod_pd_contour, type='response', plot.type = 'contour')
  dev.off()
  
  plot_gam_3d(mod_shannon_contour, main_var = pH, second_var = logSal, palette='tokyo')
  plot_gam_3d(mod_pd_contour, main_var = pH, second_var = logSal, palette='tokyo')}
##########################################################
####### 3. Tax & Tree         ############################ 
##########################################################
plot_phylot_tree <- function(phyloseq_data){
  
  rel_ab_data = phyloseq_data %>% transform_sample_counts(function(x) {x/sum(x)}) 
  
  only_bacteria = subset_taxa(rel_ab_data, Kingdom == "Bacteria")
  only_archaea = subset_taxa(rel_ab_data, Kingdom == "Archaea")
  
  melt_simple_bac <- only_bacteria %>% psmelt() %>% select(OTU, Sample, Abundance, Phylum) %>% group_by(OTU, Phylum) %>% 
    summarise(mean_ab = mean(Abundance), prevalence = sum(Abundance > 0), max_ab = max(Abundance)) %>% select(OTU, mean_ab, max_ab, prevalence, Phylum) %>% arrange(-mean_ab) %>% head(1000)
  top_ten_bac = melt_simple_bac %>% group_by(Phylum) %>% summarise(sum_ab = sum(mean_ab)) %>% top_n(10) %>% pull(Phylum)
  melt_simple_bac$Phylum = vapply(melt_simple_bac$Phylum, function(x) ifelse(x %in% top_ten_bac, x, 'Others'), FUN.VALUE = character(1))
  
  melt_simple_arc <- only_archaea %>% psmelt() %>% select(OTU, Sample, Abundance, Phylum) %>% group_by(OTU, Phylum) %>% summarise(mean_ab = mean(Abundance)) %>% select(OTU, mean_ab, Phylum) %>% filter(mean_ab > 0.00001)
  
  only_bacteria_filtered_tree <- keep.tip(phy_tree(phyloseq_data), melt_simple_bac$OTU)
  only_archaea_filtered_tree <- keep.tip(phy_tree(phyloseq_data), melt_simple_arc$OTU)

  tree_bac = ggtree(only_bacteria_filtered_tree, layout="fan", open.angle=10) 
  tree_bac_abnd <- tree_bac + geom_fruit(
      data=melt_simple_bac,
      geom=geom_bar,
      mapping=aes(y=OTU, x=mean_ab, fill=Phylum),
      size=.2,
      axis.params=list(
        axis       = "x",
        text.size  = 1.8,
        hjust      = 1,
        vjust      = 0.5,
        nbreak     = 3,
      ),
      grid.params=list(),
      pwidth=0.66, 
      orientation="y", 
      stat="identity") 
  
  tree_arc = ggtree(only_archaea_filtered_tree, layout="fan", open.angle=10)
  tree_arc_abnd <- tree_arc + geom_fruit(
    data=melt_simple_arc,
    geom=geom_bar,
    mapping=aes(y=OTU, x=mean_ab, fill=Phylum),
    size=.2,
    axis.params=list(
      axis       = "x",
      text.size  = 1.8,
      hjust      = 1,
      vjust      = 0.5,
      nbreak     = 3,
    ),
    grid.params=list(),
    pwidth=0.66, 
    orientation="y", 
    stat="identity")

    ggarrange(tree_bac_abnd, tree_arc_abnd)
    ggsave('plots/fig2.pdf', width = 20, height = 10)

}
##########################################################
####### 4. WGNCA              ############################ 
##########################################################
wgnca_analysis <- function(phyloseq_data){
  #we went up to generating the MEs and the modules with the merged dynamic modules. 
  #Next step is to correlate MEs to external data and to eachother. 
  #Module-module relationships, and module-trait relatioships
  options(scipen = 999) #disables scientific notation in R
  
  asv_table = as.data.frame(otu_table(phyloseq_data))
  tax = as.data.frame(tax_table(phyloseq_data))

  #asv <- asv_table %>% column_to_rownames("rowname")
  asv <- asv_table/colSums(asv_table)
  #tax <- read_tsv("data/filtered_tax_table.tsv")
  #tax <- tax %>% column_to_rownames("rowname")
  asv_tax <- merge(tax,asv, by = 'row.names') %>% 
    select(-Kingdom, -Phylum, -Class, -Species, -Order, -Family) %>%
    rename("rowname"="Row.names")
  
  asv_tax_genus <- asv_tax %>% select (-rowname) %>%
    aggregate(. ~ Genus, sum)
  
  asv_tax_genus <- (data.frame(t(asv_tax_genus)))
  colnames(asv_tax_genus) <- asv_tax_genus[1, ]
  asv_tax_genus = asv_tax_genus[-1,]
  
  write_csv(asv_tax_genus, "data/asv_tax_genus.csv")
  asv_tax_genus <- read_csv("data/asv_tax_genus.csv")
  rn <- row.names(asv_tax_genus)
  rownames(asv_tax_genus) <- rn
  
  ###########
  
  gsg = goodSamplesGenes(asv_tax_genus[-1], verbose = 3);
  gsg$allOK
  
  ## Cluster the samples to see if there are any obvious outliers
  sampleTree = hclust(dist(asv_tax_genus), method = "average");
  sizeGrWindow(12,9)
  par(cex = 0.6);
  par(mar = c(0,4,2,0))
  plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
       cex.axis = 1.5, cex.main = 2)
  
  CHAtraitData = read.csv("data/env_with_alpha.csv");
  dim(CHAtraitData)
  names(CHAtraitData)
  
  rnd <- (CHAtraitData$Sample)
  rownames(CHAtraitData) <- rnd
  CHAtraitData <- CHAtraitData %>% select(-Sample) %>% 
    select(-Site , -Island, -Date, -Moisture, -N, -C, -Salinity, -Latitude, -Temperature)
  
  CHAtraitData$Layer =  ifelse(CHAtraitData$Layer == 'Upper', 0, 1)
  
  ## Visualize how Environmental traits relate to clustered samples (visualized with the dendogram from above)
  ## Re-cluster samples:
  sampleTree2 = hclust(dist(asv_tax_genus), method = "average")
  ## Convert traits to a color representation: white means low, red means high, grey means missing entry:
  to_keep_vars = c('Altitude', 'Layer', 'logMoist', 'Longitude', 'logTemp', 'logSal', 'logN', 'logC', 'PD', 'Observed', 'Shannon', 'Evenness')
  traits2keepy <- CHAtraitData %>% select(all_of(to_keep_vars))
  traitColors = numbers2colors(traits2keepy, signed = FALSE);
  ## Plot the sample dendrogram and the colors underneath:
  plotDendroAndColors(sampleTree2, traitColors,
                      groupLabels = names(traits2keepy),
                      main = "Sample dendrogram and trait heatmap")
  ## White means a low value and red means a high value. Gray means missing entry.
  
  ##visualze just one trait with the dendrogram:
  traitColors = numbers2colors(CHAtraitData$logSal, signed = FALSE);
  plotDendroAndColors(sampleTree2, traitColors,
                      groupLabels = 'log Salinity',
                      main = "Sample dendrogram and trait heatmap")
  
  
  ## Save
  save(asv_tax_genus, CHAtraitData, file = "data/CHA-dataInput_genus.RData")
  
  ## Start of network analysis
  options(stringsAsFactors = FALSE);
  allowWGCNAThreads()
  ## Load the data saved in the first part
  lnames = load(file = "data/CHA-dataInput_genus.RData");
  ## The variable lnames contains the names of loaded variables.
  lnames
  ## Choose a set of soft-thresholding powers:
  powers = c(c(1:10), seq(from = 11, to=30, by=1))
  
  ## Call the network topology analysis function:
  ## Note: using a signed network because it preserves the sign of the connection (whether nodes are positively or negatively correlated); this is recommendation by authors of WGCNA:
  sft = pickSoftThreshold(asv_tax_genus, powerVector = powers, verbose = 5, networkType = "signed")
  
  ## Plot the results:
  sizeGrWindow(9, 5)
  par(mfrow = c(1,2));
  cex1 = 0.9;
  
  ## Scale-free topology fit index as a function of the soft-thresholding power:
  plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
       main = paste("Scale independence"));
  text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       labels=powers,cex=cex1,col="red");
  
  ## This line corresponds to using an R^2 cut-off of h:
  abline(h=0.8,col="red")
  
  ## Mean connectivity as a function of the soft-thresholding power:
  plot(sft$fitIndices[,1], sft$fitIndices[,5],
       xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
       main = paste("Mean connectivity"))
  text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
  
  ## Calculate the adjacencies, using the soft thresholding power of ___:
  softPower = 9;
  adjacency = adjacency(asv_tax_genus, power = softPower, type = "signed");
  
  ## Transform adjacency into Topological Overlap Matrix and calculate corresponding dissimilarity:
  ## Note: The TOM you calculate shows the topological similarity of nodes, factoring in the connection strength two nodes share with other "third party" nodes  
  ## This will minimize effects of noise and spurious associations:
  TOM = TOMsimilarity(adjacency, TOMType = "signed");
  dissTOM = 1-TOM
  ## Create a dendogram using a hierarchical clustering tree
  ## Call the hierarchical clustering function
  TaxaTree = hclust(as.dist(dissTOM), method = "average");
  ## Plot the resulting clustering tree (dendrogram)
  sizeGrWindow(12,9)
  plot(TaxaTree, xlab="", sub="", main = "Taxa clustering on TOM-based dissimilarity",
       labels = FALSE, hang = 0.04);
  ## decide the optimal module size
  minModuleSize = 30;
  ## Module identification using dynamic tree cut:
  dynamicMods = cutreeDynamic(dendro = TaxaTree, distM = dissTOM, deepSplit = 2,
                              pamRespectsDendro = FALSE, method = "tree",
                              minClusterSize = minModuleSize);
  table(dynamicMods)
  ## Convert numeric lables into colors
  dynamicColors = labels2colors(dynamicMods)
  table(dynamicColors)
  ## Plot the dendrogram with module colors underneath
  sizeGrWindow(8,6)
  plotDendroAndColors(TaxaTree, dynamicColors, "Dynamic Tree Cut",
                      dendroLabels = FALSE, hang = 0.03,
                      addGuide = TRUE, guideHang = 0.05,
                      main = "Gene dendrogram and module colors")
  
  ## Quatify co-expression similarity of the entire modules using eigengenes and cluster them based on their correlation:
  ## Note: An eigengene is 1st principal component of a module expression matrix and represents a suitably defined average OTU community
  ## Calculate eigengenes
  MEList = moduleEigengenes(asv_tax_genus, colors = dynamicColors)
  MEs = MEList$eigengenes
  MEDiss = 1-cor(MEs);## Calculate dissimilarity of module eigengenes
  METree = hclust(as.dist(MEDiss), method = "average");## Cluster module eigengenes
  sizeGrWindow(7, 6)#plot
  plot(METree, main = "Clustering of module eigengenes",
       xlab = "", sub = "")
  
  ## Chose a height cut of 0.30, corresponding to a similarity of 0.70 to merge--so nothing got merged
  MEDissThres = 0.30
  abline(h=MEDissThres, col = "red")## Plot the cut line into the dendrogram
  ## Call an automatic merging function
  merge = mergeCloseModules(asv_tax_genus, dynamicColors, cutHeight = MEDissThres, verbose = 3)
  mergedColors = merge$colors;## The merged module colors
  mergedMEs = merge$newMEs;
  
  ## If you had combined different modules then that would show in this plot:
  sizeGrWindow(12, 9)
  plotDendroAndColors(TaxaTree, cbind(dynamicColors, mergedColors),
                      c("Dynamic Tree Cut", "Merged dynamic"),
                      dendroLabels = FALSE, hang = 0.03,
                      addGuide = TRUE, guideHang = 0.05)
  
  moduleColors = mergedColors## Rename to moduleColors
  ## Construct numerical labels corresponding to the colors
  colorOrder = c("grey", standardColors(50));
  moduleLabels = match(moduleColors, colorOrder)-1;
  MEs = mergedMEs;
  ## Save module colors and labels for use in subsequent parts
  save(MEs, moduleLabels, moduleColors, TaxaTree, file = "data/CHA-networkConstruction-stepByStep_genus.RData")
  
  genus <-colnames(asv_tax_genus)
  
  g_a <-t(asv_tax_genus)%>% data.frame()
  g_a <- g_a %>% tibble::rownames_to_column() %>% rename(genus = rowname)
  
  genus_mods <- rbind(genus, moduleColors, moduleLabels) %>% data.frame() %>% t() %>% data.frame()
  
  genus_mods <- right_join(genus_mods,g_a)
  
  write_csv(genus_mods, file="data/genus_modules.csv")
  
  genera_summary_df = data.frame()
  for (module in unique(genus_mods$moduleColors)){
    genera_summary_df = rbind(genera_summary_df, data.frame(Module=module, Genera=unique(paste0(genus_mods$genus[genus_mods$moduleColors == module], collapse = ','))))
  }
  
  ## Extra Visualization--Can due after modules are assigned
  cmd1=cmdscale(as.dist(dissTOM),2)
  par(mfrow=c(1,1)) 
  plot(cmd1, col=as.character(dynamicColors),  main="
     MDS plot",xlab="Scaling 
    Dimension 1",ylab="Scaling Dimension 2",
       cex.axis=1.5,cex.lab=1.5, cex.main=1.5) 
  ## Vs. your merged modules (this will look identical to the first MDS plot since we didn't merge any modules):
  cmd1=cmdscale(as.dist(dissTOM),2)
  par(mfrow=c(1,1)) 
  plot(cmd1, col=as.character(mergedColors),  main="
    MDS plot",xlab="Scaling 
   Dimension 1",ylab="Scaling Dimension 2",
       cex.axis=1.5,cex.lab=1.5, cex.main=1.5) 
  
  #############
  
  ## Identify modules that are significantly associated with the measured Environmental traits:
  ## You already have summary profiles for each module (eigengenes), so we just have to correlate these eigengenes with Environmental traits and look for significant associations:
  ## First, define numbers of genera and samples
  #####################################################################################################################
  
  load("data/CHA-networkConstruction-stepByStep_genus.RData")
  load("data/CHA-dataInput_genus.RData")
  nTaxa = ncol(asv_tax_genus);
  nSamples = nrow(asv_tax_genus);
  
  ## Recalculate MEs (module eigengenes) with color labels
  MEs0 = moduleEigengenes(asv_tax_genus, moduleColors)$eigengenes
  MEs = orderMEs(MEs0)
  #remove diversity indexes????????????
  traits2keepy2 <- traits2keepy %>% select(-PD,-Shannon,-Evenness,-Observed)
  moduleTraitFullCor = cor(cbind(MEs, traits2keepy2), use = "p");
  moduleTraitFullPvalue = corPvalueStudent(moduleTraitFullCor, nSamples);
  moduleTraitFullPvalue = moduleTraitFullPvalue %>% as.matrix %>% as.vector %>%
    p.adjust(method='fdr') %>% matrix(ncol=ncol(moduleTraitFullCor))
  
  ## Now visualize it:
  sizeGrWindow(10,6)
  
  ## Will display correlations and their p-values
  textMatrix = paste(signif(moduleTraitFullCor, 2), "\n(",
                     signif(moduleTraitFullPvalue, 1), ")", sep = "");
  dim(textMatrix) = dim(moduleTraitFullCor)
  par(mar = c(6, 8.5, 3, 3));
  
  ## Display the correlation values within a heatmap plot
  labeledHeatmap(Matrix = moduleTraitFullCor,
                 xLabels = colnames(moduleTraitFullCor),
                 yLabels = rownames(moduleTraitFullCor),
                 ySymbols = names(moduleTraitFullCor),
                 colorLabels = FALSE,
                 colors = greenWhiteRed(50),
                 textMatrix = textMatrix,
                 setStdMargins = FALSE,
                 cex.text = 0.6,
                 zlim = c(-1,1),
                 main = paste("All Module-trait relationships"))
  
  ##################################################################################################################################
  
  ## Recalculate MEs (module eigengenes) with color labels
  MEs0 = moduleEigengenes(asv_tax_genus, moduleColors)$eigengenes
  MEs = orderMEs(MEs0)
  #CHAtraitData <- CHAtraitData %>% select(-Temperature, -C, -N, -Salinity, -Altitude, -Longitude, -Latitude, -Moisture, -PD, -Evenness, -Observed, -Shannon )
  moduleCor = cor(MEs, use = "p");
  modulePvalue = corPvalueStudent(moduleCor, nSamples);
  modulePvalue = modulePvalue %>% as.matrix %>% as.vector %>% 
    p.adjust(method='fdr') %>% matrix(ncol=ncol(moduleCor))
  
  ## Now visualize it:
  sizeGrWindow(10,6)
  
  ## Will display correlations and their p-values
  textMatrix = paste(signif(moduleCor, 2), "\n(",
                     signif(modulePvalue, 1), ")", sep = "");
  dim(textMatrix) = dim(moduleCor)
  par(mar = c(6, 8.5, 3, 3));
  
  ## Display the correlation values within a heatmap plot
  labeledHeatmap(Matrix = moduleCor,
                 xLabels = colnames(moduleCor),
                 yLabels = rownames(moduleCor),
                 ySymbols = names(moduleCor),
                 colorLabels = FALSE,
                 colors = redWhiteGreen(50),
                 textMatrix = textMatrix,
                 setStdMargins = FALSE,
                 cex.text = 0.5,
                 zlim = c(-1,1),
                 main = paste("Module relationships"))
  
  #########################################################################################################################################################################
  
  ## Recalculate MEs (module eigengenes) with color labels
  MEs0 = moduleEigengenes(asv_tax_genus, moduleColors)$eigengenes
  MEs = orderMEs(MEs0)
  #remove diversity?
  moduleTraitCor = cor(MEs, traits2keepy2, use = "p");
  moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);
  moduleTraitPvalue = moduleTraitPvalue %>% as.matrix %>% as.vector %>% p.adjust(method='fdr') %>% matrix(ncol=ncol(moduleTraitCor))
  
  ## Now visualize it:
  sizeGrWindow(10,6)
  
  ## Will display correlations and their p-values
  textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                     signif(moduleTraitPvalue, 1), ")", sep = "");
  dim(textMatrix) = dim(moduleTraitCor)
  par(mar = c(6, 8.5, 3, 3));
  
  ## Display the correlation values within a heatmap plot
  labeledHeatmap(Matrix = moduleTraitCor,
                 xLabels = colnames(moduleTraitCor),
                 yLabels = rownames(moduleTraitCor),
                 ySymbols = names(moduleTraitCor),
                 colorLabels = FALSE,
                 colors = blueWhiteRed(100),
                 textMatrix = textMatrix,
                 setStdMargins = FALSE,
                 cex.text = 0.6,
                 zlim = c(-1,1),
                 main = paste("Module-Env. variable relationships"))
  
  
  #########################################################################################################################################################################
  TraitCor = cor(traits2keepy2, use = "p");
  TraitPvalue = corPvalueStudent(TraitCor, nSamples);
  TraitPvalue = TraitPvalue %>% as.matrix %>% as.vector %>% 
    p.adjust(method='fdr') %>% matrix(ncol=ncol(TraitCor))
  
  ## Now visualize it:
  sizeGrWindow(10,6)
  
  ## Will display correlations and their p-values
  textMatrix = paste(signif(TraitCor, 2), "\n(",
                     signif(TraitPvalue, 1), ")", sep = "");
  dim(textMatrix) = dim(TraitCor)
  par(mar = c(6, 8.5, 3, 3));
  
  ## Display the correlation values within a heatmap plot
  labeledHeatmap(Matrix = TraitCor,
                 xLabels = colnames(TraitCor),
                 yLabels = rownames(TraitCor),
                 ySymbols = names(TraitCor),
                 colorLabels = FALSE,
                 colors = greenWhiteRed(50),
                 textMatrix = textMatrix,
                 setStdMargins = FALSE,
                 cex.text = 0.5,
                 zlim = c(-1,1),
                 main = paste("Environmental variable relationships"))
  
  
  
  #########################################################################################################################################################################
  
  #Eigengene clustering
  
  MET<- orderMEs(cbind(MEs, traits2keepy2))
  MET<- MET %>% select(colnames(traits2keepy2),everything())
  
  pdf("plots/WGCNA_summary.pdf", width = 22, height=22)
  
  
  plotEigengeneNetworks(MET, "Module-Module-Env-Env relationships", marHeatmap = c(3,4,1,2), 
                        xLabelsAngle=25,
                        printAdjacency = T, plotAdjacency = T, plotDendrograms = T, 
                        colorLabels = T, excludeGrey = T,
                        cex.adjacency = 1)
  
  dev.off()
}
##########################################################
####### 5. Beta diversity     ############################ 
##########################################################
beta_diversity_analyses <- function(CHA_filt, CHA_metadata){
    set.seed(19950930)
    #Hellinger transformation
    CHA_filt_rel <- transform_sample_counts(CHA_filt, function(x) x / sum(x))
    CHA_filt_hell <- CHA_filt_rel %>% microbiome::transform(transform = "hell")
    
    #Gloom phyloseq at phylum 
    CHA_filt_phylum <- tax_glom(CHA_filt, taxrank="Phylum")
    CHA_filt_otu_phylum <- data.frame(tax_table(CHA_filt_phylum)) %>% rownames_to_column("ASV") %>% dplyr::select(c(ASV, Phylum))
    bray_asv_table <- data.frame(otu_table(CHA_filt_phylum)) %>% 
      rownames_to_column("ASV") %>% 
      left_join(CHA_filt_otu_phylum) %>% 
      column_to_rownames("Phylum") %>% 
      dplyr::select(-c("ASV"))

    #Phylum abundance
    Phylum_abundance <- bray_asv_table %>% 
      mutate(Sum = rowSums(.[1:76]),
             Total = sum(Sum),
             Average = (Sum/Total)*100) %>% 
      dplyr::select(c("Average")) %>% 
      rownames_to_column("Phylum")
    
    Top_phylum <- Phylum_abundance %>% 
      filter(Average >= 0.1) 
    
    #Gloom at phylum 
    CHA_filt_hell_phylum <- tax_glom(CHA_filt_hell, taxrank="Phylum")
    CHA_filt_hell_otu_phylum <- data.frame(tax_table(CHA_filt_hell_phylum)) %>% rownames_to_column("ASV") %>% dplyr::select(c(ASV, Phylum))
    bray_asv_table_hell <- data.frame(otu_table(CHA_filt_hell_phylum)) %>% 
      rownames_to_column("ASV") %>% 
      left_join(CHA_filt_hell_otu_phylum) %>% 
      column_to_rownames("Phylum") %>% 
      dplyr::select(-c("ASV"))
    bray_asv_table_hell <- bray_asv_table %>% t() %>% data.frame()
    
    ##NMDS Weighted Unifrac 
    tree <- phy_tree(CHA_filt_hell)
    new_tree <- ape::multi2di(tree)
    phy_tree(CHA_filt_hell) <- new_tree
    
    CHA_filt_unifrac <- UniFrac(CHA_filt_hell, weighted = TRUE) %>% as.matrix()
    CHA_filt_unifrac_dist <- dist(CHA_filt_unifrac)
    nmds_CHA_unifrac <- metaMDS(CHA_filt_unifrac_dist)
    stress_unifrac <- nmds_CHA_unifrac$stress
    stress_unifrac
    scores_nmds_CHA_unifrac <- as.data.frame(nmds_CHA_unifrac$points) %>% as_tibble(rownames = "Well") %>% inner_join(., CHA_metadata, by="Well") #combined metadata
    scores_nmds_CHA_unifrac$Vegetation[scores_nmds_CHA_unifrac$Vegetation == 0] = 'none'
    scores_nmds_CHA_unifrac$Vegetation[scores_nmds_CHA_unifrac$Vegetation == 1] = 'sparse'
    scores_nmds_CHA_unifrac$Vegetation[scores_nmds_CHA_unifrac$Vegetation == 2] = 'dense'
    
    #NMDS environmental parameters Unifrac
    env_CHA_unifrac <- envfit(nmds_CHA_unifrac, env_data_with_alpha %>% select(-Well, -Site, -Date, -Sample, -Temperature, -Moisture, -N, -C, -Salinity, -Shannon, -Observed, -Evenness, -PD, -Latitude, -Longitude), permutations = 999, na.rm = TRUE) # select mostly numerical variables to display: pH and C are significant, Salinity is relevant
    fit_CHA_unifrac_vec <- data.frame(env_CHA_unifrac$vectors$arrows * sqrt(env_CHA_unifrac$vectors$r), P = env_CHA_unifrac$vectors$pvals, R = env_CHA_unifrac$vectors$r)
    fit_CHA_unifrac_vec$Condition <- rownames(fit_CHA_unifrac_vec)
    fit_CHA_unifrac_layer <- data.frame(env_CHA_unifrac$factors$centroids[c(1,2),] * sqrt(data.frame(env_CHA_unifrac$factors$r)[1,1]), P = data.frame(env_CHA_unifrac$factors$pvals)[1,1], R = env_CHA_unifrac$factors$r)
    fit_CHA_unifrac_layer <- fit_CHA_unifrac_layer %>% mutate(Condition = "Layer")
    fit_CHA_unifrac_env <- rbind(fit_CHA_unifrac_vec,fit_CHA_unifrac_layer)
    fit_CHA_unifrac_env_pv <- subset(fit_CHA_unifrac_env, P<=0.05) %>% 
      mutate(Condition = str_replace(Condition, "logSal", "Salinity"),
             Condition = str_replace(Condition, "logN", "Nitrogen"),
             Condition = str_replace(Condition, "logMoist", "Moisture"),
             Condition = str_replace(Condition, "logC", "Carbon"),
             Condition = str_replace(Condition, "logTemp", "Temperature"))
    fit_CHA_unifrac_env_pv
    
    #NMDS env Unifrac plot
    NMDS_CHA_env_unifrac <- ggplot(scores_nmds_CHA_unifrac, aes(x=MDS1, y=MDS2)) +
      scale_fill_gradient2(low = "red2", mid = "yellow", high = "darkgreen", midpoint = 5) +
      scale_shape_manual(values=c(24, 22, 21)) +
      scale_y_continuous(limits = c(-0.075, 0.05), breaks = c(-0.075, -0.05, -0.025, 0.0, 0.025, 0.05)) +
      scale_x_continuous(limits = c(-0.075, 0.05), breaks = c(-0.075, -0.05, -0.025, 0.0, 0.025, 0.05)) +
      geom_point(aes(shape = Vegetation, fill = pH), size = 3)+
      theme_classic() +
      labs(caption="Stress 0.093",
           fill = "pH") +
      theme(legend.title = element_text(size = 14, face="bold"), 
            legend.text = element_text(size = 14),
            axis.title.x = element_text(size = 16, face="bold"),
            axis.title.y = element_text(size = 16, face="bold"),
            axis.text =  element_text(size = 12),
            legend.position = "right",
            legend.justification = "left",
            strip.text = element_text(size=14, face="bold"),
            strip.background = element_blank(),
            plot.caption = element_text(size=12, hjust=1, vjust = 25),
            axis.line=element_line()) +
      geom_segment(data = fit_CHA_unifrac_env_pv, aes(x = 0, xend=NMDS1*0.08, y=0, yend=NMDS2*0.08), colour = "black", size = 0.01, arrow=arrow(length=unit(0.25, "cm")), lwd=0.5) +
      geom_text_repel(data = fit_CHA_unifrac_env_pv, aes(x=NMDS1*0.08, y=NMDS2*0.08, label = Condition), direction = "x", colour = "black", segment.size = 0.25)+
      guides(shape = guide_legend(override.aes = list(linewidth = 5)), 
             size = guide_legend(override.aes = list(shape = 5, linewidth = c(2, 3, 4, 5))))
    ggsave('plots/NMDS_CHA_env_unifrac.pdf', plot = NMDS_CHA_env_unifrac)
    
    #NMDS taxa Unifrac
    set.seed(19950930)
    sp_CHA_unifrac <- envfit(nmds_CHA_unifrac, bray_asv_table_hell, permutations = 999)
    sp.CHA_unifrac <- data.frame(sp_CHA_unifrac$vectors$arrows)
    sp.CHA_unifrac <- cbind(sp.CHA_unifrac, Phylum = rownames(sp.CHA_unifrac))
    sp.CHA_unifrac <- cbind(sp.CHA_unifrac, R2 = sp_CHA_unifrac$vectors$r)
    sp.CHA_unifrac <- cbind(sp.CHA_unifrac, pval = sp_CHA_unifrac$vectors$pvals)
    sp.CHA_unifrac$Phylum <- gsub("\\.", "-", sp.CHA_unifrac$Phylum)
    
    sp_CHA_unifrac_sign_high <- sp.CHA_unifrac %>% filter(pval <= 0.01) %>% filter(R2 > 0.15) %>% filter(Phylum %in% Top_phylum$Phylum)
    
    #NMDS taxa Unifrac plot
    NMDS_CHA_taxa_unifrac <- ggplot(scores_nmds_CHA_unifrac, aes(x=MDS1, y=MDS2)) +
      scale_fill_gradient2(low = "red2", mid = "yellow", high = "darkgreen", midpoint = 5) +
      scale_shape_manual(values=c(21, 22, 24)) +
      scale_y_continuous(limits = c(-0.075, 0.05), breaks = c(-0.075, -0.05, -0.025, 0.0, 0.025, 0.05)) +
      scale_x_continuous(limits = c(-0.075, 0.05), breaks = c(-0.075, -0.05, -0.025, 0.0, 0.025, 0.05)) +
      geom_point(aes(shape = Vegetation, fill = pH), size = 3)+
      theme_classic() +
      labs(caption="Stress 0.093",
           fill = "pH") +
      theme(legend.title = element_text(size = 14, face="bold"), 
            legend.text = element_text(size = 14),
            axis.title.x = element_text(size = 16, face="bold"),
            axis.title.y = element_text(size = 16, face="bold"),
            axis.text =  element_text(size = 12),
            legend.position = "right",
            legend.justification = "left",
            strip.text = element_text(size=14, face="bold"),
            strip.background = element_blank(),
            plot.caption = element_text(size=12, hjust=1, vjust = 25),
            axis.line=element_line()) +
      geom_segment(data = sp_CHA_unifrac_sign_high, aes(x = 0, xend=NMDS1*0.05, y=0, yend=NMDS2*0.05), colour = "darkblue", size = 0.01, arrow=arrow(length=unit(0.25, "cm")), lwd=0.5) +
      geom_text_repel(data = sp_CHA_unifrac_sign_high, aes(x=NMDS1*0.05, y=NMDS2*0.05, label = Phylum), cex = 4, direction = "x", colour = "darkblue", segment.size = 0.25)+
      guides(shape = guide_legend(override.aes = list(linewidth = 5)), 
             size = guide_legend(override.aes = list(shape = 5, linewidth = c(2, 3, 4, 5))))
    ggsave('plots/NMDS_CHA_taxa_unifrac.pdf', plot = NMDS_CHA_taxa_unifrac)
    
    nmds_unifrac_both = ggarrange(NMDS_CHA_env_unifrac, NMDS_CHA_taxa_unifrac, ncol = 2, labels = c('A', 'B'), common.legend = T)
    ggsave('plots/NMDS_unifrac.pdf', plot = nmds_unifrac_both, width = 14, height = 7)
    
    
    ##NMDS Bray Curtis 
    matrix_dist_CHA <- as.matrix(t(data.frame(otu_table(CHA_filt_hell))))
    metadata_dist_CHA <- data.frame(sample_data(CHA_filt_hell)) %>% rownames_to_column("Samples")
    dist_CHA_filt_bray <- vegdist(matrix_dist_CHA, method ="bray") %>% 
      as.matrix() %>% 
      data.frame() %>% 
      rownames_to_column("Samples")
    dist_CHA_bray <- dist_CHA_filt_bray %>% 
      dplyr::select(all_of(.[["Samples"]])) %>% 
      as.dist()
    nmds_CHA_bray <- metaMDS(dist_CHA_bray)
    stress_bray <- nmds_CHA_bray$stress
    stress_bray
    scores_nmds_CHA_bray <- as.data.frame(nmds_CHA_bray$points) %>% as_tibble(rownames = "Well") %>% inner_join(., CHA_metadata, by="Well") #combined metadata
    scores_nmds_CHA_bray$Vegetation[scores_nmds_CHA_unifrac$Vegetation == 0] = 'none'
    scores_nmds_CHA_bray$Vegetation[scores_nmds_CHA_unifrac$Vegetation == 1] = 'sparse'
    scores_nmds_CHA_bray$Vegetation[scores_nmds_CHA_unifrac$Vegetation == 2] = 'dense'
    
    #NMDS environmental parameters bray
    env_CHA_bray <- envfit(nmds_CHA_bray, env_data_with_alpha %>% select(-Well, -Site, -Date, -Sample, -Temperature, -Moisture, -N, -C, -Salinity, -Shannon, -Observed, -Evenness, -PD, -Latitude, -Longitude), permutations = 999, na.rm = TRUE) # select mostly numerical variables to display: pH and C are significant, Salinity is relevant
    fit_CHA_bray_vec <- data.frame(env_CHA_bray$vectors$arrows * sqrt(env_CHA_bray$vectors$r), P = env_CHA_bray$vectors$pvals, R = env_CHA_bray$vectors$r)
    fit_CHA_bray_vec$Condition <- rownames(fit_CHA_bray_vec)
    fit_CHA_bray_layer <- data.frame(env_CHA_bray$factors$centroids[c(1,2),] * sqrt(data.frame(env_CHA_bray$factors$r)[1,1]), P = data.frame(env_CHA_bray$factors$pvals)[1,1], R = env_CHA_bray$factors$r)
    fit_CHA_bray_layer <- fit_CHA_bray_layer %>% mutate(Condition = "Depth")
    fit_CHA_bray_env <- rbind(fit_CHA_bray_vec,fit_CHA_bray_layer)
    fit_CHA_bray_env_pv <- fit_CHA_bray_env %>% filter(P<=0.05) %>% 
      mutate(Condition = str_replace(Condition, "logTemp", "Temperature"),
             Condition = str_replace(Condition, "logSal", "Salinity"),
             Condition = str_replace(Condition, "logN", "Nitrogen"),
             Condition = str_replace(Condition, "logMoist", "Moisture"),
             Condition = str_replace(Condition, "logC", "Carbon"),
             Condition = str_replace(Condition, "VegNum", "Vegetation"))
    fit_CHA_bray_env_pv
    
    #NMDS env Bray plot
    NMDS_CHA_env_bray <- ggplot(scores_nmds_CHA_bray, aes(x=MDS1, y=MDS2)) +
      scale_fill_gradient2(low = "red2", mid = "yellow", high = "darkgreen", midpoint = 5) +
      scale_shape_manual(values=c(21, 22, 24)) +
      scale_y_continuous(limits = c(-0.35, 0.3), breaks = c(-0.3, -0.15, 0.0, 0.15, 0.3)) +
      scale_x_continuous(limits = c(-0.55, 0.6), breaks = c(-0.4, -0.2, 0.0, 0.2, 0.4, 0.6)) +
      geom_point(aes(shape = Vegetation, fill = pH), size = 3)+
      theme_classic() +
      labs(caption="Stress 0.151") +
      theme(legend.title = element_text(size = 14, face="bold"), 
            legend.text = element_text(size = 14),
            axis.title.x = element_text(size = 16, face="bold"),
            axis.title.y = element_text(size = 16, face="bold"),
            axis.text =  element_text(size = 12),
            legend.position = "right",
            legend.justification = "left",
            strip.text = element_text(size=14, face="bold"),
            strip.background = element_blank(),
            plot.caption = element_text(size=12, hjust=1, vjust = 25),
            axis.line=element_line()) +
      geom_segment(data = fit_CHA_bray_env_pv, aes(x = 0, xend=NMDS1*0.5, y=0, yend=NMDS2*0.5), colour = "black", size = 0.01, arrow=arrow(length=unit(0.25, "cm")), lwd=0.5) +
      geom_text_repel(data = fit_CHA_bray_env_pv, aes(x=NMDS1*0.5, y=NMDS2*0.5, label = Condition), cex = 4, direction = "x", colour = "black", segment.size = 0.25)+
      guides(shape = guide_legend(override.aes = list(linewidth = 5)), 
             size = guide_legend(override.aes = list(shape = 5, linewidth = c(2, 3, 4, 5))))
    ggsave('plots/NMDS_CHA_env_bray.pdf', plot = NMDS_CHA_env_bray)
    
    #NMDS taxa Bray
    set.seed(19950930)
    sp_CHA_bray <- envfit(nmds_CHA_bray, bray_asv_table_hell, permutations = 999)
    sp.CHA_bray <- as.data.frame(sp_CHA_bray$vectors$arrows)
    sp.CHA_bray <- cbind(sp.CHA_bray, Phylum = rownames(sp.CHA_bray))
    sp.CHA_bray <- cbind(sp.CHA_bray, R2 = sp_CHA_bray$vectors$r)
    sp.CHA_bray <- cbind(sp.CHA_bray, pval = sp_CHA_bray$vectors$pvals)
    sp.CHA_bray$Phylum <- gsub("\\.", "-", sp.CHA_bray$Phylum)
    
    sp_CHA_bray_sign_high <- sp.CHA_bray %>% filter(pval <= 0.01) %>% filter(R2 > 0.15) %>% filter(Phylum %in% Top_phylum$Phylum)
    
    #NMDS taxa Bray plot
    NMDS_CHA_taxa_bray <- ggplot(scores_nmds_CHA_bray, aes(x=MDS1, y=MDS2)) +
      scale_fill_gradient2(low = "red2", mid = "yellow", high = "darkgreen", midpoint = 5) +
      scale_shape_manual(values=c(21, 22, 24)) +
      geom_point(aes(shape = Vegetation, fill = pH), size = 3)+
      theme_classic() +
      labs(caption="Stress 0.151") +
      scale_y_continuous(limits = c(-0.35, 0.3), breaks = c(-0.3, -0.15, 0.0, 0.15, 0.3)) +
      scale_x_continuous(limits = c(-0.55, 0.6), breaks = c(-0.4, -0.2, 0.0, 0.2, 0.4, 0.6)) +
      theme(legend.title = element_text(size = 14, face="bold"), 
            legend.text = element_text(size = 14),
            axis.title.x = element_text(size = 16, face="bold"),
            axis.title.y = element_text(size = 16, face="bold"),
            axis.text =  element_text(size = 12),
            legend.position = "right",
            legend.justification = "left",
            strip.text = element_text(size=14, face="bold"),
            strip.background = element_blank(),
            plot.caption = element_text(size=12, hjust=1, vjust = 25),
            axis.line=element_line()) +
      geom_segment(data = sp_CHA_bray_sign_high, aes(x = 0, xend=NMDS1*0.3, y=0, yend=NMDS2*0.3), colour = "darkblue", size = 0.01, arrow=arrow(length=unit(0.25, "cm")), lwd=0.5) +
      geom_text_repel(data = sp_CHA_bray_sign_high, aes(x=NMDS1*0.3, y=NMDS2*0.3, label = Phylum), cex = 4, direction = "x", colour = "darkblue", segment.size = 0.25)+
      guides(shape = guide_legend(override.aes = list(linewidth = 5)), 
             size = guide_legend(override.aes = list(shape = 5, linewidth = c(2, 3, 4, 5))))
    ggsave('plots/NMDS_CHA_taxa_bray.pdf', plot = NMDS_CHA_taxa_bray)
    
    
    nmds_bray_both = ggarrange(NMDS_CHA_env_bray, NMDS_CHA_taxa_bray, ncol = 2, labels = c('A', 'B'), common.legend = T)
    ggsave('plots/NMDS_bray.pdf', plot = nmds_bray_both, width = 14, height = 7)


    # Distance decays    
    c_sort_collapse <- function(...){
      c(...) %>%    
        sort() %>% 
        str_c(collapse = ".")
    }
    
    CHA_filt_unifrac_df <- CHA_filt_unifrac %>% data.frame() %>% rownames_to_column() 
    CHA_filt_unifrac_long <- setNames(melt(CHA_filt_unifrac_df), c('Sample1', 'Sample2', 'UF'))
    CHA_filt_bray_long <- setNames(melt(dist_CHA_filt_bray), c('Sample1', 'Sample2', 'BC'))
    CHA_geography <- metadata_dist_CHA %>% dplyr::select(Samples, Latitude, Longitude, Island)
    CHA_distances <- CHA_filt_unifrac_long %>% 
      left_join(CHA_filt_bray_long) %>% 
      left_join(CHA_geography, by= join_by(Sample1 == Samples)) %>% 
      left_join(CHA_geography, by= join_by(Sample2 == Samples)) %>% 
      rename(Latitude1 = Latitude.x, Latitude2 = Latitude.y, Longitude1 = Longitude.x, Longitude2 = Longitude.y, Island1 = Island.x, Island2 = Island.y) %>% 
      mutate(GEO = pmap(list(a = Longitude1, 
                             b = Latitude1,
                             x = Longitude2,
                             y = Latitude2), 
                        ~ geosphere::distGeo( c(..1, ..2), c(..3, ..4))))
    CHA_distances$Sample2 <- as.character(CHA_distances$Sample2) 
    CHA_distances_all <- CHA_distances %>% 
      mutate(x_y = purrr::map2_chr(Sample1, Sample2, c_sort_collapse)) %>% 
      distinct(x_y, .keep_all = TRUE) %>% 
      dplyr::select(-x_y) %>% 
      filter(Sample1 != Sample2)
    Distances_all <- CHA_distances_all %>% pivot_longer(cols = c(BC, UF), names_to = "Distance", values_to = "Value")
    Distances_all$GEO <- as.numeric(Distances_all$GEO) 
    Distances_all <- Distances_all %>% mutate(GEO = GEO/1000) %>% data.frame()
    Distances_all$GEO <- as.integer(Distances_all$GEO)
    
    DisBC <- Distances_all %>% filter(Distance=="BC") %>% data.frame()
    DisUF <- Distances_all %>% filter(Distance=="UF") %>% data.frame()
    
    modelBC <-lm(GEO ~ Value, DisBC)
    summary(modelBC)
    modelUF <-lm(GEO ~ Value, DisUF)
    summary(modelUF)
    
    lm_eqn = function(m) {
      a <- coef(m)[2]
      b <- coef(m)[1]
      r2 <- summary(m)$r.squared
      pv <- summary(m)$coefficients[2,4]
      
      a <- as.numeric(a)
      b <- as.numeric(b)
      r2 <- as.numeric(r2)
      pv <- as.numeric(pv)
      
      a <- format(a, digits = 5)
      b <- format(b, digits = 4)
      r2 <- format(r2, digits = 1)
      pv <- format(pv, digits = 1)
      
      if (b < 0) {
        eq <- substitute(italic(y) == a ~ italic(x) ~ b ~ "," ~ italic(r)^2 ~ "=" ~ r2 ~ "," ~ italic(p) ~ "=" ~ pv,
                         list(a = a, b = b, r2 = r2, pv = pv))
      } else {
        eq <- substitute(italic(y) == a ~ italic(x) + b ~ "," ~ italic(r)^2 ~ "=" ~ r2 ~ "," ~ italic(p) ~ "=" ~ pv,
                         list(a = a, b = b, r2 = r2, pv = pv))
      }
      as.character(as.expression(eq))
    }
    
    
    p <- Distances_all %>% 
      ggplot(aes(x = GEO, y = Value, color = Distance)) + 
      geom_point() + 
      theme_classic() + 
      scale_color_manual(values = c("orange", "darkblue"), labels = c("Bray-Curtis", "wUniFrac")) + 
      xlab("Geographic distance (km)") + 
      ylab("Community dissimilarity") + 
      labs(color = "Metric") + 
      geom_smooth(method = 'lm', formula = y ~ x, se = F, size = 0.75) + 
      xlim(0, 2000) +
      scale_y_continuous(limits = c(0, 1.05), breaks = c(0.00, 0.25, 0.50, 0.75, 1.00)) +
      annotate("text", x = 1200, y = 1.05, label = lm_eqn(modelBC), color = "black", parse = TRUE, size = 5) +
      annotate("text", x = 1200, y = 0.22, label = lm_eqn(modelUF), color = "black", parse = TRUE, size = 5) +
      theme(
        legend.title = element_text(size = 14, face = "bold"), 
        legend.text = element_text(size = 14), 
        axis.title.x = element_text(size = 16, face = "bold"), 
        axis.title.y = element_text(size = 16, face = "bold"), 
        axis.text = element_text(size = 12), 
        legend.position = "right", 
        legend.justification = "left", 
        strip.text = element_text(size = 14, face = "bold"), 
        strip.background = element_blank(), 
        plot.caption = element_text(size = 12, hjust = 1, vjust = 25), 
        axis.line = element_line())
    ggsave('plots/distance_decays.pdf', plot = p, width = 8, height = 8)
    }
##########################################################
####### 7. Main function      ############################ 
##########################################################
main <- function(){
  dir.create('stats')
  dir.create('figures')
  dir.create('data')
  
  # Data loading and preprocessing
  microbial_data = load_microbial_data()
  asv_table = microbial_data$asv
  tax_table = microbial_data$tax
  phy_tree = microbial_data$tree
  
  env_data = load_environmental_data()
  
  phyloseq_data = datasets_to_phyloseq(asv_table, tax_table, phy_tree, env_data)
  env_data_with_alpha = add_alpha_to_env(phyloseq_data)
  write.csv(env_data_with_alpha, file='data/env_with_alpha.csv', quote = F, row.names = F, sep = ',')
    
  # Alpha diversity analyses --> GAMs (4x3, Observed, Shannon, PD)
  alpha_div_gams(env_data_with_alpha)
    
  # Taxonomy & Tree --> 1 figure
  plot_phylot_tree(phyloseq_data)
  
  # WGCNA --> SI plots + 1 figure
  wgnca_analysis(phyloseq_data)
    
  # Beta-diversity
  beta_diversity_analyses(phyloseq_data, env_data_with_alpha)
  }
main()

