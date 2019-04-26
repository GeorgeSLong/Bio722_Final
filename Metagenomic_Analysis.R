# Functions
Sample_Translation <- function(colname){ # Converts sample names into their treatments
	colname <- sapply(colname, function(x){ifelse(grepl("157519|15752[0-1]",x), "NK", x)})
	colname <- sapply(colname, function(x){ifelse(grepl("15752[2-4]",x), "PK", x)})
	colname <- sapply(colname, function(x){ifelse(grepl("15752[5-7]",x), "NP", x)})
	colname <- sapply(colname, function(x){ifelse(grepl("15752[8-9]|157530",x), "NPK", x)})
	colname <- sapply(colname, function(x){ifelse(grepl("15753[1-3]",x), "Low.NPK", x)})
	colname <- sapply(colname, function(x){ifelse(grepl("15753[4-6]",x), "High.NPK", x)})
	colname <- sapply(colname, function(x){ifelse(grepl("15753[7-9]",x), "Very.High.NPK", x)})
	colname <- sapply(colname, function(x){ifelse(grepl("15754[0-2]",x), "Compost", x)})
	colname <- sapply(colname, function(x){ifelse(grepl("15754[3-5]",x), "Compost.NK", x)})
	colname <- sapply(colname, function(x){ifelse(grepl("15754[6-8]",x), "Compost.PK", x)})
	colname <- sapply(colname, function(x){ifelse(grepl("157549|15755[0-1]",x), "Compost.NP", x)})
	colname <- sapply(colname, function(x){ifelse(grepl("15755[2-4]",x), "High.Compost.NPK", x)})
	return(colname)
}

plot.labels <- function(TUKEY,variable,threshold=0.05,reversed=FALSE){ # Determines Tukey's HSDs for a specific variable
      pval<-TUKEY[[variable]][,4]
      names(pval)<-row.names(TUKEY[[variable]])
      labels<-data.frame(multcompLetters(pval,threshold=threshold,reversed=FALSE)['Letters'])
      labels<-data.frame(labels,row.names(labels))
      names(labels)<-c("Letters",variable)
      return(labels)
}


############################################################################################################3

# Libraries
library(phyloseq)
library(plyr)
library(stringi)
library(dplyr)
library(ggplot2)
library(gridExtra)
library(multcompView)
library(reshape2)
library(tidyr)
library(vegan)
library(MASS)
library(igraph)
library(eulerr)
#Changing WD
setwd("/2/scratch/sam/BIO722/Project/kraken_MPA")

# First, need to know what files we're looking at
report_files <- list.files(pattern = ".*mpa")

# Getting the ID numbers - Used to ease IDing the files of interest
num <- sapply(report_files, function(x){strsplit(x, "\\.")[[1]][1]})
num <- sapply(num, function(x){strsplit(x, "_")[[1]][1]})

for(i in 1:length(report_files)){
  kraken <- read.delim(report_files[i], header = FALSE, comment.char = "#")
  colnames(kraken) <- c("Taxa", "Reads")
  
  # Filtering for requested level ID only
  kraken_Species <- kraken[grepl(x = kraken$Taxa,"g__"),] 
  # Next, we need to merge the taxa which have s__ included.
  kraken_Species$Taxa <- gsub("\\|s__.*" ,"" ,x = kraken_Species$Taxa)
  kraken_Species <- kraken_Species %>% group_by(Taxa) %>% summarize("Reads" = sum(Reads))
  
  # Fixing the names
  taxa_tmp <- kraken_Species$Taxa
  kraken_Species$Taxa <- sapply(kraken_Species$Taxa, function(x){strsplit(x, "g__")[[1]][2]})

  # Creating the overall dataframe
  if(i == 1){
    kraken_grouped <- kraken_Species
    kraken_taxa <- taxa_tmp
  }else{
    kraken_grouped <- full_join(x = kraken_grouped, kraken_Species, by = "Taxa")
    kraken_taxa <- unique(c(kraken_taxa, taxa_tmp))
  }
}

# Cleaning up the results here
rm(kraken_Species, kraken)
colnames(kraken_grouped)[2:dim(kraken_grouped)[2]] <- num #Making the Columns meaningful 
kraken_grouped[,-1] <- sapply(kraken_grouped[,-1], function(x){ifelse(is.na(x), 0, x)})
kraken_grouped  <- as.data.frame(kraken_grouped)
kraken_unfiltered <- kraken_grouped

###########################
### Filtering our Reads ###
###########################

# Next, we need to remove the reads that account for less than 1% of the total. This will account for most if not all false positives.
for(i in 2:dim(kraken_grouped)[2]){ # Skipping the names column
	tmp <- kraken_grouped[,i]
 	threshold <- ceiling(sum(tmp) * 0.01)
	tmp2 <- sapply(tmp, function(x){ifelse(x < threshold, 0, x)})
	kraken_grouped[,i] <- tmp2

}
kraken_grouped <- kraken_grouped[ifelse(rowSums(kraken_grouped[,-1]) > 0, TRUE, FALSE),] # Removing Genera that now add to 0
kraken_prop <- kraken_grouped

# Adding in the Proportions
for(i in 2:dim(kraken_prop)[2]){ # Skipping the names column
	tmp <- kraken_prop[,i]
 	tmp2 <- sapply(tmp, function(x){x/sum(tmp)})
	kraken_prop[,i] <- tmp2

}

# Transforming the results using arcsinh
kraken_asinh <- kraken_grouped
kraken_asinh[,-1] <- sapply(kraken_asinh[,-1], asinh)

# Making the OTU tables for all four results
rownames(kraken_grouped)  <- rownames(kraken_asinh) <- rownames(kraken_prop) <-  kraken_grouped$Taxa
rownames(kraken_unfiltered) <- kraken_unfiltered$Taxa
otu_prop <- otu_table(kraken_prop[,-1], TRUE)
otu_unfiltered <- otu_table(kraken_unfiltered[,-1], TRUE)
otu_grouped <- otu_table(kraken_grouped[,-1], TRUE)
otu_asinh <- otu_table(kraken_asinh[,-1], TRUE)

##########################
### Making a Taxa Table###
##########################
# Filtering the taxa to only include what was filtered out before
kraken_taxa <- kraken_taxa[sapply(kraken_taxa, function(x){strsplit(x, "g__")[[1]][2]}) %in% rownames(kraken_unfiltered)]

kraken_taxa <- strsplit(kraken_taxa, "\\|") # Making a list of the taxons

taxa_levels <- data.frame("Letter" = c("d", "k", "p", "c", "o", "f", "g"), "Level" = c("Domain", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus"))# Dictionary of taxonomy levels

kraken_taxa  <- lapply(kraken_taxa, function(x){ # Assigning the taxonomy level based its name
	letter <-  sapply(x, function(y){strsplit(y, "_")[[1]][1]}) 
	level <- taxa_levels[taxa_levels$Letter %in% letter, "Level"]
	names(x) <- level
	return(x)

})

# This section here converts our list into a dataframe, and then it'll reorder the elements in a row to make sure that the appropriate taxa is in each column.  This is because the NCBI taxonomy can be incomplete depending on the genus of interest......
kraken_taxa_df <- data.frame(t(stri_list2matrix(kraken_taxa)), stringsAsFactors = FALSE) # Converting to a dataframe
colnames(kraken_taxa_df) <- taxa_levels$Level
for(j in 1:dim(kraken_taxa_df)[2]){
	right_taxa <- grepl(paste0(taxa_levels[taxa_levels$Level %in% colnames(kraken_taxa_df)[j],"Letter"],"__"), kraken_taxa_df[,j]) # Which ones are right?
	right_taxa <- right_taxa | is.na(kraken_taxa_df[,j]) # Not much I can do about missing taxa.  Can only change those that are wrong.
	if(sum(right_taxa) != length(right_taxa)){ # If there are misassigned Taxa
		for(i in 1:length(right_taxa)){
			if(!right_taxa[i]){ # Only interested in running those that are wrong
			 	misclassified_taxa <- taxa_levels[taxa_levels$Letter == strsplit(kraken_taxa_df[i,j],"_")[[1]][1],2]
				colindex <- !is.na(kraken_taxa_df[i,-c(1:j-1)]) # Getting the length of the non NAs
				colindex <- sum(colindex) - 1 # How many?
				kraken_taxa_df[i,(dim(kraken_taxa_df)[2] - colindex):dim(kraken_taxa_df)[2]] <- kraken_taxa_df[i,j:(j + colindex)] # Assigning the taxa in their new levels
				kraken_taxa_df[i,j:(dim(kraken_taxa_df)[2] - colindex - 1)] <- NA
			}	

		}
	}
}

kraken_taxa_df <- sapply(kraken_taxa_df, function(x){gsub(".__","",x)}) # Removing the taxon indicators.  No longer needed
rownames(kraken_taxa_df) <- kraken_taxa_df[,dim(kraken_taxa_df)[2]]
kraken_tax <- tax_table(kraken_taxa_df) # Making the taxa table for phyloseq

#####################################
### Making the Experimental Frame ###
#####################################
Treatment <- Sample_Translation(colnames(kraken_grouped)[-1])

# Now, we need to separate EVERYTHING in the treatments

Treatments  <- data.frame("Treatment" = Treatment,
		"Nitrogen" = grepl("N", Treatment),
		"Phosphorus" = grepl("P", Treatment),
		"Potassium" = grepl("K", Treatment),
		"Low.Fertilizer" = grepl("Low", Treatment),
		"High.Fertilizer" = grepl("High(?!.Compost)", Treatment, perl = TRUE),
		"Very High Fertilizer" = grepl("Very.High", Treatment),
		"Compost" = grepl("Compost", Treatment),
		"High.Compost" = grepl("High.Compost", Treatment))
Treatments <- sample_data(Treatments)
####################################
### Creating the PhyloSeq Object ###
####################################

physeq_unfiltered <- phyloseq(otu_unfiltered, kraken_tax, Treatments)
physeq_grouped <- phyloseq(otu_grouped, kraken_tax, Treatments)
physeq_asinh <- phyloseq(otu_asinh, kraken_tax, Treatments)
physeq_prop <- phyloseq(otu_prop, kraken_tax, Treatments)

############################
### Analyzing the Results###
############################
# Did we get the same Phyla as the original Announcement?
merged_physeq <- merge_samples(physeq_asinh, "Treatment") # 
merged_physeq@otu_table <- merged_physeq@otu_table / 3 # Averaging the pooled results

for(i in colnames(kraken_taxa_df)){
	if(i == "Kingdom") next # Only the unfiltered reads have Eukaryotes in them
	Glommed <- tax_glom(merged_physeq, i) # Removing the useless bars for each taxa
	
	# Abundance Bar Graph
	pdf(paste0("../Graphs/Asinh_Abundances_",i,".pdf"), width = 9, height = 6)
	print(plot_bar(Glommed, fill = i) + theme(axis.text.x = element_text(angle = 45, hjust = 1)))
	dev.off()
	
	# Heatmap
	pdf(paste0("../Graphs/Asinh_Heatmap_Not_NULL",i,".pdf"), width = 9, height = 6)
	print(plot_heatmap(Glommed, method = "PCoA") + ylab(paste(i)))# Method set to Null to ensure that the results between asinh and prop are directly comparable.  Otherwise, the organisation of the two plots vary wildly
	dev.off()
}

# PCoA
# What distance do we want to use?
Distance_Testing <- list()
distances_tested <- c(distanceMethodList[[3]], distanceMethodList[[4]], distanceMethodList[[5]])
count = 1
for(i in distances_tested){
	Distance_Testing[[count]] <- plot_ordination(physeq_prop, ordinate(physeq_prop, "PCoA", i), type = "samples", color = "Treatment") + geom_point(size = 1) + theme(legend.position="none") + ggtitle(i)
	count <- count + 1
}
pdf("../Graphs/Distance_Testing_Prop.pdf", height = 40, width = 24)
	grid.arrange(grobs = Distance_Testing, ncol = 4)
dev.off()

# Actually making the plot
GP.ord <- ordinate(physeq_asinh, "PCoA", "bray")
for(i in colnames(physeq_prop@sam_data)){
	pdf(paste0("../Graphs/PCoA_",i,"_prop.pdf"), width = 6, height = 4)
		print(plot_ordination(physeq_prop, GP.ord, type = "samples", color = i) + geom_point(size = 4))
	dev.off()
}

####################
# Network Analysis #
####################
ig <- make_network(physeq_asinh, max.dist = 0.25, dist.fun="bray") #
for(i in colnames(physeq_prop@sam_data)){
	pdf(paste0("../Graphs/Network_",i,"_prop.pdf"), width = 6, height = 4)
	set.seed(1996)
	print(plot_network(ig, physeq_prop, color = i, label = NULL))
	dev.off()
}

# Making the communities via a walktrap and comparing their modularities
wc <- walktrap.community(ig)
nit  <- make_clusters(ig, membership = Treatments$Nitrogen + 1)
com  <- make_clusters(ig, membership = Treatments$Compost+ 1)
pot <- make_clusters(ig, membership = Treatments$Potassium + 1)
phos <- make_clusters(ig, membership = Treatments$Phosphorus + 1)

set.seed(1996)
layout <- layout.fruchterman.reingold(ig)

pdf("../Graphs/WalkTrap_asinh.pdf", width = 6, height = 6)
plot(com,ig, layout = layout, vertex.label = Treatments$Treatment, vertex.label.color = "black",  vertex.size = 10, edge.arrow.size = .2)
dev.off()

# Permutation Test - Compost
perm_com  <- replicate(n = 10000,{ # Making 10000 iterations of the same cluster sizes
	 mem <- sample(x = num, size = 21, rep = FALSE)
	 mem_num <- sapply(colnames(physeq_prop@otu_table), function(x){ifelse(x %in% mem, 1,2)})
	return(modularity(make_clusters(ig, membership = mem_num)))
}
)
# Graphing the Results
pdf("../Graphs/Permutation_Compost_Prop.pdf", width = 6, height = 4)
	hist(perm_com,xlim = c(-0.4, 0.4), xlab = "Modularity", main = NULL)
	abline(v = modularity(com), col = "red")
dev.off()

# P-Value of Permutation test
(sum(perm_com > modularity(com)) + 1)/ (length(perm_com) + 1)

# Permutation Test - Nitrogen
perm_nit  <- replicate(n = 10000,{
	 mem <- sample(x = num, size = 9, rep = FALSE)
	 mem_num <- sapply(colnames(physeq_prop@otu_table), function(x){ifelse(x %in% mem, 1,2)})
	return(modularity(make_clusters(ig, membership = mem_num)))
}
)

# Graphing the Results
pdf("../Graphs/Permutation_Nitrogen_Prop.pdf", width = 6, height = 4)
	hist(perm_nit,xlim = c(-0.4, 0.4), xlab = "Modularity", main = NULL)
	abline(v = modularity(nit), col = "red")
dev.off()

# P-Value of Permutation test
(sum(perm_nit > modularity(nit)) + 1)/ (length(perm_nit) + 1)

###################
# Alpha Diversity #
###################
# alpha_div  <- estimate_richness(physeq_prop, measures = c("Simpson", "Shannon"))
alpha_div  <- estimate_richness(transform_sample_counts(physeq_unfiltered, function(OTU) OTU/sum(OTU)), measures = c("Simpson", "Shannon")) #  Unfiltered Proportional
# alpha_div  <- estimate_richness(transform_sample_counts(physeq_unfiltered, function(OTU) asinh(OTU)), measures = c("Simpson", "Shannon")) #  Unfiltered asinh
alpha_div$Treatment <- Sample_Translation(rownames(alpha_div))

# Find Tukey's Groups - Shannon
shannon <- lm(data = alpha_div, formula = Shannon ~ Treatment)
shannon <- plot.labels(TukeyHSD(aov(shannon)), "Treatment")
shannon_df  <- left_join(shannon, aggregate(Shannon ~ Treatment, alpha_div, mean), by = "Treatment")

shannon_plot <- ggplot(alpha_div, aes(y = Shannon, x = Treatment)) + geom_boxplot() + geom_text(data = shannon_df, aes(label = Letters)) + ggtitle("Shannon Diversity - Proportional Transformation") + theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Find Tukey's Groups - Simpson
simpson <- lm(data = alpha_div, formula = Simpson ~ Treatment)
simpson <- plot.labels(TukeyHSD(aov(simpson)), "Treatment")
simpson_df  <- left_join(simpson, aggregate(Simpson ~ Treatment, alpha_div, mean), by = "Treatment")

simpson_plot <- ggplot(alpha_div, aes(y = Simpson, x = Treatment)) + geom_boxplot() + geom_text(data = simpson_df, aes(label = Letters))+ ggtitle("Simpson Diversity - Proportional Transformation")+ theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Let's plot these results
pdf("../Graphs/Diversity_Analysis_Prop_Unfiltered.pdf", width = 6, height = 8)
	grid.arrange(shannon_plot, simpson_plot)
dev.off()

###########
### LDA ###
###########
kraken_lda  <- data.frame(t(otu_prop)) # Reformatting
kraken_lda$Treatment <- Treatment  # Adding the Treatment
lda_results  <- lda(formula = Treatment ~., data = kraken_lda) # Making the lda

prop = lda_results$svd^2/sum(lda_results$svd^2)

plda <- predict(object = lda_results, newdata = kraken_lda)

lda_plotting  <- data.frame("Treatment" = kraken_lda$Treatment, lda = plda$x)
pdf("../Graphs/LDA_prop.pdf", width = 6, height = 4)
	ggplot(lda_plotting, aes(x = lda.LD1, y = lda.LD2, colour = Treatment)) + geom_point(size = 2.5) + ylab(paste0("LDA Dim 2 (",round(prop[2],4)*100,"%)")) + xlab(paste0("LDA Dim 1(", round(prop[1],4)*100,"%)")) # 
dev.off()

##############################################
### Eulerr Diagrams - Supplemental Figures ###
##############################################

# Nitrogen
physeq_nit <- merge_samples(physeq_asinh, "Nitrogen")
nit_euler <- list()
nit_euler$`Nitrogen Absent`  <- colnames(physeq_nit@otu_table)[physeq_nit@otu_table["FALSE",] > 0]
nit_euler$`Nitrogen Present`<- colnames(physeq_nit@otu_table)[physeq_nit@otu_table["TRUE",] > 0]
pdf("../Graphs/Nitrogen_Eulerr.pdf", width = 6, height = 4)
plot(euler(nit_euler))
dev.off()

base:::setdiff(nit_euler[[1]], nit_euler[[2]])
base:::setdiff(nit_euler[[2]], nit_euler[[1]])

#Compost
physeq_com <- merge_samples(physeq_asinh, "Compost")
com_euler <- list()
com_euler$`Compost Absent` <- colnames(physeq_com@otu_table)[physeq_com@otu_table["FALSE",] > 0]
com_euler$`Compost Present` <- colnames(physeq_com@otu_table)[physeq_com@otu_table["TRUE",] > 0]
pdf("../Graphs/Compost_Eulerr.pdf", width = 6, height = 4)
plot(euler(com_euler))
dev.off()

base:::setdiff(com_euler[[1]], com_euler[[2]])
base:::setdiff(com_euler[[2]], com_euler[[1]])
