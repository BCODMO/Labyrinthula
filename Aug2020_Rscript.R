## Laby microbe data
## August 2020
## Revision of Laby microbes

## Table of contents
# i)  load libraries/prep workspace
## I) denoising and processing ASVs.
## II) creating and modifying map file as needed.
rm(list=ls())
# i)  load libraries/prep workspace
library(phyloseq)
library(ggplot2)
library(ggthemes)
library(DESeq)
library(DESeq2)
library(car)
library(compositions)
library(coda.base)
library(dplyr)
library(vegan)
library(dada2)

theme_set(theme_bw())

## I) denoising and processing ASVs
#####
getwd()#"/Users/melissakardish/Google Drive/Davis/Research/Laby"
setwd("July2019")
system("mkdir primer_trimmed_fastqs")
system("parallel --link --jobs 4 \
  'cutadapt \
  --pair-filter any \
  --no-indels \
  --discard-untrimmed \
  -g GTGYCAGCMGCCGCGGTAA \
  -G CCGYCAATTYMTTTRAGTTT \
  -o primer_trimmed_fastqs/{1/} \
  -p primer_trimmed_fastqs/{2/} \
  {1} {2} \
  > primer_trimmed_fastqs/{1/}_cutadapt_log.txt' \
  ::: ../sequences/Kardish/*_R1_*.fastq.gz ::: ../sequences/Kardish/*_R2_*.fastq.gz ")
 
system("mkdir primerlogs")
system("mv primer_trimmed_fastqs/*log.txt primerlogs/")
path1 <- "primer_trimmed_fastqs" # CHANGE ME to the directory containing the fastq files after unzipping.
fnFs1 <- sort(list.files(path1, pattern="_R1_001.fastq.gz", full.names = TRUE))
fnRs1 <- sort(list.files(path1, pattern="_R2_001.fastq.gz", full.names = TRUE))

ii <- sample(length(fnFs1), 3)
for(i in ii) { print(plotQualityProfile(fnFs1[i]) + ggtitle("Fwd")) }
for(i in ii) { print(plotQualityProfile(fnRs1[i]) + ggtitle("Rev")) }

sample.names1 <- sapply(strsplit(basename(fnFs1), "_1"), `[`, 1)

filt_path1 <- file.path(path1, "filtered") # Place filtered files in filtered/ subdirectory


filtFs1 <- file.path(filt_path1, paste0(sample.names1, "_F_filt.fastq.gz"))
filtRs1 <- file.path(filt_path1, paste0(sample.names1, "_R_filt.fastq.gz"))



out1 <- filterAndTrim(fnFs1, filtFs1, fnRs1, filtRs1, truncLen=c(280,180),
                      maxN=1, maxEE=c(3,5), truncQ=2, rm.phix=TRUE,
                      compress=TRUE, multithread=TRUE) # On Windows set multithread=FALSE
save(out1,file="dada2progressfiles/out1.RData")
set.seed(5798)
errF1 <- learnErrors(filtFs1, multithread=TRUE)
set.seed(5798)
errR1 <- learnErrors(filtRs1, multithread=TRUE,MAX_CONSIST = 100)
save(errF1,errR1,file="dada2progressfiles/errs.RData")

derepFs1 <- derepFastq(filtFs1, verbose=TRUE)
derepRs1 <- derepFastq(filtRs1, verbose=TRUE)
names(derepFs1) <- sample.names1
names(derepRs1) <- sample.names1

save(derepFs1,derepRs1,file="dada2progressfiles/dereps.RData")


dadaFs1 <- dada(derepFs1, err=errF1, multithread=TRUE)
dadaRs1 <- dada(derepRs1, err=errR1, multithread=TRUE)

save(dadaFs1,dadaRs1,file="dada2progressfiles/dadas.RData")


mergers1 <- mergePairs(dadaFs1, derepFs1, dadaRs1, derepRs1, verbose=TRUE)
#this does kinda reasonable actually with total reads, but not with rarer things? or are these errors...
save(mergers1,file="dada2progressfiles/mergers.RData")


seqtabF1 <- makeSequenceTable(dadaFs1)
seqtab.nochim_F <- removeBimeraDenovo(seqtabF1, method="consensus", multithread=TRUE, verbose=TRUE)
seqtabmerg <- makeSequenceTable(mergers1)
seqtab.nochim_m <- removeBimeraDenovo(seqtabmerg, method="consensus", multithread=TRUE, verbose=TRUE)

  #really highchimera rates. watch out for changes when redo analysis, use forwards only because of high chimera rate

save(seqtabF1,seqtab.nochim_F,seqtabmerg,seqtab.nochim_m,file="dada2progressfiles/seqtabs.RData")

##run overnight? taking a long time
taxa_gg2F <- assignTaxonomy(seqtab.nochim_F, "~/Downloads/gg_13_8_train_set_97.fa.gz", multithread=TRUE)
taxa_silva2F <- assignTaxonomy(seqtab.nochim_F, "~/Downloads/silva_nr_v128_train_set.fa.gz", multithread=TRUE)
species_silva2F <- addSpecies(taxa_silva2F, "~/Downloads/silva_species_assignment_v128.fa.gz")
taxa_rdp2F <- assignTaxonomy(seqtab.nochim_F, "~/Downloads/rdp_train_set_16.fa.gz", multithread=FALSE)
species_rdp2F <- addSpecies(taxa_rdp2F, "~/Downloads/rdp_species_assignment_16.fa.gz")

##run overnight? taking a long time
taxa_gg2m <- assignTaxonomy(seqtab.nochim_m, "~/Downloads/gg_13_8_train_set_97.fa.gz", multithread=TRUE)
taxa_silva2m <- assignTaxonomy(seqtab.nochim_m, "~/Downloads/silva_nr_v128_train_set.fa.gz", multithread=TRUE)
species_silva2m <- addSpecies(taxa_silva2m, "~/Downloads/silva_species_assignment_v128.fa.gz")
taxa_rdp2m <- assignTaxonomy(seqtab.nochim_m, "~/Downloads/rdp_train_set_16.fa.gz", multithread=FALSE)
species_rdp2m <- addSpecies(taxa_rdp2m, "~/Downloads/rdp_species_assignment_16.fa.gz")

save(taxa_gg2F,taxa_silva2F,species_silva2F,taxa_rdp2F,species_rdp2F,file="dada2progressfiles/taxa.RData")
save(taxa_gg2m,taxa_silva2m,species_silva2m,taxa_rdp2m,species_rdp2m,file="dada2progressfiles/taxam.RData")
#####

## II) creating and modifying map file as needed.
#####
#load files 
load(file="dada2progressfiles/taxa.RData")
load("dada2progressfiles/seqtabs.RData")


##make "MAP" files
labynames<-rownames(seqtabF1)
labydata<-data.frame("sample_names"=rownames(seqtabF1),
                     "genotype"=gsub("-.*$", "", labynames),
                     "tank.pot"=gsub("^[^-]*-([^_]+).*", "\\1",labynames),
                     "tank"=gsub("-.*$","",gsub("^[^-]*-([^_]+).*", "\\1",labynames)),
                     "pot"=substr(gsub("^[^-]*-([^_]+).*", "\\1",labynames),nchar(gsub("^[^-]*-([^_]+).*", "\\1",labynames)),
                                  nchar(gsub("^[^-]*-([^_]+).*", "\\1",labynames))))
labydata$treatment<-as.factor(ifelse(labydata$tank%in%c(1,8,10,12,13),"hot","cold"))
labydata$all<-tolower(paste(labydata$genotype,labydata$tank.pot,sep=""))
labycsv<-read.csv("~/Google Drive/Davis/Research/Laby/labymap.csv",stringsAsFactors = F)
labycsv<-read.csv("~/Google Drive/Davis/Research/Laby/labyupdates_20190325.csv",stringsAsFactors = F) #spring updates
labycsv<-read.csv("~/Google Drive/Davis/Research/Laby/labyupdates_20190325_withconc.csv",stringsAsFactors = F) #spring updates
labycsv<-labycsv[,5:12] # updated 
labydata<-left_join(labydata,labycsv,by="all")
rownames(labydata)<-rownames(seqtabF1)

#####

## III) creating phyloseq object to work with
#####
ps <- phyloseq(otu_table(seqtabF1, taxa_are_rows=FALSE), 
               sample_data(labydata),
               tax_table(species_silva2F))

ps2 <-subset_taxa(ps,
                  Kingdom == "Bacteria" &
                    Family  != "Mitochondria" &
                    Class   != "Chloroplast")

laby<-prune_samples(sample_sums(ps2)>=1000,ps2)

sample_data(laby)$geno_suscep<-ifelse(sample_data(laby)$genotype%in%c("Br", "R", "O","W", "Y"),"high",
                                      ifelse(sample_data(laby)$genotype%in%c("Bl", "P"),"low",NA))
sample_data(laby)$combo<-paste(sample_data(laby)$geno_suscep,sample_data(laby)$treatment,sep="_")

#total seqs
sum(otu_table(laby))
median(rowSums(otu_table(laby)))
min(rowSums(otu_table(laby)))
#####

## IV) alpha diversity
#####
laby2<-laby
# all data
sample_data(laby2)$observed_all<-rowSums(otu_table(laby2)>0)
a<-data.frame(sample_data(laby2))
ggplot(a, aes(x=treatment,y=observed_all,colour=treatment))+
  geom_violin(aes(fill=treatment))+
  scale_color_manual(values=c("black","black"))+
  ggtitle("Observed ASVs")+scale_fill_manual(values=c("white","grey5")) + 
  theme_tufte() +xlab("Treatment") + ylab("Number observed ASVs")
t.test(a[a$treatment=="hot","observed_all"],a[a$treatment=="cold","observed_all"])
## they are different by t-test

ggplot(a, aes(x=geno_suscep,y=observed_all,colour=geno_suscep))+
  geom_violin(aes(fill=geno_suscep))+
  scale_color_manual(values=c("black","black"))+
  ggtitle("Observed ASVs")+scale_fill_manual(values=c("white","grey5")) + 
  theme_tufte() +xlab("Laby genotype susceptibility") + ylab("Number observed ASVs")
t.test(a[a$geno_suscep=="high","observed_all"],a[a$geno_suscep=="low","observed_all"])

# rarefied for richness comparisons.
plot(rowSums(otu_table(laby2)),rowSums(otu_table(laby2)>0),col=sample_data(laby2)$treatment)
summary(lm(rowSums(otu_table(laby2)>0)~rowSums(otu_table(laby2)))) 
# there is a positive relationship between sequencing depth and number of ASVs (and some of hot treatment had more reads sequenced)
# we'll rareify to test this to a more even depth
set.seed(3547)
torare<-sample(1:1000000000,200)
rareset<-lapply(torare,function(x) rarefy_even_depth(laby2,sample.size=12149,rngseed=torare,replace=F))
rareset_5000<-lapply(torare,function(x) rarefy_even_depth(laby2,sample.size=5000,rngseed=torare,replace=F))

obsrich<-sapply(rareset,function(x) rowSums(otu_table(x)>0))
a$mean_richness_rarefied<-apply(obsrich,1,mean)
obsrich5000<-sapply(rareset_5000,function(x) rowSums(otu_table(x)>0))
a$mean_richness_rarefied_5000<-apply(obsrich5000,1,mean)

chao_5000<-sapply(rareset_5000,function(x) estimate_richness(x,measures="Chao1")[,1])
a$chao_5000<-apply(chao_5000,1,mean)

a$geno_suscep<-ifelse(sample_data(laby)$genotype%in%c("Br", "R","O", "W", "Y"),"high",
                                      ifelse(sample_data(laby)$genotype%in%c("Bl", "P"),"low",NA))
a$combo<-paste(sample_data(laby)$geno_suscep,sample_data(laby)$treatment,sep="_")


# rareified data
fig_m_s1a<-ggplot(a, aes(x=treatment,y=mean_richness_rarefied,colour=treatment))+
  geom_boxplot(aes(fill=treatment))+
  scale_color_manual(values=c("black","black"))+
  ggtitle("Observed ASVs")+scale_fill_manual(values=c("white","grey25")) + 
  xlab("Treatment") + ylab("Number observed ASVs in rarefied data") + theme_bw()
t.test(a[a$treatment=="hot","mean_richness_rarefied"],a[a$treatment=="cold","mean_richness_rarefied"])
## they are different by t-test (still) p = 0.0005603
summ_data<-summarySE(a[complete.cases(a$geno_suscep),],measurevar="mean_richness_rarefied",groupvars = c("treatment"))
fig_m_s1a_i<-ggplot(data=summ_data,aes(y=mean_richness_rarefied,x=treatment,fill=treatment))+
  geom_point(shape=21,cex=2,aes(y=mean_richness_rarefied,x=treatment))+
  ggtitle("Observed ASVs")+scale_fill_manual(values=c("white","grey25")) + 
  xlab("Treatment") + ylab("Number observed ASVs in rarefied data")+
  geom_errorbar(aes(ymin=mean_richness_rarefied-se,ymax=mean_richness_rarefied+se),width=0.1)



fig_m_s1b<-ggplot(a, aes(x=geno_suscep,y=mean_richness_rarefied,colour=geno_suscep))+
  geom_boxplot(aes(fill=geno_suscep))+
  scale_color_manual(values=c("black","black"))+
  ggtitle("Observed ASVs")+scale_fill_manual(values=c("white","grey25")) + 
  theme_tufte() +xlab("Laby genotype susceptibility") + ylab("Number observed ASVs in rarefied data")
t.test(a[a$geno_suscep=="high","mean_richness_rarefied"],a[a$geno_suscep=="low","mean_richness_rarefied"]) # 0.08938


fig_m_s1c<-ggplot(a[complete.cases(a$geno_suscep),], aes(x=geno_suscep,y=mean_richness_rarefied,colour=treatment))+
  geom_boxplot(aes(fill=treatment))+
  scale_color_manual(values=c("black","black"))+
  ggtitle("Observed ASVs")+scale_fill_manual(values=c("white","grey15")) + 
  xlab("Genotype susceptibility") + ylab("Number observed ASVs in rarefied data")
loggedrich<-log(a$mean_richness_rarefied)
mod<-aov(mean_richness_rarefied~treatment*geno_suscep*Laby,data=a)
shapiro.test(mod$residuals) # fails 0.004, mod of unlogged more so still
leveneTest(mod) # passes p = 0.29
anova(mod) # only treatment significant

summ_data<-summarySE(a[complete.cases(a$geno_suscep),],measurevar="mean_richness_rarefied",groupvars = c("treatment","geno_suscep"))
summ_data$combo <- paste(summ_data$treatment, summ_data$geno_suscep)
summ_data$combo <- factor(summ_data$combo, levels = c("cold low", "hot low","cold high","hot high"))
fig_m_s1c_i<-ggplot(data=summ_data,aes(y=mean_richness_rarefied,x=combo,fill=treatment))+
  geom_point(shape=21,cex=2,aes(y=mean_richness_rarefied,x=combo))+
  ggtitle("Observed ASVs")+scale_fill_manual(values=c("white","grey15")) + 
  xlab("Genotype susceptibility") + ylab("Number observed ASVs in rarefied data")+
  geom_errorbar(aes(ymin=mean_richness_rarefied-se,ymax=mean_richness_rarefied+se),width=0.1)



# rareified data to 5000
ggplot(a, aes(x=treatment,y=mean_richness_rarefied_5000,colour=treatment))+
  geom_violin(aes(fill=treatment))+
  scale_color_manual(values=c("black","black"))+
  ggtitle("Observed ASVs")+scale_fill_manual(values=c("white","grey5")) + 
  theme_tufte() +xlab("Treatment") + ylab("Number observed ASVs")
t.test(a[a$treatment=="hot","mean_richness_rarefied_5000"],a[a$treatment=="cold","mean_richness_rarefied_5000"])
## they are different by t-test (still)

ggplot(a, aes(x=geno_suscep,y=mean_richness_rarefied_5000,colour=geno_suscep))+
  geom_violin(aes(fill=geno_suscep))+
  scale_color_manual(values=c("black","black"))+
  ggtitle("Observed ASVs")+scale_fill_manual(values=c("white","grey5")) + 
  theme_tufte() +xlab("Laby genotype susceptibility") + ylab("Number observed ASVs")
t.test(a[a$geno_suscep=="high","mean_richness_rarefied_5000"],a[a$geno_suscep=="low","mean_richness_rarefied_5000"])
## not different here either

# rareified data chao 5000
ggplot(a, aes(x=treatment,y=chao_5000,colour=treatment))+
  geom_violin(aes(fill=treatment))+
  scale_color_manual(values=c("black","black"))+
  ggtitle("Observed ASVs")+scale_fill_manual(values=c("white","grey5")) + 
  theme_tufte() +xlab("Treatment") + ylab("Chao estimate of number of ASVs")
t.test(a[a$treatment=="hot","chao_5000"],a[a$treatment=="cold","chao_5000"])
## they are different by t-test (still)

ggplot(a, aes(x=geno_suscep,y=chao_5000,colour=geno_suscep))+
  geom_violin(aes(fill=geno_suscep))+
  scale_color_manual(values=c("black","black"))+
  ggtitle("Observed ASVs")+scale_fill_manual(values=c("white","grey5")) + 
  theme_tufte() +xlab("Laby genotype susceptibility") + ylab("Chao estimate of number of ASVs")
t.test(a[a$geno_suscep=="high","chao_5000"],a[a$geno_suscep=="low","chao_5000"])
## not different here either

#####

## V) Beta Diversity
#####
# functions for this section
logplot<-function(laby,pcx,mvar.clr,variable){
  toplot<-sample_data(laby)
  toplot$pc1<-pcx$x[,1]
  toplot$pc2<-pcx$x[,2]
  toplot<-data.frame(toplot)
  p<-ggplot(data=toplot,aes(x=pc1,y=pc2,colour=toplot[,variable])) +
    geom_point(size=3) + xlab(paste("PC1: ", round(sum(pcx$sdev[1]^2)/mvar.clr, 3), sep=""))+
    ylab(paste("PC2: ", round(sum(pcx$sdev[2]^2)/mvar.clr, 3), sep=""))+theme_bw()#+scale_color_gradient(low="grey90",high="black")
  return(p)
}
# load and filter data
laby2<-laby
laby2<-subset_samples(laby2,geno_suscep%in%c("high","low")) # goes down from 61 to 55 samples
laby2_choosy<-laby2%>%filter_taxa(function(x) sum(x)>10,TRUE)%>% filter_taxa(function(x) sum(x > 3) > (0.2*length(x)),TRUE)

# normalize data
#laby2 <- transform_sample_counts(laby2, fun=function(x) log(1 + x))
#laby2_choosy <- transform_sample_counts(laby2_choosy, fun=function(x) log(1 + x))
laby2_clr<-transform_sample_counts(laby2_choosy,fun=function(x) clr(x))

adonis2(dist(otu_table(laby2_clr),method="eucli")~ treatment, 
        data = as(sample_data(laby2_clr), "data.frame")) 

adonis2(dist(otu_table(laby2_clr),method="eucli")~ treatment*geno_suscep, 
        data = as(sample_data(laby2_clr), "data.frame")) 
mod<-betadisper(dist(otu_table(laby2_clr),method="eucli"),sample_data(laby2_clr)$combo)
mod<-betadisper(dist(otu_table(laby2_clr),method="eucli"),sample_data(laby2_clr)$treatment)

anova(mod) # groups not different beta disper
permutest(mod, pairwise = TRUE, permutations = 999)
(mod.HSD <- TukeyHSD(mod))
plot(mod.HSD) # groups we're interested in not diff

toplot<-sample_data(laby2_clr)
toplot$pc1<-prcomp(otu_table(laby2_clr))$x[,1]
toplot$pc2<-prcomp(otu_table(laby2_clr))$x[,2]
toplot$pc3<-prcomp(otu_table(laby2_clr))$x[,3]
toplot$pc4<-prcomp(otu_table(laby2_clr))$x[,4]
toplot<-data.frame(toplot)
ggplot(data=toplot,aes(x=pc1,y=pc3,fill=treatment,shape=geno_suscep)) +
    geom_point(size=3) + 
    xlab(paste("PC1: ", round(sum(prcomp(otu_table(laby2_clr))$sdev[1]^2)/mvar(otu_table(laby2_clr)), 3), sep=""))+
    #ylab(paste("PC2: ", round(sum(prcomp(otu_table(laby2_clr))$sdev[2]^2)/mvar(otu_table(laby2_clr)), 3), sep=""))+
    ylab(paste("PC3: ", round(sum(prcomp(otu_table(laby2_clr))$sdev[3]^2)/mvar(otu_table(laby2_clr)), 3), sep=""))+
    scale_shape_manual(values=c(24,21))+
    scale_fill_manual(values=c("white","grey15"))+
  guides(fill = guide_legend(override.aes = list(color = c("black","grey15"),shape=c(1,21))), 
         size = guide_legend(override.aes = list(color = c("white","grey15"))))

ggplot(data=toplot,aes(x=pc1,y=pc2,fill=treatment,shape=geno_suscep)) +
  geom_point(size=3) + 
  xlab(paste("PC1: ", round(sum(prcomp(otu_table(laby2_clr))$sdev[1]^2)/mvar(otu_table(laby2_clr)), 3), sep=""))+
  ylab(paste("PC2: ", round(sum(prcomp(otu_table(laby2_clr))$sdev[2]^2)/mvar(otu_table(laby2_clr)), 3), sep=""))+
  #ylab(paste("PC3: ", round(sum(prcomp(otu_table(laby2_clr))$sdev[3]^2)/mvar(otu_table(laby2_clr)), 3), sep=""))+
  scale_shape_manual(values=c(24,21))+
  scale_fill_manual(values=c("white","grey15"))+
  guides(fill = guide_legend(override.aes = list(color = c("black","grey15"),shape=c(1,21))), 
         size = guide_legend(override.aes = list(color = c("white","grey15"))))

## to plot warmed
toplot_w<-sample_data(subset_samples(laby2_clr,treatment=="hot"))
toplot_w$pc1<-prcomp(otu_table(subset_samples(laby2_clr,treatment=="hot")))$x[,1]
toplot_w$pc2<-prcomp(otu_table(subset_samples(laby2_clr,treatment=="hot")))$x[,2]
toplot_w$pc3<-prcomp(otu_table(subset_samples(laby2_clr,treatment=="hot")))$x[,3]
toplot_w$pc4<-prcomp(otu_table(subset_samples(laby2_clr,treatment=="hot")))$x[,4]
toplot_w<-data.frame(toplot_w)
ggplot(data=toplot_w,aes(x=pc1,y=pc2,fill=treatment,shape=geno_suscep)) +
  geom_point(size=3) + 
  xlab(paste("PC1: ", round(sum(prcomp(otu_table(subset_samples(laby2_clr,treatment=="hot")))$sdev[1]^2)/mvar(otu_table(laby2_clr)), 3), sep=""))+
  ylab(paste("PC2: ", round(sum(prcomp(otu_table(subset_samples(laby2_clr,treatment=="hot")))$sdev[2]^2)/mvar(otu_table(laby2_clr)), 3), sep=""))+
  #ylab(paste("PC3: ", round(sum(prcomp(otu_table(laby2_clr))$sdev[3]^2)/mvar(otu_table(laby2_clr)), 3), sep=""))+
  scale_shape_manual(values=c(24,21))+
  scale_fill_manual(values=c("grey15"))
toplot_c<-sample_data(subset_samples(laby2_clr,treatment=="cold"))
toplot_c$pc1<-prcomp(otu_table(subset_samples(laby2_clr,treatment=="cold")))$x[,1]
toplot_c$pc2<-prcomp(otu_table(subset_samples(laby2_clr,treatment=="cold")))$x[,2]
toplot_c<-data.frame(toplot_c)
ggplot(data=toplot_c,aes(x=pc1,y=pc2,fill=treatment,shape=geno_suscep)) +
  geom_point(size=3) + 
  xlab(paste("PC1: ", round(sum(prcomp(otu_table(subset_samples(laby2_clr,treatment=="cold")))$sdev[1]^2)/mvar(otu_table(laby2_clr)), 3), sep=""))+
  ylab(paste("PC2: ", round(sum(prcomp(otu_table(subset_samples(laby2_clr,treatment=="cold")))$sdev[2]^2)/mvar(otu_table(laby2_clr)), 3), sep=""))+
  #ylab(paste("PC3: ", round(sum(prcomp(otu_table(laby2_clr))$sdev[3]^2)/mvar(otu_table(laby2_clr)), 3), sep=""))+
  scale_shape_manual(values=c(24,21))+
  scale_fill_manual(values=c("white"))


#see continued for the reset of the script