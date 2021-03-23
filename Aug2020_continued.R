## relies on libraries loaded in Aug 2020 script. 
## This is a little messier than that as it explores different combinat
#####
## VI) Taxonomic differences by DESeq2 betweent treatments and susceptibility with interaction in DESeq2
#####

# this section runs just to compare data to (if don't include susceptibiltiy)
#####
## requires treatment alone subset
laby2<-laby
laby2<-subset_samples(laby2,geno_suscep%in%c("high","low")) # goes down from 61 to 55 samples

diagdds1 = phyloseq_to_deseq2(laby2, ~ treatment)
diagdds1 = estimateSizeFactors(diagdds1)
diagdds1 = DESeq(diagdds1, test="Wald", fitType="local")
alpha=0.05
treatment_alone<-diagdds1 %>% results(cooksCutoff = FALSE,contrast=c("treatment","hot","cold")) %>% as.data.frame() %>%
  cbind(as(tax_table(laby2)[rownames(diagdds1), ], "matrix")) %>% cbind("Sequence" = rownames(diagdds1)) %>% 
  filter(padj < alpha)  # 95 these are things that vary by temp 
treatment_alone$hotcold<-ifelse(treatment_alone$log2FoldChange>0,"hot","cold")

treatment_alone2<-treatment_alone
treatment_alone2$famgen<-paste(treatment_alone2$Phylum,treatment_alone2$Class,treatment_alone2$Family,treatment_alone2$Genus)
treatment_alone2<-droplevels(treatment_alone2)
treatment_alone2<-treatment_alone2[order(-treatment_alone2$log2FoldChange),]
treatment_alone2<-treatment_alone2[order(treatment_alone2$famgen),]
treatment_alone2[!treatment_alone2$Sequence%in%a8$Sequence,"notsecond"]<-"yes"
treatment_alone2[is.na(treatment_alone2$notsecond),"notsecond"]<-"no"
treatment_alone2$notsecond<-factor(treatment_alone2$notsecond,levels=c("yes","no"))
treatment_alone2<-treatment_alone2[order(treatment_alone2$notsecond),]
treatment_alone2$Sequence<-factor(treatment_alone2$Sequence,levels=as.character(treatment_alone2$Sequence))
treatment_alone2$labels <- ifelse(!is.na(treatment_alone2$Species),paste(treatment_alone2$Genus,treatment_alone2$Species),
                                  ifelse(!is.na(treatment_alone2$Genus),paste(treatment_alone2$Genus,"sp."),
                                         ifelse(treatment_alone2$Family=="Unknown_Family", paste(treatment_alone2$Family, "in",treatment_alone2$Class), paste(treatment_alone2$Family,"spp."))
                                  ))
#####


#load data
laby2<-laby
laby2<-subset_samples(laby2,geno_suscep%in%c("high","low")) # goes down from 61 to 55 samples

# run the basic test (is something different based on combo of geno_suscep and temp)
diagdds = phyloseq_to_deseq2(laby2, ~ combo)
diagdds = estimateSizeFactors(diagdds)
diagdds = DESeq(diagdds, test="Wald", fitType="local")
alpha=.05

# do a bunch of contrasts to see what varies between treatments
a1<-diagdds %>% results(cooksCutoff = FALSE,contrast=c("combo","high_hot","high_cold")) %>% as.data.frame() %>%
  cbind(as(tax_table(laby2)[rownames(diagdds), ], "matrix")) %>% cbind("Sequence" = rownames(diagdds)) %>% 
  filter(padj < alpha) ## 69 these are things that vary by treatment in susceptable genotypes
a1.a<-diagdds %>% results(cooksCutoff = FALSE,contrast=c("combo","high_hot","high_cold")) %>% as.data.frame() %>%
  cbind(as(tax_table(laby2)[rownames(diagdds), ], "matrix")) %>% cbind("Sequence" = rownames(diagdds))## unfiltered so can pull non-sig ASVs as needed to compare

a2<-diagdds %>% results(cooksCutoff = FALSE,contrast=c("combo","high_hot","low_hot")) %>% as.data.frame() %>%
  cbind(as(tax_table(laby2)[rownames(diagdds), ], "matrix")) %>% cbind("Sequence" = rownames(diagdds)) %>% 
  filter(padj < alpha) # 12 these are things that very in susceptible genotypes when temp is hot
a2.a<-diagdds %>% results(cooksCutoff = FALSE,contrast=c("combo","high_hot","low_hot")) %>% as.data.frame() %>%
  cbind(as(tax_table(laby2)[rownames(diagdds), ], "matrix")) %>% cbind("Sequence" = rownames(diagdds)) #

a3<-diagdds %>% results(cooksCutoff = FALSE,contrast=c("combo","high_cold","low_cold")) %>% as.data.frame() %>%
  cbind(as(tax_table(laby2)[rownames(diagdds), ], "matrix")) %>% cbind("Sequence" = rownames(diagdds)) %>% 
  filter(padj < alpha)  # 8 things that vary by susceptible genotype in cold background
a3.a<-diagdds %>% results(cooksCutoff = FALSE,contrast=c("combo","high_cold","low_cold")) %>% as.data.frame() %>%
  cbind(as(tax_table(laby2)[rownames(diagdds), ], "matrix")) %>% cbind("Sequence" = rownames(diagdds)) 


a4<-diagdds %>% results(cooksCutoff = FALSE,contrast=c("combo","low_hot","low_cold")) %>% as.data.frame() %>%
  cbind(as(tax_table(laby2)[rownames(diagdds), ], "matrix")) %>% cbind("Sequence" = rownames(diagdds)) %>% 
  filter(padj < alpha)  # 37 these are things that vary by temp in the non siusceptible genotype background
a4.a<-diagdds %>% results(cooksCutoff = FALSE,contrast=c("combo","low_hot","low_cold")) %>% as.data.frame() %>%
  cbind(as(tax_table(laby2)[rownames(diagdds), ], "matrix")) %>% cbind("Sequence" = rownames(diagdds))# 47 these are things that vary by temp in the non siusceptible genotype background

# explore where there might be overlap

# ASVs that only vary between treatments in non-susceptible treatments (9)
sum(!(a4$Sequence%in%a1$Sequence)&!(a4$Sequence%in%a2$Sequence)&!(a4$Sequence%in%a3$Sequence)) 
# ASVs that vary between genotype group in hot treatments (6)
sum(!(a2$Sequence%in%a1$Sequence)&!(a2$Sequence%in%a3$Sequence)&!(a2$Sequence%in%a4$Sequence)) 
# ASVs that vary between genotype group in cold treatments (3)
sum(!(a3$Sequence%in%a1$Sequence)&!(a3$Sequence%in%a2$Sequence)&!(a3$Sequence%in%a4$Sequence)) 

# ASVs that  vary between treatments in all genotypes (26)
sum(a1$Sequence%in%a4$Sequence)
# ASVs that between genotype groups in both treatments (2) 
sum(a2$Sequence%in%a3$Sequence)
# All 3 ASVs that vary in susceptible treatments when hot also vary in by treatment in susceptible genotypes
sum(a1$Sequence%in%a2$Sequence) # 4
# All 3 ASVs that vary in susceptible treatments when cold also vary in by treatment in susceptible genotypes
sum(a1$Sequence%in%a3$Sequence) # 3
# one ASV that varies between suspetible genotypes when temp is hot also varies in temperture when temp in non susceptible genotypes
sum(a2$Sequence%in%a4$Sequence) #3

dat<-c(as.character(a1[a1$Sequence%in%a2$Sequence,"Sequence"]),as.character(a1[a1$Sequence%in%a3$Sequence,"Sequence"]),as.character(a2[a2$Sequence%in%a4$Sequence,"Sequence"]))
data <- plotCounts(diagdds, gene=as.character(dat[6]), intgroup=c("treatment","geno_suscep"), returnData=TRUE)

length(a1[(a1$Sequence%in%a2$Sequence&a1$Sequence%in%a3$Sequence&a1$Sequence%in%a4$Sequence),"test"]) # 0 
length(a1[(a1$Sequence%in%a2$Sequence&a1$Sequence%in%a3$Sequence&!a1$Sequence%in%a4$Sequence),"test"]) # 0
length(a1[(a1$Sequence%in%a2$Sequence&!a1$Sequence%in%a3$Sequence&!a1$Sequence%in%a4$Sequence),"test"]) # 0

length(a1[(a1$Sequence%in%a2$Sequence&a1$Sequence%in%a3$Sequence),"test"]) # 0

#adding identifiers, there's some extra text but this is just to set up identifiers
a1[(a1$Sequence%in%a2$Sequence),"test"]  
## 
a1[(a1$Sequence%in%a2$Sequence),"test"]<-"s_hvc and h_svn"
a2[(a2$Sequence%in%a1$Sequence),"test"]<-"s_hvc and h_svn"

# All 3 ASVs that vary in susceptible treatments when cold also vary in by treatment in susceptible treatments
a1[a1$Sequence%in%a3$Sequence,"test"]<-"s_hvc and c_svn"
a3[a3$Sequence%in%a1$Sequence,"test"]<-"s_hvc and c_svn"

a1[a1$Sequence%in%a4$Sequence,"test"]<-"n_hvc and s_hvc"
a4[a4$Sequence%in%a1$Sequence,"test"]<-"n_hvc and s_hvc"


a2[a2$Sequence%in%a3$Sequence,"test"]<-"c_svn and h_svn"
a3[a3$Sequence%in%a2$Sequence,"test"]<-"c_svn and h_svn"

# one ASV that varies between suspetible genotypes when temp is hot also varies in temperture when temp in non susceptible genotypes
a2[a2$Sequence%in%a4$Sequence,"test"]<-"n_hvc and h_svn"
a4[a4$Sequence%in%a2$Sequence,"test"]<-"n_hvc and h_svn"

a3[a3$Sequence%in%a4$Sequence,"test"]<-"n_hvc and c_svn"
a4[a4$Sequence%in%a3$Sequence,"test"]<-"n_hvc and c_svn"

#genovary
a5<-rbind(a2,a3[!a3$Sequence%in%a2$Sequence,])
a6<-rbind(a1,a4)
a5$highlow<-ifelse(a5$log2FoldChange>0,"high","low")
a5$famgen<-paste(a5$Phylum,a5$Class,a5$Family,a5$Genus)
a5<-droplevels(a5)
a5<-a5[order(-a5$log2FoldChange),]
a5<-a5[order(substr(a5$test,nchar(a5$test)-5,nchar(a5$test))),]
a5$Sequence<-factor(a5$Sequence,levels=as.character(a5$Sequence))
a5$Sequence<-factor(a5$Sequence,levels=c(as.character(a5$Sequence)[c(1:5,7)],as.character(a5$Sequence[6])))
a5$labels <- ifelse(!is.na(a5$Species),paste(a5$Genus,a5$Species),
                    ifelse(!is.na(a5$Genus),paste(a5$Genus,"sp."),
                           ifelse(a5$Family=="Unknown_Family", paste(a5$Family, "in",a5$Class), paste(a5$Family,"spp."))
                    ))
a6.a<-a6[!is.na(a6$test),]
a6.a$Sequence<-factor(a6.a$Sequence,levels=c(as.character(a5$Sequence)[c(1:5,7)],as.character(a5$Sequence[6])))
a6$combo<-ifelse(a6$log2FoldChange>0,"hot treatment","cold treatment")


#treatmentvary
a1[(a1$Sequence%in%a4$Sequence),"test"]<-"temp_always"
a4[(a4$Sequence%in%a1$Sequence),"test"]<-"temp_always"
a1$which<-"a1"
a4$which<-"a4"
a7<-rbind(a1,a4)
a8<-rbind(filter(a7,is.na(test)),filter(a7,test!="temp_always"),filter(a7,test=="temp_always",which=="a1"))
a8$highlow<-ifelse(a8$log2FoldChange>0,"high","low")
a8$famgen<-paste(a8$Phylum,a8$Class,a8$Family,a8$Genus)
a8<-droplevels(a8)
a8<-a8[order(-a8$log2FoldChange),]
a8$sorting<-paste(a8$test,a8$which)
a8[a8$Sequence%in%treatment_alone2$Sequence,"first"]<-"yes"
a8[is.na(a8[,"first"]),"first"]<-"no"
a8<-a8[order(a8$first),]
a8$sorting<-factor(a8$sorting,levels=c("s_hvc and c_svn a1","s_hvc and h_svn a1",
                                       "NA a1","n_hvc and h_svn a4","n_hvc and c_svn a4", "NA a4","temp_always a1"))
a8<-a8[order(a8$sorting),]

a8$Sequence<-factor(a8$Sequence,levels=as.character(a8$Sequence))
a8$labels <- ifelse(!is.na(a8$Species),paste(a8$Genus,a8$Species),
                    ifelse(!is.na(a8$Genus),paste(a8$Genus,"sp."),
                           ifelse(a8$Family=="Unknown_Family", paste(a8$Family, "in",a8$Class), paste(a8$Family,"spp."))
                    ))
a9<-rbind(filter(a7,is.na(test)),filter(a7,test!="temp_always"),filter(a7,test=="temp_always",which=="a4"))
a9$Sequence<-factor(a9$Sequence,levels=as.character(a9$Sequence))
a9$highlow<-ifelse(a9$log2FoldChange>0,"high","low")
a9$famgen<-paste(a9$Phylum,a9$Class,a9$Family,a9$Genus)
a9<-droplevels(a9)
a9<-a9[order(-a8$log2FoldChange),]
a9$sorting<-paste(a9$test,a9$which)
a9[a9$Sequence%in%treatment_alone2$Sequence,"first"]<-"yes"
a9[is.na(a9[,"first"]),"first"]<-"no"
a9<-a9[order(a9$first),]
a9[a9$sorting=="temp_always a4","sorting"]<-"temp_always"
a9$sorting<-factor(a9$sorting,levels=c("s_hvc and c_svn a1","s_hvc and h_svn a1",
                                       "NA a1","n_hvc and h_svn a4", "NA a4","temp_always"   ))
a9<-a9[order(a9$sorting),]
a9$Sequence<-factor(a9$Sequence,levels=as.character(a8[,"Sequence"]))
a9<-a9[order(a9$Sequence),]
a9$labels <- ifelse(!is.na(a9$Species),paste(a9$Genus,a9$Species),
                    ifelse(!is.na(a9$Genus),paste(a9$Genus,"sp."),
                           ifelse(a9$Family=="Unknown_Family", paste(a9$Family, "in",a9$Class), paste(a9$Family,"spp."))
                    ))

test<-ggplot(a8, aes(x=Sequence, y=log2FoldChange,fill=sorting,shape=first)) + 
  geom_hline(yintercept = 0.0, color = "gray", size = 0.5) +  
  geom_point(aes(x=a9$Sequence,y=a9$log2FoldChange),size=2) + 
  geom_point(size=2) +  
  geom_point(data=a9[55:80,],aes(x=Sequence,y=log2FoldChange),size=2,fill="blue") + 
  geom_point(data=a8[55:80,],aes(x=Sequence,y=log2FoldChange),size=2,fill="firebrick") + 
  scale_shape_manual(values=c(24,21))+
  geom_errorbar(aes(ymin=a9$log2FoldChange-a9$lfcSE,ymax=a9$log2FoldChange+a9$lfcSE),width=0.5) +
  geom_errorbar(aes(ymin=log2FoldChange-lfcSE,ymax=log2FoldChange+lfcSE),width=0.5) +
  coord_flip()+
  geom_vline(xintercept=54.5) +
  geom_vline(xintercept=43.5) +
  scale_x_discrete(breaks=a8$Sequence,labels=a8$labels) +
  scale_fill_manual(values=c("firebrick","firebrick","firebrick","blue","blue","blue","black"))  +
  scale_y_continuous(expand=c(0,0))+
  ylab("log2 fold change")+
  theme(panel.grid.major.y = element_blank(),
        axis.title.y=element_blank(),legend.position = "none",
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent", color = NA)) +
  ylim(c(-50,50)) +  scale_y_continuous(expand=c(0,0)) 
test
ggsave("~/Desktop/temp.png",test) # this is supp figure for the MS
# so this shows the log2fold change for temperature treatment with groups including also varies in treatment alone comparisons and varying by genotype susceptibility

# comparing genotype group log2fold changes controlling for temp diffs
sum(a1$Sequence%in%a2$Sequence) #4
sum(a1$Sequence%in%a3$Sequence) # 3
sum(a1$Sequence%in%a4$Sequence) # 26
sum(a2$Sequence%in%a1$Sequence) # 4
sum(a2$Sequence%in%a3$Sequence) # 2
sum(a2$Sequence%in%a4$Sequence) # 3
sum(a3$Sequence%in%a1$Sequence) #3
sum(a3$Sequence%in%a2$Sequence) #2
sum(a3$Sequence%in%a4$Sequence) # 3
sum(a4$Sequence%in%a1$Sequence) #26
sum(a4$Sequence%in%a2$Sequence) #3
sum(a4$Sequence%in%a3$Sequence) #3
sum(a1$Sequence%in%a2$Sequence & a1$Sequence%in%a3$Sequence) #1
sum(a1$Sequence%in%a2$Sequence & a1$Sequence%in%a4$Sequence) #2
sum(a1$Sequence%in%a3$Sequence & a1$Sequence%in%a4$Sequence) #2
sum(a2$Sequence%in%a1$Sequence & a2$Sequence%in%a3$Sequence) #1
sum(a2$Sequence%in%a1$Sequence & a2$Sequence%in%a4$Sequence) #2
sum(a2$Sequence%in%a3$Sequence & a2$Sequence%in%a4$Sequence) #1
sum(a3$Sequence%in%a1$Sequence & a3$Sequence%in%a2$Sequence) #1
sum(a3$Sequence%in%a1$Sequence & a3$Sequence%in%a4$Sequence) #2
sum(a3$Sequence%in%a2$Sequence & a3$Sequence%in%a4$Sequence) #1
sum(a4$Sequence%in%a1$Sequence & a4$Sequence%in%a2$Sequence) #2
sum(a4$Sequence%in%a1$Sequence & a4$Sequence%in%a3$Sequence) #2
sum(a4$Sequence%in%a2$Sequence & a4$Sequence%in%a3$Sequence) #1
sum(a4$Sequence%in%a2$Sequence & a4$Sequence%in%a3$Sequence & a4$Sequence%in%a1$Sequence) #1
# some overlap among three groups (1 in everything)
a4[a4$Sequence%in%a2$Sequence & a4$Sequence%in%a3$Sequence & a4$Sequence%in%a1$Sequence,] # maribacter



treatment_alone<-diagdds1 %>% results(cooksCutoff = FALSE,contrast=c("treatment","hot","cold")) %>% as.data.frame() %>%
  cbind(as(tax_table(laby2)[rownames(diagdds1), ], "matrix")) %>% cbind("Sequence" = rownames(diagdds1)) %>% 
  filter(padj < alpha)  # 54 these are things that vary by temp in the non siusceptible genotype background
# write.csv(treatment_alone,file="treatonly.csv") 

sum(a1$Sequence%in%a2$Sequence) #3
sum(a1$Sequence%in%a3$Sequence) # 4
sum(a1$Sequence%in%a4$Sequence) # 26
sum(a2$Sequence%in%a1$Sequence) #4
sum(a2$Sequence%in%a3$Sequence) # 2
sum(a2$Sequence%in%a4$Sequence) # 3
sum(a3$Sequence%in%a1$Sequence) #3
sum(a3$Sequence%in%a2$Sequence) #2
sum(a3$Sequence%in%a4$Sequence) # 3
sum(a4$Sequence%in%a1$Sequence) #26
sum(a4$Sequence%in%a2$Sequence) #3
sum(a4$Sequence%in%a3$Sequence) #3
treatment_alone$Sequence%in%rbind(a1,a2,a3,a4)$Sequence
unique(rbind(a1,a2,a3,a4)$Sequence)%in%treatment_alone$Sequence

#write.csv(rbind(a1,a2,a3,a4),file="~/Desktop/labytax2.csv")

## more less vary:
(a6[!is.na(a6$test),"Sequence"]) %in% a5$Sequence
a2[(a2$Sequence%in%a3$Sequence),"test"]<-"effect_always"
a3[(a3$Sequence%in%a2$Sequence),"test"]<-"effect_always"
a2$which<-"a2"
a3$which<-"a3"
a7<-rbind(a2,a3)
a8<-rbind(filter(a7,is.na(test)),filter(a7,test!="effect_always"),filter(a7,test=="effect_always",which=="a2"))
a8$highlow<-ifelse(a8$log2FoldChange>0,"high","low")
a8$famgen<-paste(a8$Phylum,a8$Class,a8$Family,a8$Genus)
a8<-droplevels(a8)
a8<-a8[order(-a8$log2FoldChange),]
a8$sorting<-paste(a8$test,a8$which)
a8[a8$Sequence%in%treatment_alone2$Sequence,"first"]<-"yes"
a8[is.na(a8[,"first"]),"first"]<-"no"
a8<-a8[order(a8$first),]
#a9[a9$sorting=="effect_always a2","sorting"]<-"effect_always"
a8$sorting<-factor(a8$sorting,levels=c("s_hvc and c_svn a3", "n_hvc and c_svn a3", "NA a3","s_hvc and h_svn a2",
                                       "n_hvc and h_svn a2", "NA a2",
                                       "effect_always a3", "effect_always a2"))
a8<-a8[order(a8$sorting),]

a8$Sequence<-factor(a8$Sequence,levels=as.character(a8$Sequence))
a8$labels <- ifelse(!is.na(a8$Species),paste(a8$Genus,a8$Species),
                    ifelse(!is.na(a8$Genus),paste(a8$Genus,"sp."),
                           ifelse(a8$Family=="Unknown_Family", paste(a8$Family, "in",a8$Class), paste(a8$Family,"spp."))
                    ))
a9<-rbind(filter(a7,is.na(test)),filter(a7,test!="effect_always"),filter(a7,test=="effect_always",which=="a3"))
a9$Sequence<-factor(a9$Sequence,levels=as.character(a9$Sequence))
a9$highlow<-ifelse(a9$log2FoldChange>0,"high","low")
a9$famgen<-paste(a9$Phylum,a9$Class,a9$Family,a9$Genus)
a9<-droplevels(a9)
a9<-a9[order(-a8$log2FoldChange),]
a9$sorting<-paste(a9$test,a9$which)
a9[a9$Sequence%in%treatment_alone2$Sequence,"first"]<-"yes"
a9[is.na(a9[,"first"]),"first"]<-"no"
a9<-a9[order(a9$first),]
a9[a9$sorting=="effect_always a3","sorting"]<-"effect_always"
a9$sorting<-factor(a9$sorting,levels=c("NA a2", "effect_always a3","NA NA","s_hvc and h_svn a2", "NA a3",
                                       "n_hvc and h_svn a2", "n_hvc and c_svn a3", "s_hvc and c_svn a3" ))
a9<-a9[order(a9$sorting),]
a9$Sequence<-factor(a9$Sequence,levels=as.character(a8[,"Sequence"]))
a9<-a9[order(a9$Sequence),]
a9$labels <- ifelse(!is.na(a9$Species),paste(a9$Genus,a9$Species),
                    ifelse(!is.na(a9$Genus),paste(a9$Genus,"sp."),
                           ifelse(a9$Family=="Unknown_Family", paste(a9$Family, "in",a9$Class), paste(a9$Family,"spp."))
                    ))

test<-ggplot(a8, aes(x=Sequence, y=log2FoldChange,fill=sorting,shape=sorting,alpha=sorting)) + 
  geom_hline(yintercept = 0.0, color = "gray", size = 0.5) +  
  geom_errorbar(aes(ymin=a9$log2FoldChange-a9$lfcSE,ymax=a9$log2FoldChange+a9$lfcSE),width=0.5,alpha=1) +
  geom_errorbar(aes(ymin=log2FoldChange-lfcSE,ymax=log2FoldChange+lfcSE),width=0.5,alpha=1) +
  #annotate("rect",xmin=0,xmax=6.5,ymin=-40,ymax=30,fill="red",alpha=0.15)+
  #annotate("rect",xmin=6.5,xmax=16.5,ymin=-40,ymax=30,fill="blue",alpha=0.15)+  
  geom_point(aes(x=a9$Sequence,y=a9$log2FoldChange),size=2) + 
  geom_point(size=2) +  
  geom_point(aes(x=a9$Sequence[17],y=a9$log2FoldChange[17]),size=2,fill="white") + 
  geom_point(aes(x=a9$Sequence[18],y=a9$log2FoldChange[18]),size=2,fill="white") + 
  geom_point(aes(x=Sequence[17],y=log2FoldChange[17]),size=2,fill="black") + 
  geom_point(aes(x=Sequence[18],y=log2FoldChange[18]),size=2,fill="black") + 
  scale_alpha_manual(values=c(1,1,1,1,1,1,1))+
  scale_shape_manual(values=c(21,21,21,21,21,21,21))+
  coord_flip()+
  scale_x_discrete(breaks=a8$Sequence,labels=a8$labels) +
  scale_fill_manual(values=c("grey15","grey15","grey15","white","white","white","white","black"))  +
  geom_vline(xintercept=16.5) +
  geom_vline(xintercept=6.5) +
  scale_y_continuous(expand=c(0,0))+
  ylab("log2 fold change")+
  theme(panel.grid.major.y = element_blank(),
        axis.title.y=element_blank(),legend.position = "none",
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent", color = NA)) +
  ylim(c(-40,30)) #+  scale_y_continuous(expand=c(0,0)) 
test  
# this is the main treatment comparison fig




