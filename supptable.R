## Here is the code for the supplemental table for microbiomes. This relies on the other two scripts being run first

treatment_alone$which<-"treatment_alone"
treatment_alone$test<-"treatment_alone"
a1$which<-"a1"
a2$which<-"a2"
a3$which<-"a3"
a4$which<-"a4"
l<-rbind(a1,a2,a3,a4,treatment_alone)

colnames(l)
l<-l[,c(1,2,3,6,7,8,9,10,11,12,13,14,16)]
l
write.csv(l,file="alllaby.csv")
unique(l$Sequence)

l$Sequence
a1t<-diagdds %>% results(cooksCutoff = FALSE,contrast=c("combo","high_hot","high_cold")) %>% as.data.frame() %>%
  cbind(as(tax_table(laby2)[rownames(diagdds), ], "matrix")) %>% cbind("Sequence" = rownames(diagdds))
a1t<-a1t[a1t$Sequence %in% l$Sequence,]

a2t<-diagdds %>% results(cooksCutoff = FALSE,contrast=c("combo","high_hot","low_hot")) %>% as.data.frame() %>%
  cbind(as(tax_table(laby2)[rownames(diagdds), ], "matrix")) %>% cbind("Sequence" = rownames(diagdds)) 
a2t<-a2t[a2t$Sequence %in% l$Sequence,]

a3t<-diagdds %>% results(cooksCutoff = FALSE,contrast=c("combo","high_cold","low_cold")) %>% as.data.frame() %>%
  cbind(as(tax_table(laby2)[rownames(diagdds), ], "matrix")) %>% cbind("Sequence" = rownames(diagdds)) 
a3t<-a3t[a3t$Sequence %in% l$Sequence,]

a4t<-diagdds %>% results(cooksCutoff = FALSE,contrast=c("combo","low_hot","low_cold")) %>% as.data.frame() %>%
  cbind(as(tax_table(laby2)[rownames(diagdds), ], "matrix")) %>% cbind("Sequence" = rownames(diagdds)) 
a4t<-a4t[a4t$Sequence %in% l$Sequence,]

treatment_alonet<-diagdds1 %>% results(cooksCutoff = FALSE,contrast=c("treatment","hot","cold")) %>% as.data.frame() %>%
  cbind(as(tax_table(laby2)[rownames(diagdds1), ], "matrix")) %>% cbind("Sequence" = rownames(diagdds1))

treatment_alonet<-treatment_alonet[treatment_alonet$Sequence %in% l$Sequence,]
colnames(a1t)
l2<-cbind(a1t[,c(7:13,1:3,6)],a2t[,c(2,3,6)],a3t[,c(2,3,6)],a4t[,c(2,3,6)],treatment_alonet[,c(1,2,3,6,13)])
write.csv(l2,file="alllaby2.csv")
