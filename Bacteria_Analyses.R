# Fire coral microbiome dataset

# Authors: Caroline Dubé
# Copyright (C) 2021 Caroline Dubé
# License: GPL

### Dependencies ===============================================================
library(vegan)
library(ggplot2)
library(plyr)
library(dplyr)
library(reshape2)
library(tidyverse)
library(expm)
library(indicspecies)
library(pairwiseAdonis)

### Open Files =================================================================
Design_Files <- read.csv2("~/Desktop/Research/Microbiome/Review/qiime2/Design_Files.csv")
Design_Files$Group<-paste(Design_Files$Ind,ifelse(Design_Files$LOC=="M",substring(Design_Files$Ind, 1, 4),substring(Design_Files$Ind, 1, 2)),sep="_")
row.names(Design_Files)<-Design_Files$Ind

ASVTable_Rar <- read.csv("~/Desktop/Research/Microbiome/Review/qiime2/ASVTable_Rar.csv", sep=";")
row.names(ASVTable_Rar)<-Design_Files$Ind

Mill<-merge(Design_Files,ASVTable_Rar,by="row.names",all.y=T)
row.names(Mill)<-Mill[,1]

FamBar <- read.csv("~/Desktop/Research/Microbiome/Review/qiime2/taxonomy_file.csv", sep=";")

size<-read.csv2("~/Desktop/Research/Microbiome/Review/qiime2/ColonySize_Analyses.csv", header=T, sep=";")

### Select 6 genotypes ("46","49","113","174","268","370") =====================
Select<-subset(Mill,MLL%in%c("46","49","113","174","268","370"))
Select<-data.frame(t(Select[13:42719]))
Gen<-data.frame(apply(Select == 0, 1, sum))
Gen$sumrow<-apply(Select, 1, sum)
Gen$Presence<-135-Gen[,1]
Gen<-subset(Gen,!Presence==0)
Gen<-merge(Gen,Select,by="row.names",all.x=T)
row.names(Gen)<-Gen$Row.names
clean.1<-Gen[,c(5:139)]
row.names(clean.1)<-Gen[,1]
clean.2<-data.frame(t(clean.1))
clean.2<-merge(Design_Files,clean.2,by="row.names",all.y=T)
row.names(clean.2)<-clean.2[,1]

### Core Microbiome ============================================================
## >80% of occurrence = core members
clean.2<-subset(clean.2[c(13:21804),])
clean.2$Feature.ID<-row.names(clean.1)
clean.3<-merge(clean.2,FamBar,by.x="row.names",by.y="OTU", all.x=T)
row.names(clean.3)<-clean.3$Row.names
clean.3<-subset(clean.3[,c(2:136)])
clean.3<-data.frame(t(clean.3))
Core<-data.frame(ifelse(clean.3[,c(13:21804)]>1,1,0))
Core$Ind<-row.names(Core)
Core.2<-melt(Core,id="Ind")
Core.2<-ddply(Core.2,~variable,function(x){c(sum(sum(x$value)))})
Core.2$pourcent<-Core.2$V1/135
Core.3<-subset(Core.2,pourcent>=0.8)
Core.3<-merge(Core.3,FamBar,by.x="variable",by.y="OTU")
write.csv2(Core.3,"80pourcent")

### Bacteria community composition analysis ASV level ==========================
## Two-way Permanova
MLL.1<-clean.2[,c(13:21804)]
MLL.1<-log(MLL.1+1)
adonis(MLL.1~as.factor(MLL)*Hab, data=clean.2, permutations = 9999, method = "bray")
adonis(MLL.1~as.factor(LifSta), data=clean.4bis, permutations = 9999, method = "bray")

## Betadisper analysis
dis<-vegdist(MLL.1)
groups_genhab<- factor(c(rep(1,38),rep(2,16),rep(3,15),rep(4,32),rep(5,13),rep(6,21),rep(7,19),rep(8,90),rep(9,26)), labels = c("46","49","113","174","268","370","BR", "UP", "MD")) 
betagenhab<-betadisper(dis, groups_genhab, type=c("median"), sqrt.dist=FALSE, bias.adjust=TRUE)

plot(betagenhab, ellipse=TRUE, hull=FALSE)
boxplot(betagen, ellipse=TRUE, hull=FALSE)

anova(betagen)
permutest(betagen, pairwise=TRUE, permutations=9999)

groups_gen<- factor(c(rep(1,38), rep(2,16), rep(3,15), rep(4,32), rep(5,13), rep(6,21)), labels = c("46","49","113","174","268","370")) 
betagen<-betadisper(dis, groups_gen, type=c("median"), sqrt.dist=FALSE, bias.adjust=TRUE)

plot(betagen, ellipse=TRUE, hull=FALSE)
boxplot(betagen, ellipse=TRUE, hull=FALSE)

anova(betagen)
permutest(betagen, pairwise=TRUE, permutations=9999)

groups_hab<- factor(c(rep(7,19), rep(8,90), rep(9,26)), labels = c("BR", "UP", "MD"))
betahab<-betadisper(dis, groups_hab, type=c("median"), sqrt.dist=FALSE, bias.adjust=TRUE)

plot(betahab, ellipse=TRUE, hull=FALSE)
boxplot(betahab, ellipse=TRUE, hull=FALSE)

anova(betahab)
permutest(betagen, pairwise=TRUE, permutations=9999)

## Anosim
veg<- vegdist(MLL.1)
simHAB<- with(clean.2, anosim(veg, Hab, permutations = 9999))
t<-summary(simHAB)
plot(simHAB)

veg<- vegdist(MLL.1)
simGEN<- with(clean.2, anosim(veg, MLL, permutations = 9999))
t<-summary(simGEN)
plot(simGEN)

## Analyse log per habitat: "Bray curtis" + one-way Permanova + nMDS
x=c("BR","MD","UP")
for(i in x){
  MLL.2<-data.frame(subset(clean.2,Hab==i))
  MLL.1<-MLL.2[,c(13:21804)]
  MLL.1<-log(MLL.1+1)
  y<-adonis(MLL.1~MLL, data=MLL.2, permutations = 9999, method = "bray")
  x.MDS<-metaMDS(x)
  Mill.scores <- as.data.frame(scores(x.MDS)) 
  Mill.scores<-merge(Mill.scores,clean.2[,c(2:11)],by="row.names",all.x=T)
  
  environment_basis<-ggplot(Mill.scores)+ 
    geom_point(aes(x=NMDS1,y=NMDS2,shape=factor(Mill.scores$MLL),colour=factor(Mill.scores$MLL),fill=factor(Mill.scores$MLL)),size=3)+
    ggtitle(y$aov.tab$Pr[1])+
    stat_ellipse(geom = "polygon",aes(x=NMDS1,y=NMDS2,color=factor(Mill.scores$MLL),fill =factor(Mill.scores$MLL)),alpha=0.2,level=0.7)+
    scale_fill_manual(name="Genotypes",breaks = c("46","49","113","174","268","370"),labels=c("G1","G2","G3","G4","G5","G6"),values=c("chocolate1","seagreen1","deeppink", "deepskyblue", "darkorchid3","olivedrab"))+ 
    scale_color_manual(name="Genotypes",breaks = c("46","49","113","174","268","370"),labels=c("G1", "G2","G3","G4","G5","G6"),values=c("chocolate1","seagreen1","deeppink", "deepskyblue", "darkorchid3","olivedrab"))+
    theme_bw()+theme(panel.grid.major=element_line(colour="white"),panel.grid.minor=element_line(colour="white"),strip.background = element_rect(fill="white"),legend.position="none")
  environment_basis
  ggsave(environment_basis, height = 4.5, width = 6.5, filename=paste(i,"Gen",".pdf", sep=""), path = "~/Desktop/Research/Microbiome/Review/R_analyses/Gen/Final2")
  
}

## Analysis log per genotype: "Bray curtis" + one-way Permanova + nMDS
x=c("46","49","113","174","268","370")
for(i in x){
  MLL.2<-data.frame(subset(clean.2,MLL==i))
  MLL.1<-MLL.2[,c(13:21804)]
  MLL.1<-log(MLL.1+1)
  y<-adonis(MLL.1~Hab, data=MLL.2, permutations = 9999, method = "bray")
  x.MDS<-metaMDS(x)
  Mill.scores <- as.data.frame(scores(x.MDS)) 
  Mill.scores<-merge(Mill.scores,clean.2[,c(2:11)],by="row.names",all.x=T)
  
  genetic_basis<-ggplot(Mill.scores)+ 
    geom_point(aes(x=NMDS1,y=NMDS2,shape=factor(Mill.scores$Hab),colour=factor(Mill.scores$Hab),fill=factor(Mill.scores$Hab)),size=3)+
    ggtitle(y$aov.tab$Pr[1])+
    stat_ellipse(geom = "polygon",aes(x=NMDS1,y=NMDS2,color=factor(Mill.scores$Hab),fill =factor(Mill.scores$Hab)),alpha=0.2,level=0.7)+
    scale_fill_manual(name="Habitat",values=c("MD"="dodgerblue4","UP"="gold3", "BR"="firebrick1"))+
    scale_colour_manual(name="Habitat",values=c("MD"="dodgerblue4","UP"="gold3", "BR"="firebrick1"))+
    scale_shape_manual(name="Habitat",values=c("MD"=16,"UP"=17, "BR"=15))+ 
    theme_bw()+theme(panel.grid.major=element_line(colour="white"),panel.grid.minor=element_line(colour="white"),strip.background = element_rect(fill="white"),legend.position="none")
  genetic_basis  
  ggsave(genetic_basis,height = 4.5, width = 6.5, filename=paste(i,"Hab",".pdf",sep=""), path = "~/Desktop/Research/Microbiome/Review/R_analyses/Hab/Final2")
  
}

## Adonis + Pairwise comparisons
# Mid slope habitat
MLL.3<-data.frame(subset(clean.2,Hab=="MD"))
MLL.1<-MLL.3[,c(13:21804)]
MLL.1<-log(MLL.1+1)
adonis(MLL.1~as.factor(MLL), data=MLL.3, permutations = 9999, method = "bray")

MLL.3$MLL <- factor(MLL.3$MLL, levels =c("46","49","113","370"))
pairwise.adonis(MLL.1, MLL.3$MLL)

# Upper slope habitat
MLL.3<-data.frame(subset(clean.2,Hab=="UP"))
MLL.1<-MLL.3[,c(13:21804)]
MLL.1<-log(MLL.1+1)
adonis(MLL.1~as.factor(MLL), data=MLL.3, permutations = 9999, method = "bray")

MLL.3$MLL <- factor(MLL.3$MLL, levels =c("46","49","113","174","268","370"))
pairwise.adonis(MLL.1, MLL.3$MLL)

# Back reef habitat
MLL.3<-data.frame(subset(clean.2,Hab=="BR"))
MLL.1<-MLL.3[,c(13:21804)]
MLL.1<-log(MLL.1+1)
adonis(MLL.1~as.factor(MLL), data=MLL.3, permutations = 9999, method = "bray")

# Genotype G1
MLL.3<-data.frame(subset(clean.2,MLL=="46"))
MLL.1<-MLL.3[,c(13:21804)]
MLL.1<-log(MLL.1+1)
adonis(MLL.1~Hab, data=MLL.3, permutations = 9999, method = "bray")

MLL.3$Hab <- factor(MLL.3$Hab, levels =c("MD","UP","BR"))
pairwise.adonis(MLL.1, MLL.3$Hab)

# Genotype G2
MLL.3<-data.frame(subset(clean.2,MLL=="49"))
MLL.1<-MLL.3[,c(13:21804)]
MLL.1<-log(MLL.1+1)
adonis(MLL.1~Hab, data=MLL.3, permutations = 9999, method = "bray")

MLL.3$Hab <- factor(MLL.3$Hab, levels =c("MD","UP","BR"))
pairwise.adonis(MLL.1, MLL.3$Hab)

# Genotype G3
MLL.3<-data.frame(subset(clean.4bis,MLL=="113"))
MLL.1<-MLL.3[,c(13:21804)]
MLL.1<-log(MLL.1+1)
adonis(MLL.1~Hab, data=MLL.3, permutations = 9999, method = "bray")

MLL.3$Hab <- factor(MLL.3$Hab, levels =c("MD","UP","BR"))
pairwise.adonis(MLL.1, MLL.3$Hab)

# Genotype G4
MLL.3<-data.frame(subset(clean.2,MLL=="174"))
MLL.1<-MLL.3[,c(13:21804)]
MLL.1<-log(MLL.1+1)
adonis(MLL.1~Hab, data=MLL.3, permutations = 9999, method = "bray")

MLL.3$Hab <- factor(MLL.3$Hab, levels =c("MD","UP","BR"))
pairwise.adonis(MLL.1, MLL.3$Hab)

# Genotype G5
MLL.3<-data.frame(subset(clean.2,MLL=="268"))
MLL.1<-MLL.3[,c(13:21804)]
MLL.1<-log(MLL.1+1)
adonis(MLL.1~Hab, data=MLL.3, permutations = 9999, method = "bray")

MLL.3$Hab <- factor(MLL.3$Hab, levels =c("MD","UP","BR"))
pairwise.adonis(MLL.1, MLL.3$Hab)

# Genotype G6
MLL.3<-data.frame(subset(clean.2,MLL=="370"))
MLL.1<-MLL.3[,c(13:21804)]
MLL.1<-log(MLL.1+1)
adonis(MLL.1~Hab, data=MLL.3, permutations = 9999, method = "bray")

MLL.3$Hab <- factor(MLL.3$Hab, levels =c("MD","UP","BR"))
pairwise.adonis(MLL.1, MLL.3$Hab)

### Barplot Bacterial families =================================================
## Create Table data ASV --> family
taxon<-melt(clean.2,id=1:13)
FamBar <-FamBar[,c(1,7)]
taxon.1<-merge(FamBar,taxon,by.x="OTU",by.y="variable")
taxon.1<-ddply(taxon.1,~Ind+Morpho+Hab+MLL+LifSta+family,function(x){c(value=sum(x$value))})
taxon.2<-ddply(taxon.1,~Ind,function(x){c(value=sum(x$value))})
taxon.1

## Selection 20 most abundant families
taxon.3<-ddply(taxon.1,~family,function(x){c(Sumread=sum(x$value))})
taxon.3$family2<-ifelse(taxon.3$Sumread>sort(taxon.3$Sumread, decreasing=TRUE)[21],as.character(taxon.3$family),"Other")
taxon.1<-merge(taxon.1,taxon.3,by="family")
taxon.1

## Rank most abundant families
Taxon<-ddply(taxon.1,~family2,function(x){c(readtaxon=sum(x$value))})
Taxon <- Taxon[order(Taxon[,2],decreasing=T),1]

## Family order
Taxon$family2 <- factor(Taxon$family2, levels =c("c__Gammaproteobacteria_unclassified","p__Firmicutes_unclassified","f__Spirochaetaceae","c__Alphaproteobacteria_unclassified","o__Thalassobaculales_unclassified","c__Cyanobacteriia_unclassified","f__Rhodobacteraceae","f__Cyclobacteriaceae","f__Kiloniellaceae","f__Amoebophilaceae","f__Rhizobiaceae","f__Flavobacteriaceae","f__Fodinicurvataceae","f__Thalassospiraceae","d__Bacteria_unclassified","f__Unknown_Family","f__Vibrionaceae","f__Sandaracinaceae","f__Microtrichaceae","o__Rickettsiales_unclassified","Other"))

## Barplot
# Barplot par samples
taxon.2<-data.frame(taxon.1[,c(9,2,7)])
taxon.2<-data.frame(acast(taxon.2,Ind~family2,sum))
taxon.2$Ind<-row.names(taxon.2)
taxon.2<-merge(Design_Files,taxon.2,by="Ind",all.y=T,all.x=F)
taxon.2$sumrow<-apply(taxon.2[,c(12:32)], 1, sum)
taxon.3<-(taxon.2[,c(12:32)]/taxon.2$sumrow)*100
taxon.3$Ind<-taxon.2$Ind
row.names(taxon.3)<-taxon.2$Ind
plot<-melt(taxon.3,id=22)
plot<-merge(Design_Files,plot,by="Ind",all.y=T,all.x=F)
Barplot.2<-plot
plot$Ind <- factor(plot$Ind, levels = plot$Ind[order(plot[,4],plot[,6])])
Barplot.2<-plot
Barplot.2$variable<- factor(Barplot.2$variable, levels =c("c__Gammaproteobacteria_unclassified","p__Firmicutes_unclassified","f__Spirochaetaceae","c__Alphaproteobacteria_unclassified","o__Thalassobaculales_unclassified","c__Cyanobacteriia_unclassified","f__Rhodobacteraceae","f__Cyclobacteriaceae","f__Kiloniellaceae","f__Amoebophilaceae","f__Rhizobiaceae","f__Flavobacteriaceae","f__Fodinicurvataceae","f__Thalassospiraceae","d__Bacteria_unclassified","f__Unknown_Family","f__Vibrionaceae","f__Sandaracinaceae","f__Microtrichaceae","o__Rickettsiales_unclassified","Other"))
cbPalette <- c("c__Gammaproteobacteria_unclassified"="#75A5F8","p__Firmicutes_unclassified"="#FE2E64","f__Spirochaetaceae"="#E69F00","c__Alphaproteobacteria_unclassified"="#812599","o__Thalassobaculales_unclassified"="#009E73","c__Cyanobacteriia_unclassified"="#2773F6","f__Rhodobacteraceae"="#E23434","f__Cyclobacteriaceae"="#EE8041","f__Kiloniellaceae"="#999999","f__Amoebophilaceae"="#F9FF33","f__Rhizobiaceae"="#000000","f__Flavobacteriaceae"="cyan","f__Fodinicurvataceae"="#CC79A7","f__Thalassospiraceae"="#258499","d__Bacteria_unclassified"="#82EE41","f__Unknown_Family"="dodgerblue4","f__Vibrionaceae"="#A663D1","f__Sandaracinaceae"="#FF5FBB","f__Microtrichaceae"="#173B0B","o__Rickettsiales_unclassified"="#900C3F","Other"="#E6E6E6")

# Barplot per sample --> mean
plot.2<-ddply(plot,~MLL+Hab+variable,function(x){c(value=mean(x$value),n=nrow(x))})
Barplot.2<-plot.2
Barplot.2$Ind<-paste(Barplot.2$MLL,Barplot.2$Hab)
Barplot.2$variable<- factor(Barplot.2$variable, levels =c("c__Gammaproteobacteria_unclassified","p__Firmicutes_unclassified","f__Spirochaetaceae","c__Alphaproteobacteria_unclassified","o__Thalassobaculales_unclassified","c__Cyanobacteriia_unclassified","f__Rhodobacteraceae","f__Cyclobacteriaceae","f__Kiloniellaceae","f__Amoebophilaceae","f__Rhizobiaceae","f__Flavobacteriaceae","f__Fodinicurvataceae","f__Thalassospiraceae","d__Bacteria_unclassified","f__Unknown_Family","f__Vibrionaceae","f__Sandaracinaceae","f__Microtrichaceae","o__Rickettsiales_unclassified","Other"))
Barplot.2$Ind<- factor(Barplot.2$Ind, levels =c("46 MD","46 UP","46 BR","49 MD","49 UP","49 BR","113 MD","113 UP","113 BR","174 UP","174 BR","268 UP","268 BR","370 MD","370 UP"))
barplot_family<-ggplot(Barplot.2)+geom_histogram(aes(x=Ind,y=value,fill=variable,order=variable,width=.75),stat="identity",position = position_stack(reverse = TRUE))+theme_bw() + guides(fill = guide_legend(keywidth = 1.8, keyheight = 1.2))+scale_fill_manual(values=cbPalette)+geom_text(data=Barplot.2,aes(x=Ind,y=105,label=paste("n=",n)))
barplot_family
ggsave(barplot_family, height = 8, width = 14, filename=paste("sample","bacterie3",".pdf",sep=""),path = "~/Desktop/Research/Microbiome/Review/R_analyses/Family2")

### Bacterial genera ===========================================================
## Create Table data ASV --> genus
taxon<-melt(clean.2,id=1:13)
FamBar.g <- read.csv("~/Desktop/Research/Microbiome/Review/qiime2/taxonomy_file.csv", sep=";")
FamBar.g <-FamBar.g[,c(1,8)]
taxon.1g<-merge(FamBar.g,taxon,by.x="OTU",by.y="variable")
taxon.1g<-ddply(taxon.1g,~Ind+Morpho+Hab+MLL+LifSta+genus,function(x){c(value=sum(x$value))})
taxon.2g<-ddply(taxon.1g,~Ind,function(x){c(value=sum(x$value))})
taxon.1g

## Selection 20 most abundant genera
taxon.3g<-ddply(taxon.1g,~genus,function(x){c(Sumread=sum(x$value))})
taxon.3g$genus2<-ifelse(taxon.3g$Sumread>sort(taxon.3g$Sumread, decreasing=TRUE)[21],as.character(taxon.3g$genus),"Other")
taxon.1g<-merge(taxon.1g,taxon.3g,by="genus")
taxon.1g

## Rank most abundant genera
Taxon.g<-ddply(taxon.1g,~genus2,function(x){c(readtaxon=sum(x$value))})
Taxon.g <- Taxon.g[order(Taxon.g[,2],decreasing=T),1]

### Bacterial phylum ===========================================================
## Create Table data ASV --> phyla
taxon<-melt(clean.2,id=1:13)
FamBar.p <- read.csv("~/Desktop/Research/Microbiome/Review/qiime2/taxonomy_file.csv", sep=";")
FamBar.p <-FamBar.g[,c(1,4)]
taxon.1p<-merge(FamBar.p,taxon,by.x="OTU",by.y="variable")
taxon.1p<-ddply(taxon.1p,~Ind+Morpho+Hab+MLL+LifSta+phylum,function(x){c(value=sum(x$value))})
taxon.2p<-ddply(taxon.1p,~Ind,function(x){c(value=sum(x$value))})
taxon.1p

## Selection 20 must abundant phyla
taxon.3p<-ddply(taxon.1p,~phylum,function(x){c(Sumread=sum(x$value))})
taxon.3p$phylum2<-ifelse(taxon.3p$Sumread>sort(taxon.3p$Sumread, decreasing=TRUE)[21],as.character(taxon.3p$phylum),"Other")
taxon.1p<-merge(taxon.1p,taxon.3p,by="phylum")
taxon.1g

## Rank most abundant phyla
Taxon.p<-ddply(taxon.1p,~phylum2,function(x){c(readtaxon=sum(x$value))})
Taxon.g <- Taxon.g[order(Taxon.g[,2],decreasing=T),1]

Phyl<-ddply(Taxon.p,~phylum2,function(x){c(sum=sum(x$readtaxon))})
Phyl<-ddply(Phyl,~phylum2,function(x){c(pourcent=(x$readtaxon))})

### Bacterial community composition analysis family level ======================
## One-way Permanova all families
# By genotype
taxon.a<-data.frame(taxon.1)
taxon.a<-data.frame(acast(taxon.a,Ind~family,sum))
taxon.a2<-log(taxon.a+1)
taxon.a$sumrow<-apply(taxon.a[,c(1:505)], 1, sum)
taxon.b<-data.frame((taxon.a[,c(1:505)]/taxon.a$sumrow)*100)
taxon.b$Ind<-row.names(taxon.b)
taxon.c<-merge(Design_Files,taxon.b,by="Ind",all.y=T,all.x=F)
row.names(taxon.c)<-taxon.b$Ind
# on percentage 
adonis(taxon.b~MLL, data=taxon.c, permutations = 9999, method = "bray")

# Pairwise comparison most abundant families
pairwise.adonis(taxon.b[,c(1:505)],taxon.c$MLL)
pairwise.adonis(taxon.a2,taxon.c$MLL)

# By habitat
# on percentage 
adonis(taxon.b~Hab, data=taxon.c, permutations = 9999, method = "bray")

# Pairwise comparison most abundant families
pairwise.adonis(taxon.b[,c(1:505)],taxon.c$Hab)
pairwise.adonis(taxon.a2,taxon.c$Hab)

## Betadisper families
taxon.b2<-data.frame(sapply(taxon.b,function(x) as.numeric(x)))
dis<-vegdist(taxon.b2, na.rm = TRUE)

groups_gen<- factor(c(rep(1,38), rep(2,16), rep(3,15), rep(4,32), rep(5,13), rep(6,21)), labels = c("46","49","113","174","268","370")) 
betagen<-betadisper(dis, groups_gen, type=c("median"), sqrt.dist=FALSE, bias.adjust=TRUE)

plot(betagen, ellipse=TRUE, hull=FALSE)
boxplot(betagen, ellipse=TRUE, hull=FALSE)

anova(betagen)
permutest(betagen, pairwise=TRUE, permutations=9999)

groups_hab<- factor(c(rep(7,19), rep(8,90), rep(9,26)), labels = c("BR", "UP", "MD"))
betahab<-betadisper(dis, groups_hab, type=c("median"), sqrt.dist=FALSE, bias.adjust=TRUE)

plot(betahab, ellipse=TRUE, hull=FALSE)
boxplot(betahab, ellipse=TRUE, hull=FALSE)

anova(betahab)
permutest(betagen, pairwise=TRUE, permutations=9999)

## Simper tests 
# Genotype within BR 
taxon.cBR<-data.frame(subset(taxon.c,Hab=="BR"))
taxon.cBR$MLL<-factor(taxon.cBR$MLL, levels =c("46","49","113","174","268"))
taxon.bBR<-data.frame(subset(taxon.cBR[,c(12:516)]))
adonis(taxon.bBR~as.factor(MLL), data=taxon.cBR, permutations = 9999, method = "bray")

sim <- with(taxon.cBR, simper(taxon.cBR[,c(12:516)], MLL, permutations = 9999))
options(max.print=1000000)
s<-summary(sim)
s
capture.output(s, file = "summarysimperfamilies%_GEN-BR.txt")

# Genotype within MD
taxon.cMD<-data.frame(subset(taxon.c,Hab=="MD"))
taxon.cMD$MLL<-factor(taxon.cMD$MLL, levels =c("46","49","113","370"))
taxon.bMD<-data.frame(subset(taxon.cMD[,c(12:516)]))
adonis(taxon.bMD~as.factor(MLL), data=taxon.cMD, permutations = 9999, method = "bray")

taxon.cMD<-subset(taxon.cMD, MLL%in%c("49","370"))
sim <- with(taxon.cMD, simper(taxon.cMD[,c(12:516)], MLL, permutations = 9999))
options(max.print=1000000)
s<-summary(sim)
s
capture.output(s, file = "summarysimperfamilies%_GEN-MD.txt")

kruskal.test(f__Roseiflexaceae ~ MLL, data=taxon.cMD)

# Genotype within UP
taxon.cUP<-data.frame(subset(taxon.c,Hab=="UP"))
taxon.cUP$MLL<-factor(taxon.cUP$MLL, levels =c("46","49","113","174","268","370"))
taxon.bUP<-data.frame(subset(taxon.cUP[,c(12:516)]))
adonis(taxon.bUP~as.factor(MLL), data=taxon.cUP, permutations = 9999, method = "bray")
pairwise.adonis(taxon.bUP,taxon.cUP$MLL)

simpcom<-simper(comm =taxon.cUP[,c(12:516)], group = taxon.cUP$MLL)
simpcom
summary(simpcom, ordered=TRUE, digits=3)
lapply(simpcom, FUN=function(x){x$overall})

sim <- with(taxon.cUP, simper(taxon.cUP[,c(12:516)], MLL, permutations = 9999))
options(max.print=1000000)
s<-summary(sim)
s
capture.output(s, file = "summarysimperfamilies%_GEN-UP.txt")

taxon.cUP<-subset(taxon.cUP, MLL%in%c("46","113")) 
taxon.cUP<-subset(taxon.cUP, MLL%in%c("46","174")) 
taxon.cUP<-subset(taxon.cUP, MLL%in%c("46","370")) 
taxon.cUP<-subset(taxon.cUP, MLL%in%c("49","113")) 
taxon.cUP<-subset(taxon.cUP, MLL%in%c("49","174")) 
taxon.cUP<-subset(taxon.cUP, MLL%in%c("49","370")) 
taxon.cUP<-subset(taxon.cUP, MLL%in%c("113","174")) 
taxon.cUP<-subset(taxon.cUP, MLL%in%c("113","268")) 
taxon.cUP<-subset(taxon.cUP, MLL%in%c("113","370")) 
taxon.cUP<-subset(taxon.cUP, MLL%in%c("174","370")) 
taxon.cUP<-subset(taxon.cUP, MLL%in%c("268","370")) 

kruskal.test(f__Lineage_IIb~ MLL, data=taxon.cUP)

# Habitat per genotype
# Hab for Genotype G1
taxon.c46<-data.frame(subset(taxon.c,MLL=="46"))
taxon.c46$Hab<-factor(taxon.cUP$Hab, levels =c("BR","MD","UP"))
taxon.b46<-data.frame(subset(taxon.c46[,c(12:516)]))
adonis(taxon.b46~as.factor(Hab), data=taxon.c46, permutations = 9999, method = "bray")
pairwise.adonis(taxon.b46,taxon.c46$Hab)

sim <- with(taxon.c46, simper(taxon.c46[,c(12:516)], Hab, permutations = 9999))
lapply(sim, FUN=function(x){x$overall})
options(max.print=1000000)
s<-summary(sim)
s
capture.output(s, file = "summarysimperfamilies%_HAB-46.txt")

# Hab for Genotype G2
taxon.c49<-data.frame(subset(taxon.c,MLL=="49"))
taxon.c49$Hab<-factor(taxon.c49$Hab, levels =c("BR","MD","UP"))
taxon.b49<-data.frame(subset(taxon.c49[,c(12:516)]))
adonis(taxon.b49~as.factor(Hab), data=taxon.c49, permutations = 9999, method = "bray")
pairwise.adonis(taxon.b49,taxon.c49$Hab)

sim <- with(taxon.c49, simper(taxon.c49[,c(12:516)], Hab, permutations = 999))
lapply(sim, FUN=function(x){x$overall})
options(max.print=1000000)
s<-summary(sim)
s
capture.output(s, file = "summarysimperfamilies%_HAB-49.txt")

# Hab for Genotype G3
taxon.c113<-data.frame(subset(taxon.c,MLL=="113"))
taxon.c113$Hab<-factor(taxon.c113$Hab, levels =c("BR","MD","UP"))
taxon.b113<-data.frame(subset(taxon.c113[,c(12:516)]))
adonis(taxon.b113~as.factor(Hab), data=taxon.c113, permutations = 9999, method = "bray")
pairwise.adonis(taxon.b113,taxon.c113$Hab)

# Hab for Genotype G4
taxon.c174<-data.frame(subset(taxon.c,MLL=="174"))
taxon.c174$Hab<-factor(taxon.c174$Hab, levels =c("BR","MD","UP"))
taxon.b174<-data.frame(subset(taxon.c174[,c(12:516)]))
adonis(taxon.b174~as.factor(Hab), data=taxon.c174, permutations = 9999, method = "bray")
pairwise.adonis(taxon.b174,taxon.c174$Hab)

# Hab for Genotype G5
taxon.c268<-data.frame(subset(taxon.c,MLL=="268"))
taxon.c268$Hab<-factor(taxon.c268$Hab, levels =c("BR","MD","UP"))
taxon.b268<-data.frame(subset(taxon.c268[,c(12:516)]))
adonis(taxon.b268~as.factor(Hab), data=taxon.c268, permutations = 9999, method = "bray")
pairwise.adonis(taxon.b268,taxon.c268$Hab)

# Hab for Genotype G6
taxon.c370<-data.frame(subset(taxon.c,MLL=="370"))
taxon.c370$Hab<-factor(taxon.c370$Hab, levels =c("BR","MD","UP"))
taxon.b370<-data.frame(subset(taxon.c370[,c(12:516)]))
adonis(taxon.b370~as.factor(Hab), data=taxon.c370, permutations = 9999, method = "bray")
pairwise.adonis(taxon.b370,taxon.c370$Hab)

sim <- with(taxon.c370, simper(taxon.c370[,c(12:516)], Hab, permutations = 999))
lapply(sim, FUN=function(x){x$overall})
options(max.print=1000000)
s<-summary(sim)
s
capture.output(s, file = "summarysimperfamilies%_HAB-370.txt")

### Heatmap for bacterial indicator taxa =======================================
## Create files for ASV alpha 0.01
# Habitat
Multi.1a<-clean.2
Multi.3a<-data.frame(Multi.1a[,c(13:21804)])
# Selection ASV specific
cluster<-factor(Multi.1a$Hab)
Multi.2a<-multipatt(Multi.3a,cluster)
capture.output(Multi.2a, file = "ASVSpe_HAB.txt")
summary(Multi.2a,alpha=0.01)
Multi.4a<-data.frame(Multi.2a$sign)
Multi.4a<-subset(Multi.4a,s.MD+s.UP+s.BR==1)
Multi.4a<-subset(Multi.4a,p.value<=0.01)
write.csv2(Multi.4a,"ASVSpe_HABsig.csv")
Multi.4a$OTU<-rownames(Multi.4a)
Multi.4a<-melt(Multi.4a[,c(7,1:3)],id=1)
Multi.4a<-subset(Multi.4a,value==1)
names(Multi.4a)[2]<-"Habspe"
OTU<-factor(Multi.4a$OTU)
write.csv2(Multi.4a,"OTUSpecific_Habitat.csv")

# Creation heatmap
OTUHab <- read.csv("~/Desktop/Research/Microbiome/Review/For figures/Figure 3/OTUSpecific_Habitat.csv", header=TRUE, sep=";")
OTU<-OTUHab[,c(1)]
clean.3<-data.frame(clean.2[,c(13:21804)])
clean.3<-data.frame(scale(clean.3, center = TRUE, scale = TRUE))
clean.3$Ind<-clean.2[,c(2)]
heat<-melt(clean.3,id="Ind")
heat<-subset(heat,variable%in%OTU)
OTU.2<-factor(OTUHab$OTU)
heat$OTU2 <- factor(heat$variable, levels = OTU.2)
heat<-merge(heat,Design_Files,by="Ind",all.x=T)
heat<-merge(heat,OTUHab,by.x="variable",by.y="OTU",all.x=T)
x<-heat[order(heat[,9],decreasing=F),]
MLL.2<-factor(x$Ind)
heat$Ind2 <- factor(heat$Ind, levels = unique(MLL.2))
ggplot(heat, aes(x = Ind2, y = OTU2)) +geom_tile(aes(fill = value))+scale_fill_continuous(low="white",high="red",na.value = "white")
ggplot(heat, aes(x = Ind2, y = OTU2)) +geom_tile(aes(fill = zscore))+scale_fill_continuous(low="black",high="#red",na.value = "black")

# Genotype
Multi.1<-clean.2
Multi.3<-data.frame(Multi.1[,c(13:21804)])
# selection ASV specific
cluster<-factor(Multi.1$MLL)
Multi.2<-multipatt(Multi.3,cluster)
write.csv2(Multi.2,"OTUSpecific_Genotype.csv")
summary(Multi.2,alpha=0.01)
Multi.4<-data.frame(Multi.2$sign)
Multi.4<-subset(Multi.4,s.46+s.49+s.113+s.174+s.268+s.370==1)
Multi.4<-subset(Multi.4,p.value<=0.01)
write.csv2(Multi.4,"ASVSpe_GENsig.csv")
Multi.4$OTU<-rownames(Multi.4)
Multi.4<-melt(Multi.4[,c(10,1:6)],id=1)
Multi.4<-subset(Multi.4,value==1)
names(Multi.4)[2]<-"Genspe"
OTU<-factor(Multi.4$OTU)
write.csv2(Multi.4,"OTUSpecific_Genotype.csv")

# Create heatmap
OTUGen <- read.csv("~/Desktop/Research/Microbiome/Review/For figures/Figure 3/OTUSpecific_Genotype.csv", header=TRUE, sep=";")
OTU<-OTUGen[,c(1)]
clean.3<-data.frame(clean.2[,c(13:21804)])
clean.3<-data.frame(scale(clean.3, center = TRUE, scale = TRUE))
clean.3$Ind<-clean.2[,c(2)]
heat<-melt(clean.3,id="Ind")
heat<-subset(heat,variable%in%OTU)
OTU.2<-factor(OTUGen$OTU)
heat$OTU2 <- factor(heat$variable, levels = OTU.2)
heat<-merge(heat,Design_Files,by="Ind",all.x=T)
heat<-merge(heat,OTUGen,by.x="variable",by.y="OTU",all.x=T)
x<-heat[order(heat[,7],decreasing=F),]
MLL.2<-factor(x$Ind)
heat$Ind2 <- factor(heat$Ind, levels = unique(MLL.2))
ggplot(heat, aes(x = Ind2, y = OTU2)) +geom_tile(aes(fill = value))+scale_fill_continuous(low="white",high="red",na.value = "white")
ggplot(heat, aes(x = Ind2, y = OTU2)) +geom_tile(aes(fill = zscore))+scale_fill_continuous(low="white",high="red",na.value = "white")

### Families associated to bacterial indicator taxa ============================
##Select par hab
Select<-clean.2
Select<-data.frame(Select[,c(7,13:21804)])
Select<-data.frame(melt(Select,id=1))
Select<-ddply(Select,~Hab+variable,function(x){c(value=sum(x$value))})
Selecttot<-ddply(Select,~Hab,function(x){c(value=sum(x$value))})
Select<-merge(Select,Selecttot,by="Hab")

OTUHab2<-OTUHab[,c(1,5)]
names(OTUHab2)[1]<-"OTU"
names(OTUHab2)[2]<-"Habspe"

OTUHab2<-merge(OTUHab2,FamBarSpe,by.x="OTU",by.y="OTU",x.all=T)
OTUHab2<-merge(OTUHab2,Select,by.x=c("OTU","Habspe"),by.y=c("variable","Hab"),x.all=T)
OTUHab2$NbrOTU<-1
write.csv2(OTUHab2,"Habspe_TAXO.csv")
OTUHab3<-ddply(OTUHab2,~Habspe+family,function(x){c(sumOTU=sum(x$NbrOTU),sumread=sum(x$value.x))})
OTUHab3<-merge(OTUHab3,Selecttot,by.x="Habspe",by.y="Hab",all.x=T)
OTUHab3$percent<-OTUHab3$sumread/OTUHab3$value*100
OTUHabtot<-ddply(OTUHab3,~Habspe,function(x){c(sumreadspe=sum(x$sumread))})
OTUHab3<-merge(OTUHab3,OTUHabtot,by="Habspe",all.x=T)
OTUHab3$percentspe<-OTUHab3$sumread/OTUHab3$sumreadspe*100
write.csv2(OTUHab3,"OTUheatmapHabspe.csv")

OTUHab3$tax<-ifelse(OTUHab3$percentspe>2,as.character(OTUHab3[,2]),"Otherspe")

OTUHabvf<-ddply(OTUHab3,~Habspe+tax,function(x){c(sumOTU=sum(x$sumOTU),sumread=sum(x$sumread),percent=sum(x$percent),percentspe=sum(x$percentspe))})
OTUHabvf$tax <- factor(OTUHabvf$tax)
x<-ddply(OTUHabvf,~tax,function(x)c(sum=sum(x$percent)))
x <- x[order(x[,2],decreasing=F),] 
TAX<-factor(x$tax)
OTUHabvf$tax <- factor(OTUHabvf$tax ,levels=as.character(TAX))

OTUHabvf$tax <- factor(OTUHabvf$tax , levels =c("f__Rhodobacteraceae","c__Alphaproteobacteria_unclassified","f__Flavobacteriaceae","o__Thalassobaculales_unclassified","f__Sandaracinaceae","c__Cyanobacteriia_unclassified","c__Gammaproteobacteria_unclassified","f__Cyclobacteriaceae","f__Halieaceae","f__Parvularculaceae","f__Kiloniellaceae","f__Saprospiraceae","f__Rhizobiaceae","c__Bacteroidia_unclassified","f__Defluviicoccaceae","f__Cryomorphaceae","Otherspe"))
cbPalette3 <- c("f__Cyclobacteriaceae"="#EE8041","f__Rhizobiaceae"="#000000","o__Thalassobaculales_unclassified"="#009E73","f__Rhodobacteraceae"="#E23434","f__Sandaracinaceae"="#FF5FBB","f__Cryomorphaceae"="thistle1","c__Alphaproteobacteria_unclassified"="#812599","f__Kiloniellaceae"="#999999","c__Bacteroidia_unclassified"="gold1","f__Flavobacteriaceae"="cyan","c__Cyanobacteriia_unclassified"="#2773F6","c__Gammaproteobacteria_unclassified"="#75A5F8","f__Halieaceae"="#00CCCC","f__Parvularculaceae"="violet","f__Saprospiraceae"="red4","f__Defluviicoccaceae"="darkseagreen1","Otherspe"="#E6E6E6")

OTUHabvf$Habspe<-factor(OTUHabvf$Habspe,levels=c("BR","MD","UP")) 

OTUHabFamily<-ggplot(OTUHabvf)+
  geom_bar(aes(x=Habspe,y=percent, fill=tax),stat="identity",position = position_stack(reverse = TRUE))+
  scale_fill_manual(values=cbPalette3)+theme_bw()+theme(panel.grid.major=element_line(colour="white"),panel.grid.minor=element_line(colour="white"),strip.background = element_rect(fill="white")) + guides(fill = guide_legend(keywidth = 1.8, keyheight = 1.2)) 
OTUHabFamily
ggsave(OTUHabFamily,filename=paste("HabSpe2",".pdf",sep=""),path = "~/Desktop/Research/Microbiome/Review/For figures/Figure 3")

## Selection by genotype
Select<-clean.2
Select<-data.frame(Select[,c(5,13:21804)])
Select<-data.frame(melt(Select,id=1))
Select<-ddply(Select,~MLL+variable,function(x){c(value=sum(x$value))})
Selecttot<-ddply(Select,~MLL,function(x){c(value=sum(x$value))})
Select<-merge(Select,Selecttot,by="MLL")

OTUGen2<-OTUGen[,c(1,5)]
names(OTUGen2)[1]<-"OTU"
names(OTUGen2)[2]<-"Genspe"

FamBarSpe <- read.csv("~/Desktop/Research/Microbiome/Review/qiime2/taxonomy_file.csv", sep=";")
FamBarSpe <-FamBarSpe[,c(1,5:8)]

OTUGen2<-merge(OTUGen2,FamBarSpe,by.x="OTU",by.y="OTU",x.all=T)
OTUGen2<-merge(OTUGen2,Select,by.x=c("OTU","Genspe"),by.y=c("variable","MLL"),x.all=T)
OTUGen2$NbrOTU<-1
OTUGen3<-ddply(OTUGen2,~Genspe+family,function(x){c(sumOTU=sum(x$NbrOTU),sumread=sum(x$value.x))})
OTUGen3<-merge(OTUGen3,Selecttot,by.x="Genspe",by.y="MLL",all.x=T)
OTUGen3$percent<-OTUGen3$sumread/OTUGen3$value*100
OTUGentot<-ddply(OTUGen3,~Genspe,function(x){c(sumreadspe=sum(x$sumread))})
OTUGen3<-merge(OTUGen3,OTUGentot,by="Genspe",all.x=T)
OTUGen3$percentspe<-OTUGen3$sumread/OTUGen3$sumreadspe*100
write.csv2(OTUGen3,"OTUheatmapGenspefamille.csv")

cbPalette2 <- c("c__Gammaproteobacteria_unclassified"="#75A5F8","f__Spirochaetaceae"="#E69F00","c__Alphaproteobacteria_unclassified"="#812599","o__Thalassobaculales_unclassified"="#009E73","f__Rhodobacteraceae"="#E23434","f__Cyclobacteriaceae"="#EE8041","f__Kiloniellaceae"="#999999","f__Flavobacteriaceae"="cyan","f__Unknown_Family"="dodgerblue4","f__Sandaracinaceae"="#FF5FBB","o__Rickettsiales_unclassified"="#900C3F","f__Balneolaceae"="#79B123","c__Bacteroidia_unclassified"="gold1","f__Woeseiaceae"="#A663D1","f__Alteromonadaceae"="#1313E3","f__Arcobacteraceae"="#86CDEC","f__Saccharospirillaceae"="sienna3","f__Brevibacteriaceae"="#CBFFB0","f__Candidatus_Kaiserbacteria"="darkslateblue","f__Dermabacteraceae"="#C6BD09","p__Planctomycetota_unclassified"="lightslategray","f__Dietziaceae"="purple4","f__Halieaceae"="#00CCCC","f__Simkaniaceae"="#E31384","f__NB1-j"="lightgoldenrodyellow","f__Marinomonadaceae"="#F9D6E0","o__Defluviicoccales_unclassified"="#CC9999","Otherspe"="#E6E6E6")
OTUGen3$tax<-ifelse(OTUGen3$percentspe>2,as.character(OTUGen3[,2]),"Otherspe")

OTUGenvf<-ddply(OTUGen3,~Genspe+tax,function(x){c(sumOTU=sum(x$sumOTU),sumread=sum(x$sumread),percent=sum(x$percent),percentspe=sum(x$percentspe))})
OTUGenvf$tax <- factor(OTUGenvf$tax)
x<-ddply(OTUGenvf,~tax,function(x)c(sum=sum(x$percent)))
x <- x[order(x[,2],decreasing=F),] 
TAX<-factor(x$tax)
OTUGenvf$tax <- factor(OTUGenvf$tax ,levels=as.character(TAX))
OTUGenvf$tax <- factor(OTUGenvf$tax , levels =c("f__Brevibacteriaceae","f__Halieaceae","f__Sandaracinaceae","f__Spirochaetaceae","c__Alphaproteobacteria_unclassified","f__Kiloniellaceae","o__Rickettsiales_unclassified","f__Alteromonadaceae","c__Gammaproteobacteria_unclassified","f__Rhodobacteraceae","f__Simkaniaceae","o__Defluviicoccales_unclassified","f__NB1-j","o__Thalassobaculales_unclassified","f__Unknown_Family","f__Cyclobacteriaceae","f__Dietziaceae","f__Flavobacteriaceae","f__Dermabacteraceae","f__Woeseiaceae","c__Bacteroidia_unclassified","f__Marinomonadaceae","p__Planctomycetota_unclassified","f__Candidatus_Kaiserbacteria","f__Arcobacteraceae","f__Saccharospirillaceae","f__Balneolaceae","Otherspe"))
OTUGenvf$Genspe<-factor(OTUGenvf$Genspe,levels=c("46","49","113","174","268","370"))

OTUGenFamily<-ggplot(OTUGenvf)+
  geom_bar(aes(x=Genspe,y=percent, fill=tax),stat="identity",position = position_stack(reverse = TRUE))+
  scale_fill_manual(values=cbPalette2)+theme_bw()+theme(panel.grid.major=element_line(colour="white"),panel.grid.minor=element_line(colour="white"),strip.background = element_rect(fill="white")) + guides(fill = guide_legend(keywidth = 1.8, keyheight = 1.2)) 
OTUGenFamily
ggsave(OTUGenFamily, height = 8, width = 13, filename=paste("GenSpe2",".pdf",sep=""),path = "~/Desktop/Research/Microbiome/Review/For figures/Figure 3")

### Colony size analysis =======================================================
## Two-way Permanova 
# all genotypes
sizeALL<-subset(clean.2)
sizeALL$SizeClass<-cut(sizeALL$Area, breaks=c(0.01,3.14,12.56,28.26,50.24,78.5,113.04,153.86,200.96,254.34,314,379.94,452.16,530.66,615.44,706.5,803.84,907.46,1017.36,1133.54,1256,2826,5024,14000))
sizeALL$SizeClass_f = factor(sizeALL$SizeClass, levels=c('(3.14,12.6]','(12.6,28.3]','(28.3,50.2]','(50.2,78.5]','(78.5,113]','(113,154]','(154,201]','(201,254]','(254,314]','(314,380]','(380,452]','(452,531]','(531,615]','(615,706]','(706,804]','(804,907]','(907,1.02e+03]','(1.02e+03,1.13e+03]','(1.13e+03,1.26e+03]','(1.26e+03,2.83e+03]','(2.83e+03,5.02e+03]','(5.02e+03,1.4e+04]'))
sizeALL$sumrow<-apply(sizeALL[,c(13:21804)], 1, sum)
sizeALLclean<-sizeALL[,c(13:21804)]/sizeALL[,21808]*100
adonis(sizeALLclean~SizeClass_f, data=sizeALL, permutations = 9999, method = "bray")
pairwise.adonis(sizeALLclean,sizeALL$SizeClass_f)

# Test G1 Adult BR/UP
size46<-subset(clean.2,MLL==46)
size46$SizeClass<-cut(size46$Area, breaks=c(0.01,3.14,12.56,28.26,50.24,153.86,200.96,254.34,314,379.94,452.16,530.66,615.44,706.5,803.84,907.46,1017.36,1133.54,1256,2826,5024,14000))
size46$SizeClass_f = factor(size46$SizeClass, levels=c('(3.14,12.6]','(12.6,28.3]','(28.3,50.2]','(50.2,154]','(154,201]','(201,254]','(254,314]','(314,380]','(380,452]','(452,531]','(531,615]','(615,706]','(706,804]','(804,907]','(907,1.02e+03]','(1.02e+03,1.13e+03]','(1.13e+03,1.26e+03]','(1.26e+03,2.83e+03]','(2.83e+03,5.02e+03]','(5.02e+03,1.4e+04]'))
size46$sumrow<-apply(size46[,c(13:21804)], 1, sum)
size46clean<-size46[,c(13:21804)]/size46[,21808]*100
adonis(size46clean~SizeClass_f, data=size46, permutations = 9999, method = "bray")
pairwise.adonis(size46clean,size46$SizeClass_f)

# Test G2 Adult MD/UP
size49<-subset(clean.2,MLL==49)
size49$SizeClass<-cut(size49$Area, breaks=c(0.01,3.14,12.56,28.26,50.24,78.5,113.04,153.86,200.96,254.34,314,379.94,452.16,530.66,615.44,706.5,803.84,907.46,1017.36,1133.54,1256,2826,5024,14000))
size49$SizeClass_f = factor(size49$SizeClass, levels=c('(3.14,12.6]','(12.6,28.3]','(28.3,50.2]','(50.2,78.5]','(78.5,113]','(113,154]','(154,201]','(201,254]','(254,314]','(314,380]','(380,452]','(452,531]','(531,615]','(615,706]','(706,804]','(804,907]','(907,1.02e+03]','(1.02e+03,1.13e+03]','(1.13e+03,1.26e+03]','(1.26e+03,2.83e+03]','(2.83e+03,5.02e+03]','(5.02e+03,1.4e+04]'))
size49$sumrow<-apply(size49[,c(13:21804)], 1, sum)
size49clean<-size49[,c(13:21804)]/size49[,21808]*100
adonis(size49clean~SizeClass_f, data=size49, permutations = 9999, method = "bray")
pairwise.adonis(size49clean,size49$SizeClass_f)

# Test G3 MD/UP Life Stage
size113<-subset(clean.2,MLL==113)
size113$SizeClass<-cut(size113$Area, breaks=c(0.01,3.14,12.56,28.26,50.24,78.5,113.04,153.86,200.96,254.34,314,379.94,452.16,530.66,615.44,706.5,803.84,907.46,1017.36,1133.54,1256,2826,5024,14000))
size113$SizeClass_f = factor(size113$SizeClass, levels=c('(3.14,12.6]','(12.6,28.3]','(28.3,50.2]','(50.2,78.5]','(78.5,113]','(113,154]','(154,201]','(201,254]','(254,314]','(314,380]','(380,452]','(452,531]','(531,615]','(615,706]','(706,804]','(804,907]','(907,1.02e+03]','(1.02e+03,1.13e+03]','(1.13e+03,1.26e+03]','(1.26e+03,2.83e+03]','(2.83e+03,5.02e+03]','(5.02e+03,1.4e+04]'))
size113$Group<-paste(size113$Hab,size113$SizeClass)
size113$sumrow<-apply(size113[,c(13:21804)], 1, sum)
size113clean<-size113[,c(13:21804)]/size113[,21808]*100
adonis(size113clean~SizeClass_f, data=size113, permutations = 9999, method = "bray")
pairwise.adonis(size113clean,size113$SizeClass_f)

# Test G4 BR/UP Life Stage
size174<-subset(clean.2,MLL==174)
size174$SizeClass<-cut(size174$Area, breaks=c(0.01,3.14,12.56,28.26,50,78.5,113.04,153.86,200.96,254.34,314,379.94,452.16,530.66,615.44,706.5,803.84,907.46,1017.36,1133.54,1256,2826,5024,14000))
size174$SizeClass_f = factor(size174$SizeClass, levels=c('(3.14,12.6]','(12.6,28.3]','(28.3,50]','(50,78.5]','(78.5,113]','(113,154]','(154,201]','(201,254]','(254,314]','(314,380]','(380,452]','(452,531]','(531,615]','(615,706]','(706,804]','(804,907]','(907,1.02e+03]','(1.02e+03,1.13e+03]','(1.13e+03,1.26e+03]','(1.26e+03,2.83e+03]','(2.83e+03,5.02e+03]','(5.02e+03,1.4e+04]'))
size174$sumrow<-apply(size174[,c(13:21804)], 1, sum)
size174clean<-size174[,c(13:21804)]/size174[,21808]*100
adonis(size174clean~SizeClass_f, data=size174, permutations = 9999, method = "bray")
pairwise.adonis(size174clean,size174$SizeClass_f)

# Test G5 Adult MD/UP
size268<-subset(clean.4bis,MLL==268)
size268$SizeClass<-cut(size268$Area, breaks=c(0.01,3.14,12.56,28.26,50.24,78.5,113.04,153.86,200.96,254.34,314,379.94,452.16,530.66,615.44,706.5,803.84,907.46,1017.36,1133.54,1256,2826,5024,14000))
size268$SizeClass_f = factor(size268$SizeClass, levels=c('(3.14,12.6]','(12.6,28.3]','(28.3,50.2]','(50.2,78.5]','(78.5,113]','(113,154]','(154,201]','(201,254]','(254,314]','(314,380]','(380,452]','(452,531]','(531,615]','(615,706]','(706,804]','(804,907]','(907,1.02e+03]','(1.02e+03,1.13e+03]','(1.13e+03,1.26e+03]','(1.26e+03,2.83e+03]','(2.83e+03,5.02e+03]','(5.02e+03,1.4e+04]'))
size268$sumrow<-apply(size268[,c(13:21804)], 1, sum)
size268clean<-size268[,c(13:21804)]/size268[,21808]*100
adonis(size268clean~SizeClass_f, data=size268, permutations = 9999, method = "bray")
pairwise.adonis(size268clean,size268$SizeClass_f)

# Test G6 UP Life Stage
size370<-subset(clean.4bis,MLL==370)
size370$SizeClass<-cut(size370$Area, breaks=c(0.01,3.14,12.56,28.26,50.24,78.5,113.04,153.86,200.96,254.34,314,379.94,452.16,530.66,615.44,706.5,803.84,907.46,1017.36,1133.54,1256,2826,5024,14000))
size370$SizeClass_f = factor(size370$SizeClass, levels=c('(3.14,12.6]','(12.6,28.3]','(28.3,50.2]','(50.2,78.5]','(78.5,113]','(113,154]','(154,201]','(201,254]','(254,314]','(314,380]','(380,452]','(452,531]','(531,615]','(615,706]','(706,804]','(804,907]','(907,1.02e+03]','(1.02e+03,1.13e+03]','(1.13e+03,1.26e+03]','(1.26e+03,2.83e+03]','(2.83e+03,5.02e+03]','(5.02e+03,1.4e+04]'))
size370$sumrow<-apply(size370[,c(13:21804)], 1, sum)
size370clean<-size370[,c(13:21804)]/size370[,21808]*100
adonis(size370clean~SizeClass_f, data=size370, permutations = 9999, method = "bray")
pairwise.adonis(size370clean,size370$SizeClass_f)

## Correlation number ASV and colony size
CorrSize <- read.csv("~/Desktop/Research/Microbiome/Review/CorrSize.csv", sep=";")
CorrSize2<-merge(CorrSize,size,by.x="Ind",by.y="Ind",x.all=T)
cor.test(CorrSize2$NumASV, CorrSize2$Area)

### Supplementary figures for Predicted functions Habitats ==========================================
## MetaCyc function analysis
PathwayHAB <- read.csv2("~/Desktop/Research/Microbiome/Review/Picrust/PathwayHAB.csv",header=T, sep=";")
kegg2<-PathwayHAB

order <- kegg2[order(kegg2[,4],decreasing=T),]
kegg2$Kegg<- factor(kegg2$Kegg, levels =c("Superpathway of taurine degradation","Heterolactic fermentation","Pyridoxal 5'-phosphate biosynthesis I","Superpathway of thiamine diphosphate biosynthesis II","Superpathway of ubiquinol-8 biosynthesis (early decarboxylation)","Ubiquinol-8 biosynthesis (early decarboxylation)","Ubiquinol-10 biosynthesis (early decarboxylation)","Ubiquinol-9 biosynthesis (early decarboxylation)","Ubiquinol-7 biosynthesis (early decarboxylation)","Superpathway of GDP-mannose-derived O-antigen building blocks biosynthesis","Tricarboxylic acid (TCA) cycle ","L-lysine biosynthesis III","Inosine-5'-phosphate biosynthesis III","Inosine 5'-phosphate degradation","Superpathway of L-isoleucine biosynthesis I","Guanosine nucleotides degradation III","Queuosine biosynthesis I (de novo)","Pyruvate fermentation to propanoate I","2-amino-3-carboxymuconate semialdehyde degradation to 2-hydroxypentadienoate","L-histidine degradation I","Pyrimidine deoxyribonucleotides de novo biosynthesis II","S-adenosyl-L-methionine salvage (SAM) I","Superpathway of histidine, purine, and pyrimidine biosynthesis","Homolactic fermentation"))
kegg2$Habitat<- factor(kegg2$Habitat, levels =c("MD", "UP", "BR"))
kegg2$LDA<-as.numeric(kegg2$LDA)
cbPalette4bis <- c("MD"="dodgerblue4","UP"="gold3","BR"="firebrick1")

KeggHAB<-ggplot(kegg2)+
  geom_bar(aes(x=Kegg,y=LDA,fill=Habitat),stat="identity",position = position_stack(reverse = TRUE))+ ylim(0,3.2) +
  theme_bw()+theme(panel.grid.major=element_line(colour="white"),panel.grid.minor=element_line(colour="white"),strip.background = element_rect(fill="white")) + guides(fill = guide_legend(keywidth = 1, keyheight = 1))+scale_fill_manual(values=cbPalette4bis)+coord_flip()
KeggHAB

ggsave(KeggHAB,height = 5, width = 13,filename=paste("KeggHab2",".pdf",sep=""),path = "~/Desktop/Research/Microbiome/Review/Picrust")

## KEGG level 3 analysis
PathwayKeggHAB <- read.csv2("~/Desktop/Research/Microbiome/Review/Picrust/PathwayKeggHAB.csv",header=T, sep=";")
kegg3<-PathwayKeggHAB

order <- kegg3[order(kegg3[,4],decreasing=T),]
kegg3$Kegg<- factor(kegg3$Kegg, levels =c("Metabolism.Metabolism of Terpenoids and Polyketides.Terpenoid backbone biosynthesis","Genetic Information Processing.Folding sorting and degradation.Proteasome","Environmental Information Processing.Membrane transport.Bacterial secretion system","Metabolism.Lipid Metabolism.Biosynthesis of unsaturated fatty acids","Cellular Processes.Transportand catabolism.Lysosome","Metabolism.Glycan biosynthesis and metabolism.Other glycan degradation","Metabolism.Carbohydrate Metabolism.Galactose metabolism","Unclassified.metabolism.Carbohydrate metabolism","Organismal Systems.Aging.Longevity regulating pathway_multiple species","Unclassified.metabolism.Energy metabolism"))
kegg3$Habitat<- factor(kegg3$Habitat, levels =c("MD", "UP", "BR"))
kegg3$LDA<-as.numeric(kegg3$LDA)
cbPalette4bis <- c("MD"="dodgerblue4","UP"="gold3","BR"="firebrick1")

KeggFunHAB<-ggplot(kegg3)+
  geom_bar(aes(x=Kegg,y=LDA,fill=Habitat),stat="identity",position = position_stack(reverse = TRUE))+ ylim(0,3.5) +
  theme_bw()+theme(panel.grid.major=element_line(colour="white"),panel.grid.minor=element_line(colour="white"),strip.background = element_rect(fill="white")) + guides(fill = guide_legend(keywidth = 1, keyheight = 1))+scale_fill_manual(values=cbPalette4bis)+coord_flip()
KeggFunHAB

ggsave(KeggFunHAB,height = 2.5, width = 13,filename=paste("KeggFunHab2",".pdf",sep=""),path = "~/Desktop/Research/Microbiome/Review/Picrust")