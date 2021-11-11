##################
#read in faststructure output
##################

# install dependencies and devtools
#install.packages(c("Cairo","ggplot2","gridExtra","gtable","tidyr","devtools"),dependencies=T)
# install pophelper package from GitHub
#devtools::install_github('royfrancis/pophelper')
library(pophelper)
library(dplyr)
library(data.table)

setwd("~/commongarden_merged_missing_faststruct")
# collect aligned files. Default prefix is pop
sfiles <- list.files(pattern = "*2.meanQ")
slist <- sortQ(readQ(sfiles))
slist_1 <- alignK(slist)
slist_final<- mergeQ(slist_1)

#setwd("~/structure_CG_run1/")
# collect aligned files. Default prefix is pop
sfiles <- list.files(pattern = "*3.meanQ")
slist <- sortQ(readQ(sfiles))
slist_1 <- alignK(slist)
slist_final[2] <- mergeQ(slist_1)

slist_final<-as.qlist(slist_final)

###################
#structure plot, fig 2C
###################

meta_ordered<-read.table("commongarden_and_PNAS_metadata_ancestry.txt",header = T)
inds<-as.data.frame(meta_ordered$Ind)
names(inds)<-"ind"
if(length(unique(sapply(slist_final,nrow)))==1) slist_final <- lapply(slist_final,"rownames<-",inds$ind) #label individuals in your structure matrix for plotting

group<-as.data.frame(meta_ordered[,c(4,5,8,3)]) #ignore these depending on your labelling needs - here im just subsetting columns that I'll use for labelling, and have a broad region label as well as a population label
names(group)<-c("Pop","Env","Pair","Long") 
group$Long<-round(as.numeric(group$Long),digits = 3)
group$Long<-group$Long*-1
group$Pop<-as.character(group$Pop)

#plotting
clist <- list(
  "shiny"=c("#1D72F5","#DF0101","#77CE61", "#FF9326","#A945FF","#0089B2","#FDF060","#FFA6B2","#BFF217","#60D5FD","#CC1577","#F2B950","#7FB21D","#EC496F","#326397","#B26314","#027368","#A4A4A4","#610B5E"),
  "strong"=c("#11A4C8","#63C2C5","#1D4F9F","#0C516D","#2A2771","#396D35","#80C342","#725DA8","#B62025","#ED2224","#ED1943","#ED3995","#7E277C","#F7EC16","#F8941E","#8C2A1C","#808080"),
  "oceanfive"=c("#00A0B0", "#6A4A3C", "#CC333F", "#EB6841", "#EDC951"),
  "keeled"=c("#48B098", "#91CB62", "#FFEE3B", "#FB9013", "#FF3C28"),
  "vintage"=c("#400F13", "#027368", "#A3BF3F", "#F2B950", "#D93A2B"),
  "muted"=c("#46BDDD","#82DDCE","#F5F06A","#F5CC6A","#F57E6A"),
  "teal"=c("#CFF09E","#A8DBA8","#79BD9A","#3B8686","#0B486B"),
  "merry"=c("#5BC0EB","#FDE74C","#9BC53D","#E55934","#FA7921"),
  "funky"=c("#A6CEE3", "#3F8EAA", "#79C360", "#E52829", "#FDB762","#ED8F47","#9471B4"),
  "retro"=c("#01948E","#A9C4E2","#E23560","#01A7B3","#FDA963","#323665","#EC687D"),
  "cb_paired"=c("#A6CEE3","#1F78B4","#B2DF8A","#33A02C","#FB9A99","#E31A1C","#FDBF6F","#FF7F00","#CAB2D6","#6A3D9A","#FFFF99","#B15928"),
  "cb_set3"=c("#8DD3C7","#FFFFB3","#BEBADA","#FB8072","#80B1D3","#FDB462","#B3DE69","#FCCDE5","#D9D9D9","#BC80BD","#CCEBC5","#FFED6F"),
  "morris"=c("#4D94CC","#34648A","#8B658A","#9ACD32","#CC95CC","#9ACD32","#8B3A39","#CD6601","#CC5C5B","#8A4500"),
  "wong"=c("#000000","#E69F00","#56B4E9","#009E73","#F0E442","#006699","#D55E00","#CC79A7"),
  "krzywinski"=c("#006E82","#8214A0","#005AC8","#00A0FA","#FA78FA","#14D2DC","#AA0A3C","#FA7850","#0AB45A","#F0F032","#A0FA82","#FAE6BE"))


names(slist_final)<-c("K2","K3")
kelly <- c("#222222","#e6e7e8","grey60")


plotQ(slist_final[1:2],imgoutput = "join",grplab=group, sppos = "left",
      ordergrp=T,showlegend=F,useindlab=T,grplabsize = 1.75, grplabangle = 65,
      showtitle=F,showsubtitle=F,divsize = .5,splabsize = 5,
      height=2,width=30,indlabheight=0.4,indlabspacer=5, grplabpos = 1,
      splab=c("K=2","K=3"),divgrp="Pop",selgrp = "Long", 
      barbordercolour="white",barbordersize=0,outputfilename="faststruct_commongarden_mergedmissing_K2K3",imgtype="pdf",
      clustercol = kelly,  sharedindlab=F,showindlab=F)


###########################
#ancestry by pair, fig2D
###########################

some<-read.table("~/ancestry_byaccession.txt",header=T)
long_some<-melt(some,id.vars=c(1:5,8))
std <- function(x) sd(x)/sqrt(length(x))
library(forcats)
long_some$Pair<-as.factor(long_some$Pair)
long_some$Pair<-as.factor(long_some$Pair)
#install.packages("ggtext")
library(ggtext)  # remotes::install_github("clauswilke/ggtext")


bypop<-long_some %>% filter(variable== "K1") %>% 
  group_by(Region) %>% dplyr::summarise(Kmean=mean(value), diff_lower= quantile(value,probs =  0.05), diff_upper= quantile(value,probs =  0.95), Env=first(Env), Pair=first(Pair), Long=first(Long.x))

#anxestry by pair, sorted by longitude
long_some %>% filter(variable== "K1") %>% 
  group_by(Pair,Env) %>% dplyr::summarise(Kmean=mean(value), Ksd= std(value), Long=first(Long.x)) %>%
  ggplot(aes(x=reorder(Pair,Long), Kmean,color=Env))  +
  geom_errorbar(aes(ymin=Kmean-Ksd, ymax=Kmean+Ksd), width=.2) +
  #geom_line(aes(group = fct_reorder(Pair,Long))) +
  geom_point(cex=2) +
  theme_bw() +
  ylab("Proportion of \nvar. rudis Ancestry") +
  xlab("Population Pair") +
  scale_color_manual(values=c(bay[5],bay[1])) +
  ylim(0,1)


#regression of ancestry by variables, get lsmeans
long_some_K1<-long_some[long_some$variable == "K1",]
k1byenv<-lm(data=long_some_K1, value ~  Pair + Env + Long.x + Lat.x)
Anova(k1byenv, type=2)
summary(k1byenv)
library(lsmeans)
env_means<-as.data.frame(lsmeans(k1byenv, "Env"))
#lsmip(k1byenv, ~ Env)
head(env_means)

#plot of least squares mean rudis ancestry by pair
ggplot(data=env_means,aes(Env, lsmean, color=Env)) +
  geom_point(cex=2) +
  geom_errorbar(aes(ymin=lsmean-SE, ymax=lsmean+SE),width=.4) +
  scale_color_manual(values=c(bay[5],bay[1])) +
  #ylab("Least Squares Mean \nProportion of \nvar. rudis Ancestry") +
  theme_bw() +
  ylim(0,1)
  

ag <- bypop [bypop$Env == "Ag",] 
nat <- bypop [bypop$Env == "Nat",] 

bothenv<-inner_join(ag,nat,by="Pair")
bothenv$diff<-bothenv$Kmean.x - bothenv$Kmean.y
median(bothenv[bothenv$diff > 0,]$diff) #percent excess rudis ancestry for pairs where ag > nat


############################
#PCA plotting and analyses, fig 2A-B
############################

#analyses on both common garden and kreiner 2019 (PNAS) samples
eigen<-fread("~/commongarden_and_PNAS_snps_missing2.eigenvec",header = F)
CG_pca<-read.table("commongarden_and_PNAS_snps.eigenvec")

CG_key<-read.table("~/3waymerged_sampleinfo.txt", na.strings = c(""),sep = "\t",header=T)
CG_key<-CG_key[-c(188:295),]
CG_pca<-CG_pca[,-1]
names(CG_pca)[1]<-"sample"
merged<-inner_join(CG_pca,CG_key,by="sample")


library(ggplot2)
library(car)
merged$Region <- merged$state
merged$Region[merged$year == 2019] <- "Common Garden"

merged$Region<-as.factor(merged$Region)
library(PNWColors)
bay<-pnw_palette("Bay",6,type="continuous")

#plot of joint PCA
ggplot(data=merged, aes(V3,V4,color=Region)) +
  geom_point(size=2,alpha=.8) +
  labs(x="PC1 (18%)",y="PC2 (7%)") +
  theme_bw() +
  scale_color_manual(values=c(bay[2],bay[4],bay[3],bay[6],"grey30","grey60")) +
  theme(legend.title = element_blank())


#analyses on just common garden samples

CG_pca<-read.table("commongarden_mergedSNPS_missing2_justCGinds.eigenvec")
CG_key<-read.table("contemp_wpair.txt",header=F)
names(CG_key)<-c("sample","env","sex","pair")
CG_key2<-read.table("3waymerged_sampleinfo.txt",na.strings = c("","NA"),sep = "\t",header=T)
CG_key$sample<-as.character(as.integer(CG_key$sample))
both<-inner_join(CG_key,CG_key2,by="sample")                   

CG_pca<-CG_pca[,-1]
names(CG_pca)[1]<-"sample"
CG_pca$sample<-as.character(as.integer(CG_pca$sample))
merged<-inner_join(CG_pca,both,by="sample")
test<-merged[c(1,2,3,22,23,24,27,28,30)]
head(test)

library(reshape2)
names(test)[c(2,3)]<-c("PC1","PC2")
longpc<-melt(test, id.vars = c("sample","env.x","sex.x","pair","lat","long","sex.x","state"))

#PCA regressions
lm1<-lm(data=merged, V3 ~ env.x  + long +lat + pair + sex.x)
Anova(lm1, type = 3)
lm1<-lm(data=merged, V4 ~ env.x  + long +lat + pair + sex.x)
Anova(lm1, type = 3)

#plot fig 1B
names(merged)[c(22,28)]<-c("Environment","Longitude")
ggplot(data=merged, aes(-V3,V4,color=Longitude, shape=Environment)) +
  geom_point(size=2) +
  labs(x="PC1 (18%)",y="PC2 (15%)") +
  theme_bw() +
  scale_colour_gradient(low = "white", high = "black")

#grouped regression of PCs
herb_contemp<-lm(data=longpc, value ~ env.x*variable + long*variable + pair*variable )
Anova(herb_contemp,type = 3)
