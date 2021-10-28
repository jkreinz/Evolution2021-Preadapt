replicated_bytreat<-read.table("Treatmentspecificmeans_byfamily_bysex.txt",header = T)
replicated_bytreat$Treatment<-factor(replicated_bytreat$Treatment,levels = c("Water","Control","Soy"))
replicated_bytreat_males<-replicated_bytreat[replicated_bytreat$Sex == "M",]
replicated_bytreat_females<-replicated_bytreat[replicated_bytreat$Sex == "F",]


#standardize predictors

replicated_bytreat_females$stnd_biomass<-scale(replicated_bytreat_females$BMmeans,center = T)
replicated_bytreat_females$stnd_StemWidth<-scale(replicated_bytreat_females$SWmeans,center = T)
replicated_bytreat_females$stnd_NodeNum<-scale(replicated_bytreat_females$NNmeans,center = T)
replicated_bytreat_females$stnd_PlantHeight<-scale(replicated_bytreat_females$PHmeans,center = T)
replicated_bytreat_females$stnd_DaystoFlower<-scale(replicated_bytreat_females$FTmeans,center = T)
replicated_bytreat_females$stnd_Germ_JD<-scale(replicated_bytreat_females$Germmeans,center = T)
replicated_bytreat_females$stnd_Stem.Colour<-scale(replicated_bytreat_females$SCmeans,center = T)
#merged_wFlower_heteige$stnd_Hypocotyl.mm.<-scale(replicated_bytreat_females$,center = T)
replicated_bytreat_females$stnd_RateofCY<-scale(replicated_bytreat_females$COTYmeans,center = T)
replicated_bytreat_females$stnd_RateofLN<-scale(replicated_bytreat_females$LNmeans,center = T)
replicated_bytreat_females$stnd_Flower.Color<-scale(replicated_bytreat_females$FCmeans,center = T)
replicated_bytreat_females$stnd_EPSPS<-scale(replicated_bytreat_females$EPSPS_mapped_cov_mean,center = T)
replicated_bytreat_females$stnd_GenomeSize<-scale(replicated_bytreat_females$Genome.Size..Mb.,center = T)
replicated_bytreat_females$stnd_Het<-scale(replicated_bytreat_females$Het,center = T)
replicated_bytreat_females$stnd_hypo<-scale(replicated_bytreat_females$HYPOmeans,center = T)

replicated_bytreat_males$stnd_biomass<-scale(replicated_bytreat_males$BMmeans,center = T)
replicated_bytreat_males$stnd_StemWidth<-scale(replicated_bytreat_males$SWmeans,center = T)
replicated_bytreat_males$stnd_NodeNum<-scale(replicated_bytreat_males$NNmeans,center = T)
replicated_bytreat_males$stnd_PlantHeight<-scale(replicated_bytreat_males$PHmeans,center = T)
replicated_bytreat_males$stnd_DaystoFlower<-scale(replicated_bytreat_males$FTmeans,center = T)
replicated_bytreat_males$stnd_Germ_JD<-scale(replicated_bytreat_males$Germmeans,center = T)
replicated_bytreat_males$stnd_Stem.Colour<-scale(replicated_bytreat_males$SCmeans,center = T)
#merged_wFlower_heteige$stnd_Hypocotyl.mm.<-scale(replicated_bytreat_males$,center = T)
replicated_bytreat_males$stnd_RateofCY<-scale(replicated_bytreat_males$COTYmeans,center = T)
replicated_bytreat_males$stnd_RateofLN<-scale(replicated_bytreat_males$LNmeans,center = T)
replicated_bytreat_males$stnd_Flower.Color<-scale(replicated_bytreat_males$FCmeans,center = T)
replicated_bytreat_males$stnd_EPSPS<-scale(replicated_bytreat_males$EPSPS_mapped_cov_mean,center = T)
replicated_bytreat_males$stnd_GenomeSize<-scale(replicated_bytreat_males$Genome.Size..Mb.,center = T)
replicated_bytreat_males$stnd_Het<-scale(replicated_bytreat_males$Het,center = T)
replicated_bytreat_males$stnd_long<-scale(replicated_bytreat_males$Long,center = T)
replicated_bytreat_males$stnd_lat<-scale(replicated_bytreat_males$Lat,center = T)
replicated_bytreat_males$stnd_sex<-scale(replicated_bytreat_males$Sex.x,center = T)
replicated_bytreat_males$stnd_hypo<-scale(replicated_bytreat_males$HYPOmeans,center = T)


#selection gradients

library(car)
levels(replicated_bytreat_males$Treatment)

male_selgrad_quad<-lm(BMmeans ~ stnd_Germ_JD + stnd_RateofLN + stnd_RateofCY + stnd_hypo +
                        stnd_PlantHeight + stnd_DaystoFlower + stnd_StemWidth + stnd_NodeNum + 
                        + stnd_Stem.Colour  + stnd_Flower.Color +
                        Lat.x  + Long.x + Env + K1*Treatment + I(K1^2) , data=replicated_bytreat_males)


male_selgrad_quad_ft<-lm(FTmeans ~ stnd_Germ_JD + stnd_RateofLN + stnd_RateofCY + stnd_hypo +
                           stnd_PlantHeight + stnd_biomass + stnd_StemWidth + stnd_NodeNum + 
                           + stnd_Stem.Colour  + stnd_Flower.Color +
                           Lat.x  + Long.x + Env + K1*Treatment + I(K1^2), data=replicated_bytreat_males)

Anova(male_selgrad_quad, type=3)
Anova(male_selgrad_quad_ft, type=3)

summary(male_selgrad_quad)
summary(male_selgrad_quad_ft)


female_selgrad_quad<-lm(BMmeans ~ stnd_Germ_JD + stnd_RateofLN + stnd_RateofCY + stnd_hypo +
                          stnd_PlantHeight + stnd_DaystoFlower + stnd_StemWidth + stnd_NodeNum + 
                          + stnd_Stem.Colour  + stnd_Flower.Color +
                          Lat.x+ Long.x + Env + K1 *Treatment + I(K1^2) , data=replicated_bytreat_females)

female_selgrad_quad_ft<-lm(FTmeans ~ stnd_Germ_JD + stnd_RateofLN + stnd_RateofCY + stnd_hypo +
                             stnd_PlantHeight + stnd_biomass + stnd_StemWidth + stnd_NodeNum + 
                             + stnd_Stem.Colour  + stnd_Flower.Color +
                             Lat.x+ Long.x + Env + K1 *Treatment + I(K1^2) , data=replicated_bytreat_females)


Anova(female_selgrad_quad,type=3)
Anova(female_selgrad_quad_ft,type=3)
summary(female_selgrad_quad)



library(lsmeans)
lsmeans(male_selgrad_quad, specs=c("Treatment","Env"))
test1<-lsmip(male_selgrad_quad, Treatment ~  K1 + Env, ylab = "Biomass/Day",
             at=list(K1=c(0, 0.2, 0.4, 0.6,0.8,1)),
             xlab="Proportion var. rudis ancestry",
             type="response",
             plotit=F)

male_selgrad_noquad<-lm(BMmeans ~ stnd_Germ_JD + stnd_RateofLN + stnd_RateofCY + stnd_hypo +
                          stnd_PlantHeight + stnd_DaystoFlower + stnd_StemWidth + stnd_NodeNum + 
                          + stnd_Stem.Colour  + stnd_Flower.Color +  
                          Lat.x+ Long.x + Env + K1 * Treatment  , data=replicated_bytreat_males)

test2<-lsmip(male_selgrad_noquad, Treatment ~  K1 + Env, ylab = "Biomass/Day",
             at=list(K1=c(0, 0.2, 0.4, 0.6,0.8,1)),
             xlab="Proportion var. rudis ancestry",
             type="response",
             plotit=F)

AIC(male_selgrad_noquad,male_selgrad_quad)

test1$Treatment <- factor(test1$Treatment, levels = c("Soy", "Control", "Water"))
p1<- ggplot(data=test1, aes(K1,yvar, color=Treatment)) +
  geom_ribbon(aes(ymin=test1$yvar-test1$SE, ymax=test1$yvar+test1$SE, fill=Env), alpha=.15, linetype = 0) +
  # geom_point() +
  geom_smooth() +
  #geom_smooth(data=test2, lty="dashed",aes(alpha=.8)) +
  ylim(0.03, 0.085) +
  facet_wrap(~Env) +
  theme_bw() +
  xlab("Proportion var. rudis ancestry") +
  ylab("Biomass/Day (g)") +
  theme(legend.position = "none") +
  scale_color_manual(values=c("black","grey50","grey70")) +
  scale_fill_manual(values=c(bay[5],bay[1]))

test<-lsmip(male_selgrad_quad_ft, Treatment ~ K1 , ylab = "Days to Flowering",
            at=list(K1=c(0, 0.2, 0.4, 0.6,0.8,1)),
            xlab="Proportion var. rudis ancestry",
            type="response",
            plotit=F)

test$Treatment <- factor(test$Treatment, levels = c("Soy", "Control", "Water"))
p2<- ggplot(data=test, aes(K1,yvar, color=Treatment, group=Treatment)) +
  #geom_point() +
  geom_ribbon(aes(ymin=test$yvar-test$SE, ymax=test$yvar+test$SE), alpha=.1, linetype = 0) +
  geom_smooth() +
  ylim(40, 50.5) +
  theme_bw() +
  xlab("Proportion var. rudis ancestry") +
  ylab("Days until Flowering") +
  scale_color_manual(values=c("black","grey50","grey70"))

