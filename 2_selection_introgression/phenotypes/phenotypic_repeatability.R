library(tidyverse)
theme_set(theme_classic())

#not including lon at 40 because Alaska has low sample size & Washington is just north of 40
migrateVars<-c("fall_30","spring_30",
               "wintering_long","wintering_lat",
               "distance_fall","distance_spring",
               "fall_day","spring_day")

morphoVars<-c("tarsus.length","tail.length","wing.cord","kipps","distal")

specVars<-c("u","s","m","l")

pheno<-read.csv("C:/Users/Steph/GitHub/thrush_hybrids/genomic_clines/phenos/SWTH_trait_data.csv")

pheno.sp<-pheno%>%mutate_at(vars(all_of(c(migrateVars,morphoVars,specVars))),as.numeric)%>%
  filter(species!="hybrid")%>%
  pivot_longer(cols=all_of(c(migrateVars,morphoVars,specVars)),names_to="trait",values_to="traitVal")


pheno.l<-pheno%>%mutate_at(vars(all_of(c(migrateVars,morphoVars,specVars))),as.numeric)%>%
  filter(species=="hybrid")%>%
  pivot_longer(cols=all_of(c(migrateVars,morphoVars,specVars)),names_to="trait",values_to="traitVal")

sampleSizes<-pheno.l%>%select(trait,reference,release_site,traitVal)%>%
  distinct()%>%
  drop_na()%>%
  group_by(trait,release_site)%>%
  summarise(n())%>%
  pivot_wider(values_from=`n()`,names_from=release_site)

write.csv(sampleSizes,'C:/Users/Steph/GitHub/thrush_hybrids/genomic_clines/phenos/sampleSizes_byTrait.csv',
          row.names=F)


ggplot(pheno.l,aes(x=aims_ancestry,y=traitVal,colour=release_site))+
  geom_point(alpha=0.5)+geom_smooth(method="lm",se=F)+
  facet_wrap(vars(trait),scales="free")

pheno.l%>%group_by(trait,release_site)%>%
  filter(!is.na(traitVal)&release_site=="Washington")%>%
  summarise(ancestry_range=(max(aims_ancestry)-min(aims_ancestry)),
            sample_size=n())
pheno.sp%>%filter(!is.na(traitVal))%>%pull(trait)%>%unique()

phenoOut<-data.frame(); phenoNew<-data.frame()
for(tr1 in c(migrateVars,morphoVars,specVars)){
  for(pop1 in unique(pheno.l$release_site)){
    df1<-pheno.l%>%filter(trait==tr1&release_site==pop1&!is.na(traitVal))%>%
      mutate(ynorm=(traitVal-min(traitVal))/(max(traitVal)-min(traitVal)))
    
    sp_diff<-pheno.sp%>%filter(trait==tr1)%>%
      summarise(sp_diff=mean(traitVal[species=="inland"],na.rm=T)-mean(traitVal[species=="coastal"],na.rm=T))%>%
      pull(sp_diff)
    #if(sp_diff<0){df1<-df1%>%mutate(ynorm=ynorm*(-1))}
    
    df2<-df1%>%summarise(ancestry_range=(max(aims_ancestry)-min(aims_ancestry)),
                         sample_size=n())
    #if(df2$ancestry_range>0.4&df2$sample_size>7){
      lm1<-lm(traitVal~aims_ancestry,df1)
      lm2<-lm(ynorm~aims_ancestry,df1)
      phenoOut<-rbind(phenoOut,
                      data.frame(release_site=pop1,trait=tr1,
                                 slope=as.vector(lm1$coefficients["aims_ancestry"]),
                                 slopeNorm=as.vector(lm2$coefficients["aims_ancestry"])))
      phenoNew<-rbind(phenoNew,df1)
    #}
  }
}






phenoRepeat<-data.frame()
for(tr1 in c(migrateVars,morphoVars,specVars)){
  df1<-pheno.l%>%filter(trait==tr1)
  lm1<-lm(traitVal~aims_ancestry*release_site,df1)
  aovlm1<-anova(lm1)
  S=1-aovlm1$`Sum Sq`[3]/sum(aovlm1$`Sum Sq`[c(1,3)])
  phenoRepeat<-rbind(phenoRepeat,data.frame(trait=tr1,repeatability=S))
}
phenoRepeat<-phenoRepeat%>%mutate(category=case_when(trait%in%migrateVars~"migratory",
                                        trait%in%morphoVars~"morphology",
                                        trait%in%specVars~"colour"))

library(ggpubr)
library(viridis)

a<-ggplot(phenoRepeat,aes(x=category,y=repeatability))+
  geom_point(size=3,alpha=0.7,colour="thistle4")+
  xlab("trait category")+
  ylim(0,1)

b<-ggplot(phenoNew%>%filter(trait%in%c("spring_30","wing.cord","l"))%>%
            mutate(trait=case_when(trait=="l"~"short wavelengths",
                   trait=="spring_30"~"spring longitude",
                   trait=="wing.cord"~"wing chord")),
       aes(x=aims_ancestry,y=traitVal,colour=release_site))+
  geom_point(alpha=0.4)+geom_smooth(method="lm",se=F)+
  scale_colour_viridis(discrete=T,option="mako",end=0.9)+
  xlab("ancestry")+ylab("trait value")+
  facet_grid(vars(trait),scales="free")



ggarrange(a,b,widths=c(2,2.5))




ggplot(phenoNew,
       aes(x=aims_ancestry,y=traitVal,colour=release_site))+
  geom_point(alpha=0.4)+geom_smooth(method="lm",se=F)+
  scale_colour_viridis(discrete=T,option="mako",end=0.9)+
  facet_wrap(vars(trait),scales="free")

#get correlation between population pairs
phenoOut<-phenoOut%>%mutate(category=case_when(trait%in%migrateVars~"migratory",
                                     trait%in%morphoVars~"morphology",
                                     trait%in%specVars~"colour"))




phenoCorr<-data.frame()
popCombos<-t(combn(sort(unique(phenoOut$release_site)),2))
for(i in 1:nrow(popCombos)){
  for(cat1 in unique(phenoOut$category)){
    
    pop1=popCombos[i,1]
    pop2=popCombos[i,2]
    
    tr1<-phenoOut%>%filter(release_site==pop1&category==cat1)%>%arrange(trait)
    tr2<-phenoOut%>%filter(release_site==pop2&category==cat1)%>%arrange(trait)
    if(all(tr1$trait==tr2$trait)){
    r1<-cor(tr1%>%pull(slopeNorm),
            tr2%>%pull(slopeNorm))
    phenoCorr<-rbind(phenoCorr,
                     data.frame(populations=paste(pop1,pop2,sep=" x "),
                                trait=cat1,corCoef=r1))
    }else{print(paste("check trait matching",pop1,pop2))}
    }}

clines<-read.csv("C:/Users/Steph/GitHub_data/bgchm_clines/gradients_withZ.csv")
colnames(clines)



for(i in 1:nrow(popCombos)){
    
    pop1=popCombos[i,1]
    pop2=popCombos[i,2]
    
    tr1<-clines%>%filter(population==pop1)%>%arrange(locus)
    tr2<-clines%>%filter(population==pop2)%>%arrange(locus)
    if(all(tr1$trait==tr2$trait)){
      r1<-cor(tr1%>%pull(median),
              tr2%>%pull(median))
      phenoCorr<-rbind(phenoCorr,
                       data.frame(populations=paste(pop1,pop2,sep=" x "),
                                  trait="genomic",corCoef=r1))
    }else{print(paste("check trait matching",pop1,pop2))}
  }

phenoCorr$trait<-factor(phenoCorr$trait,levels=c("colour","migratory","morphology","genomic"))

pop.cols<-c("#202D26","#49654A","#849D74","#ADB593")
a<-ggplot(phenoCorr,aes(x=trait,y=corCoef,colour=populations))+
  geom_hline(yintercept=0,lty=2,colour="grey60")+
  geom_point(size=4,alpha=0.7)+
  xlab("trait category")+ylab("correlation")+
  scale_colour_viridis(discrete=T,option="mako",end=0.8)

b<-ggplot(phenoNew%>%filter(trait%in%c("spring_30","wing.cord","l"))%>%
            mutate(trait=case_when(trait=="l"~"short wavelengths",
                                   trait=="spring_30"~"spring longitude",
                                   trait=="wing.cord"~"wing chord")),
          aes(x=aims_ancestry,y=ynorm,colour=release_site))+
  geom_point(alpha=0.4)+geom_smooth(method="lm",se=F)+
  #scale_colour_viridis(discrete=T,option="mako",end=0.9,name="population")+
  scale_colour_manual(values=pop.cols,guide="none")+
  xlab("ancestry")+ylab("trait value")+
  facet_grid(cols=vars(trait))+
  scale_x_continuous(breaks=c(0,0.5,1))+
  scale_y_continuous(breaks=c(0,0.5,1))
b

c<-ggplot(phenoOut,aes(x=category,y=abs(slopeNorm)))+
  geom_boxplot(fill=NA)+
  geom_jitter(aes(colour=release_site),size=4,alpha=0.9,height=0,width=0.1)+
  scale_colour_manual(values=pop.cols,name="population")+
  ylab("| slope |")


p1<-ggarrange(c,b,nrow=2,heights=c(3,2),labels=c("A","B"))
fig7<-ggarrange(p1,a,widths=c(2,2.3),labels=c("","C"))

ggsave("C:/Users/Steph/GitHub/thrush_hybrids/genomic_clines/phenos/Fig7-revise.png",
      plot=fig7,width=10,height=6,units="in",bg="white")

phenoCorr%>%group_by(trait)%>%
  summarise(min(corCoef),max(corCoef))

phenoOut%>%group_by(category)%>%
  summarise(min(abs(slopeNorm)),max(abs(slopeNorm)))

library(lme4); library(car)
Anova(lmer(abs(slopeNorm)~category+(1|trait)+(1|release_site),phenoOut))
Anova(lmer(abs(corCoef)~trait+(1|populations),phenoCorr))



lmMigrate<-lm(slope~trait*release_site,data=phenoOut%>%filter(trait%in%migrateVars))
aovMigrate<-anova(lmMigrate)
1-aovMigrate$`Sum Sq`[3]/sum(aovMigrate$`Sum Sq`[c(1,3)])

lmMorpho<-lm(slope~trait*release_site,data=phenoOut%>%filter(trait%in%morphoVars))
aovMorpho<-anova(lmMorpho)
1-aovMorpho$`Sum Sq`[3]/sum(aovMorpho$`Sum Sq`[c(1,3)])

lmSpec<-lm(slope~trait*release_site,data=phenoOut%>%filter(trait%in%specVars))
aovSpec<-anova(lmSpec)
1-aovSpec$`Sum Sq`[3]/sum(aovSpec$`Sum Sq`[c(1,3)])

