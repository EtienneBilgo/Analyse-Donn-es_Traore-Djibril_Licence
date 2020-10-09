### Analyses des donnnées Licence Djibril Traore ###
# Data importation
Virulence_Dj <- read.csv("~/Desktop/Virulence_Dj.csv")
View(Virulence_Dj)
levels(Virulence_Dj$Treatment)
##Traitement des doneées
colnames(Virulence_Dj)
colnames(Virulence_Dj)[2]="Replicate"
colnames(Virulence_Dj)
library(plyr)
df2=ddply(Virulence_Dj, .(Replicate,Treatment), transform,
          nval=sum(Motartility, na.rm=TRUE))
View(df2)
df3=ddply(df2, .(Replicate,Treatment), transform,
          Dead=cumsum(Motartility))
View(df3) 
df3$Survival=(1-df3$Dead/df3$nval) 
View(df3) 
colnames(df3)[3]="Treatments" 
View(df3)
df4=subset(df3, Day!="Still Alives") 
df4$Day=as.numeric(as.character(df4$Day)) 
View(df4)
df5=ddply(df4, .(Day,Treatments), summarize,
          mean=mean(Survival), replicates=length(Survival),
          se=sd(Survival)/sqrt(length(Survival)))
View(df5)
library(ggplot2)
library(scales)
limits=aes(ymax=mean+se, ymin=mean-se)
theme = theme_bw()+theme(text = element_text(size=25),
                         axis.title.x = element_text(size=25), axis.title.y = element_text(size=25), title = element_text(size=25), legend.title = element_text(size=25), legend.text = element_text(size=20))
levels(df5$Treatments) 
cbPalette <-c("#1D07F3","#A807F3","#12F307","#FDD7DB","#E16572","#780410")
colnames(df5)[2]="Traitements"
levels(df5$Traitements)[1]="Contrôle"
Plt=ggplot(df5, aes(Day, mean,color=Traitements))+geom_line(size=2)+geom_errorbar(limits, width=.1, size=1)+theme+scale_colour_manual(values=cbPalette)+ylab("Survie")+xlab("Jours post-infection")+scale_y_continuous(labels=percent)
summary(Plt)
#Le graphique
Plt
library(plotly)
ggplotly(Plt)
Plt
#Calcul des LT50s
####Load packages####
library(tidyverse)
library(reshape2)
library(plyr)
library(scales)
library(survival)
library(MASS)
View(df4)
LTdat=df4
View(LTdat)
LTdat$Alives=LTdat$nval-LTdat$Dead
View(LTdat)
attach(LTdat)
surv.per=0.50
colnames(LTdat)
LTdat2=ddply(LTdat, .(Treatments, Replicate), summarize,LT=as.numeric(dose.p(glm(cbind(Alives,Dead)~Day,binomial),p=surv.per)))
View(LTdat2)
LTdat2[LTdat2$LT>14|LTdat2$LT<0,]$LT=NA
LT50.Error=ddply(LTdat2, .(Treatments), summarize,"LT50 Mean"=mean(LT,na.rm=T), se=sd(LT,na.rm=T)/sqrt(length(LT[!is.na(LT)])),Replicates=length(LT[!is.na(LT)])) 
View(LT50.Error)
### Comparisons des LT50s
FF_dat <- read.csv("~/Desktop/FF_dat.csv")
View(FF_dat)
FF_dat <- read.csv("~/Desktop/FV_dat.csv")
View(FF_dat)
colnames(FF_dat)
library(ggplot2)
library(scales)
FF_dat$Fecondite=FF_dat$Oeufs/FF_dat$Total.Moustiques
View(FF_dat)
colnames(FF_dat)[3]="Traitements"
boxplot(FF_dat$Fecondite ~Traitements,data= FF_dat,outpch=NA) 
View(FF_dat)
FF_dat$Viability=FF_dat$Larves/FF_dat$Oeufs
View(FF_dat)
library(ggplot2)
library(scales)
colnames(FF_dat)
levels(FF_dat$Traitements)
cbPalette <- c("#C217EC","#1721EC","#17EC2E","#E93163 ")
theme = theme_bw()+theme(text = element_text(size=20),axis.title.x=element_text(size=30),axis.text.x=element_text(angle=0,hjust =.5,vjust=.5, size=25), axis.text.y=element_text(size=25), title=element_text(size=35),legend.title=element_text(size =25),legend.text=element_text(size =20)) 
Plt<- ggplot(FF_dat,aes(Traitements,Viability,fill=Traitements))+geom_boxplot()+theme+scale_colour_manual(values=cbPalette)+xlab("Traitements")+ylab("Pourcentage de Viabilité")+scale_y_continuous(breaks=pretty_breaks(n = 10), labels=percent)
Plt+scale_fill_manual(values=c("#C217EC","#1721EC","#17EC2E","#E93163")) 

### Analyses statistiques
### Stats pour Fecondite##
View(FF_dat)
hist(FF_dat$Fecondite)
pairwise.t.test(FF_dat$Fecondite,FF_dat$Traitements)
## Resultats comparaisons t-test pour la Fecondite##
## Pairwise comparisons using t tests with pooled SD 
## ata:  FF_dat$Fecondite and FF_dat$Traitements 

#          C. violaceum Controle Met_S26
#  Controle      1            -        -      
#  Met_S26       1            1        -      
#  Met_S26_CV    1            1        1    
#  P value adjustment method: holm 

data:  FF_dat$Viability and FF_dat$Traitements 

## Stats pour viabilite##
View(FF_dat)
hist(FF_dat$Viability)
pairwise.t.test(FF_dat$Viability,FF_dat$Traitements)
## Resultats comparaisons t-test pour la Viabilité##
## Pairwise comparisons using t tests with pooled SD 
#           C. violaceum Controle  Met_S26
# Controle     3.6e-05      -        -      
# Met_S26     0.01715      0.00058   -      
# Met_S26_CV  0.63751      4.2e-05  0.02344
# P value adjustment method: holm 
