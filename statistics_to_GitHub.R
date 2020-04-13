library(vegan)
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggExtra)

###Occupancy-frequency bimodal distribution test
##This script was cordially provided by Markus Lindh
##Reference DOI: 10.1111/1462-2920.13650

setwd ("C:/Recuperar/Dropbox/Erick_Hugo/Artigos/Bimodalidade/R")
setwd ("C:/Recuperar/Dropbox/Erick_Hugo/Bioinformática/Biota_ASV/")
bimodal=read.table(file = "asv_rarefied.txt", header=TRUE, sep="\t", dec=".", na.strings="NA")
bimodal=bimodal[,-1]

x=t(bimodal)
x[which(x > 0)] = 1

print(ggplot(data=data.frame(apply(x, 1, sum)[which(apply(x, 2, sum) >0)]),aes(x=data.frame(apply(x, 2, sum)[which(apply(x, 2, sum) >0)])[,1]))+
        geom_histogram(fill="grey80",colour="black",binwidth=1)+
        ylim(c(0,ncol(x)))+
        ylab(paste("Frequency","(Number of OTUs)",sep="\n"))+
        xlab(paste("Occupancy","(Number of sites occupied)",sep="\n"))+
        ggtitle("")+
        coord_cartesian(ylim = c(0,120))+
        theme_bw())



table.out=as.data.frame(table(apply(x, 2, sum)[which(apply(x, 2, sum) >0)]))
colnames(table.out)= c("Occupancy","Frequency")
print(table.out)
#(ratio=(table.out[nrow(table.out),2])/(table.out[1,2]))
write.table(table.out,file=paste("C:/Recuperar/Dropbox/Erick_Hugo/Artigos/Bimodalidade/R/novas figuras/bimodal_otu.txt",sep="_"),sep="\t",row.names=F)


test.input=as.data.frame(table(apply(x, 2, sum)[which(apply(x, 2, sum) >0)]))
MOS.test=MOStest(as.numeric(table.out[,1]), 
                 as.numeric(table.out[,2]),
                 family=gaussian(link = "identity"), 
                 maxit = 200)
print(MOS.test)
capture.output(MOS.test,file=paste("MOStest_otu.txt",sep="_")) #Saves table with stats

#Fig. 2
library(vegan)
library(ggplot2)
library(dplyr)
library(tidyr)
library(scales)

setwd ("C:/Recuperar/Dropbox/Erick_Hugo/Artigos/Bimodalidade/R")
otu=read.table(file = "count.txt", header=TRUE, row.names = 1, sep="\t", dec=".", na.strings="NA")
count=as.data.frame(ifelse(otu>0,1,0))
count=decostand(count[1:60],"total",2)
count$Occ.type=otu$Occ.type

a=gather(otu,key = "OTUId", value = "abund", -Occ.type)
b=gather(count,key = "OTUId", value = "count", -Occ.type)

ggplot(b, aes(x=OTUId,y=count,fill=Occ.type,colour=Occ.type))+
  geom_col(width = 1)+
  #stat_count(aes(y=abund))+
  ylab(paste("Frequency (Number of OTUs)",sep="\n"))+
  xlab(paste("Site",sep="\n"))+
  scale_fill_manual(name="", labels=c("Core","Non-Core","Satellite"),values = c("red","blue","yellow"))+
  scale_color_manual(name="",labels=c("Core","Non-Core","Satellite"),values = c("red","blue","yellow"))+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, hjust = TRUE))

ggplot(a, aes(x=OTUId,y=abund,fill=Occ.type,colour=Occ.type))+
  geom_col(width = 1)+
  ylab(paste("Relative Abundance",sep="\n"))+
  xlab(paste("Site",sep="\n"))+
  scale_fill_manual(name="", labels=c("Core","Non-Core","Satellite"),values = c("red","blue","yellow"))+
  scale_color_manual(name="",labels=c("Core","Non-Core","Satellite"),values = c("red","blue","yellow"))+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, hjust = TRUE))


write.table(a,"abund.txt")  
write.table(b,"pa.txt")

###SpADs - Spatial Abundance Distributions
##This script was cordially provided by Juan Pablo Niño-Garcia and Clara Ruíz-Gonzáléz
##Reference DOI: 10.1111/ele.12704
##Reference DOI: 10.1111/mec.15026

setwd ("C:/Recuperar/Dropbox/Erick_Hugo/Artigos/Bimodalidade/R")

############modeling of species abundance distributions

library(diptest)
library(fitdistrplus)
library(scales)
library(plyr)
library(tidyr)
library(dplyr)
library(ggplot2)


otu_total=read.table(file = "count.txt", header=TRUE, row.names = 1, sep="\t", dec=".", na.strings="NA")
otu_core=filter(otu_total,Occ.type=="Core")
#otu_n.core=filter(otu_total,Occ.type!="Core")
otu_n.core=filter(otu_total,Occ.type=="n.core")
otu_satellite=filter(otu_total,Occ.type=="Satellite")
otu_total=otu_total[,1:60]
otu_core=otu_core[,1:60]
row.names(otu_core)=row.names(otu_total[1:15,])
#otu_n.core=otu_n.core[,1:60]
#row.names(otu_n.core)=row.names(otu_total[16:3946,])
otu_n.core=otu_n.core[,1:60]
row.names(otu_n.core)=row.names(otu_total[16:1780,])
otu_satellite=otu_satellite[,1:60]
row.names(otu_satellite)=row.names(otu_total[1781:3946,])



otusl=otu_n.core#trocar aqui pelas outras tabelas pra não causar confusão nas análises subsequentes
otusl=t(otusl)

###OTUSl IS MY RAREFIED OTU TABLE
# Create a dataframe with statistics per OTU
{
  stats_sk=matrix(0,10,ncol(otusl))
  colnames(stats_sk)=colnames(otusl)
  stats_sk<-as.data.frame(stats_sk)
  
  distmod<-matrix(0,6)
  newnames<-c("normal","weibull","gamma", "lnormal", "logistic", "cauchy")
  rownames(distmod)<-newnames
  row.names(stats_sk)[1]<-"min"
  row.names(stats_sk)[2]<-"max"
  row.names(stats_sk)[3]<-"median"
  row.names(stats_sk)[4]<-"mean"
  row.names(stats_sk)[5]<-"sd"
  row.names(stats_sk)[6]<-"skew"
  row.names(stats_sk)[7]<-"kurt"
  row.names(stats_sk)[8]<-"dist"
  row.names(stats_sk)[9]<-"meandist"
  row.names(stats_sk)[10]<-"sddist"
}

for(i in 1:ncol(otusl)){
  (colnames(otusl)[i])
  m=1
  j=1
  print(paste("i:",i))
  # dist=1
  t=j+1
  
  #print(paste("j:",j))
  
  m=m+1
  
  #print(paste("k:",k))
  dist=descdist(log10(otusl[,i]+5), boot = 1000, graph = FALSE)
  diptest<-dip.test(log10(otusl[,i]+5), simulate.p.value = FALSE, B = 2000)
  normal<-try(fitdist(log10(otusl[,i]+5), "norm"))
  if(inherits( normal, "try-error"))
  {
    #error handling code, maybe just skip this iteration using
    normal$aic=NA
  }
  weibull <-try(fitdist(log10(otusl[,i]+5), "weibull"))
  if(inherits( weibull, "try-error"))
  {
    #error handling code, maybe just skip this iteration using
    weibull$aic=NA
  }
  gamma <- try(fitdist(log10(otusl[,i]+5), "gamma"))
  if(inherits( gamma, "try-error"))
  {
    #error handling code, maybe just skip this iteration using
    gamma$aic=NA
  }
  lnormal <- try(fitdist(log10(otusl[,i]+5), "lnorm"))
  if(inherits( lnormal, "try-error"))
  {
    #error handling code, maybe just skip this iteration using
    lnormal$aic=NA
  }
  
  logistic <- try(fitdist(log10(otusl[,i]+5), "logis"))
  if(inherits( logistic, "try-error"))
  {
    #error handling code, maybe just skip this iteration using
    logistic$aic=NA
  }
  cauchy <- try(fitdist(log10(otusl[,i]+5), "cauchy"))
  if(inherits( cauchy, "try-error"))
  {
    #error handling code, maybe just skip this iteration using
    cauchy$aic=NA
  }
  distmod[1,]= normal$aic
  distmod[2,]= weibull$aic
  distmod[3,]= gamma$aic
  distmod[4,]= lnormal$aic
  distmod[5,]= logistic$aic
  distmod[6,]= cauchy$aic
  AICmin = rownames(which(distmod == min(distmod, na.rm=TRUE), arr.ind=TRUE))
  parameters<-paste(as.name(AICmin),"$estimate", sep="")
  eval(parse(text=parameters))
  parameters<-as.matrix(eval(parse(text=parameters)))
  
  m=m+1
  #otuss.relsa[1,i]=mean(dist)
  stats_sk[1,i]=dist$min
  stats_sk[2,i]=dist$max
  stats_sk[3,i]=dist$median
  stats_sk[4,i]=dist$mean
  stats_sk[5,i]=dist$sd
  stats_sk[6,i]=dist$skew
  stats_sk[7,i]=dist$kurt
  
  if (diptest$p.value < 0.01){
    stats_sk[8,i]= print("bimodal")
  } else {
    stats_sk[8,i]=print(AICmin)
  }
  stats_sk[9,i]=parameters[1,1]
  stats_sk[10,i]=parameters[2,1]
  
} 

######MAKE DIFFERENT OTU TABLES FOR EACH SPAD CATEGORY
str(stats_sk)
rownames(stats_sk$min)

aa<-data.frame(t(stats_sk))
str(aa)

aa$dist
all(rownames(aa)==colnames(otusl))
#aa$dist1<-aa$dist
#aa[aa$dist=="normal"|aa$dist=="weibull"|aa$dist=="cauchy"|aa$dist=="gamma",]$dist1<-"normal"
#aa<-droplevels(aa)

#Se o de cima não der certo, fazer assim
#normal<-otusl[aa$dist=="normal"|aa$dist=="weibull"|aa$dist=="cauchy"|aa$dist=="gamma"]#these distributions are normal-like
normal<-otusl[,aa$dist=="normal"]
weibull<-otusl[,aa$dist=="weibull"]
cauchy<-as.data.frame(otusl[,aa$dist=="cauchy"])
colnames(cauchy)= "OTU_3"
gamma<-as.data.frame(otusl[,aa$dist=="gamma"])
#normal.like<- cbind(normal,weibull,cauchy,gamma)


bimodal<-as.data.frame(otusl[,aa$dist=="bimodal"])

logistic<-as.data.frame(otusl[,aa$dist=="logistic"])
colnames(logistic)="OTU_6"

lognormal<-as.data.frame(otusl[,aa$dist=="lnormal"])

other<- cbind(cauchy,logistic,lognormal)

#########making a graph with ggplot
a=gather(bimodal, key = "otu", value = "bimodal")
a1=gather(a,key = "pattern", value = "abund", -otu)

b=gather(gamma, key = "otu", value = "gamma")
b1=gather(b,key = "pattern", value = "abund", -otu)

c=gather(other, key = "otu", value = "other")
c1=gather(c,key = "pattern", value = "abund", -otu)

d=rbind(a1,b1,c1)
write.table(d,"C:/Recuperar/Dropbox/Erick_Hugo/Artigos/Bimodalidade/R/novas figuras/d_total.txt")

#write.table(a1,"a1.txt")
#write.table(b1,"b1.txt")
#write.table(c1,"c1.txt")

abund=read.table(file = "C:/Recuperar/Dropbox/Erick_Hugo/Artigos/Bimodalidade/R/novas figuras/d_total.txt", header=TRUE, sep="\t", dec=".", na.strings="NA")
abund=d
abund$log10=log10(abund$abund)

ggplot(abund,aes(x=pattern,y=log10, group=pattern, fill=pattern))+
  geom_boxplot()+
  theme_bw()+
  scale_fill_manual(name="",labels=c("Bimodal","Gamma", "Other"),
                    values = c("green","lightblue","grey"))+
  scale_x_discrete(labels=c("Bimodal","Gamma","Other"))+
  xlab("")+
  scale_y_continuous(name= "Mean Abundances", limits=c(-4.5,0))


#presence-absence
#normalpa<-ifelse(normal.like>0,1,0)
#cauchypa<-ifelse(cauchy>0,1,0)
gammapa<-ifelse(gamma>0,1,0)
gammapa=as.data.frame(gammapa)
bimodalpa<-ifelse(bimodal>0,1,0)
bimodalpa=as.data.frame(bimodalpa)
#logisticpa<-ifelse(logistic>0,1,0)
#lognormalpa<-ifelse(lognormal>0,1,0)
otherpa<-ifelse(other>0,1,0)
otherpa=as.data.frame(otherpa)


a=gather(bimodalpa, key = "otu", value = "bimodal")
a1=gather(a,key = "pattern", value = "pa", -otu)

b=gather(gammapa, key = "otu", value = "gamma")
b1=gather(b,key = "pattern", value = "pa", -otu)

c=gather(otherpa, key = "otu", value = "other")
c1=gather(c,key = "pattern", value = "pa", -otu)

e=rbind(a1,b1,c1)
write.table(e,"C:/Recuperar/Dropbox/Erick_Hugo/Artigos/Bimodalidade/R/novas figuras/e_total.txt")


#write.table(a1,"a1.txt")
#write.table(b1,"b1.txt")
#write.table(c1,"c1.txt")

pa=read.table(file = "C:/Recuperar/Dropbox/Erick_Hugo/Artigos/Bimodalidade/R/novas figuras/e_total.txt", header=TRUE, sep="\t", dec=".", na.strings="NA") %>% 
  group_by(pattern,otu) %>%
  summarise(sum(pa == 1))
names(pa)[names(pa) == "sum(pa == 1)"] <- "pa" %>% as.matrix

ggplot(pa,aes(x=pattern,y=pa, group=pattern, fill=pattern))+
  geom_violin()+
  theme_bw()+
  scale_fill_manual(name="",labels=c("Bimodal","Gamma", "Other"),
                    values = c("green","lightblue","grey"))+
  scale_x_discrete(labels=c("Bimodal","Gamma","Other"))+
  xlab("")+
  ylab("Occurrence")+
  scale_y_continuous(name= "Occurrence", limits=c(0,61))


#dados corretamente categorizados
#n<-data.frame(cbind(colnames(normal.like),"normal"))
#c<-data.frame(cbind(colnames(cauchy),"cauchy"))
g<-data.frame(cbind(colnames(gamma),"gamma"))
b<-data.frame(cbind(colnames(bimodal),"bimodal"))
#lg<-data.frame(cbind(colnames(logistic),"logistic"))#or you can call these two 'rare'
#ln<-data.frame(cbind(colnames(lognormal),"lognormal"))#or you can call these two 'rare'
o<-data.frame(cbind(colnames(other),"other"))

cat1<-data.frame(rbind(b,g,o))
head(cat1)
colnames(cat1)<-c("OTUID","cat")
dim(cat1)
pos<-match(colnames(otusl),cat1$OTUID)
cat2<-cat1[pos,]
all(cat2$OTUID==colnames(otusl))
cat<-cat2
head(cat1)
cat=as.data.frame(cat)

#write.table(cat,"C:/Recuperar/Dropbox/Erick_Hugo/Artigos/Bimodalidade/R/novas figuras/core-categories.txt")
#write.table(otusl, "core-abundancies.txt")


count=count(cat, vars = cat)
head(count)


#Gráfico de pizza
ggplot(count, aes(x="", y=n, fill=vars))+
  geom_bar(width = 1, stat = "identity")+
  coord_polar("y", start=0)+
  scale_fill_manual(name="",labels=c("Bimodal","Gamma", "Other"),
                    values = c("green","lightblue","grey"))+
  theme_minimal()+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank(),panel.grid=element_blank(),axis.ticks = element_blank())+
  theme(axis.text.x=element_blank())

#write.table(count,"other_p4 - satellite.txt")

#################AQUI####################
#Manually check the abundance distribution of individual OTUs
ggplot(bimodal, aes(x=log10(OTU_9+5)))+
  geom_density(color="darkgreen",fill="green",alpha=0.4)+
  theme_bw()+
  theme(axis.text = element_blank(),axis.ticks = element_blank())+
  xlab("") + 
  ylab("")

ggplot(cauchy, aes(x=log10(OTU_3+5)))+
  geom_density(color="red4",fill="red",alpha=0.4)+
  theme_bw()+
  theme(axis.text = element_blank(),axis.ticks = element_blank())+
  xlab("") + 
  ylab("")

ggplot(gamma, aes(x=log10(OTU_19+5)))+
  geom_density(color="blue",fill="lightblue",alpha=0.4)+
  theme_bw()+
  theme(axis.text = element_blank(),axis.ticks = element_blank())+
  xlab("") + 
  ylab("")

ggplot(logistic, aes(x=log10(OTU_6+5)))+
  geom_density(color="brown3",fill="brown4",alpha=0.4)+
  theme_bw()+
  theme(axis.text = element_blank(),axis.ticks = element_blank())+
  xlab("") + 
  ylab("")

ggplot(lognormal, aes(x=log10(OTU_7+5)))+
  geom_density(color="purple4",fill="purple",alpha=0.4)+
  theme_bw()+
  theme(axis.text = element_blank(),axis.ticks = element_blank())+
  xlab("") + 
  ylab("")

########INDIVIDUAL DYNAMICS OF e.g. NORMAL-LIKE OTUS 

#I check their individual abundance distribution:
library(dplyr)
library(tidyr)
library(ggplot2)

nn=gather(as.data.frame(otusl),key = "otuID", value = "abund") %>% 
  transform(otuID=factor(otuID,levels=c("OTU_1","OTU_3","OTU_5","OTU_2","OTU_4","OTU_10","OTU_8","OTU_7","OTU_9","OTU_19","OTU_32","OTU_28","OTU_33","OTU_98","OTU_640","OTU_6","OTU_20","OTU_117","OTU_54","OTU_11","OTU_14","OTU_12","OTU_21","OTU_39")))


head(nn)

########INDIVIDUAL DYNAMICS OF NORMAL OTUs
ggplot(nn, aes(log10(abund+5)))+
  geom_histogram(color="grey50",fill="grey50",alpha=0.4)+
  geom_density(aes(y = (..scaled..)/0.1))+
  facet_wrap(~otuID,scales="free",ncol=5)+
  theme(strip.background=element_blank())+
  theme(panel.grid.major.y=element_blank(),panel.grid.minor.y=element_blank(),
        panel.grid.major.x=element_blank(),panel.grid.minor.x=element_blank(),
        axis.ticks=element_line(size=0.2),aspect.ratio=2/(1+sqrt(2)),
        text=element_text(size=8),axis.text.y=element_text(size=8),
        axis.title.y=element_text(margin=margin(0,20,0,0)),
        axis.title.x=element_text(margin=margin(20,0,0,0)),
        axis.text.x=element_text(size=8),
        plot.title=element_text(size=9, face='bold', family="Times",
                                margin=margin(0,0,20,0)))+labs(x="Log10 Abundance",y="Density")

###########################
#Environmental vs abundance tests
#mantel test
library(vegan)
library(fields)
library(dplyr)
library(tidyr)
library(ggplot2)

setwd ("")
otu=read.table(file = "otu_table.txt", header=TRUE, row.names = 1, sep="\t", dec=".", na.strings="NA")
otu=t(otu)

amb<-read.table(file="environmental_table.txt", header=T, row.names=1, sep="\t", dec=".", na.strings="NA")
denv=decostand(amb, na.rm=T, method="standardize")
geo<-read.table(file="geographic_location_table.txt", header=T, row.names=1, sep="\t", dec=".", na.strings="NA")

#Distance matrices for the abundance otu table
bac.dist <- vegdist(otu, "bray") #para diversidade melhor usar Bray

#Geographical Distance matrices
geo.dist <- rdist.earth(geo, geo, miles = FALSE, R = NULL)
rownames(geo.dist)=rownames(geo)
colnames(geo.dist)=rownames(geo)
rownames(geo.dist)==rownames(geo)
colnames(geo.dist)==rownames(geo)

geo.dist[upper.tri(geo.dist, diag=TRUE)]<-""
geo.dist<-as.dist(geo.dist)

#Mantel test Abundance vs Geography
mantel(bac.dist, geo.dist, method="spear") 

#Environmental Distance matrices
#Mantel test: Abundance vs Environment 
env.dist <- vegdist(denv[,1], "euclidean") #aqui eu escolho uma variável ambiental da minha tabela para fazer a matriz de dissimilaridade
mantel(bac.dist, env.dist, method="spear") #mantel abundancias vs altitude

env.dist <- vegdist(denv[,2], "euclidean") #aqui eu escolho uma variável ambiental da minha tabela para fazer a matriz de dissimilaridade
mantel(bac.dist, env.dist, method="spear") #mantel abundancias vs pH

env.dist <- vegdist(denv[,3], "euclidean") #aqui eu escolho uma variável ambiental da minha tabela para fazer a matriz de dissimilaridade
mantel(bac.dist, env.dist, method="spear") #mantel abundancias vs Temperatura

env.dist <- vegdist(denv[,4], "euclidean") #aqui eu escolho uma variável ambiental da minha tabela para fazer a matriz de dissimilaridade
mantel(bac.dist, env.dist, method="spear") #mantel abundancias vs Condutividade

env.dist <- vegdist(denv[,5], "euclidean") #aqui eu escolho uma variável ambiental da minha tabela para fazer a matriz de dissimilaridade
mantel(bac.dist, env.dist, method="spear") #mantel abundancias vs DOC

env.dist <- vegdist(denv[,6], "euclidean") #aqui eu escolho uma variável ambiental da minha tabela para fazer a matriz de dissimilaridade
mantel(bac.dist, env.dist, method="spear") #mantel abundancias vs DIC

env.dist <- vegdist(denv[,7], "euclidean") #aqui eu escolho uma variável ambiental da minha tabela para fazer a matriz de dissimilaridade
mantel(bac.dist, env.dist, method="spear") #mantel abundancias vs DN

env.dist <- vegdist(denv[,8], "euclidean") #aqui eu escolho uma variável ambiental da minha tabela para fazer a matriz de dissimilaridade
mantel(bac.dist, env.dist, method="spear") #mantel abundancias vs DIN

env.dist <- vegdist(denv[,9], "euclidean") #aqui eu escolho uma variável ambiental da minha tabela para fazer a matriz de dissimilaridade
mantel(bac.dist, env.dist, method="spear") #mantel abundancias vs Fluorescência

env.dist <- vegdist(denv[,10], "euclidean") #aqui eu escolho uma variável ambiental da minha tabela para fazer a matriz de dissimilaridade
mantel(bac.dist, env.dist, method="spear") #mantel abundancias vs Chlorofila a


#Make the Heatmap
setwd ("C:/Recuperar/Dropbox/Erick_Hugo/Artigos/Bimodalidade/R")
mante=read.table(file = "my/documents/mantel.txt", header=TRUE, sep="\t", dec=".", na.strings="NA")

mantel %>% 
  gather(key = "type", value = "valor", -variable) %>% 
  #transform(type=factor(type,levels=c("Core","Non.core","Satellite"))) %>% 
  transform(type=factor(type,levels=c("bimodal","gamma","other"))) %>% 
  transform(variable=factor(variable,levels=c("Geography","FC","condut","DIC","DIN","pH"))) %>% 
  ggplot(aes(x=type,y=variable)) +
  geom_tile(aes(fill = valor))+
  theme_bw()+
  #scale_x_discrete(labels=c("Core","Non-core","Satellite"))+
  scale_x_discrete(labels=c("Bimodal","Gamma","Other"))+
  scale_fill_gradient(low = "lightblue", high = "Navy")+
  xlab("")+
  ylab("")

#dbRDA
library(vegan)

setwd ("C:/Recuperar/Dropbox/Erick_Hugo/Bimodalidade/R")
otu=read.table(file = "otu_table.txt", header=TRUE, row.names = 1, sep="\t", dec=".", na.strings="NA")
otu = log(otu+1)
otu=t(otu)
amb<-read.table(file="amb_corrigido.txt", header=T, row.names=1, sep="\t", dec=".", na.strings="NA")
denv=decostand(amb, na.rm=T, method="standardize")

dbRDA=capscale(otu ~ alt+pH+temp+condut+DOC+DIC+DIN+FC+chla, denv, dist="bray")
plot(dbRDA) 
anova(dbRDA)
plot(dbRDA, type="n", xlim = c(-2, 2), ylim = c(-3, 3), ylab="CAP2 (28.14%)", xlab="CAP1 (39.24%)")
points(dbRDA, pch=16, col="black", bg="black", cex=0.8)
text(dbRDA, dis="cn", col="red3", cex=1, font = 2)
text(dbRDA,"sites", col="black", cex=0.4, pos=1)
text("p=0.001", x=-2, y=2)

summary(dbRDA)
screeplot(dbRDA)