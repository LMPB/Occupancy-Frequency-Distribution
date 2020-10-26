###Occupancy-frequency bimodal distribution test
##This script was cordially provided by Markus Lindh
##Reference DOI: 10.1111/1462-2920.13650

library(vegan)
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggExtra)

bimodal=read.table(file = "C:/The/way/to/your/data/otu_table_cleared_rarefied.txt", header=TRUE, sep="\t", dec=".", na.strings="NA")
bimodal=bimodal[,-1]

x=t(bimodal)
x[which(x > 0)] = 1

a=
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
write.table(table.out,file=paste("bimodal.txt",sep="_"),sep="\t",row.names=F)


#test.input=as.data.frame(table(apply(x, 2, sum)[which(apply(x, 2, sum) >0)]))
MOS.test=MOStest(as.numeric(table.out[,1]), 
                 as.numeric(table.out[,2]),
                 family=gaussian(link = "log"), 
                 maxit = 100)
print(MOS.test)
capture.output(MOS.test,file=paste("MOStest.txt",sep="_")) #Saves table with stats


otus=read.table(file = "C:/The/way/to/your/data/otu_table_cleared_rarefied.txt", header=TRUE,row.names = 1, sep="\t", dec=".", na.strings="NA")
otus.ra <- otus/colSums(otus)
otus.ra[,61] <- rowMeans(otus.ra)
otus.ra[,62] <- rowSums(otus.ra[1:60] > 0)
otus.ra <- otus.ra[61:62]
colnames(otus.ra)= c("ab","occ")

b=
ggplot(otus.ra,aes(x=occ, y=ab))+
  geom_jitter()+
  ylab(paste("Mean Abundance (%)"))+
  xlab(paste("Occupancy","(Number of sites occupied)",sep="\n"))+
  ggtitle("")+
  theme_bw()

x = arrangeGrob(a, b,
                ncol = 2, nrow = 1,
                layout_matrix = rbind(c(1,2)))
as_ggplot(x) +
  draw_plot_label(label = c ("A","B", "*"), size = 15,
                  x = c(0, 0.5, 0.4655), y = c(1, 1, 0.355))


###SpADs - Spatial Abundance Distributions
##This script was cordially provided by Juan Pablo Niño-Garcia and Clara Ruíz-Gonzáléz
##Reference DOI: 10.1111/ele.12704
##Reference DOI: 10.1111/mec.15026


#################################
#Environment vs composition tests
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
setwd (" ")
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

setwd (" ")
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
