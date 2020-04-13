###Occupancy-frequency bimodal distribution test
##This script was cordially provided by Markus Lindh
##Reference DOI: 10.1111/1462-2920.13650


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
