#Moisture analyses from 2015

#phyloseq object: datITS315a (rarefied), datITS15a (not rarefied)
#otu table: datITS315aotu, datITS15aotu (not rarefied)


datITS315aotu

hist(datITS315aotu$moisture)

head(datITS315aotu)[,1:40]

mois<-ifelse(datITS315aotu$moisture<7.5,1,0)
ind<-which(datITS315aotu$moisture>=7.5&datITS315aotu$moisture<11)
length(ind)
mois[ind]<-2
ind<-which(datITS315aotu$moisture>=11&datITS315aotu$moisture<16)
length(ind)
mois[ind]<-3
ind<-which(datITS315aotu$moisture>=16&datITS315aotu$moisture<100)
length(ind)
mois[ind]<-4

#In moist meadow plots, soil mosture in july ranged from 22-170, mean = 65

env<-data.frame(cbind(datITS315aotu[,1:32],mois=as.factor(mois)))
head(env)

spe<-datITS315aotu[,33:4575] #rarefied
spe<-datITS15aotu[,33:11230] #not rarefied
head(spe)

#take out doubletons and singletons
ind<-which(colSums(spe>0)>2)
spe2<-spe[,ind]


#Ordination

m1<-metaMDS(spe2, distance="bray", k=2, autotransform=F,trymax=200)
plot(scores(m1),col=env$mois)

m2<-capscale(spe2~mois, distance="bray",data=env,na.action = na.omit)
m2
plot(scores(m2)$sites,col=env$mois,bg=env$mois,pch=21)
text(scores(m2)$centroids,labels=c("Vlow","Low","Med","High"),col=1:4)




#Calculating overlap in composition

spevlo<-spe[which(env$mois==1),]
spelo<-spe[which(env$mois==2),]
spemed<-spe[which(env$mois==3),]
spehi<-spe[which(env$mois==4),]

otuvlo<-colnames(spevlo)[(which(colSums(spevlo)>0))]
otulo<-colnames(spelo)[(which(colSums(spelo)>0))]
otume<-colnames(spemed)[(which(colSums(spemed)>0))]
otuhi<-colnames(spehi)[(which(colSums(spehi)>0))]

#jaccard
length(intersect(otuvlo,otulo))/length(union(otuvlo,otulo))
length(intersect(otulo,otume))/length(union(otulo,otume))
length(intersect(otume,otuhi))/length(union(otume,otuhi))

length(intersect(otuvlo,otuhi))/length(union(otuvlo,otuhi))

#sorenson, probably more common
(2*length(intersect(otuvlo,otulo)))/(length(intersect(otuvlo,otulo))+length(union(otuvlo,otulo)))
(2*length(intersect(otulo,otume)))/(length(intersect(otulo,otume))+length(union(otulo,otume)))
(2*length(intersect(otume,otuhi)))/(length(intersect(otume,otuhi))+length(union(otume,otuhi)))

(2*length(intersect(otuvlo,otuhi)))/(length(intersect(otuvlo,otuhi))+length(union(otuvlo,otuhi)))

aggregate.data.frame(env$moisture,by=list(env$mois),mean)





