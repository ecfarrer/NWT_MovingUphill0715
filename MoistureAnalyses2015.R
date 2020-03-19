#Moisture analyses from 2015

#phyloseq object: datITS315a (rarefied, relabun), datITS15a (not rarefied, counts), datITS15b (not rarefied, relabun)
#otu table: datITS315aotu (rarefied, relabun), datITS15aotu (not rarefied,counts), datITS15botu (not rarefied, relabun)


hist(datITS315aotu$moisture)

head(datITS315aotu)[,1:40]


# Moisture analysis

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
spe<-datITS15botu[,33:11230] #not rarefied
head(spe)


#take out doubletons and singletons
ind<-which(colSums(spe>0)>2)
spe2<-spe[,ind]



##### Ordination #####

#NMDS
m1<-metaMDS(spe2, distance="bray", k=2, autotransform=F,trymax=200)
plot(scores(m1),col=env$mois)

#dbRDA
#moisture and snowdepth are not correlated, and if anything are inversely correlated
m2<-capscale(spe2~snowdepth+VascPlant_Dens, distance="bray",data=env,na.action = na.omit)#+VascPlant_Dens+snowdepth
m2
plot(m2)

anova(m2,by="margin")

plot(scores(m2)$sites,col=env$mois,bg=env$mois,pch=21)
text(scores(m2)$centroids,labels=c("Vlow","Low","Med","High"),col=1:4)

col<-ifelse(env$mois==1,"#9350a1",ifelse(env$mois==2,"#62ad64",ifelse(env$mois==3,"#697cd4","#b8475f")))
#vlow purple, low green, med blue, high red

plot(scores(m2)$sites,col=col,bg=col,pch=21,cex=2)
text(scores(m2)$centroids,labels=c("Vlow","Low","Med","High"),col=c("#9350a1","#62ad64","#697cd4","#b8475f"),cex=2)
text(scores(m2)$sites,labels=env$Sample_name)


#plot with snowdepth and vasc plants for Liz Kimborough

env<-cbind(env,snow)
write.csv(spe2,"/Users/farrer/Desktop/spe2.csv")
write.csv(env,"/Users/farrer/Desktop/env.csv")

spe2<-read.csv("/Users/farrer/Desktop/spe2.csv",row.names=1)
env<-read.csv("/Users/farrer/Desktop/env.csv",row.names=1)

m2<-capscale(spe2~snowdepth+VascPlant_Dens, distance="bray",data=env,na.action = na.omit)#+VascPlant_Dens+snowdepth
m2
plot(m2)

col<-ifelse(env$mois==1,"#9350a1",ifelse(env$mois==2,"#62ad64",ifelse(env$mois==3,"#697cd4","#b8475f")))
#colored by soil moisture: vlow purple, low green, med blue, high red

par(mai=c(1,1,.42,.42))
plot(scores(m2)$sites,col=col,bg=col,pch=21,xlab=paste("Axis 1 (",sprintf("%.1f",m2$CCA$eig["CAP1"]/m2$tot.chi*100,3),"%)",sep=""),ylab=paste("Axis 2 (",sprintf("%.1f",m2$CCA$eig["CAP2"]/m2$tot.chi*100,3),"%)",sep=""),cex.lab=1.4,ylim=c(-2.5,2),xlim=c(-2.5,2))
abline(h=0,col=gray(.70))
abline(v=0,col=gray(.70))
text(m2$CCA$biplot*2-.25,labels=c("Snowdepth","PlantDensity"),col=1) #I multiplied the centroid scores by 2 so the graph would look prettier, as long as they are multiplied by a constant, it doesn't change how you would read the result. the (-.25) just makes the name slighly offset from the head of the arrow created below
arrows(0,0,m2$CCA$biplot[,1]*2,m2$CCA$biplot[,2]*2,length=.05,col=1,lwd=2) #these scores have to be multiplid by 2 as well





##### Calculating overlap in composition #####

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





# Snow analysis

#Group by snow depth

snow<-ifelse(datITS315aotu$snowdepth<155,1,0)
sum(snow)
ind<-which(datITS315aotu$snowdepth>=155&datITS315aotu$snowdepth<262)
length(ind)
snow[ind]<-2
ind<-which(datITS315aotu$snowdepth>=262&datITS315aotu$snowdepth<500)
length(ind)
snow[ind]<-3

env<-data.frame(cbind(datITS315aotu[,1:32],snow=as.factor(snow)))
head(env)

spe<-datITS315aotu[,33:4575] #rarefied
#spe<-datITS15botu[,33:11230] #not rarefied
head(spe)

env$snow<-as.factor(snow)

#take out doubletons and singletons
ind<-which(colSums(spe>0)>2)
spe2<-spe[,ind]


#NMDS

m1<-metaMDS(spe2, distance="bray", k=3, autotransform=F,trymax=200)

col<-ifelse(env$snow==1,"#ab3b57",NA) #low snow is red
ind<-which(env$snow==2) #medium snow is blue
col[ind]<-"#5268cb"
ind<-which(env$snow==3) #high snow is green
col[ind]<-"#639e51"

#env$yearsnow<-as.factor(paste(env$year,env$snow,sep="."))

#I just changed the width for the 2020 NSF EAGER version
pdf("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Proposals/NSF/NSF2020/Figs/fungisnowordination2015.pdf",height=5,width=4.3)
plot(scores(m1)[,1:2],col=col,pch=21,bg=col,cex=1.5)
ordiellipse(m1,env$snow,col=c("#ab3b57","#5268cb","#639e51"),conf=.99999,kind="se",lwd=2.5,lty=c(1,1,1))#
legend("topleft",c("Low snow","Med snow","High snow"),col=c("#ab3b57","#5268cb","#639e51"),pch=21,pt.bg=c("#ab3b57","#5268cb","#639e51"),lty=1,bty="n",y.intersp=1)
dev.off()

#text(-.15,.74,"Shift over time")
#segments(-.4,.6,-.2,.6,lwd=2)
#arrows(-.2,.6,.1,.6,lwd=2,lty=2)

adonis(spe2~snow, data=env, method="bray", by="terms")




m2<-capscale(spe2~snow, distance="bray",data=env,na.action = na.omit)#snowdepth+moisture+WHC
m2
plot(m2)

