MoistureAnalysis0715

datasets: datITS9otu (rarified, relaabun)

spe<-datITS9otu[33:5524]

#datITS9otu %>% 
  #select(moisture,year,Sample_name) %>%
  #group_by(year,Sample_name) %>%
  #summarise(mean=mean(moisture))

datITS9otu$Sample_name<-as.factor(datITS9otu$Sample_name)

ggplot(datITS9otu,aes(x=year,y=moisture,col=Sample_name)) +
  geom_point(stat="identity") +
  geom_line()

env<-datITS9otu[,1:32]
head(env)
env2<-env %>%
  select(Sample_name,year,moisture)%>%
  spread(year,moisture)
colnames(env2)<-c("Sample_name","X2007","X2015")

ind<-which(env2$X2015<env2$X2007)
env3<-env2[ind,]
env2[-ind,]

dim(env2)

ind<-which(env$Sample_name%in%env3$Sample_name)
spe2<-spe[ind,]
env4<-env[ind,]
env4$mois<-ifelse(env4$Sample_name%in%c(2,3,25,26,30,47,65),"lo","hi")

#Ordination
m1<-metaMDS(spe2, distance="bray", k=2, autotransform=F,trymax=200)
plot(scores(m1),col=ifelse(env4$year==2015,"blue","red"))
text(scores(m1),labels=env4$Sample_name)

m2<-capscale(spe2~year, distance="bray",data=env4,na.action = na.omit)
m2
plot(m2)
adonis(spe2~year,data=env4)
anova(m2,by="margin")

plot(scores(m2)$sites,col=env$mois,bg=env$mois,pch=21)
text(scores(m2)$centroids,labels=c("Vlow","Low","Med","High"),col=1:4)
