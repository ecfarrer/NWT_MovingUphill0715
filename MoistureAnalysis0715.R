

datasets: datITS9otu (rarified, relaabun)

spe<-datITS9otu[33:5524]
datITS9otu$Sample_name<-as.factor(datITS9otu$Sample_name)
env<-datITS9otu[,1:32]


##### Moisture Analysis 07-15 ####

#datITS9otu %>% 
  #select(moisture,year,Sample_name) %>%
  #group_by(year,Sample_name) %>%
  #summarise(mean=mean(moisture))

ggplot(datITS9otu,aes(x=year,y=moisture,col=Sample_name)) +
  geom_point(stat="identity") +
  geom_line()

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



##### Snow Analysis 07-15: did plots with different snow depth respond differently over time ####

#group by snow depth
head(env)
hist(env$snowdepth)
sort(env$snowdepth)

snow<-ifelse(env$snowdepth<120,1,0)
sum(snow)
ind<-which(env$snowdepth>=120&env$snowdepth<210)
length(ind)
snow[ind]<-2
ind<-which(env$snowdepth>=210&env$snowdepth<500)
length(ind)
snow[ind]<-3

env$snow<-as.factor(snow)

#take doubletons and singletons out
ind<-which(colSums(spe>0)>2)
spe2<-spe[,ind]



##### Ordination ####

#NMDS
m1<-metaMDS(spe2, distance="bray", k=2, autotransform=F,trymax=200)

col<-ifelse(env$snow==1&env$year==2007,"blue",NA)
ind<-which(env$snow==1&env$year==2015)
col[ind]<-"lightblue"
ind<-which(env$snow==2&env$year==2007)
col[ind]<-"darkgreen"
ind<-which(env$snow==2&env$year==2015)
col[ind]<-"lightgreen"
ind<-which(env$snow==3&env$year==2007)
col[ind]<-"darkorange4"
ind<-which(env$snow==3&env$year==2015)
col[ind]<-"tan1"

env$yearsnow<-as.factor(paste(env$year,env$snow,sep="."))

plot(scores(m1),col=col,pch=21,bg=col)
ordiellipse(m1,env$yearsnow,col=c("blue","darkgreen","darkorange4","lightblue","lightgreen","tan1"),conf=.99999,kind="se",lwd=2)#
#ordiellipse(m1,env$yearsnow,col=c("blue","darkgreen","darkorange4","lightblue","lightgreen","tan1"))#
#ordiellipse(m1,env$yearsnow,col=c("blue","darkgreen","darkorange4","lightblue","lightgreen","tan1"),kind="ehull")#
#ordihull(m1,env$yearsnow,col=c("blue","darkgreen","darkorange4","lightblue","lightgreen","tan1"))
text(scores(m1),labels=env$Sample_name)

adonis(spe2~year+snow+year*snow, data=env, method="bray", by="terms")


#dbRDA
m2<-capscale(spe2~year+snow+year*snow+Condition(Sample_name), distance="bray",data=env,na.action = na.omit)

anova(m2,by="terms")

plot(m2)

col<-ifelse(env$snow==1,"blue",ifelse(env$snow==2,"green","orange"))
col<-ifelse(env$year==2007,"red","blue")

plot(scores(m2)$sites,col=col,bg=col,pch=21)
text(scores(m2)$centroids,labels=c("Vlow","Low","Med","High"),col=1:4)



