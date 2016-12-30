

##### Get data together #####
#datasets: datITS9otu (rarified, relaabun), datITS6otu (not rarified, relabun)

#rarifieid
spe<-datITS9otu[33:5524]
datITS9otu$Sample_name<-as.factor(datITS9otu$Sample_name)
env<-datITS9otu[,1:32]

#not rarified
spe<-datITS6otu[33:12650]
datITS6otu$Sample_name<-as.factor(datITS9otu$Sample_name)
env<-datITS6otu[,1:32]



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

#select plots that only got drier
ind<-which(env2$X2015<env2$X2007)
env3<-env2[ind,]
env2[-ind,]

dim(env2)

ind<-which(env$Sample_name%in%env3$Sample_name)



#Did the plots that got drier all have the same response in microbial communtiy (answer: no)

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



#Compare plots that got drier to plots that got wetter

env3<-env
env3$gotdrier<-NA
env3$gotdrier[ind]<-"yes"
env3$gotdrier[-ind]<-"no"
env3$yeargotdrier<-paste(env3$year,env3$gotdrier,sep=".")

#take doubletons and singletons out
ind<-which(colSums(spe>0)>2)
spe5<-spe[,ind]

m1<-metaMDS(spe5, distance="bray", k=2, autotransform=F,trymax=200)
plot(scores(m1),col=as.numeric(as.factor(env3$yeargotdrier)))#ifelse(env3$gotdrier=="yes","red","blue")
#text(scores(m1),labels=env5$Sample_name)
ordiellipse(m1,env3$yeargotdrier,col=1:4,conf=.99999,kind="se",lwd=2)#

#yes 2007-2015: red -> blue
#no 2007-2015: black -> green

adonis(spe5~year*gotdrier,data=env3)
#there is an affect of year and gotdrier (which is odd) but no sig interaction. the year effect is the same across plots (2015 is different regardless of whether plots got wetter or drier). Thus I think this indicates that our moisture measurement is not very useful.









##### Snow Analysis 07-15: did plots with different snow depth respond differently over time ####

#Group by snow depth
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

#write.csv(envsnow,"")



##### Ordination ####

#Take doubletons and singletons out
ind<-which(colSums(spe>0)>2)
spe2<-spe[,ind]


#old colors
# col<-ifelse(env$snow==1&env$year==2007,"red3",NA)
# ind<-which(env$snow==1&env$year==2015)
# col[ind]<-"indianred1"
# ind<-which(env$snow==2&env$year==2007)
# col[ind]<-"darkorange4"
# ind<-which(env$snow==2&env$year==2015)
# col[ind]<-"tan1"
# ind<-which(env$snow==3&env$year==2007)
# col[ind]<-"blue"
# ind<-which(env$snow==3&env$year==2015)
# col[ind]<-"lightblue"
# 

http://www.hexcolortool.com/#ab3b57


#NMDS
m1<-metaMDS(spe2, distance="bray", k=2, autotransform=F,trymax=200)

col<-ifelse(env$snow==1&env$year==2007,"#ab3b57",NA)
ind<-which(env$snow==1&env$year==2015)
col[ind]<-"#db95a6"
ind<-which(env$snow==2&env$year==2007)
col[ind]<-"#5268cb"
ind<-which(env$snow==2&env$year==2015)
col[ind]<-"#a9b4e5"
ind<-which(env$snow==3&env$year==2007)
col[ind]<-"#639e51"
ind<-which(env$snow==3&env$year==2015)
col[ind]<-"#abcea1"

env$yearsnow<-as.factor(paste(env$year,env$snow,sep="."))

plot(scores(m1),col=col,pch=21,bg=col)
ordiellipse(m1,env$yearsnow,col=c("#ab3b57","#5268cb","#639e51","#db95a6","#a9b4e5","#abcea1"),conf=.99999,kind="se",lwd=2,lty=c(1,1,1,2,2,2))#
legend("topleft",c("Long GSL","Med GSL","Short GSL"),col=c("#ab3b57","#5268cb","#639e51"),pch=21,pt.bg=c("#ab3b57","#5268cb","#639e51"),lty=1,bty="n")
text(-.15,.74,"Shift over time")
segments(-.4,.6,-.2,.6,lwd=2)
arrows(-.2,.6,.1,.6,lwd=2,lty=2)
#ordiellipse(m1,env$yearsnow,col=c("red3","darkorange4","blue","indianred1","tan1","lightblue"),conf=.99999,kind="se",lwd=2)#
#ordiellipse(m1,env$yearsnow,col=c("blue","darkgreen","darkorange4","lightblue","lightgreen","tan1"),kind="ehull")#
#ordihull(m1,env$yearsnow,col=c("blue","darkgreen","darkorange4","lightblue","lightgreen","tan1"))
text(scores(m1),labels=env$Sample_name)


adonis(spe2~year+snow+year*snow, data=env, method="bray", by="terms")


#Dissimilarities
dis0715<-data.frame(Sample_name=rep(NA,69),dist=rep(NA,69))
for(i in 1:69){
  plots<-sort(unique(env$Sample_name))
  tempplot<-plots[i]
  tempind<-which(env$Sample_name==tempplot)
  tempspe<-spe2[tempind,]
  dis0715[i,1]<-as.character(tempplot)
  dis0715[i,2]<-vegdist(tempspe)
}

env2<-merge(env,dis0715,"Sample_name")
env3<-env2%>%
  filter(year==2015)%>%
  group_by(snow)%>%
  summarise(meandist=mean(dist),sedist=std.error(dist))
env3

pdf("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Proposals/NSFpreproposal2017/Figs/microbe0715dissimilaritysnowdepth.pdf",width=4,height=4)
ggplot(env3,aes(x=snow,y=meandist,col=snow)) +
  theme(legend.position="none",text = element_text(size=20)) +
  labs(x ="Snow depth",y="Dissimilarity 2007-2015") +
  geom_point(size=3) +
  geom_errorbar(aes(ymin=meandist-sedist,ymax=meandist+sedist),size=1,width=.5,stat="identity")
dev.off()

#This is from the network code, using Guem and SB dominants cutoffs for community type
dmplots<-c(0,104,106,113,120,122,132,133,135,18,36,66,7,74,79,9,98)
sbplots<-c( 107,109,111,119,123,124,125,152,154,155,157,30,31,75 ,83,84,89)
dmsbplots<-data.frame(cbind(Sample_name=c(dmplots,sbplots),com=rep(c("DM","SB"),each=17)))
env5<-merge(env2,dmsbplots)

env6<-env5%>%
  filter(year==2015)%>%
  group_by(com)%>%
  summarise(meandist=mean(dist),sedist=std.error(dist))
env6$com<-factor(env6$com,levels=c("SB","DM"))

pdf("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Proposals/NSFpreproposal2017/Figs/microbe0715dissimilaritydmsb.pdf",width=4,height=4)
ggplot(env6,aes(x=com,y=meandist,col=com)) +
  theme(legend.position="none",text = element_text(size=20)) +
  labs(x ="Community type",y="Dissimilarity 2007-2015") +
  geom_point(size=3) +
  geom_errorbar(aes(ymin=meandist-sedist,ymax=meandist+sedist),size=1,width=.5,stat="identity")
dev.off()

env4<-env2%>%
  filter(year==2015)

ggplot(env4,aes(x=snowdepth,y=dist)) +
  theme(legend.position="none",text = element_text(size=20)) +
  labs(x ="Snow depth",y="Disimilarity between 2007, 2015") +
  geom_point(size=3) +
  geom_smooth(method="lm",se=F)




#dbRDA
m2<-capscale(spe2~year+year*snow+Condition(Sample_name), distance="bray",data=env,na.action = na.omit)

anova(m2,by="margin")

anova(m2,by="margin",permutations=how(blocks=env$Sample_name,nperm=1000))
anova(m2,by="terms",permutations=how(blocks=env$Sample_name,nperm=1000))

plot(m2)

plot(scores(m2)$sites,pch=21,bg=col,col=col)
ordiellipse(m2,env$yearsnow,col=c("red3","darkorange4","blue","indianred1","tan1","lightblue"),conf=.99999,kind="se",lwd=2)#

col<-ifelse(env$snow==1,"blue",ifelse(env$snow==2,"green","orange"))
col<-ifelse(env$year==2007,"red","blue")

plot(scores(m2)$sites,col=col,bg=col,pch=21)
text(scores(m2)$centroids,labels=c("Vlow","Low","Med","High"),col=1:4)



##### Overlap and Richness #####

#Richness
env$rich<-rowSums(spe>0)
aggregate.data.frame(env$rich,by=list(env$year,env$snow),mean)

env2<-env%>%
  group_by(Sample_name)

ggplot(env2,aes(x=year,y=rich,col=snow)) +
  geom_point(stat="identity") +
  geom_line(aes(group=Sample_name))

env2<-env%>%
  group_by(snow,year)%>%
  summarise(mean=mean(rich),se=std.error(rich))

ggplot(env2,aes(x=year,y=mean,col=snow)) +
  geom_point(stat="identity") +
  geom_line()

env2<-env%>%
  group_by(year)

ggplot(env2,aes(x=snowdepth,y=rich,col=year)) +
  geom_point(stat="identity") +
  geom_line(stat="smooth",method = "lm",size=.8,aes(group=year))


#Overlap
spe2007lo<-spe[which(env$snow==1&env$year==2007),]
spe20015lo<-spe[which(env$snow==1&env$year==2015),]
spe2007me<-spe[which(env$snow==2&env$year==2007),]
spe2015me<-spe[which(env$snow==2&env$year==2015),]
spe2007hi<-spe[which(env$snow==3&env$year==2007),]
spe2015hi<-spe[which(env$snow==3&env$year==2015),]

otu2007lo<-colnames(spe2007lo)[(which(colSums(spe2007lo)>0))]
otu2015lo<-colnames(spe20015lo)[(which(colSums(spe20015lo)>0))]
otu2007me<-colnames(spe2007me)[(which(colSums(spe2007me)>0))]
otu2015me<-colnames(spe2015me)[(which(colSums(spe2015me)>0))]
otu2007hi<-colnames(spe2007hi)[(which(colSums(spe2007hi)>0))]
otu2015hi<-colnames(spe2015hi)[(which(colSums(spe2015hi)>0))]

#sorenson, probably more common
(2*length(intersect(otu2007lo,otu2015lo)))/(length(intersect(otu2007lo,otu2015lo))+length(union(otu2007lo,otu2015lo)))
(2*length(intersect(otu2007me,otu2015me)))/(length(intersect(otu2007me,otu2015me))+length(union(otu2007me,otu2015me)))
(2*length(intersect(otu2007hi,otu2015hi)))/(length(intersect(otu2007hi,otu2015hi))+length(union(otu2007hi,otu2015hi)))

length(setdiff(otu2007lo,otu2015lo))#unique taxa to 2007 lo
length(setdiff(otu2015lo,otu2007lo))
length(setdiff(otu2007me,otu2015me))
length(setdiff(otu2015me,otu2007me))
length(setdiff(otu2007hi,otu2015hi))
length(setdiff(otu2015hi,otu2007hi))

aggregate.data.frame(env$snowdepth,by=list(env$snow),mean)


#jaccard
length(intersect(otuvlo,otulo))/length(union(otuvlo,otulo))
length(intersect(otulo,otume))/length(union(otulo,otume))
length(intersect(otume,otuhi))/length(union(otume,otuhi))


#Calculating species overlap on a pairwise plot basis

plots<-sort(unique(env$Sample_name))
ov<-rep(NA,69)

for (i in 1:69){
  curplot<-plots[i]
  temp07<-subset(spe,env$Sample_name==curplot&env$year==2007)
  temp15<-subset(spe,env$Sample_name==curplot&env$year==2015)
  sp07<-colnames(temp07)[(which(colSums(temp07)>0))]
  sp15<-colnames(temp15)[(which(colSums(temp15)>0))]
  ov[i]<-(2*length(intersect(sp07,sp15)))/(length(intersect(sp07,sp15))+length(union(sp07,sp15)))
  }

pairwise<-data.frame(plots,ov)
pairwise$Sample_name<-as.factor(pairwise$plots)
envsub<-subset(env,year==2007)
pairwise2<-merge(pairwise,envsub,by="Sample_name")
plot(pairwise2$snowdepth,pairwise2$ov)
abline(lm(pairwise2$ov~pairwise2$snowdepth))
aggregate.data.frame(pairwise2$ov,by=list(pairwise2$snow),mean)
aggregate.data.frame(pairwise2$ov,by=list(pairwise2$snow),std.error)

pairwise3<-pairwise2%>%
  group_by(snow)
means<-pairwise3%>%
  summarise(mean=mean(ov),se=std.error(ov))
ggplot(means,aes(x=snow,y=mean)) +
  geom_bar(stat="identity") +
  geom_errorbar(aes(ymax=mean+se,ymin=mean-se),width=.25)

anova(gls(ov~snowdepth,data=pairwise2))


