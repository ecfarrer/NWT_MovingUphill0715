Taxa 2007-2015

#rarified
datITS9

#spe<-datITS9otu[33:5524]
#datITS9otu$Sample_name<-as.factor(datITS9otu$Sample_name)
#env<-datITS9otu[,1:32]



#Pathogens -  also have opposite pattern, they increase in snowbeds but stay the same in DM

datITS9b<-subset_taxa(datITS9, Rank6=="g__Fusarium"|Rank6=="g__Alternaria"|Rank6=="g__Botrytis"|Rank6=="g__Aspergillus"|Rank6=="g__Gibberella"|Rank6=="g__Olpidium"|Rank6=="g__Phoma"|Rank6=="g__Acremonium")

datITS9b
datITS9c<-tax_glom(datITS9b,"Rank6")

tax_table(datITS9c)
sample_data(datITS9c)

datITS9d<-subset_samples(datITS9c,Sample_name%in%dmsbplots$Sample_name)

spe<-t(otu_table(datITS9d))
env<-sample_data(datITS9d)
speenv<-cbind(env,spe)

speenv2<-merge(speenv,dmsbplots)

speenv3<-speenv2%>%
  gather(species,abun,denovo5818:denovo35175)%>%
  group_by(com,year)%>%
  summarize(meanabun=mean(abun),sebun=std.error(abun))
speenv3$com<-factor(speenv3$com,levels=c("SB","DM"))

ggplot(speenv3,aes(x=com,y=meanabun,color=year)) +
  #theme(legend.position="none",text = element_text(size=20)) +
  labs(x ="Community type",y="Pathogen relative abundance") +
  geom_point(size=3) +
  geom_errorbar(aes(ymin=meanabun-sebun,ymax=meanabun+sebun),size=1,width=.5,stat="identity")




#AMF- results do not mesh with predictions, AMF increase much more in the DM plots

datITS9b<-subset_taxa(datITS9,Rank2=="p__Glomeromycota")
datITS9b
datITS9c<-tax_glom(datITS9b,"Rank6")

tax_table(datITS9c)
sample_data(datITS9c)


datITS9d<-subset_samples(datITS9c,Sample_name%in%dmsbplots$Sample_name)

spe<-t(otu_table(datITS9d))
env<-sample_data(datITS9d)
speenv<-cbind(env,spe)

speenv2<-merge(speenv,dmsbplots)

speenv3<-speenv2%>%
  gather(species,abun,denovo446:denovo41340)%>%
  group_by(com,year)%>%
  summarize(meanabun=mean(abun),sebun=std.error(abun))
speenv3$com<-factor(speenv3$com,levels=c("SB","DM"))

ggplot(speenv3,aes(x=com,y=meanabun,color=year)) +
  #theme(legend.position="none",text = element_text(size=20)) +
  labs(x ="Community type",y="AMF relative abundance") +
  geom_point(size=3) +
  geom_errorbar(aes(ymin=meanabun-sebun,ymax=meanabun+sebun),size=1,width=.5,stat="identity")





#Ericoid, snowbed has more ericoid
datITS9b<-subset_taxa(datITS9,Rank6%in%c("g__Rhizoscyphus","g__Meliniomyces","g__Cairneyella","g__Oidiodendron","g__Capronia","g__Sabacina","g__Hymenoscyphus")) 

datITS9b
datITS9c<-tax_glom(datITS9b,"Rank6")

tax_table(datITS9c)

datITS9d<-subset_samples(datITS9c,Sample_name%in%dmsbplots$Sample_name)

spe<-t(otu_table(datITS9d))
env<-sample_data(datITS9d)
speenv<-cbind(env,spe)

speenv2<-merge(speenv,dmsbplots)

speenv3<-speenv2%>%
  gather(species,abun,denovo6246:denovo24500)%>%
  group_by(com,year)%>%
  summarize(meanabun=mean(abun),sebun=std.error(abun))
speenv3$com<-factor(speenv3$com,levels=c("SB","DM"))

ggplot(speenv3,aes(x=com,y=meanabun,color=year)) +
  #theme(legend.position="none",text = element_text(size=20)) +
  labs(x ="Community type",y="Ericoid relative abundance") +
  geom_point(size=3) +
  geom_errorbar(aes(ymin=meanabun-sebun,ymax=meanabun+sebun),size=1,width=.5,stat="identity")








#DSE
datITS9b<-subset_taxa(datITS9, Rank6%in%c("g__Cadophora","g__Rhizopycnis","g__Periconia","g__Curvularia","g__Microdochium","g__Phialophora","g__Cladosporium","g__Phialocephala","g__Leptodontidium","g__Lecythophora"))  #first five taken from knapp et al 2012 dark side not fastidious, Phialophora taken from Schadt work in Ranunculus adoneus, Cladosporium taken from Terri's paper, (Microdochium and Periconia from Mandyham dissertation too), phialocephala might also be a pathogen from Reininger and Sieber 2012, leptodontidium and lecythophora from Day dn Curah 2011
datITS9b
datITS9c<-tax_glom(datITS9b,"Rank6")

tax_table(datITS9c)
sample_data(datITS9c)

datITS9d<-subset_samples(datITS9c,Sample_name%in%dmsbplots$Sample_name)

spe<-t(otu_table(datITS9d))
env<-sample_data(datITS9d)
speenv<-cbind(env,spe)

speenv2<-merge(speenv,dmsbplots)

speenv3<-speenv2%>%
  gather(species,abun,denovo3284:denovo41077)%>%
  group_by(com,year)%>%
  summarize(meanabun=mean(abun),sebun=std.error(abun))
speenv3$com<-factor(speenv3$com,levels=c("SB","DM"))

ggplot(speenv3,aes(x=com,y=meanabun,color=year)) +
  #theme(legend.position="none",text = element_text(size=20)) +
  labs(x ="Community type",y="DSE relative abundance") +
  geom_point(size=3) +
  geom_errorbar(aes(ymin=meanabun-sebun,ymax=meanabun+sebun),size=1,width=.5,stat="identity")



