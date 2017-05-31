#Looking into patterns among individual taxa

#rarified from 2015
datITS315a



#Pathogens 

datITS315b<-subset_taxa(datITS315a, Rank6%in%c("g__Fusarium","g__Alternaria","g__Botrytis","g__Aspergillus","g__Gibberella","g__Olpidium","g__Phoma","g__Acremonium"))#added phoma |Rank6=="g__Phoma"  adding phoma makes the pattern a little weaker, Acremonium could be added (from Mandyan dissertation)

datITS315c<-tax_glom(datITS315b,"Rank6")

otu_table(datITS315c)
tax_table(datITS315c)
sample_data(datITS315b)
taxa_names(datITS315b)

#library("data.table")

#Took this from MoistureAnalyses2015
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

spe<-datITS315aotu[,33:4575]

spepath<-spe[,which(colnames(spe)%in%taxa_names(datITS315b))]
plot(env$snowdepth,spepath$denovo6247)



spe<-t(otu_table(datITS315c))
env<-sample_data(datITS315c)
speenv<-cbind(env,spe)

speenv2<-speenv%>%
  gather(species,abun,denovo5818:denovo38212)

ggplot(speenv2,aes(x=snowdepth,y=sqrt(abun),color=species)) +
  geom_point(stat="identity") +
  geom_line(stat="smooth",method = "lm",size=.8)+
  facet_wrap(~species,scales="free")


#looking for pathogens only in DM and SB plots
dmsbplots
datITS315d<-subset_samples(datITS315c,Sample_name%in%dmsbplots$Sample_name)


spe<-t(otu_table(datITS315d))
env<-sample_data(datITS315d)
speenv<-cbind(env,spe)

speenv2<-merge(speenv,dmsbplots)

speenv3<-speenv2%>%
  gather(species,abun,denovo5818:denovo38212)%>%
  group_by(com)%>%
  summarize(meanabun=mean(abun),sebun=std.error(abun))
speenv3$com<-factor(speenv3$com,levels=c("SB","DM"))

pdf("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Proposals/NSFpreproposal2017/Figs/pathogenabundmsb.pdf",width=4,height=4)
ggplot(speenv3,aes(x=com,y=meanabun,color=com)) +
  theme(legend.position="none",text = element_text(size=20)) +
  labs(x ="Community type",y="Pathogen relative abundance") +
  geom_point(size=3) +
  geom_errorbar(aes(ymin=meanabun-sebun,ymax=meanabun+sebun),size=1,width=.5,stat="identity")
dev.off()





#mycorrhizae - all, totally overlab, DM larger variance than SB
datITS315b<-subset_taxa(datITS315a, Rank6=="g__Glomus"|Rank6=="g__Acaulospora"|Rank6=="g__Entrophosphora"|Rank6=="g__Gigaspora"|Rank6=="g__Scutelospora"|Rank6=="g__Rhizoscyphus"|Rank6=="g__Meliniomyces"|Rank6=="g__Cairneyella"|Rank6=="g__Oidiodendron"|Rank6=="g__Suillus"|Rank6=="g__Rhizopogon"|Rank6=="g__Gomphidius") #rhizoscyphum to cairneyella are ericoid, last are ectos

datITS315c<-tax_glom(datITS315b,"Rank6")

tax_table(datITS315c)

env<-data.frame(cbind(datITS315aotu[,1:32],mois=as.factor(mois)))
head(env)

spe<-datITS315aotu[,33:4575]

spemut<-spe[,which(colnames(spe)%in%taxa_names(datITS315b))]
plot(env$snowdepth,spemut$denovo815)



spe<-t(otu_table(datITS315c))
env<-sample_data(datITS315c)
speenv<-cbind(env,spe)

speenv2<-speenv%>%
  gather(species,abun,denovo6246:denovo39838)

ggplot(speenv2,aes(x=snowdepth,y=sqrt(abun),color=species)) +
  geom_point(stat="identity") +
  geom_line(stat="smooth",method = "lm",size=.8)

#looking for pathogens only in DM and SB plots
dmsbplots
datITS315d<-subset_samples(datITS315c,Sample_name%in%dmsbplots$Sample_name)

spe<-t(otu_table(datITS315d))
env<-sample_data(datITS315d)
speenv<-cbind(env,spe)

speenv2<-merge(speenv,dmsbplots)

speenv3<-speenv2%>%
  gather(species,abun,denovo6246:denovo39838)%>%
  group_by(com)%>%
  summarize(meanabun=mean(abun),sebun=std.error(abun))
speenv3$com<-factor(speenv3$com,levels=c("SB","DM"))

pdf("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Proposals/NSFpreproposal2017/Figs/mutualistabundmsb.pdf",width=4,height=4)
ggplot(speenv3,aes(x=com,y=meanabun,color=com)) +
  theme(legend.position="none",text = element_text(size=20)) +
  labs(x ="Community type",y="Pathogen relative abundance") +
  geom_point(size=3) +
  geom_errorbar(aes(ymin=meanabun-sebun,ymax=meanabun+sebun),size=1,width=.5,stat="identity")
dev.off()



#Endomycorrhizae - DM has more, but huge error bars
#datITS315b<-subset_taxa(datITS315a, Rank6=="g__Glomus"|Rank6=="g__Acaulospora"|Rank6=="g__Entrophosphora"|Rank6=="g__Gigaspora"|Rank6=="g__Scutelospora") #rhizoscyphum to cairneyella are ericoid, last are ectos

datITS315b<-subset_taxa(datITS315a, Rank2=="p__Glomeromycota")

datITS315c<-tax_glom(datITS315b,"Rank6")

tax_table(datITS315c)

env<-data.frame(cbind(datITS315aotu[,1:32],mois=as.factor(mois)))
head(env)

spe<-datITS315aotu[,33:4575]

spe<-t(otu_table(datITS315c))
env<-sample_data(datITS315c)
speenv<-cbind(env,spe)

speenv2<-speenv%>%
  gather(species,abun,denovo446:denovo39838)

ggplot(speenv2,aes(x=snowdepth,y=sqrt(abun),color=species)) +
  geom_point(stat="identity") +
  geom_line(stat="smooth",method = "lm",size=.8)

#dmsbplots
datITS315d<-subset_samples(datITS315c,Sample_name%in%dmsbplots$Sample_name)

spe<-t(otu_table(datITS315d))
env<-sample_data(datITS315d)
speenv<-cbind(env,spe)

speenv2<-merge(speenv,dmsbplots)

speenv3<-speenv2%>%
  gather(species,abun,denovo446:denovo39838)%>%
  group_by(com)%>%
  summarize(meanabun=mean(abun),sebun=std.error(abun))
speenv3$com<-factor(speenv3$com,levels=c("SB","DM"))

pdf("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Proposals/NSFpreproposal2017/Figs/AMFabundmsb.pdf",width=4,height=4)
ggplot(speenv3,aes(x=com,y=meanabun,color=com)) +
  theme(legend.position="none",text = element_text(size=20)) +
  labs(x ="Community type",y="Pathogen relative abundance") +
  geom_point(size=3) +
  geom_errorbar(aes(ymin=meanabun-sebun,ymax=meanabun+sebun),size=1,width=.5,stat="identity")
dev.off()


#ericoid, snowbed has more ericoid. pattern dissappears when you add Capronia, then overall the DM has more ericoid
datITS315b<-subset_taxa(datITS315a, Rank6%in%c("g__Rhizoscyphus","g__Meliniomyces","g__Cairneyella","g__Oidiodendron","g__Capronia","g__Sabacina","g__Hymenoscyphus")) #capronia to hymenoscyphus from Allen et al 2003
datITS315c<-tax_glom(datITS315b,"Rank6")

tax_table(datITS315c)

env<-data.frame(cbind(datITS315aotu[,1:32],mois=as.factor(mois)))
head(env)

spe<-datITS315aotu[,33:4575]

spe<-t(otu_table(datITS315c))
env<-sample_data(datITS315c)
speenv<-cbind(env,spe)

speenv2<-speenv%>%
  gather(species,abun,denovo6246:denovo36871)

ggplot(speenv2,aes(x=snowdepth,y=sqrt(abun),color=species)) +
  geom_point(stat="identity") +
  geom_line(stat="smooth",method = "lm",size=.8)

dmsbplots
datITS315d<-subset_samples(datITS315c,Sample_name%in%dmsbplots$Sample_name)

spe<-t(otu_table(datITS315d))
env<-sample_data(datITS315d)
speenv<-cbind(env,spe)

speenv2<-merge(speenv,dmsbplots)

speenv3<-speenv2%>%
  gather(species,abun,denovo6246:denovo36871)%>%
  group_by(com)%>%
  summarize(meanabun=mean(abun),sebun=std.error(abun))
speenv3$com<-factor(speenv3$com,levels=c("SB","DM"))

pdf("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Proposals/NSFpreproposal2017/Figs/ericoidabundmsb.pdf",width=4,height=4)
ggplot(speenv3,aes(x=com,y=meanabun,color=com)) +
  theme(legend.position="none",text = element_text(size=20)) +
  labs(x ="Community type",y="Ericoid relative abundance") +
  geom_point(size=3) +
  geom_errorbar(aes(ymin=meanabun-sebun,ymax=meanabun+sebun),size=1,width=.5,stat="identity")
dev.off()



#ectos - not sure if I want to go there b/c these are ectos, not endophytes???
datITS315b<-subset_taxa(datITS315a, Rank6=="g__Suillus"|Rank6=="g__Rhizopogon"|Rank6=="g__Gomphidius"|Rank6=="g__Amanita"|Rank6=="g__Laccaria"|Rank6=="g__Lactarius"|Rank6=="g__Pisolithus"|Rank6=="g__Alpova"|Rank6=="g__Sebacina"|Rank6=="g__Descolea"|Rank6=="g__Dicymbe") #rhizoscyphum to oidio are ericoid, last are ectos. Suillus is the only one that is really abundant in high snow plots, it is associated with pines. this was a random list from wikipedia or something

datITS315b<-subset_taxa(datITS315a, Rank6%in%c("g__Cenococcum","g__Phialocephala","g__Cadophora","g__Cortinarius","g__Tomentella","g__Caloplaca","g__Hymenoscyphus","g__Leohumicola","g__Russula","g__Amanita","g__Boletus","g__Sebacina","g__Inocybe","g__Hebeloma","g__Laccaria","g__Lactarius","g__Hymenoscyphus","g__Lachnum","g__Acephala","g__Phialocephala","g__Leptodontidium","g__Pleiochaeta","g__Cistella","g__Naeviopsis","g__Microglossum","g__Hyalacrotes","g__Pseudeurotium"))  #taken from a paper about Dryas and a paper about Kobresia

datITS315c<-tax_glom(datITS315b,"Rank6")

tax_table(datITS315c)

env<-data.frame(cbind(datITS315aotu[,1:32],mois=as.factor(mois)))
head(env)

spe<-datITS315aotu[,33:4575]

spe<-t(otu_table(datITS315c))
env<-sample_data(datITS315c)
speenv<-cbind(env,spe)

speenv2<-speenv%>%
  gather(species,abun,denovo5109:denovo28824)

ggplot(speenv2,aes(x=snowdepth,y=sqrt(abun),color=species)) +
  geom_point(stat="identity") +
  geom_line(stat="smooth",method = "lm",size=.8) +
  facet_wrap(~species)

#dmsbplots
datITS315d<-subset_samples(datITS315c,Sample_name%in%dmsbplots$Sample_name)

spe<-t(otu_table(datITS315d))
env<-sample_data(datITS315d)
speenv<-cbind(env,spe)

speenv2<-merge(speenv,dmsbplots)

speenv3<-speenv2%>%
  gather(species,abun,denovo1242:denovo44425)%>%
  group_by(com)%>%
  summarize(meanabun=mean(abun),sebun=std.error(abun))
speenv3$com<-factor(speenv3$com,levels=c("SB","DM"))

pdf("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Proposals/NSFpreproposal2017/Figs/mutualistabundmsb.pdf",width=4,height=4)
ggplot(speenv3,aes(x=com,y=meanabun,color=com)) +
  theme(legend.position="none",text = element_text(size=20)) +
  labs(x ="Community type",y="Pathogen relative abundance") +
  geom_point(size=3) +
  geom_errorbar(aes(ymin=meanabun-sebun,ymax=meanabun+sebun),size=1,width=.5,stat="identity")
dev.off()




#some DSE taxa, only two present

datITS315b<-subset_taxa(datITS315a, Rank6%in%c("g__Cadophora","g__Rhizopycnis","g__Periconia","g__Curvularia","g__Microdochium","g__Phialophora","g__Cladosporium","g__Phialocephala","g__Leptodontidium","g__Lecythophora"))  #first five taken from knapp et al 2012 dark side not fastidious, Phialophora taken from Schadt work in Ranunculus adoneus, Cladosporium taken from Terri's paper, (Microdochium and Periconia from Mandyham dissertation too), phialocephala might also be a pathogen from Reininger and Sieber 2012, leptodontidium and lecythophora from Day dn Curah 2011

datITS315c<-tax_glom(datITS315b,"Rank6")

tax_table(datITS315c)

env<-data.frame(cbind(datITS315aotu[,1:32],mois=as.factor(mois)))
head(env)

spe<-datITS315aotu[,33:4575]

spe<-t(otu_table(datITS315c))
env<-sample_data(datITS315c)
speenv<-cbind(env,spe)

speenv2<-speenv%>%
  gather(species,abun,denovo3284:denovo44096)

ggplot(speenv2,aes(x=snowdepth,y=sqrt(abun),color=species)) +
  geom_point(stat="identity") +
  geom_line(stat="smooth",method = "lm",size=.8) +
  facet_wrap(~species)

#dmsbplots
datITS315d<-subset_samples(datITS315c,Sample_name%in%dmsbplots$Sample_name)

spe<-t(otu_table(datITS315d))
env<-sample_data(datITS315d)
speenv<-cbind(env,spe)

speenv2<-merge(speenv,dmsbplots)

speenv3<-speenv2%>%
  gather(species,abun,denovo3284:denovo44096)%>%
  group_by(com)%>%
  summarize(meanabun=mean(abun),sebun=std.error(abun))
speenv3$com<-factor(speenv3$com,levels=c("SB","DM"))

#pdf("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Proposals/NSFpreproposal2017/Figs/mutualistabundmsb.pdf",width=4,height=4)
ggplot(speenv3,aes(x=com,y=meanabun,color=com)) +
  theme(legend.position="none",text = element_text(size=20)) +
  labs(x ="Community type",y="DSE relative abundance") +
  geom_point(size=3) +
  geom_errorbar(aes(ymin=meanabun-sebun,ymax=meanabun+sebun),size=1,width=.5,stat="identity")
#dev.off()


