#Looking into patterns among individual taxa

#rarified from 2015
datITS315a

datITS315b<-subset_taxa(datITS315a, Rank6=="g__Fusarium"|Rank6=="g__Alternaria"|Rank6=="g__Botrytis"|Rank6=="g__Aspergillus"|Rank6=="g__Gibberella"|Rank6=="g__Olpidium")

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
  gather(species,abun,denovo14211:denovo38212)

ggplot(speenv2,aes(x=snowdepth,y=sqrt(abun),color=species)) +
  geom_point(stat="identity") +
  geom_line(stat="smooth",method = "lm",size=.8)
  #facet_wrap(~species,scales="free")


#looking for pathogens only in DM and SB plots
dmsbplots
datITS315d<-subset_samples(datITS315c,Sample_name%in%dmsbplots$Sample_name)


spe<-t(otu_table(datITS315d))
env<-sample_data(datITS315d)
speenv<-cbind(env,spe)

speenv2<-merge(speenv,dmsbplots)

speenv3<-speenv2%>%
  gather(species,abun,denovo14211:denovo35175)%>%
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



