#Data cleaning for working with 2007 and 2015 soil data from the King plots


##### Read in ITS files, filtered with singletons removed #####
otufileITS <-import_biom("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/Figures&Stats/kingdata/ITS/ITS_ALL_97_OTU_tablefiltsingnonchimericFunSoilsing.biom")
head(tax_table(otufileITS))
head(otu_table(otufileITS))

#Import mapping and tree file
mapITS<-import_qiime_sample_data("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/Figures&Stats/kingdata/MappingFiles/ITS_Niwot_20072015_All_MapFilenewlomehisnow.txt")

datITS<-merge_phyloseq(otufileITS,mapITS)

#Figure out which plots we have data for in both 2007 and 2015
temp<-sample_data(datITS)
temp07<-subset(temp,year==2007)
temp15<-subset(temp,year==2015)

inboth<-sort(as.numeric(intersect(temp07$Sample_name,temp15$Sample_name)))

#setdiff(temp07$Sample_name,temp15$Sample_name)
#setdiff(temp15$Sample_name,temp07$Sample_name)



#****
#take out other samples for final reads and OTU counts
#datITSfincount<-subset_samples(datITS,Sample_name!=81&Sample_name!=61&Sample_name!=33&Sample_name!=56&Sample_name!=78&Sample_name!=126&Sample_name!=5&Sample_name!=34)
#datITSfincount2<-filter_taxa(datITSfincount, function(x) sum(x) > (0), prune=T)
#sum(otu_table(datITSfincount2))
#****


#unrarified data
datITS15<-subset_samples(datITS,year==2015)
datITS15a<-filter_taxa(datITS15, function(x) sum(x) > (0), prune=T)

datITS15aotu<-cbind(sample_data(datITS15a),t(otu_table(datITS15a)))
datITS15aotu$Sample_name<-as.numeric(as.character(datITS15aotu$Sample_name))


#rarefy and transform to relative abundance
min(sample_sums(datITS))#rarefy to 1000
datITS2<-rarefy_even_depth(datITS,sample.size=min(sample_sums(datITS)),rngseed=10,replace=F) #7181 OTUs were removed because they are no longer present in any sample after random subsampling
datITS3 = transform_sample_counts(datITS2, function(x) x/sum(x) )

head(tax_table(datITS3))
unique(tax_table(datITS3)[,"Rank2"])

#remove the p__ with substring
#labelsITS<-substring(tax_table(datITS3)[,"Rank2"],4)

#colnames(labelsITS)<-"labels"

#replace tax table, datITS3 (relativized) datITS2 (not relativized)
#tax_table(datITS3)<-cbind(tax_table(datITS3),labelsITS)
#tax_table(datITS2)<-cbind(tax_table(datITS2),labelsITS)


#make dataset of only 2015 data, there are many more samples in 2015 compared to 2007, so single year analyses should be done on that data

datITS315<-subset_samples(datITS3,year==2015)
datITS315a<-filter_taxa(datITS315, function(x) sum(x) > (0), prune=T)

datITS315aotu<-cbind(sample_data(datITS315a),t(otu_table(datITS315a)))
datITS315aotu$Sample_name<-as.numeric(as.character(datITS315aotu$Sample_name))



