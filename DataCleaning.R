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



##### Unrarified data #####

#2015 only
datITS15<-subset_samples(datITS,year==2015)
datITS15a<-filter_taxa(datITS15, function(x) sum(x) > (0), prune=T)

datITS15b = transform_sample_counts(datITS15a, function(x) x/sum(x) )

datITS15aotu<-cbind(sample_data(datITS15a),t(otu_table(datITS15a)))
datITS15aotu$Sample_name<-as.numeric(as.character(datITS15aotu$Sample_name))

datITS15botu<-cbind(sample_data(datITS15b),t(otu_table(datITS15b)))
datITS15botu$Sample_name<-as.numeric(as.character(datITS15botu$Sample_name))

#07 and 15
datITS4<-subset_samples(datITS,Sample_name%in%inboth)
datITS5<-filter_taxa(datITS4, function(x) sum(x) > (0), prune=T)

datITS6 = transform_sample_counts(datITS5, function(x) x/sum(x) )

datITS5otu<-cbind(sample_data(datITS5),t(otu_table(datITS5)))
datITS5otu$Sample_name<-as.numeric(as.character(datITS5otu$Sample_name))

datITS6otu<-cbind(sample_data(datITS6),t(otu_table(datITS6)))
datITS6otu$Sample_name<-as.numeric(as.character(datITS6otu$Sample_name))


##### Rarefied data ####

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


#make dataset of only 2015 data, there are many more samples in 2015 compared to 2007, so single year analyses should be done on 2015

datITS315<-subset_samples(datITS3,year==2015)
datITS315a<-filter_taxa(datITS315, function(x) sum(x) > (0), prune=T)

datITS315aotu<-cbind(sample_data(datITS315a),t(otu_table(datITS315a)))
datITS315aotu$Sample_name<-as.numeric(as.character(datITS315aotu$Sample_name))


#2007 and 2015
datITS7<-subset_samples(datITS3,Sample_name%in%inboth)
datITS8<-filter_taxa(datITS7, function(x) sum(x) > (0), prune=T)

datITS9 = transform_sample_counts(datITS8, function(x) x/sum(x) )

datITS9otu<-cbind(sample_data(datITS9),t(otu_table(datITS9)))
datITS9otu$Sample_name<-as.numeric(as.character(datITS9otu$Sample_name))






