#Data cleaning for working with 2007 and 2015 soil data from the King plots

##### FUNGI #####

##### Read in ITS files, filtered with singletons removed #####
otufileITS <-import_biom("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/Figures&Stats/kingdata/ITS/ITS_ALL_97_OTU_tablefiltsingnonchimericFunSoilsing.biom") #throws an error but it still looks like it imports fine
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

#Rarefy everything here then subset years below

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



##### For Katie and Marko AGU talk #####
# I"m not fitering for plots that have data in both years, I'm just giving them everything.
datITS3
datITS3otu<-cbind(sample_data(datITS3),t(otu_table(datITS3)))
datITS3otu[1:10,30:40]
head(tax_table(datITS3))

datITS3otu2<-datITS3otu[,-c(2:4,6:7,9:32)]

write.csv(datITS3otu2,"/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/Figures&Stats/kingdata/ITS/datITS3otu2.csv",row.names = F)
write.csv(tax_table(datITS3),"/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/Figures&Stats/kingdata/ITS/datITS3otu_taxtable.csv")


#aggregating at order level, maybe for future analyses
datITS3order<-tax_glom(datITS3,taxrank=rank_names(datITS3)[4])



##### Back to my stuff #####
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









##### BACTERIA #####

otufile16S <-import_biom("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/Figures&Stats/kingdata/Bact/16S_ALL_97_OTU_tablefiltsingnonchimerickeepbactarcfiltmitchlSoilreadsing.biom") 
head(tax_table(otufile16S))
head(otu_table(otufile16S))

#Import mapping and tree file #need to put snow in the mapping file if I want to work on that
map16S<-import_qiime_sample_data("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/Figures&Stats/kingdata/MappingFiles/515BC_Niwot_20072015_All_MapFilenewlomehi.txt")
dat16S<-merge_phyloseq(otufile16S,map16S)

#Rarefy everything here then subset years below

min(sample_sums(dat16S))#rarefy to 2699
dat16S2<-rarefy_even_depth(dat16S,sample.size=min(sample_sums(dat16S)),rngseed=10,replace=F) #22837 OTUs were removed because they are no longer present in any sample after random subsampling
dat16S3 = transform_sample_counts(dat16S2, function(x) x/sum(x) )

head(tax_table(dat16S3))
unique(tax_table(dat16S3)[,"Rank2"])

#remove the p__ with substring
#labelsITS<-substring(tax_table(datITS3)[,"Rank2"],4)
#colnames(labelsITS)<-"labels"

#replace tax table, datITS3 (relativized) datITS2 (not relativized)
#tax_table(datITS3)<-cbind(tax_table(datITS3),labelsITS)
#tax_table(datITS2)<-cbind(tax_table(datITS2),labelsITS)



##### For Katie and Marko AGU talk #####
# I"m not fitering for plots that have data in both years, I'm just giving them everything.
dat16S3
dat16S3otu<-cbind(sample_data(dat16S3),t(otu_table(dat16S3))) #takes about 10-15 min
dat16S3otu[1:10,30:40]
head(tax_table(dat16S3))

dat16S3otu2<-dat16S3otu[,-c(2:4,6:7,9:31)]
dat16S3otu2[order(dat16S3otu2)]

dat16S3otu3<-dat16S3otu2[order( as.numeric(as.character(dat16S3otu2[,2])), dat16S3otu2[,3] ),]

write.csv(dat16S3otu3,"/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/Figures&Stats/kingdata/Bact/dat16S3otu3.csv",row.names = F)
write.csv(tax_table(dat16S3),"/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/Figures&Stats/kingdata/Bact/dat16S3otu_taxtable.csv")






##### EUKS SOIL #####

Euk_ALL_97_S111_OTU_tablefiltsingnonchimericbactarcplantEukSoilnometazoanfungising.biom
otufile18S <-import_biom("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/Figures&Stats/kingdata/Euks/Euk_ALL_97_S111_OTU_tablefiltsingnonchimericbactarcplantEukSoilnometazoanfungireadsing.biom") #throw an error but it looks like it imports fine
head(tax_table(otufile18S))
head(otu_table(otufile18S))

#Import mapping and tree file #need to put snow in the mapping file if I want to work on that
map18S<-import_qiime_sample_data("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/Figures&Stats/kingdata/MappingFiles/EukBr_Niwot_20072015_All_MapFilenewlomehi.txt")
dat18S<-merge_phyloseq(otufile18S,map18S)

#take out sample with 490 reads
#GP.chl = prune_samples(sample_sums(GP.chl)>=20, GP.chl)
dat18Sb <- prune_samples(sample_names(dat18S)!="S.61.2015", dat18S)


#Rarefy everything here then subset years below

min(sample_sums(dat18Sb))#rarefy to 810
dat18S2<-rarefy_even_depth(dat18Sb,sample.size=min(sample_sums(dat18Sb)),rngseed=10,replace=F) #3049 OTUs were removed because they are no longer present in any sample after random subsampling
dat18S3 = transform_sample_counts(dat18S2, function(x) x/sum(x) )

head(tax_table(dat18S3))

#remove the p__ with substring
#labelsITS<-substring(tax_table(datITS3)[,"Rank2"],4)
#colnames(labelsITS)<-"labels"

#replace tax table, datITS3 (relativized) datITS2 (not relativized)
#tax_table(datITS3)<-cbind(tax_table(datITS3),labelsITS)
#tax_table(datITS2)<-cbind(tax_table(datITS2),labelsITS)



##### For Katie and Marko AGU talk #####
# I"m not fitering for plots that have data in both years, I'm just giving them everything.
dat18S3
dat18S3otu<-cbind(sample_data(dat18S3),t(otu_table(dat18S3))) #takes a few minutes
dat18S3otu[1:10,30:40]
head(tax_table(dat18S3))

dat18S3otu2<-dat18S3otu[,-c(2:4,6:7,9:31)]

dat18S3otu3<-dat18S3otu2[order( as.numeric(as.character(dat18S3otu2[,2])), dat18S3otu2[,3] ),]

write.csv(dat18S3otu3,"/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/Figures&Stats/kingdata/Euks/dat18S3otu3.csv",row.names = F)
write.csv(tax_table(dat18S3),"/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/Figures&Stats/kingdata/Euks/dat18S3otu3_taxtable.csv")




##### EUKS NEMATODE #####


otufile18N <-import_biom("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/Figures&Stats/kingdata/Euks/Euk_ALL_97_S111_OTU_tablefiltsingnonchimericbactarcplantEukNematodemetazoancransing.biom") #throw an error but it looks like it imports fine
head(tax_table(otufile18N))
head(otu_table(otufile18N))

#Import mapping and tree file #need to put snow in the mapping file if I want to work on that
map18N<-import_qiime_sample_data("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/Figures&Stats/kingdata/MappingFiles/EukBr_Niwot_20072015_All_MapFilenewlomehi.txt")
dat18N<-merge_phyloseq(otufile18N,map18N)

#take out samples 2A and 78. 2A is a replicate and only has 3 taxa in it which is odd. 78 only has 164 reads
dat18Nb <- prune_samples(sample_names(dat18N)!="N.2A.2015", dat18N)
dat18Nc <- prune_samples(sample_names(dat18Nb)!="N.78.2015", dat18Nb)




#Rarefy everything here then subset years below

min(sample_sums(dat18Nc))#rarefy to 692
dat18N2<-rarefy_even_depth(dat18Nc,sample.size=min(sample_sums(dat18Nc)),rngseed=10,replace=F) #2601 OTUs were removed because they are no longer present in any sample after random subsampling
dat18N3 = transform_sample_counts(dat18N2, function(x) x/sum(x) )

head(tax_table(dat18N3))

#remove the p__ with substring
#labelsITS<-substring(tax_table(datITS3)[,"Rank2"],4)
#colnames(labelsITS)<-"labels"

#replace tax table, datITS3 (relativized) datITS2 (not relativized)
#tax_table(datITS3)<-cbind(tax_table(datITS3),labelsITS)
#tax_table(datITS2)<-cbind(tax_table(datITS2),labelsITS)



##### For Katie and Marko AGU talk #####
# I"m not fitering for plots that have data in both years, I'm just giving them everything.
dat18N3
dat18N3otu<-cbind(sample_data(dat18N3),t(otu_table(dat18N3))) #takes a few minutes
dat18N3otu[1:10,30:40]
head(tax_table(dat18N3))

dat18N3otu2<-dat18N3otu[,-c(2:4,6:7,9:31)]

dat18N3otu3<-dat18N3otu2[order( as.numeric(as.character(dat18N3otu2[,2])), dat18N3otu2[,3] ),]

write.csv(dat18N3otu3,"/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/Figures&Stats/kingdata/Euks/dat18N3otu3.csv",row.names = F)
write.csv(tax_table(dat18N3),"/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/Figures&Stats/kingdata/Euks/dat18N3otu3_taxtable.csv")





