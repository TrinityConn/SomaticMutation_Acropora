

#R script  for organizing exported allele frequency, mutation and location id information for all called mutations in all combinations.

#load in required packages 
library(ggplot2)
library(tidyverse)
library(data.table)
library(cowplot)
library(vegan)
library(TSA)
library(wesanderson)

#loading in finalized VAF file. These files were exported with the example column names "# [1]CHROM", "[2]POS" "[3]145_1:AF",
# "[4]145_1:DP", "[5]145_3:AF", "[6]145_3:DP"

#set working directory 

setwd("/Users/trini/Documents/PSU/Acropora_Manuscript/VAF_All")


#read in files in directory 
filelist=list.files(pattern="*_VAF*")

datalist=lapply(filelist, function(x)read.table(x, header=FALSE))

#add in column for filename 
for (i in 1:length(datalist)){datalist[[i]]<-cbind(datalist[[i]],filelist[i])}

#create function to label columns based on if column 5 sum is greater than column 3 
tumor2<-function(x) {if (sum(x$V5) > sum(x$V3)) {
  colnames(x)[3]="Normal_AF"
  colnames(x)[4]="Normal_DP"
  colnames(x)[5]="Tumor_AF"
  colnames(x)[6]="Tumor_DP"
  x
}
}

#create function to label columns based on if column 3 is greater than column 5 
tumor3<-function(x) {if (sum(x$V3) > sum(x$V5)) {
    colnames(x)[3]="Tumor_AF"
    colnames(x)[4]="Tumor_DP"
    colnames(x)[5]="Normal_AF"
    colnames(x)[6]="Normal_DP"
    x
}
}

#apply function tumor3 to data list
data2=lapply(datalist, tumor3) 


#apply function tumor2 to data list 
data3=lapply(datalist, tumor2)

#combine the two lists 
data4<-c(data2,data3)

#remove duplicate 0 length files 
data5<-data4[lapply(data4, length)>0]


#combine rows of all imported files
data_old=do.call("rbind", data5)

#add in column names 

colnames(data_old)[1]="CHRO"
colnames(data_old)[2]="POS"


#add in normal v tumor classifications 
#load in metadata

data_names<-as.data.frame(unique(data_old$`filelist[i]`))
colnames(data_names)[1]="FileNames"
write.csv(data_names, file="Data_Names.csv")

metadata<-read.csv("./Mutect_Metadata.csv")

#merge metadata with data 
colnames(metadata)[1]="filelist[i]"
mutations<-merge(data_old, metadata, by=c("filelist[i]"))

write.csv(mutations, file="Mutations_Merged.csv")

#remove samples that had error -- these are A41,A43,A44,A45,A46,A47,A48,A50,A51,A52
#A54,A55,57,A59,A62,A66,A69,A70,A71,A72,A79

problem_samples<-c('A41', 'A43', 'A44', 'A45', 'A46', 'A47', 'A48', 'A50', 'A51', 'A52', 'A53',
                   'A54', 'A55', 'A57', 'A59', 'A62', 'A66', 'A69', 'A70', 'A71', 'A72', 'A79',
                   'A75', 'A58', 'Y22','Y23', 'Y24', 'Y25', 'Y26', 'Y27', 'Y28')

problem<-as.data.frame(problem_samples)
colnames(problem)[1]="Comparison"

mutations2<-mutations[ !mutations$Comparison %in% problem_samples, ]

#calculating mutation burden per sample 
#for each 'tumor' sample, how many rows that contain unique CHRO and POS 
#create new dataframe with tumor sample and the mutation burden 

#separate cpg and whole genome mutations 
whole<-mutations2%>%
  filter(FileType=="all")

cpg<-mutations2%>%
  filter(FileType=="cpg")


############Identifying Fixed Mutations########################

#filter out mutations so that it is just normal AF less than or equal to 0.01, 
#tumor VAF greater than or equal to 0.5, tumor depth greater than or equal 
#to 23 and normal depth greater than or equal to 12
fixed<-whole%>%
  filter(Normal_AF<=0.01)%>%
  filter(Tumor_AF>=0.5)%>%
  filter(Normal_DP>=12)%>%
  filter(Tumor_DP>=23)


#filter mutations by depth 

whole_filter<-whole%>%
  filter(Normal_DP>=12)%>%
  filter(Tumor_DP>=60)

cpg_filter<-cpg%>%
  filter(Normal_DP>=12)%>%
  filter(Tumor_DP>=23)


#mutation burden per sample for whole dataset 

burdenwhole1<-whole_filter%>%
  group_by(Tumor)%>%
  summarise(count=n_distinct(CHRO,POS))

burdenwholeu<-unique(whole_filter)

burdenwhole2<-burdenwholeu%>%
  group_by(Tumor)%>%
  summarise(count=n_distinct(CHRO,POS))


#mutation burden for sites in cpg heavy regions 

burdencpg1<-cpg_filter%>%
  group_by(Tumor)%>%
  summarise(count=n_distinct(CHRO,POS))

burdencpg2<-cpg_filter%>%
  group_by(Comparison)%>%
  summarise(count=n_distinct(CHRO,POS))


##merge in age and colony information 

age_metadata<-mutations2[, c('Age', 'Colony', 'Comparison', 'Tumor', 'Normal')]

burden_whole<-merge(age_metadata, burdenwhole2, by=c("Comparison"))
burden_whole<-unique(burden_whole)

burden_cpg<-merge(age_metadata, burdencpg2, by=c("Comparison"))
burden_cpg<-unique(burden_cpg)

colnames(burden_whole)[6]="Mutations Called"
colnames(burden_cpg)[6]="Mutations Called"
ggplot(burden_whole, aes(x=Age, y=`Mutations Called`, fill=Age))+geom_boxplot()+
  theme_bw()+scale_fill_manual(values=wes_palette("GrandBudapest2", n=4))

ggplot(burden_cpg, aes(x=Age, y=`Mutations Called`, fill=Age))+geom_boxplot()+geom_jitter()+theme_bw()+scale_fill_grey()

burden_whole$Age<- factor(burden_whole$Age, levels=c("6", "8", "10", ">100"))
burden_whole$Age[is.na(burden_whole$Age)]<-">100"
burden_cpg$Age<- factor(burden_cpg$Age, levels=c("6", "8", "10", ">50"))

write.csv(burden_whole, file="burden_whole.csv")
write.csv(burden_cpg, file="burden_cpg.csv")
write.csv(mutations2, file="Mutations2.csv")

#summary_tables

VAF<-whole_filter%>%
  group_by(Age)%>%
  summarise(mean=mean(Tumor_AF), sd=sd(Tumor_AF))

Mutations<-burden_whole%>%
  group_by(Age)%>%
  summarise(mean=mean(`Mutations Called`), sd=sd(`Mutations Called`))

Mutations

###plotting Variant Allele Frequencies 

VAF_unique<-unique(whole_filter$POS)
V2<-unique(whole_filter)

V2_6<-V4%>%filter(Age=="6")
V2_8<-V4%>%filter(Age=="8")
V2_10<-V4%>%filter(Age=="10")
V2_100<-V4%>%filter(Colony=="P2194")

ggplot(V2_6, aes(x=Tumor_AF, ))+geom_histogram(bins=50)+xlim(0,1.0)+theme_bw()
ggplot(V2_8, aes(x=Tumor_AF))+geom_histogram(bins=50)+xlim(0,1.0)+theme_bw()


ggplot(V2_10, aes(x=Tumor_AF))+geom_histogram(bins=50)+xlim(0,1.0)+theme_bw()
ggplot(V2_100, aes(x=Tumor_AF))+geom_histogram(bins=50)+xlim(0,1.0)+theme_bw()

ks.test(V2_10$Tumor_AF, V2_100$Tumor_AF)
ks<-ks.boot(V2_10$Tumor_AF, V2_100$Tumor_AF, nboots=10000)
 
minMax<-seq(min(V2_6$Tumor_AF, V2_10$Tumor_AF), max(V2_6$Tumor_AF, V2_10$Tumor_AF), length.out=length(V2_6$Tumor_AF))
x0<-minMax[which( abs(cdf1(minMax)-cdf2(minMax)) == max(abs(cdf1(minMax) - cdf2(minMax))) )]
y0<-cdf1(x0)
y1<-cdf2(x0)

plot(cdf1, verticals=TRUE, do.points=FALSE, col="blue") 
plot(cdf2, verticals=TRUE, do.points=FALSE, col="green", add=TRUE) 

## alternatine, use standard R plot of ecdf 
#plot(f.a, col="blue") 
#lines(f.b, col="green") 

points(c(x0, x0), c(y0, y1), pch=16, col="red") 
segments(x0, y0, x0, y1, col="red", lty="dotted") 

cdf1<-ecdf(V2_6$Tumor_AF)
cdf2<-ecdf(V2_10$Tumor_AF)

V3<-V2%>%
  filter(Tumor_AF>=0.1)%>%
  filter(Tumor_AF<=0.49)

VAF_unique2<-VAF_unique%>%
  filter(Normal_AF<=0.01)

V4<-V2%>%
  filter(Normal_DP>=12)%>%
  filter(Tumor_DP>=60)
V4[is.na(V4)]="100"

#calculating best bin-width 

#
bw<-2 * IQR(V2$Tumor_AF/length(V2$Tumor_AF)^1/3)
ggplot(V4, aes(Tumor_AF, fill=Age))+geom_histogram(aes(y=..density..), color='gray50',
                                                             position="identity", bins=30)+geom_density(alpha=0.2)+
  theme_bw()+scale_fill_manual(values=wes_palette("GrandBudapest2", n=4))+facet_wrap(~Age)+xlab("VAF")+xlim(0, 1.0)

ggplot(Keri_Kinship, aes(dist,))+geom_histogram(aes(y=..density..), color='gray50', position="identity")+geom_density(alpha=0.2)+
  theme_bw()+xlim(-0.5, .5)+xlab("Kinship Coefficient")

V4

ggplot(V4, aes(x=Age, y=Tumor_AF, fill=Age))+geom_violin(alpha=0.2)+theme_bw()+
  geom_boxplot(width=0.1, color='black')

V4<-V4%>%
  filter(Tumor_AF>0.2)%>%
  filter(Tumor_AF<0.5)
V4$Age<- factor(V4$Age, levels=c("6", "8", "10", ">100"))
V4$Colony<-factor(V4$Colony, levels=c("M132", "M145", "M148", "S19", "S15", "S18", "S37" ,"P2194"))
V2$Age<-factor(V2$Age, levels=c("6", "8", "10", ">100"))
V4$Age[is.na(V4$Age)]<-">100"

write.csv(V4, file="Final_AlleleFrequencies_Acropora.csv")
##cpg as a proportion of all mutations 

# per comparison, what proportion of all mutations are cpg mutations 

##merge cpg sum and all sum 

TotalSum<-merge(burden_cpg, burden_whole, by=c("Comparison", "Age", "Colony", "Normal", "Tumor"))
colnames(TotalSum)[6]="cpg"
colnames(TotalSum)[7]="all"

TotalSum$Proportion_CpG<-TotalSum$cpg/TotalSum$all

Median<-TotalSum%>%
  group_by(Colony)%>%
  summarise(mean=mean(Proportion_CpG), sd=sd(Proportion_CpG))

Median

prop_colony<-ggplot(TotalSum, aes(x=Colony, y=Proportion_CpG, fill=Colony))+geom_boxplot()+
  theme_bw()+scale_fill_manual(values=wes_palette("GrandBudapest2", n=8, type=c("continuous")))

prop_age<-ggplot(TotalSum, aes(x=Age, y=Proportion_CpG, fill=Age))+geom_boxplot()+theme_bw()+
  scale_fill_manual(values=wes_palette("GrandBudapest2",n=4))+theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.4))
  
TotalSum$Colony<-factor(TotalSum$Colony, levels=c("M132", "M145", "M148", "S19", "S15", "S18", "S37" ,"P2194"))
TotalSum$Age<-factor(TotalSum$Age, levels=c("6", "8", "10", ">100"))


plot_grid(prop_colony, prop_age, labels=c("A", "B"))





#diversity statistics 

#first pivot table from long to wide 



whole_filter$mutid<-paste(whole_filter$CHRO, whole_filter$POS)

M145<-whole_filter%>%filter(Colony=="M145")
M132<-whole_filter%>%filter(Colony=="M132")
M148<-whole_filter%>%filter(Colony=="M148")
S15<-whole_filter%>%filter(Colony=="S15")
S18<-whole_filter%>%filter(Colony=="S18")
S19<-whole_filter%>%filter(Colony=="S19")
S37<-whole_filter%>%filter(Colony=="S37")
P2194<-whole_filter%>%filter(Colony=="P2194")



#make a column that is the sum of the number of times that mutations is found in a tumor site 
#for each unique muteid and tumor sample combination, how many rows int eh datafram contain that combination and then make a 
#column that describes that column 

#calculating alpha diversity 

S18_n<-S18%>%
  count(mutid, Tumor)

M145_n<-M145%>%
  count(mutid, Tumor)

M132_n<-M132%>%
  count(mutid, Tumor)

M148_n<-M148%>%
  count(mutid, Tumor)

S15_n<-S15%>%
  count(mutid, Tumor)

S37_n<-S37%>%
  count(mutid, Tumor)

S19_n<-S19%>%
  count(mutid, Tumor)

P2194_n<-P2194%>%
  count(mutid, Tumor)

P2194_alpha<-P2194_n%>%
  group_by(Tumor)%>%
  summarize(sobs=specnumber(n), shannon=diversity(n, index="shannon"), simpson=diversity(n, index="simpson"), invsimpson=diversity(n, index="invsimpson"), 
           Pielou=diversity(n, index="shannon")/log(specnumber(n)))



S15_alpha<-S15_n%>%
  group_by(Tumor)%>%
  summarize(sobs=specnumber(n), shannon=diversity(n, index="shannon"), simpson=diversity(n, index="simpson"), invsimpson=diversity(n, index="invsimpson"), 
            Pielou=diversity(n, index="shannon")/log(specnumber(n)))

S19_alpha<-S19_n%>%
  group_by(Tumor)%>%
  summarize(sobs=specnumber(n), shannon=diversity(n, index="shannon"), simpson=diversity(n, index="simpson"), invsimpson=diversity(n, index="invsimpson"),
            Pielou=diversity(n, index="shannon")/log(specnumber(n)))

S18_alpha<-S18_n%>%
  group_by(Tumor)%>%
  summarize(sobs=specnumber(n), shannon=diversity(n, index="shannon"), simpson=diversity(n, index="simpson"), invsimpson=diversity(n, index="invsimpson"),
            Pielou=diversity(n, index="shannon")/log(specnumber(n)))

S37_alpha<-S37_n%>%
  group_by(Tumor)%>%
  summarize(sobs=specnumber(n), shannon=diversity(n, index="shannon"), simpson=diversity(n, index="simpson"), invsimpson=diversity(n, index="invsimpson"),
            Pielou=diversity(n, index="shannon")/log(specnumber(n)))

M145_alpha<-M145_n%>%
  group_by(Tumor)%>%
  summarize(sobs=specnumber(n), shannon=diversity(n, index="shannon"), simpson=diversity(n, index="simpson"), invsimpson=diversity(n, index="invsimpson"),
            Pielou=diversity(n, index="shannon")/log(specnumber(n)))

M148_alpha<-M148_n%>%
  group_by(Tumor)%>%
  summarize(sobs=specnumber(n), shannon=diversity(n, index="shannon"), simpson=diversity(n, index="simpson"), invsimpson=diversity(n, index="invsimpson"),
            Pielou=diversity(n, index="shannon")/log(specnumber(n)))

M132_alpha<-M132_n%>%
  group_by(Tumor)%>%
  summarize(sobs=specnumber(n), shannon=diversity(n, index='shannon'), simpson=diversity(n, index="simpson"), invsimpson=diversity(n, index="invsimpson"),
            Pielou=diversity(n, index="shannon")/log(specnumber(n)))

AlphaDiversity<-rbind(M132_alpha, M148_alpha, M145_alpha, S19_alpha, S18_alpha, S37_alpha, S15_alpha, P2194_alpha)

#plot alpha diversity 
AlphaDiversity$Colony<-factor(AlphaDiversity$Colony, levels=c("M132", "M145", "M148", "S19", "S15", "S18", "S37" ,"P2194"))

plot_shan<-ggplot(AlphaDiversity, aes(x=Colony, y=Shannon_Diversity, fill=Age))+
  geom_boxplot(alpha=0.8)+geom_point()+
  scale_fill_manual(values=wes_palette("GrandBudapest2", n=4))+
  ylab("Shannon's H' ")+
  xlab("Colony")+
  theme_bw()+
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.4))

plot_rich<-ggplot(AlphaDiversity, aes(x=Colony, y=Richness, fill=Age))+
  geom_boxplot(alpha=0.8)+ 
  ylab("Richness")+
  xlab("Colony")+
  scale_fill_manual(values=wes_palette("GrandBudapest2", n=4))+
  theme_bw()+
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.4))

plot_simp<-ggplot(AlphaDiversity, aes(x=Colony, y=Simpson_Index, fill=Age))+
  geom_boxplot(alpha=0.8)+
  ylab("Simpson's Index")+
  xlab("Colony")+
  scale_fill_manual(values=wes_palette("GrandBudapest2", n=4))+
  theme_bw()+theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.4))

plot_invsimp<-ggplot(AlphaDiversity, aes(x=Colony, y=AlphaDiversity$Inverse.Simpson.Index, fill=Age))+
  geom_boxplot(alpha=0.8)+
  ylab("Inverse Simpson's Index")+
  xlab("Colony")+
  scale_fill_manual(values=wes_palette("GrandBudapest2", n=4))+
  theme_bw()+theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.4))

plot_evenness<-ggplot(AlphaDiversity, aes(x=Colony, y=Pielou_Evenness, fill=Age))+
  geom_boxplot(alpha=0.8)+geom_jitter()+
  ylab("Pielou's Evenness")+
  xlab("Colony")+scale_fill_manual(values=wes_palette("GrandBudapest2", n=4))+
  theme_bw()+theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.4))

plot_grid(plot_rich, plot_shan, plot_simp, plot_invsimp, plot_evenness)
#beta diversity calculations 
M132_beta<-M132_n%>%
  group_by(Tumor)
 
M132_beta2<-betadiver(M132_beta$n, "w")
betadiver(help)
P2194_beta2<-betadiver()
plobeta_dist<-vegdist(M132_beta$n, index="bray")
mds<-metaMDS(beta_dist)
mds_data<-as.data.frame(mds$points)
ggplot(mds_data, aes(x=MDS1, y=MDS2))+geom_point()

M145_beta<-M145_n%>%
  group_by(Tumor)
M145_beta2<-as.matrix.data.frame(M145_beta)
rownames(M145_beta3)<-M145_n$Tumor
M145_beta3<-M145_beta2[, 3]
M145_beta<-
M145_bray<-vegdist(M145_beta$n, method="bray")
beta_dist1<-vegdist(M145_beta$n, index="bray")
rownames(beta_dist1)<-M145_n$Tumor

mds2<-metaMDS(beta_dist1)
 mds_data2<-as.data.frame(mds2$points)
ggplot(mds_data2, aes(x=MDS1, y=MDS2, colour=tumor))+geom_point()+
mds_data2$tumor<-M145_beta$Tumor
M148_beta<-M148_n%>%
  group_by(Tumor)

beta_dist2<-vegdist(M148_beta$n, index="bray")
mds3<-metaMDS(beta_dist2)
mds_data3<-as.data.frame(mds3$points)
ggplot(mds_data3, aes(x=MDS1, y=MDS2))+geom_point()

S18_beta<-S18_n%>%
  group_by(Tumor)
beta_dist3<-vegdist(S18_beta$n, index="bray")
mds4<-metaMDS(beta_dist3)
mds_data4<-as.data.frame(mds4$points)
df<-melt(as.matrix(beta_dist2), varnames=c("row", "col"))
ggplot(mds_data4, aes(x=MDS1, y=MDS2))+geom_point()

S15_beta<-S15_n%>%
  group_by(Tumor)
beta_dist4<-vegdist(S15_beta$n, index="bray")
mds5<-metaMDS(beta_dist4)
mds_data5<-as.data.frame(mds5$points)
ggplot(mds_data5, aes(x=MDS1, y=MDS2))+geom_point()

S37_beta<-S37_n%>%
  group_by(Tumor)
beta_dist5<-vegdist(S37_beta$n, index="bray")
mds6<-metaMDS(beta_dist5)
mds_data6<-as.data.frame(mds6$points)
ggplot(mds_data6, aes(x=MDS1, y=MDS2))+geom_point()

S19_beta<-S19_n%>%
  group_by(Tumor)
beta_dist6<-vegdist(S19_beta$n, index="bray")
mds7<-metaMDS(beta_dist6)
mds_data7<-as.data.frame(mds7$points)

P2194_beta<-P2194_n%>%
  group_by(Tumor)
beta_dist7<-vegdist(P2194_beta$n, index="bray")
mds8<-metaMDS(beta_dist7)
mds_data8<-as.data.frame(mds8$points)

#rarefy data 
s2_alpha<-

#combine column with rest of data frame 
S18_merge<-merge(S18, S18_n, by=c("mutid", "Tumor"))
S15_merge<-merge(S15, S15_n, by=c("mutid", "Tumor"))
S19_merge<-merge(S19, S19_n, by=c("mutid", "Tumor"))
S37_merge<-merge(S37, S37_n, by=c("mutid", "Tumor"))
M145_merge<-merge(M145, M145_n, by=c("mutid", "Tumor"))
M148_merge<-merge(M148, M148_n, by=c("mutid", "Tumor"))
M132_merge<-merge(M132, M132_n, by=c("mutid", "Tumor"))
P2194_merge<-merge(P2194, P2194_n, by=c("mutid", "Tumor"))

#remove duplicates 
S18_dis<-S18_merge%>%
  distinct(mutid, .keep_all=TRUE)
S15_dis<-S15_merge%>%
  distinct(mutid, .keep_all=TRUE)
S19_dis<-S19_merge%>%
  distinct(mutid, .keep_all=TRUE)
S37_dis<-S37_merge%>%
  distinct(mutid, .keep_all = TRUE)
M145_dis<-M145_merge%>%
  distinct(mutid, .keep_all=TRUE)
M148_dis<-M148_merge%>%
  distinct(mutid, .keep_all=TRUE)
M132_dis<-M132_merge%>%
  distinct(mutid, .keep_all=TRUE)
P2194_dis<-P2194_merge%>%
  distinct(mutid, .keep_all=TRUE)
#pivot wider 
S18_wide<-S18_dis%>%
   (names_from = mutid, values_from = n)
S15_wide<-S15_dis%>%
  pivot_wider(names_from = mutid, values_from = n)
S19_wide<-S19_dis%>%
  pivot_wider(names_from = mutid, values_from = n)
S37_wide<-S37_dis%>%
  pivot_wider(names_from = mutid, values_from = n)
M148_wide<-M148_dis%>%
  pivot_wider(names_from = mutid, values_from = n)
M132_wide<-M132_dis%>%
  pivot_wider(names_from = mutid, values_from = n)
M145_wide<-M145_dis%>%
  pivot_wider(names_from = mutid, values_from = n)
P2194_wide<-P2194_dis%>%
  pivot_wider(names_from = mutid, values_from = n)

#P2194 too big, split it into smaller dataframes to perform the pivot and then recombine 
#split 
s<-20000
P2194_1<-P2194_dis[row.names(P2194_dis) %in% 1:s, ]
P2194_2<-P2194_dis[row.names(P2194_dis) %in% (s+1):nrow(P2194_dis), ]

P2194_wid1<-P2194_1%>%
  pivot_wider(names_from=mutid, values_from=n)

P2194_wid2<-P2194_2%>%
  pivot_wider(names_from = mutid, values_from=n)


P2194_wide<-bind_rows(P2194_wid1, P2194_wid2)

#make NAs 0s 

P2194_wide[is.na(P2194_wide)] <- 0 
S15_wide[is.na(S15_wide)]<- 0
S18_wide[is.na(S18_wide)]<- 0
S37_wide[is.na(S37_wide)]<- 0
S19_wide[is.na(S19_wide)]<- 0
M145_wide[is.na(M145_wide)]<- 0 
M148_wide[is.na(M148_wide)]<- 0 
M132_wide[is.na(M132_wide)]<-0


##make upset plot 
drop<-c("filelist[i]", "CHRO", "POS", "Tumor_AF", "Tumor_DP", "Normal_AF", "Normal_DP", 
                                       "Comparison", "Normal", "Age", "Colony", "FileType", "X", "X.1", "X.2")
M132_wide2<-M132_wide[,!(names(M132_wide)%in% drop)]
drop2<-c("Tumor")
M132_set<-M132_wide2[,!(names(M132_wide2)%in% drop2)]
set_vars<-colnames(M132_set)

upset(M132_wide2, order="freq")

M132_diff2<-M132_wide%>%
  row_to_names(row_number=2)

M132_diff2%>%
  rowise()%>%
  mutate(Sample_4=sum())

M132_diff3<-M132_diff2%>>%
  full_join(M132_diff)


M132_melted<-melt(M132_diff2, id.vars=4, measure.vars=patterns(c("132_2", "132_3", "132_4", "132_5"),
                        
 M132_2                                                                                          value.name=c("132_2", "132_3", "132_4", "132_5")))
##alpha diversity measurements 






##structural variant detection 

setwd("/Users/trini/Documents/PSU/Acropora_Manuscript/Structural")


#read in files in directory 
filelist=list.files(pattern="*_SVinfo.txt")

datalist=lapply(filelist, function(x)read.table(x, header=FALSE, sep="\t"))

for (i in 1:length(datalist)){datalist[[i]]<-cbind(datalist[[i]],filelist[i])}

sv<-do.call("rbind", datalist)

#merge with metadata

sv2<-merge(sv, Manta_Info, by=c("filelist[i]"))


#get rid of problem samples 


problem_samples<-c('manta185', 'manta187', 'manta188', 'manta189', 'manta190', 'manta191', 'manta192',
                   'manta194', 'manta195', 'manta196', 'manta198', 'manta199', 'manta201', 'manta203',
                   'manta206', 'manta210', 'manta213', 'manta214', 'manta215', 'manta216', 'manta223')

problem<-as.data.frame(problem_samples)
colnames(problem)[1]="Comparison"

sv3<-sv2[ !sv2$Comparison %in% problem_samples, ]

delsv<-sv3%>%
  filter(V3=="DEL")

insv<-sv3%>%
  filter(V3=="INS")
insv$V4<-as.numeric(insv$V4)
insv$insertionlength

View(delsv)
delsv$V4<-as.numeric(delsv$V4)
delsv$deletionlength<-abs(delsv$V4)
del<-ggplot(delsv, aes(x=Colony, y=deletionlength, colour=Age))+geom_boxplot()+scale_colour_viridis_d(option="mako")+
  ylab("Deletion Length")+
  xlab("Colony")+
  theme_bw()+
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.4))

ins<-ggplot(insv, aes(x=Colony, y=V4, colour=Age))+geom_boxplot()+ scale_colour_viridis_d(option="mako")+
  ylab("Richness")+
  xlab("Colony")+
  theme_bw()+
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.4))+
  ylab("Insertion Length")+
  xlab("Colony")+
  theme_bw()+
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.4))
plot_grid(ins, del)

#shannon means

a10<-AlphaDiversity%>%
  filter(Age=="10")
mean(a10$Pielou_Evenness)
            
#plotting medians 
Age<-c("6", "8", "10", ">100")
Shannon<-c("6.447563", "6.938804", "6.222099", "8.547895")
Simpsons<-c(".9981969", "0.9987337", "0.9977614", "0.999734")
InvSimpsons<-c("554.6002", "789.7897", "447.2088", "3758.938")
alphamedians$InvSimpsons<-as.numeric(alphamedians$InvSimpsons)
alphamedians$Simpsons<-as.numeric(alphamedians$Simpsons)
alphamedians$Richness<-as.numeric(alphamedians$Richness)
alphamedians$Shannon<-as.numeric(alphamedians$Shannon)
alphamedians$Evenness<-as.numeric(alphamedians$Evenness)
Richness<-c("748.5", "1308", "536", "6879")
Evenness<-c("0.977118", "0.9677522", "0.9867883", "0.9673692")
alphamedians<-data.frame(Age, Shannon, Simpsons, InvSimpsons, Richness, Evenness)
medians<-read.csv("./Distance_mosaic.csv", header=T)
       
alphamedians$Age<-factor(alphamedians$Age, levels=c("6", "8", "10", ">100"))
plot_evenness_med<-ggplot(alphamedians, aes(x=Age, y=Evenness, color=Age))+
  geom_point(size=3)+
  ylab("Pielou's Evenness")+
  xlab("Age")+scale_color_manual(values=wes_palette("GrandBudapest2", n=4))+
  theme_bw()+theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.4))

plot_shan_med<-ggplot(alphamedians, aes(x=Age, y=Shannon, color=Age))+
  geom_point(size=3)+
  ylab("Shannon's H'")+
  xlab("Age")+scale_color_manual(values=wes_palette("GrandBudapest2", n=4))+
  theme_bw()+theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.4))

plot_simp_med<-ggplot(alphamedians, aes(x=Age, y=Simpsons, color=Age))+
  geom_point(size=3)+
  ylab("Simpson's Index")+
  xlab("Age")+scale_color_manual(values=wes_palette("GrandBudapest2", n=4))+
  theme_bw()+theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.4))

plot_invsimp_med<-ggplot(alphamedians, aes(x=Age, y=InvSimpsons, color=Age))+
  geom_point(size=3)+
  ylab("Inverse Simpson's Index")+
  xlab("Age")+scale_color_manual(values=wes_palette("GrandBudapest2", n=4))+
  theme_bw()+theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.4))


plot_rich_med<-ggplot(alphamedians, aes(x=Age, y=Richness, color=Age))+
  geom_point(size=3)+
  ylab("Richness")+
  xlab("Age")+scale_color_manual(values=wes_palette("GrandBudapest2", n=4))+
  theme_bw()+theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.4))

plot_grid(plot_evenness, plot_evenness_med)


#Plotting D-values from ks 

#using correlattion matrix 
corrplot(Distance_mosaic, is.corr=FALSE, method="color", number.digits = 4, type="lower")
Distance_mosaic<-as.numeric(Distance_mosaic)
Distance_mosaic$X6<-as.numeric(Distance_mosaic$X6)
Distance_mosaic$X8<-as.numeric(Distance_mosaic$X8)
Distance_mosaic$X10<-as.numeric(Distance_mosaic$X10)
Distance_mosaic$X.100<-as.numeric(Distance_mosaic$X.100)
Distance_mosaic<-as.matrix(Distance_mosaic)

colnames(Distance_mosaic)[1]="6"
colnames(Distance_mosaic)[2]="8"
colnames(Distance_mosaic)[3]="10"
colnames(Distance_mosaic)[4]=">100"




##plotting Upset plot 

#starting with M145
M145_1<-M145%>%
  filter(Tumor=="145_1")

M145_1_n<-M145_1%>%
  count(mutid, Tumor)
M145_1_counts<-pivot_wider(M145_1_n, names_from = Tumor, values_from = n)


M145_2<-M145%>%
  filter(Tumor=="145_2")
M145_2_n<-M145_2%>%
  count(mutid, Tumor)
M145_2_counts<-pivot_wider(M145_2_n, names_from = Tumor, values_from = n)

M145_3<-M145%>%
  filter(Tumor=="145_3")
M145_3_n<-M145_3%>%
  count(mutid, Tumor)
M145_3_counts<-pivot_wider(M145_3_n, names_from = Tumor, values_from = n)

M145_4<-M145%>%
  filter(Tumor=="145_4")

M145_4_n<-M145_4%>%
  count(mutid, Tumor)
M145_4_counts<-pivot_wider(M145_4_n, names_from = Tumor, values_from = n)

M145_5<-M145%>%
  filter(Tumor=="145_5")
M145_5_n<-M145_5%>%
  count(mutid, Tumor)
M145_5_counts<-pivot_wider(M145_5_n, names_from = Tumor, values_from = n)


M145_upset<-full_join(M145_1_counts, M145_2_counts)
M145_upset2<-full_join(M145_upset, M145_3_counts)
M145_upset3<-full_join(M145_upset2, M145_4_counts)
M145_upset4<-full_join(M145_upset3, M145_5_counts)


#replace NAs with 0 
M145_upset4[is.na(M145_upset4)] <- 0
M145_upset4<-as.data.frame(M145_upset4)
#rename columsn to remove numeric constant 

colnames(M145_upset4)[2]="M145_1"
colnames(M145_upset4)[3]="M145_2"
colnames(M145_upset4)[4]="M145_3"
colnames(M145_upset4)[5]="M145_4"
colnames(M145_upset4)[6]="M145_5"
#convert all numbers greater than 1 to 1 
M145_upset4$M145_1<-ifelse(M145_upset4$M145_1,1,M145_upset4$M145_1)
M145_upset4$M145_2<-ifelse(M145_upset4$M145_2,1,M145_upset4$M145_2)
M145_upset4$M145_3<-ifelse(M145_upset4$M145_3,1,M145_upset4$M145_3)
M145_upset4$M145_4<-ifelse(M145_upset4$M145_4,1,M145_upset4$M145_4)
M145_upset4$M145_5<-ifelse(M145_upset4$M145_5,1,M145_upset4$M145_5)
#load upset package
set_vars<-c("M145_1", "M145_2", "M145_3M", "M145_4", "M145_5")
M145_list<-upset(M145_upset4, keep.order=T, sets=c("M145_1", "M145_2", "M145_3", "M145_4", "M145_5"))

#extract sets 
Data145<-as.data.frame(M145_list$New_data)
Data145<-Data145%>%
  mutate(ID=row_number())

#mutations where the row sum is 1 
private_M145<-rowSums(Data145[,2:6])==1
private_M145<-as.data.frame(private_M145)
private_M145<-private_M145%>%
  mutate(ID=row_number())

M145_Identi<-merge(Data145, private_M145, by=c("ID"))
M145_private<-M145_Identi%>%
  filter(private_M145=="TRUE")
View(M145_private)

#add in colony and mutation type 

M145_private<-M145_private%>%
  add_column(Intersection="Private")

#merge with original dataframe
private_mutations<-merge(whole_filter, M145_private, by=c("mutid"))
private_mutations<-private_mutations%>%
  select(mutid, Tumor_AF, Tumor_DP, Colony, Tumor, Intersection)

ggplot(private_mutations, aes(x=Tumor_AF))+geom_histogram()+theme_bw()+xlim(0,1.0)

#extract shared with two 
#mutations where the row sum is 2
private_M1452<-rowSums(Data145[,2:6])==2
private_M1452<-as.data.frame(private_M1452)
private_M1452<-private_M1452%>%
  mutate(ID=row_number())
M145_Identi2<-merge(Data145, private_M1452, by=c("ID"))
M145_private2<-M145_Identi2%>%
  filter(private_M1452=="TRUE")

#add in colony and mutation type 


M145_private2<-M145_private2%>%
  add_column(Intersection="Two")
View(M145_private2)

twoshare<-merge(whole_filter, M145_private2, by=c("mutid"))
twoshare<-twoshare%>%
  select(mutid, Tumor_AF, Tumor_DP, Colony, Tumor, Intersection)+theme_bw()
ggplot(twoshare, aes(x=Tumor_AF))+geom_histogram()+theme_bw()+xlim(0,1.0)
M145_Upset_AF<-rbind(private_mutations, twoshare)



#extract shared with three
#mutations where the row sum is 3
private_M1453<-rowSums(Data145[,2:6])==3
private_M1453<-as.data.frame(private_M1453)
private_M1453<-private_M1453%>%
  mutate(ID=row_number())

View(M145_private3)
M145_Identi3<-merge(Data145, private_M1453, by=c("ID"))
M145_private3<-M145_Identi3%>%
  filter(private_M1453=="TRUE")
M145_private3<-M145_private3%>%
  add_column(Intersection="Three")
threeshare<-threeshare%>%
  select(mutid, Tumor_AF, Tumor_DP, Colony, Tumor, Intersection)
threeshare<-merge(whole_filter, M145_private3, by=c("mutid"))
ggplot(threeshare, aes(x=Tumor_AF))+geom_histogram()+theme_bw()+xlim(0,1.0)
M145_Upset_AF<-rbind(M145_Upset_AF, fourshare)

#extract shared with four 
#mutations where the row sum is 4
private_M1454<-rowSums(Data145[,2:6])==4
private_M1454<-as.data.frame(private_M1454)
private_M1454<-private_M1454%>%
  mutate(ID=row_number())
M145_Identi4<-merge(Data145, private_M1454, by=c("ID"))
M145_private4<-M145_Identi4%>%
  filter(private_M1454=="TRUE")
M145_private4<-M145_private4%>%
  add_column(Intersection="Four")


fourshare<-merge(whole_filter, M145_private4, by=c("mutid"))
fourshare<-fourshare%>%
  select(mutid, Tumor_AF, Tumor_DP, Colony, Tumor, Intersection)

M145_Upset_AF1<-M145_Upset_AF%>%
  filter(Intersection=="Private")
M145_Upset_AF4<-M145_Upset_AF_filter%>%
  filter(Intersection=="Four")

ggplot(M145_Upset_AF1, aes(x=Tumor_AF))+geom_histogram(bins=100
                                                      )+xlim(0,1.0)+theme_bw()
M145_Upset_AF_filter<-M145_Upset_AF%>%
  filter(Tumor_DP>=60)

#M148
M148_1<-M148%>%
  filter(Tumor=="148_1")

M148_1_n<-M148_1%>%
  count(mutid, Tumor)
M148_1_counts<-pivot_wider(M148_1_n, names_from = Tumor, values_from = n)


M148_2<-M148%>%
  filter(Tumor=="148_2")
M148_2_n<-M148_2%>%
  count(mutid, Tumor)
M148_2_counts<-pivot_wider(M148_2_n, names_from = Tumor, values_from = n)

M148_3<-M148%>%
  filter(Tumor=="148_3")
M148_3_n<-M148_3%>%
  count(mutid, Tumor)
M148_3_counts<-pivot_wider(M148_3_n, names_from = Tumor, values_from = n)

M148_4<-M148%>%
  filter(Tumor=="148_4")

M148_4_n<-M148_4%>%
  count(mutid, Tumor)
M148_4_counts<-pivot_wider(M148_4_n, names_from = Tumor, values_from = n)

M148_5<-M148%>%
  filter(Tumor=="148_5")
M148_5_n<-M148_5%>%
  count(mutid, Tumor)
M148_5_counts<-pivot_wider(M148_5_n, names_from = Tumor, values_from = n)


M148_upset<-full_join(M148_1_counts, M148_2_counts)
M148_upset2<-full_join(M148_upset, M148_3_counts)
M148_upset3<-full_join(M148_upset2, M148_4_counts)
M148_upset4<-full_join(M148_upset3, M148_5_counts)


#replace NAs with 0 
M148_upset4[is.na(M148_upset4)] <- 0
M148_upset4<-as.data.frame(M148_upset4)
#rename columsn to remove numeric constant 

colnames(M148_upset4)[2]="M148_1"
colnames(M148_upset4)[3]="M148_2"
colnames(M148_upset4)[4]="M148_3"
colnames(M148_upset4)[5]="M148_4"
colnames(M148_upset4)[6]="M148_5"

#convert all numbers greater than 1 to 1 
M148_upset4$M148_1<-ifelse(M148_upset4$M148_1,1,M148_upset4$M148_1)
M148_upset4$M148_2<-ifelse(M148_upset4$M148_2,1,M148_upset4$M148_2)
M148_upset4$M148_3<-ifelse(M148_upset4$M148_3,1,M148_upset4$M148_3)
M148_upset4$M148_4<-ifelse(M148_upset4$M148_4,1,M148_upset4$M148_4)
M148_upset4$M148_5<-ifelse(M148_upset4$M148_5,1,M148_upset4$M148_5)

#load upset package
set_vars<-c("M145_1", "M145_2", "M145_3", "M145_4", "M145_5")
M148_List<-upset(M148_upset4, keep.order=T, sets=c("M148_1", "M148_2", "M148_3", "M148_4", "M148_5"))

#extract sets 
Data148<-as.data.frame(M148_List$New_data)
Data148<-Data148%>%
  mutate(ID=row_number())

#mutations where the row sum is 1 
private_M148<-rowSums(Data148[,2:6])==1
private_M148<-as.data.frame(private_M148)
private_M148<-private_M148%>%
  mutate(ID=row_number())

M148_Identi<-merge(Data148, private_M148, by=c("ID"))
M148_private<-M148_Identi%>%
  filter(private_M148=="TRUE")

#add in colony and mutation type 

M148_private<-M148_private%>%
  add_column(Intersection="Private")


#merge with original dataframe
private_mutations148<-merge(whole_filter, M148_private, by=c("mutid"))

ggplot(private_mutations148, aes(x=Tumor_AF))+geom_histogram()+xlim(0,1.0)+theme_bw()

#extract shared with two 
#mutations where the row sum is 2
private_M1482<-rowSums(Data148[,2:6])==2

private_M1482<-as.data.frame(private_M1482)
private_M1482<-private_M1482%>%
  mutate(ID=row_number())
M148_Identi2<-merge(Data148, private_M1482, by=c("ID"))
M148_private2<-M148_Identi2%>%
  filter(private_M1482=="TRUE")
twoshare6<-full_join(twoshare, twoshare148)
#add in colony and mutation type 


M148_private<-M148_private%>%
  add_column(Intersection="Private")

twoshare148<-merge(whole_filter, M148_private2, by=c("mutid"))

ggplot(twoshare148, aes(x=Tumor_AF))+geom_histogram()+xlim(0,1.0)+theme_bw()
#extract shared with three
#mutations where the row sum is 3
private_M1483<-rowSums(Data148[,2:6])==3
private_M1483<-as.data.frame(private_M1483)
private_M1483<-private_M1483%>%
  mutate(ID=row_number())
M148_Identi3<-merge(Data148, private_M1483, by=c("ID"))
M148_private3<-M148_Identi3%>%
  filter(private_M1483=="TRUE")

#add in colony and mutation type 


M148_private<-M148_private%>%
  add_column(Intersection="Private")


threeshare148<-merge(whole_filter, M148_private3, by=c("mutid"))

#extract shared with four 
#mutations where the row sum is 4
private_M1484<-rowSums(Data148[,2:6])==4
private_M1484<-as.data.frame(private_M1484)
private_M1484<-private_M1484%>%
  mutate(ID=row_number())
M148_Identi4<-merge(Data148, private_M1484, by=c("ID"))
M148_private4<-M148_Identi4%>%
  filter(private_M1484=="TRUE")



fourshare148<-merge(whole_filter, M148_private4, by=c("mutid"))

threeshare148<-threeshare148%>%
  select(c("mutid", "Tumor_AF"))
fourshare148<-fourshare148%>%
  select(c("mutid", "Tumor_AF"))

multi<-rbind(fourshare148, threeshare148)
ggplot(multi, aes(x=Tumor_AF))+geom_histogram()+xlim(0,1.0)+theme_bw()
#merge all sets together 



#M132
M132_1<-M132%>%
  filter(Tumor=="132_1")

M132_1_n<-M132_1%>%
  count(mutid, Tumor)
M132_1_counts<-pivot_wider(M132_1_n, names_from = Tumor, values_from = n)


M132_2<-M132%>%
  filter(Tumor=="132_2")
M132_2_n<-M132_2%>%
  count(mutid, Tumor)
M132_2_counts<-pivot_wider(M132_2_n, names_from = Tumor, values_from = n)

M132_3<-M132%>%
  filter(Tumor=="132_3")
M132_3_n<-M132_3%>%
  count(mutid, Tumor)
M132_3_counts<-pivot_wider(M132_3_n, names_from = Tumor, values_from = n)

M132_4<-M132%>%
  filter(Tumor=="132_4")
M132_4_n<-M132_4%>%
  count(mutid, Tumor)
M132_4_counts<-pivot_wider(M132_4_n, names_from = Tumor, values_from = n)

M132_5<-M132%>%
  filter(Tumor=="132_5")
M132_5_n<-M132_5%>%
  count(mutid, Tumor)
M132_5_counts<-pivot_wider(M132_5_n, names_from = Tumor, values_from = n)


M132_upset<-full_join(M132_1_counts, M132_2_counts)
M132_upset2<-full_join(M132_upset, M132_3_counts)
M132_upset3<-full_join(M132_upset2, M132_4_counts)
M132_upset4<-full_join(M132_upset3, M132_5_counts)


#replace NAs with 0 
M132_upset4[is.na(M132_upset4)] <- 0
M132_upset4<-as.data.frame(M132_upset4)
#rename columsn to remove numeric constant 

colnames(M132_upset4)[2]="M132_1"
colnames(M132_upset4)[3]="M132_2"
colnames(M132_upset4)[4]="M132_3"
colnames(M132_upset4)[5]="M132_4"


#convert all numbers greater than 1 to 1 
M132_upset4$M132_1<-ifelse(M132_upset4$M132_1,1,M132_upset4$M132_1)
M132_upset4$M132_2<-ifelse(M132_upset4$M132_2,1,M132_upset4$M132_2)
M132_upset4$M132_3<-ifelse(M132_upset4$M132_3,1,M132_upset4$M132_3)
M132_upset4$M132_4<-ifelse(M132_upset4$M132_4,1,M132_upset4$M132_4)

#load upset package

M132_list<-upset(M132_upset4, keep.order=T, sets=c("M132_1", "M132_2", "M132_3", "M132_4"))

#extract sets 
Data132<-as.data.frame(M132_list$New_data)
Data132<-Data132%>%
  mutate(ID=row_number())

#mutations where the row sum is 1 
private_M132<-rowSums(Data132[,2:5])==1
private_M132<-as.data.frame(private_M132)
private_M132<-private_M132%>%
  mutate(ID=row_number())

M132_Identi<-merge(Data132, private_M132, by=c("ID"))
M132_private<-M132_Identi%>%
  filter(private_M132=="TRUE")

#add in colony and mutation type 

M132_private<-M132_private%>%
  add_column(Intersection="Private")


#merge with original dataframe
private_mutations132<-merge(whole_filter, M132_private, by=c("mutid"))

ggplot(private_mutations132, aes(x=Tumor_AF))+geom_histogram()+xlim(0,1.0)+theme_bw()

#mutations where the row sum is 2
private_M1322<-rowSums(Data132[,2:5])==2

private_M1322<-as.data.frame(private_M1322)
private_M1322<-private_M1322%>%
  mutate(ID=row_number())
M132_Identi2<-merge(Data132, private_M1322, by=c("ID"))
M132_private2<-M132_Identi2%>%
  filter(private_M1322=="TRUE")

#add in colony and mutation type 


M132_private2<-M132_private2%>%
  add_column(Intersection="Two")

twoshare132<-merge(whole_filter, M132_private2, by=c("mutid"))

ggplot(twoshare132, aes(x=Tumor_AF))+geom_histogram()+xlim(0,1.0)+theme_bw()
#extract shared with three
private_M1323<-rowSums(Data132[,2:5])==3

private_M1323<-as.data.frame(private_M1323)
private_M1323<-private_M1323%>%
  mutate(ID=row_number())
M132_Identi3<-merge(Data132, private_M1323, by=c("ID"))
M132_private3<-M132_Identi3%>%
  filter(private_M1323=="TRUE")

#add in colony and mutation type 


M132_private3<-M132_private3%>%
  add_column(Intersection="Three")

twoshare133<-merge(whole_filter, M132_private3, by=c("mutid"))

ggplot(twoshare133, aes(x=Tumor_AF))+geom_histogram()+xlim(0,1.0)+theme_bw()
#S15
S15_1<-S15%>%
  filter(Tumor=="15_1")

S15_1_n<-S15_1%>%
  count(mutid, Tumor)
S15_1_counts<-pivot_wider(S15_1_n, names_from = Tumor, values_from = n)


S15_2<-S15%>%
  filter(Tumor=="15_2")
S15_2_n<-S15_2%>%
  count(mutid, Tumor)
S15_2_counts<-pivot_wider(S15_2_n, names_from = Tumor, values_from = n)

S15_3<-S15%>%
  filter(Tumor=="15_3")
S15_3_n<-S15_3%>%
  count(mutid, Tumor)
S15_3_counts<-pivot_wider(S15_3_n, names_from = Tumor, values_from = n)

S15_4<-S15%>%
  filter(Tumor=="15_4")

S15_4_n<-S15_4%>%
  count(mutid, Tumor)
S15_4_counts<-pivot_wider(S15_4_n, names_from = Tumor, values_from = n)

S15_5<-S15%>%
  filter(Tumor=="15_5")
S15_5_n<-S15_5%>%
  count(mutid, Tumor)
S15_5_counts<-pivot_wider(S15_5_n, names_from = Tumor, values_from = n)


S15_upset<-full_join(S15_1_counts, S15_2_counts)
S15_upset2<-full_join(S15_upset, S15_3_counts)
S15_upset3<-full_join(S15_upset2, S15_4_counts)
S15_upset4<-full_join(S15_upset3, S15_5_counts)


#replace NAs with 0 
S15_upset4[is.na(S15_upset4)] <- 0
S15_upset4<-as.data.frame(S15_upset4)
#rename columsn to remove numeric constant 

colnames(S15_upset4)[2]="S15_1"

colnames(S15_upset4)[3]="S15_3"
colnames(S15_upset4)[4]="S15_4"
colnames(S15_upset4)[5]="S15_5"

#convert all numbers greater than 1 to 1 
S15_upset4$S15_1<-ifelse(S15_upset4$S15_1,1,S15_upset4$S15_1)
M148_upset4$M148_2<-ifelse(M148_upset4$M148_2,1,M148_upset4$M148_2)
S15_upset4$S15_3<-ifelse(S15_upset4$S15_3,1,S15_upset4$S15_3)
S15_upset4$S15_4<-ifelse(S15_upset4$S15_4,1,S15_upset4$S15_4)
S15_upset4$S15_5<-ifelse(S15_upset4$S15_5,1,S15_upset4$S15_5)
#load upset package
set_vars<-c("M145_1", "M145_2", "M145_3", "M145_4", "M145_5")
upset(S15_upset4, keep.order=T, sets=c("S15_1",  "S15_3", "S15_4", "S15_5"))

#S18
S18_1<-S18%>%
  filter(Tumor=="18_1")

S18_1_n<-S18_1%>%
  count(mutid, Tumor)
S18_1_counts<-pivot_wider(S18_1_n, names_from = Tumor, values_from = n)


S18_2<-S18%>%
  filter(Tumor=="18_2")
S18_2_n<-S18_2%>%
  count(mutid, Tumor)
S18_2_counts<-pivot_wider(S18_2_n, names_from = Tumor, values_from = n)

S18_3<-S18%>%
  filter(Tumor=="18_3")
S18_3_n<-S18_3%>%
  count(mutid, Tumor)
S18_3_counts<-pivot_wider(S18_3_n, names_from = Tumor, values_from = n)

S18_4<-S18%>%
  filter(Tumor=="18_4")
S18_4_n<-S18_4%>%
  count(mutid, Tumor)
S18_4_counts<-pivot_wider(S18_4_n, names_from = Tumor, values_from = n)

S18_5<-S18%>%
  filter(Tumor=="18_5")
S18_5_n<-S18_5%>%
  count(mutid, Tumor)
S18_5_counts<-pivot_wider(S18_5_n, names_from = Tumor, values_from = n)


S18_upset<-full_join(S18_1_counts, S18_2_counts)
S18_upset2<-full_join(S18_upset, S18_3_counts)
S18_upset3<-full_join(S18_upset2, S18_4_counts)
S18_upset4<-full_join(S18_upset3, S18_5_counts)


#replace NAs with 0 
S18_upset4[is.na(S18_upset4)] <- 0
S18_upset4<-as.data.frame(S18_upset4)
#rename columsn to remove numeric constant 

colnames(S18_upset4)[2]="S18_1"
colnames(S18_upset4)[3]="S18_2"
colnames(S18_upset4)[4]="S18_3"
colnames(S18_upset4)[4]="S18_4"
colnames(S18_upset4)[5]="S18_5"

#convert all numbers greater than 1 to 1 
S18_upset4$S18_1<-ifelse(S18_upset4$S18_1,1,S18_upset4$S18_1)
S18_upset4$S18_2<-ifelse(S18_upset4$S18_2,1,S18_upset4$S18_2)
S18_upset4$S18_3<-ifelse(S18_upset4$S18_3,1,S18_upset4$S18_3)
S18_upset4$S18_4<-ifelse(S18_upset4$S18_4,1,S18_upset4$S18_4)
S18_upset4$S18_5<-ifelse(S18_upset4$S18_5,1,S18_upset4$S18_5)
#load upset package
set_vars<-c("S18_1", "S18_2", "S18_4", "S18_5")
S18_List<-upset(S18_upset4, keep.order=T, sets=c("S18_1", "S18_2", "S18_4", "S18_5"))

#extract sets 
Data18<-as.data.frame(S18_List$New_data)
Data18<-Data18%>%
  mutate(ID=row_number())

#mutations where the row sum is 1 
private_S18<-rowSums(Data18[,2:5])==1
private_S18<-as.data.frame(private_S18)
private_S18<-private_S18%>%
  mutate(ID=row_number())

S18_Identi<-merge(Data18, private_S18, by=c("ID"))
S18_private<-S18_Identi%>%
  filter(private_S18=="TRUE")

#add in colony and mutation type 

S18_private<-S18_private%>%
  add_column(Intersection="Private")


#merge with original dataframe
private_mutations18<-merge(whole_filter, S18_private, by=c("mutid"))

ggplot(private_mutations18, aes(x=Tumor_AF))+geom_histogram()+xlim(0,1.0)+theme_bw()

#mutations where the row sum is 2
private_S182<-rowSums(Data18[,2:5])==2

private_S182<-as.data.frame(private_S182)
private_S182<-private_S182%>%
  mutate(ID=row_number())
S18_Identi2<-merge(Data18, private_S182, by=c("ID"))
S18_private2<-S18_Identi2%>%
  filter(private_S182=="TRUE")

#add in colony and mutation type 


S18_private2<-S18_private2%>%
  add_column(Intersection="Two")

twoshare18<-merge(whole_filter, S18_private2, by=c("mutid"))
twoshare18<-twoshare18%>%
  filter(Tumor_AF<=0.5)

ggplot(twoshare18, aes(x=Tumor_AF))+geom_histogram()+xlim(0,1.0)+theme_bw()
#extract shared with three
private_S183<-rowSums(Data18[,2:5])==3

private_S183<-as.data.frame(private_S183)
private_S183<-private_S183%>%
  mutate(ID=row_number())
S18_Identi3<-merge(Data18, private_S183, by=c("ID"))
S18_private3<-S18_Identi3%>%
  filter(private_S183=="TRUE")

#add in colony and mutation type 


S18_private3<-S18_private3%>%
  add_column(Intersection="Three")

threeshare18<-merge(whole_filter, S18_private3, by=c("mutid"))

ggplot(threeshare18, aes(x=Tumor_AF))+geom_histogram()+xlim(0,1.0)+theme_bw()
#S37

S37_1<-S37%>%
filter(Tumor=="37_1")

S37_1_n<-S37_1%>%
  count(mutid, Tumor)

S37_1_counts<-pivot_wider(S37_1_n, names_from = Tumor, values_from = n)


S37_2<-S37%>%
  filter(Tumor=="37_2")
S37_2_n<-S37_2%>%
  count(mutid, Tumor)
S37_2_counts<-pivot_wider(S37_2_n, names_from = Tumor, values_from = n)

S37_3<-S37%>%
  filter(Tumor=="37_3")
S37_3_n<-S37_3%>%
  count(mutid, Tumor)
S37_3_counts<-pivot_wider(S37_3_n, names_from = Tumor, values_from = n)

S37_4<-S37%>%
  filter(Tumor=="37_4")
S37_4_n<-S37_4%>%
  count(mutid, Tumor)
S37_4_counts<-pivot_wider(S37_4_n, names_from = Tumor, values_from = n)




S37_upset<-full_join(S37_1_counts, S37_2_counts)
S37_upset2<-full_join(S37_upset, S37_3_counts)
S37_upset3<-full_join(S37_upset2, S37_4_counts)
S37_upset4<-S37_upset3


#replace NAs with 0 
S37_upset4[is.na(S37_upset4)] <- 0
S37_upset4<-as.data.frame(S37_upset4)
#rename columsn to remove numeric constant 

colnames(S37_upset4)[2]="S37_1"
colnames(S37_upset4)[3]="S37_2"
colnames(S37_upset4)[4]="S37_3"
colnames(S37_upset4)[5]="S37_4"


#convert all numbers greater than 1 to 1 
S37_upset4$S37_1<-ifelse(S37_upset4$S37_1,1,S37_upset4$S37_1)
S37_upset4$S37_2<-ifelse(S37_upset4$S37_2,1,S37_upset4$S37_2)
S37_upset4$S37_3<-ifelse(S37_upset4$S37_3,1,S37_upset4$S37_3)
S37_upset4$S37_4<-ifelse(S37_upset4$S37_4,1,S37_upset4$S37_4)

#load upset package
set_vars<-c("S37_1", "S37_2", "S37_3", "S_5")
S37_List<-upset(S37_upset4, keep.order=T, sets=c("S37_1", "S37_2", "S37_3", "S37_4"))
#extract sets 
Data37<-as.data.frame(S37_List$New_data)
Data37<-Data37%>%
  mutate(ID=row_number())

#mutations where the row sum is 1 
private_S37<-rowSums(Data37[,2:5])==1
private_S37<-as.data.frame(private_S37)
private_S37<-private_S37%>%
  mutate(ID=row_number())

S37_Identi<-merge(Data37, private_S37, by=c("ID"))
S37_private<-S37_Identi%>%
  filter(private_S37=="TRUE")

#add in colony and mutation type 

S37_private<-S37_private%>%
  add_column(Intersection="Private")


#merge with original dataframe
private_mutations37<-merge(whole_filter, S37_private, by=c("mutid"))

ggplot(private_mutations37, aes(x=Tumor_AF))+geom_histogram()+xlim(0,1.0)+theme_bw()

#mutations where the row sum is 2
private_S372<-rowSums(Data37[,2:5])==2

private_S372<-as.data.frame(private_S372)
private_S372<-private_S372%>%
  mutate(ID=row_number())
S37_Identi2<-merge(Data37, private_S372, by=c("ID"))
S37_private2<-S37_Identi2%>%
  filter(private_S372=="TRUE")

#add in colony and mutation type 


S37_private2<-S37_private2%>%
  add_column(Intersection="Two")

twoshare37<-merge(whole_filter, S37_private2, by=c("mutid"))

ggplot(twoshare37, aes(x=Tumor_AF))+geom_histogram()+xlim(0,1.0)+theme_bw()
#extract shared with three
private_S373<-rowSums(Data37[,2:5])==3

private_S373<-as.data.frame(private_S373)
private_S373<-private_S373%>%
  mutate(ID=row_number())
S37_Identi3<-merge(Data37, private_S373, by=c("ID"))
S37_private3<-S37_Identi3%>%
  filter(private_S373=="TRUE")

#add in colony and mutation type 


S37_private3<-S37_private3%>%
  add_column(Intersection="Three")

threeshare37<-merge(whole_filter, S37_private3, by=c("mutid"))

ggplot(threeshare37, aes(x=Tumor_AF))+geom_histogram()+xlim(0,1.0)+theme_bw()
#S19
S37_1<-S37%>%
  filter(Tumor=="37_1")

S37_1_n<-S37_1%>%
  count(mutid, Tumor)

S37_1_counts<-pivot_wider(S37_1_n, names_from = Tumor, values_from = n)


S37_2<-S37%>%
  filter(Tumor=="37_2")
S37_2_n<-S37_2%>%
  count(mutid, Tumor)
S37_2_counts<-pivot_wider(S37_2_n, names_from = Tumor, values_from = n)

S37_3<-S37%>%
  filter(Tumor=="37_3")
S37_3_n<-S37_3%>%
  count(mutid, Tumor)
S37_3_counts<-pivot_wider(S37_3_n, names_from = Tumor, values_from = n)

S37_4<-S37%>%
  filter(Tumor=="37_4")
S37_4_n<-S37_4%>%
  count(mutid, Tumor)
S37_4_counts<-pivot_wider(S37_4_n, names_from = Tumor, values_from = n)




S37_upset<-full_join(S37_1_counts, S37_2_counts)
S37_upset2<-full_join(S37_upset, S37_3_counts)
S37_upset3<-full_join(S37_upset2, S37_4_counts)
S37_upset4<-S37_upset3


#replace NAs with 0 
S37_upset4[is.na(S37_upset4)] <- 0
S37_upset4<-as.data.frame(S37_upset4)
#rename columsn to remove numeric constant 

colnames(S37_upset4)[2]="S37_1"
colnames(S37_upset4)[3]="S37_2"
colnames(S37_upset4)[4]="S37_3"
colnames(S37_upset4)[5]="S37_4"


#convert all numbers greater than 1 to 1 
S37_upset4$S37_1<-ifelse(S37_upset4$S37_1,1,S37_upset4$S37_1)
S37_upset4$S37_2<-ifelse(S37_upset4$S37_2,1,S37_upset4$S37_2)
S37_upset4$S37_3<-ifelse(S37_upset4$S37_3,1,S37_upset4$S37_3)
S37_upset4$S37_4<-ifelse(S37_upset4$S37_4,1,S37_upset4$S37_4)

#load upset package
set_vars<-c("S37_1", "S37_2", "S37_4", "S18_5")
upset(S37_upset4, keep.order=T, sets=c("S37_1", "S37_2", "S37_3", "S37_4"))

#P2194

S37_1<-S37%>%
  filter(Tumor=="37_1")

S37_1_n<-S37_1%>%
  count(mutid, Tumor)

S37_1_counts<-pivot_wider(S37_1_n, names_from = Tumor, values_from = n)


S37_2<-S37%>%
  filter(Tumor=="37_2")
S37_2_n<-S37_2%>%
  count(mutid, Tumor)
S37_2_counts<-pivot_wider(S37_2_n, names_from = Tumor, values_from = n)

S37_3<-S37%>%
  filter(Tumor=="37_3")
S37_3_n<-S37_3%>%
  count(mutid, Tumor)
S37_3_counts<-pivot_wider(S37_3_n, names_from = Tumor, values_from = n)

S37_4<-S37%>%
  filter(Tumor=="37_4")
S37_4_n<-S37_4%>%
  count(mutid, Tumor)
S37_4_counts<-pivot_wider(S37_4_n, names_from = Tumor, values_from = n)




S37_upset<-full_join(S37_1_counts, S37_2_counts)
S37_upset2<-full_join(S37_upset, S37_3_counts)
S37_upset3<-full_join(S37_upset2, S37_4_counts)
S37_upset4<-S37_upset3


#replace NAs with 0 
S37_upset4[is.na(S37_upset4)] <- 0
S37_upset4<-as.data.frame(S37_upset4)
#rename columsn to remove numeric constant 

colnames(S37_upset4)[2]="S37_1"
colnames(S37_upset4)[3]="S37_2"
colnames(S37_upset4)[4]="S37_3"
colnames(S37_upset4)[5]="S37_4"


#convert all numbers greater than 1 to 1 
S37_upset4$S37_1<-ifelse(S37_upset4$S37_1,1,S37_upset4$S37_1)
S37_upset4$S37_2<-ifelse(S37_upset4$S37_2,1,S37_upset4$S37_2)
S37_upset4$S37_3<-ifelse(S37_upset4$S37_3,1,S37_upset4$S37_3)
S37_upset4$S37_4<-ifelse(S37_upset4$S37_4,1,S37_upset4$S37_4)

#load upset package
set_vars<-c("S37_1", "S37_2", "S37_4", "S18_5")
S37_list<-upset(S37_upset4, keep.order=T, sets=c("S37_1", "S37_2", "S37_3", "S37_4"))

#S37

S37_1<-S37%>%
  filter(Tumor=="37_1")

S37_1_n<-S37_1%>%
  count(mutid, Tumor)

S37_1_counts<-pivot_wider(S37_1_n, names_from = Tumor, values_from = n)


S37_2<-S37%>%
  filter(Tumor=="37_2")
S37_2_n<-S37_2%>%
  count(mutid, Tumor)
S37_2_counts<-pivot_wider(S37_2_n, names_from = Tumor, values_from = n)

S37_3<-S37%>%
  filter(Tumor=="37_3")
S37_3_n<-S37_3%>%
  count(mutid, Tumor)
S37_3_counts<-pivot_wider(S37_3_n, names_from = Tumor, values_from = n)

S37_4<-S37%>%
  filter(Tumor=="37_4")
S37_4_n<-S37_4%>%
  count(mutid, Tumor)
S37_4_counts<-pivot_wider(S37_4_n, names_from = Tumor, values_from = n)




S37_upset<-full_join(S37_1_counts, S37_2_counts)
S37_upset2<-full_join(S37_upset, S37_3_counts)
S37_upset3<-full_join(S37_upset2, S37_4_counts)
S37_upset4<-S37_upset3


#replace NAs with 0 
S37_upset4[is.na(S37_upset4)] <- 0
S37_upset4<-as.data.frame(S37_upset4)
#rename columsn to remove numeric constant 

colnames(S37_upset4)[2]="S37_1"
colnames(S37_upset4)[3]="S37_2"
colnames(S37_upset4)[4]="S37_3"
colnames(S37_upset4)[5]="S37_4"


#convert all numbers greater than 1 to 1 
S37_upset4$S37_1<-ifelse(S37_upset4$S37_1,1,S37_upset4$S37_1)
S37_upset4$S37_2<-ifelse(S37_upset4$S37_2,1,S37_upset4$S37_2)
S37_upset4$S37_3<-ifelse(S37_upset4$S37_3,1,S37_upset4$S37_3)
S37_upset4$S37_4<-ifelse(S37_upset4$S37_4,1,S37_upset4$S37_4)

#load upset package
set_vars<-c("S37_1", "S37_2", "S37_4", "S18_5")
upset(S37_upset4, keep.order=T, sets=c("S37_1", "S37_2", "S37_3", "S37_4"))

#S19
S19_1<-S19%>%
  filter(Tumor=="16621")

S19_1_n<-S19_1%>%
  count(mutid, Tumor)

S19_1_counts<-pivot_wider(S19_1_n, names_from = Tumor, values_from = n)


S19_2<-S19%>%
  filter(Tumor=="16622")
S19_2_n<-S19_2%>%
  count(mutid, Tumor)
S19_2_counts<-pivot_wider(S19_2_n, names_from = Tumor, values_from = n)

S19_3<-S19%>%
  filter(Tumor=="16624")
S19_3_n<-S19_3%>%
  count(mutid, Tumor)
S19_3_counts<-pivot_wider(S19_3_n, names_from = Tumor, values_from = n)

S19_4<-S19%>%
  filter(Tumor=="16625")
S19_4_n<-S19_4%>%
  count(mutid, Tumor)
S19_4_counts<-pivot_wider(S19_4_n, names_from = Tumor, values_from = n)

S19_5<-S19%>%
  filter(Tumor=="16627")
S19_5_n<-S19_5%>%
  count(mutid, Tumor)
S19_5_counts<-pivot_wider(S19_5_n, names_from=Tumor, values_from=n)

S19_6<-S19%>%
  filter(Tumor=='16628')
S19_6_n<-S19_6%>%
  count(mutid, Tumor)
S19_6_counts<-pivot_wider(S19_6_n, names_from=Tumor, values_from=n )

S19_7<-S19%>%
  filter(Tumor=="16630")
S19_7_n<-S19_7%>%
  count(mutid, Tumor)
S19_7_counts<-pivot_wider(S19_7_n, names_from=Tumor, values_from=n)

S19_8<-S19%>%
  filter(Tumor=="16631") 
S19_8_n<-S19_8%>%
  count(mutid, Tumor)
S19_8_counts<-pivot_wider(S19_8_n, names_from=Tumor, values_from=n)


S19_upset<-full_join(S19_1_counts, S19_2_counts)
S19_upset2<-full_join(S19_upset, S19_3_counts)
S19_upset3<-full_join(S19_upset2, S19_4_counts)
S19_upset4<-full_join(S19_upset3, S19_5_counts)
S19_upset5<-full_join(S19_upset4, S19_6_counts )
S19_upset6<-full_join(S19_upset5, S19_7_counts)
S19_upset7<-full_join(S19_upset6, S19_8_counts)

#replace NAs with 0 
S19_upset7[is.na(S19_upset7)] <- 0
S19_upset7<-as.data.frame(S19_upset7)

#rename column to remove numeric constant 
colnames(S19_upset7)[2]="S19_1"
colnames(S19_upset7)[3]="S19_2"
colnames(S19_upset7)[4]="S19_3"
colnames(S19_upset7)[5]="S19_4"
colnames(S19_upset7)[6]="S19_5"
colnames(S19_upset7)[7]="S19_6"
colnames(S19_upset7)[8]="S19_7"
colnames(S19_upset7)[9]="S19_8"

#convert all numbers greater than 1 to 1 
S19_upset7$S19_1<-ifelse(S19_upset7$S19_1,1,S19_upset7$S19_1)
S19_upset7$S19_2<-ifelse(S19_upset7$S19_2,1,S19_upset7$S19_2)
S19_upset7$S19_3<-ifelse(S19_upset7$S19_3,1,S19_upset7$S19_3)
S19_upset7$S19_4<-ifelse(S19_upset7$S19_4,1,S19_upset7$S19_4)
S19_upset7$S19_5<-ifelse(S19_upset7$S19_5,1,S19_upset7$S19_5)
S19_upset7$S19_6<-ifelse(S19_upset7$S19_6,1, S19_upset7$S19_6)
S19_upset7$S19_7<-ifelse(S19_upset7$S19_7,1, S19_upset7$S19_7)
S19_upset7$S19_8<-ifelse(S19_upset7$S19_8,1, S19_upset7$S19_8)

#load upset package
#library(UpSetR)
set_vars<-c("S19_1", "S19_2", "S19_4", "S18_5")
S19_List<-upset(S19_upset7, keep.order=T, sets=c("S19_1","S19_2", "S19_3", "S19_4", "S19_5", "S19_6", "S19_7", "S19_8"))

#extract sets 
Data19<-as.data.frame(S19_List$New_data)
Data19<-Data19%>%
  mutate(ID=row_number())

#mutations where the row sum is 1 
private_S19<-rowSums(Data19[,2:8])==1
private_S19<-as.data.frame(private_S19)
private_S19<-private_S19%>%
  mutate(ID=row_number())


S19_Identi<-merge(Data19, private_S19, by=c("ID"))
S19_private<-S19_Identi%>%
  filter(private_S19=="TRUE")

#add in colony and mutation type 

S19_private<-S19_private%>%
  add_column(Intersection="Private")


#merge with original dataframe
private_mutationsS19<-merge(whole_filter, S19_private, by=c("mutid"))

ggplot(private_mutationsS19, aes(x=Tumor_AF))+geom_histogram()+xlim(0,1.0)+theme_bw()

#mutations where the row sum is 2
private_S192<-rowSums(Data19[,2:8])==2

private_S192<-as.data.frame(private_S192)
private_S192<-private_S192%>%
  mutate(ID=row_number())
S19_Identi2<-merge(Data19, private_S192, by=c("ID"))
S19_private2<-S19_Identi2%>%
  filter(private_S192=="TRUE")




#add in colony and mutation type 


S19_private2<-S19_private2%>%
  add_column(Intersection="Two")

twoshare19<-merge(whole_filter, S19_private2, by=c("mutid"))

ggplot(twoshare19, aes(x=Tumor_AF))+geom_histogram()+xlim(0,1.0)+theme_bw()

#three share
private_S193<-rowSums(Data19[,2:8])==3

private_S193<-as.data.frame(private_S193)
private_S193<-private_S193%>%
  mutate(ID=row_number())
S19_Identi3<-merge(Data19, private_S193, by=c("ID"))
S19_private3<-S19_Identi3%>%
  filter(private_S193=="TRUE")




#add in colony and mutation type 


S19_private3<-S19_private3%>%
  add_column(Intersection="Three")

threeshare19<-merge(whole_filter, S19_private3, by=c("mutid"))

ggplot(threeshare19, aes(x=Tumor_AF))+geom_histogram()+xlim(0,1.0)+theme_bw()


#P2194

P2194_1<-P2194%>%
  filter(Tumor=="5838")

P2194_1_n<-P2194_1%>%
  count(mutid, Tumor)

P2194_1_counts<-pivot_wider(P2194_1_n, names_from = Tumor, values_from = n)


P2194_2<-P2194%>%
  filter(Tumor=="5839")
P2194_2_n<-P2194_2%>%
  count(mutid, Tumor)
P2194_2_counts<-pivot_wider(P2194_2_n, names_from = Tumor, values_from = n)

P2194_3<-P2194%>%
  filter(Tumor=="5840")
P2194_3_n<-P2194_3%>%
  count(mutid, Tumor)
P2194_3_counts<-pivot_wider(P2194_3_n, names_from = Tumor, values_from = n)

P2194_4<-P2194%>%
  filter(Tumor=="5841")
P2194_4_n<-P2194_4%>%
  count(mutid, Tumor)
P2194_4_counts<-pivot_wider(P2194_4_n, names_from = Tumor, values_from = n)

P2194_5<-P2194%>%
  filter(Tumor=="5843")
P2194_5_n<-P2194_5%>%
  count(mutid,Tumor)
P2194_5_counts<-pivot_wider(P2194_5_n, names_from=Tumor, values_from=n)

P2194_6<-P2194%>%
  filter(Tumor=="5845")
P2194_6_n<-P2194_6%>%
  count(mutid, Tumor)
P2194_6_counts<-pivot_wider(P2194_6_n, names_from=Tumor, values_from=n)

P2194_7<-P2194%>%
  filter(Tumor=="5847")
P2194_7_n<-P2194_7%>%
  count(mutid, Tumor)
P2194_7_counts<-pivot_wider(P2194_7_n, names_from=Tumor, values_from=n)

P2194_8<-P2194%>%
  filter(Tumor=="5848")
P2194_8_n<-P2194_8%>%
  count(mutid, Tumor)
P2194_8_counts<-pivot_wider(P2194_8_n, names_from=Tumor, values_from = n)

P2194_9<-P2194%>%
  filter(Tumor=="5849")
P2194_9_n<-P2194_9%>%
  count(mutid, Tumor)
P2194_9_counts<-pivot_wider(P2194_9_n,names_from = Tumor, values_from = n)

P2194_upset<-full_join(P2194_1_counts, P2194_2_counts)
P2194_upset2<-full_join(P2194_upset, P2194_3_counts)
P2194_upset3<-full_join(P2194_upset2, P2194_4_counts)
P2194_upset4<-full_join(P2194_upset3, P2194_5_counts)
P2194_upset5<-full_join(P2194_upset4, P2194_6_counts)
P2194_upset6<-full_join(P2194_upset5, P2194_7_counts)
P2194_upset7<-full_join(P2194_upset6, P2194_8_counts)
P2194_upset8<-full_join(P2194_upset7, P2194_9_counts)

#replace NAs with 0 
P2194_upset8[is.na(P2194_upset8)] <- 0
P2194_upset8<-as.data.frame(P2194_upset8)
#rename columsn to remove numeric constant 

colnames(P2194_upset8)[2]="P2194_1"
colnames(P2194_upset8)[3]="P2194_2"
colnames(P2194_upset8)[4]="P2194_3"
colnames(P2194_upset8)[5]="P2194_4"
colnames(P2194_upset8)[6]="P2194_5"
colnames(P2194_upset8)[7]="P2194_6"
colnames(P2194_upset8)[8]="P2194_7"
colnames(P2194_upset8)[9]="P2194_8"
colnames(P2194_upset8)[10]="P2194_9"

#convert all numbers greater than 1 to 1 
P2194_upset8$P2194_1<-ifelse(P2194_upset8$P2194_1,1,P2194_upset8$P2194_1)
P2194_upset8$P2194_2<-ifelse(P2194_upset8$P2194_2,1,P2194_upset8$P2194_2)
P2194_upset8$P2194_3<-ifelse(P2194_upset8$P2194_3,1,P2194_upset8$P2194_3)
P2194_upset8$P2194_4<-ifelse(P2194_upset8$P2194_4,1,P2194_upset8$P2194_4)
P2194_upset8$P2194_5<-ifelse(P2194_upset8$P2194_5,1,P2194_upset8$P2194_5)
P2194_upset8$P2194_6<-ifelse(P2194_upset8$P2194_6,1, P2194_upset8$P2194_6)
P2194_upset8$P2194_7<-ifelse(P2194_upset8$P2194_7,1, P2194_upset8$P2194_7)
P2194_upset8$P2194_8<-ifelse(P2194_upset8$P2194_8,1, P2194_upset8$P2194_8)
P2194_upset8$P2194_9<-ifelse(P2194_upset8$P2194_9,1, P2194_upset8$P2194_9)


#load upset package
 
P2194_list<-upset(P2194_upset8, keep.order=T, sets=c("P2194_1","P2194_2", "P2194_3", "P2194_4", "P2194_5", "P2194_6", "P2194_7", "P2194_8", "P2194_9"))

#extract sets 
Dataold<-as.data.frame(P2194_list$New_data)
Dataold<-Dataold%>%
  mutate(ID=row_number())
private_old<-rowSums(Dataold[,2:8])==1

private_old<-as.data.frame(private_old)
private_old<-private_old%>%
  mutate(ID=row_number())
old_Identi<-merge(Dataold, private_old, by=c("ID"))
old_private<-old_Identi%>%
  filter(private_old=="TRUE")


#add in colony and mutation type 


old_private<-old_private%>%
  add_column(Intersection="Private")

oldprivate<-merge(whole_filter, old_private, by=c("mutid"))

ggplot(oldprivate, aes(x=Tumor_AF))+geom_histogram()+xlim(0,1.0)+theme_bw()

#twoshare

private_old2<-rowSums(Dataold[,2:8])==2

private_old2<-as.data.frame(private_old2)
private_old2<-private_old2%>%
  mutate(ID=row_number())
old_Identi2<-merge(Dataold, private_old2, by=c("ID"))
old_private2<-old_Identi2%>%
  filter(private_old2=="TRUE")


#add in colony and mutation type 


old_private2<-old_private2%>%
  add_column(Intersection="Private")

oldprivate2<-merge(whole_filter, old_private2, by=c("mutid"))

ggplot(oldprivate2, aes(x=Tumor_AF))+geom_histogram()+xlim(0,1.0)+theme_bw()



#calculating the number in genic regions 
#load in data 

genic<-Mutation_Annotation_Acro

#counting how many unique mutations in a genic region each tumor sample has 
genic2<-genic%>%
  group_by(sampleID)%>%
  count(sampleID)

colnames(genic2)[1]<-"Tumor"

#counting how many total unique mutations each tumor sample has 
 
all<-burden_whole
all<-burdenwhole2


genic3<-merge(genic2, all, by=c("Tumor"))
colnames(genic3)[2]<-"Total_Genic_Mutations"
colnames(genic3)[3]<-"Total_Unique_Mutations"

genic3$propgen<-genic3$Total_Genic_Mutations/genic3$Total_Unique_Mutations

write.csv(genic3, file="Proportion_Genic_Sample.csv")
gen<-Proportion_Genic_Sample
ggplot(gen, aes(x=Age, y=propgen))+geom_point()+ylim(0,.2)+theme_bw()+geom_hline(yintercept=0.13,
                                                                                 linetype="dashed", size=1)



gen$Age<-factor(gen$Age, levels=c("6", "8", "10", ">100"))
chisq.test(gen$propgen)

#mutation annotation 
af<-AF
annotation<-Mutation_Annotation_Acro

colnames(af)[3]<-"chr"
colnames(af)[4]<-"pos"
colnames(af)[11]<-"sampleID"

af_annotation<-merge(af,annotation, by=c("sampleID","chr","pos"))

af_annotation<-af_annotation[c("sampleID", "chr", "pos", "Tumor_AF", "impact", "Colony", "Age")]

anno_6<-af_annotation%>%
  filter(Age=="6")
anno_8<-af_annotation%>%
  filter(Age=="8")
anno_10<-af_annotation%>%
  filter(Age=="10")
anno_old<-af_annotation%>%
  filter(Age==">100")


#plot frequency of different impacts

ggplot(anno_6, aes(x=Tumor_AF, fill=impact, group=impact))+geom_density(alpha=0.6
                                                                        )+xlim(0,1.0)+theme_bw()+
  scale_fill_manual(values=wes_palette("GrandBudapest2", n=3))
)))

ggplot(anno_8,  aes(x=Tumor_AF, fill=impact, group=impact))+geom_density(alpha=0.6)+xlim(0,1.0)+theme_bw()+
  scale_fill_manual(values=wes_palette("GrandBudapest2", n=3))

ggplot(anno_10,  aes(x=Tumor_AF, fill=impact, group=impact))+geom_density(alpha=0.6)+xlim(0,1.0)+theme_bw()+
  scale_fill_manual(values=wes_palette("GrandBudapest2", n=3))

ggplot(anno_old,  aes(x=Tumor_AF, fill=impact, group=impact))+geom_density(alpha=0.6)+xlim(0,1.0)+theme_bw()+
  scale_fill_manual(values=wes_palette("GrandBudapest2", n=3))


##export dataframe

write.csv(af_annotation, file="AlleleFreq_withAnnotation.csv")


##plotting individual variant allele frequencies 

#filter by depth
af<-V4

af_6<-af%>%
  filter(Age=="6")
af_8<-af%>%
  filter(Age=="8")
af_10<-af%>%
  filter(Age=="10")
af_old<-af%>%
  filter(Age==">100")
ggplot(af_10, aes(x=Tumor_AF))+geom_histogram(bins=40)+geom_density()+theme_bw()+xlim(0,1.0)



##
ggplot(Neutrality_Metrics, aes(x=Age, y=R2))+geom_point(size=3)+theme_bw()

Neutrality_Metrics$Age<-factor(Neutrality_Metrics$Age, levels=c("6", "8", "10", ">100"))
