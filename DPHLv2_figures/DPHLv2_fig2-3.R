rm(list = ls())
options(stringsAsFactors = F)
library(ggplot2)
library(plyr)
library(reshape2)
library(tidyverse)
library(grid)
library(futile.logger)
library(VennDiagram)
library(readr)
library(data.table)
library(UpSetR)
require(gridExtra)
library(circlize)
total<-read.csv("20210701_summary_of 4lib.csv",stringsAsFactors = F,header = T)
total1<-melt(total,c("library"),na.rm=T,factorsAsStrings=TRUE)
# total1<-data.frame(t(total))
# total1$name<-row.names(total1)
a<-ggplot(total1,aes(x=library,y=value,fill=variable))+geom_bar(stat="identity",position = 'stack',colour="black")+
  coord_polar(direction=-1)+scale_fill_manual(values=c("#F7ACA9","#BCA4EA","#A0C19A"))+
  geom_text(aes(label=value),size=4,position = position_stack())+guides(fill=guide_legend(reverse=F))
a

ggsave("20210421_DPHLv2_summary_4lib_compair.pdf", plot=a,device = NULL, width = 12, height = 8)

total2<- total %>% mutate(library = fct_reorder(library,proteins ))
a<-ggplot(total2[,c(1,2)],aes(x=library,y=proteins),size=8)+geom_bar(stat="identity",fill="gray",colour="black",width=0.4)+
  geom_text(aes(label=proteins),size=5,vjust=-0.2)+theme(axis.text = element_text(size=rel(1.2)))+theme_bw()+
  theme(axis.title.x =element_text(size=14), axis.title.y=element_text(size=14))+
  theme(axis.text.x = element_text(size = 14,color="black"),axis.text.y = element_text(size = 14,color="black")) 
a
ggsave("20210702_DPHLv2_proteins_4lib_compair.pdf", plot=a,device = NULL, width = 8, height = 8)

total2<- total %>% mutate(library = fct_reorder(library,peptides ))
a<-ggplot(total2[,c(1,3)],aes(x=library,y=peptides),size=8)+geom_bar(stat="identity",fill="gray",colour="black",width=0.4)+
  geom_text(aes(label=peptides),size=5,vjust=-0.2)+theme(axis.text = element_text(size=rel(1.2)))+theme_bw()+
  theme(axis.title.x =element_text(size=14), axis.title.y=element_text(size=14))+
  theme(axis.text.x = element_text(size = 14,color="black"),axis.text.y = element_text(size = 14,color="black")) 
a
ggsave("20210702_DPHLv2_peptides_4lib_compair.pdf", plot=a,device = NULL, width = 8, height = 8)

total2<- total %>% mutate(library = fct_reorder(library,precursors ))
a<-ggplot(total2[,c(1,4)],aes(x=library,y=precursors),size=8)+geom_bar(stat="identity",fill="gray",colour="black",width=0.4)+
  geom_text(aes(label=precursors),size=5,vjust=-0.2)+theme(axis.text = element_text(size=rel(1.2)))+theme_bw() +
  theme(axis.title.x =element_text(size=14), axis.title.y=element_text(size=14))+
  theme(axis.text.x = element_text(size = 14,color="black"),axis.text.y = element_text(size = 14,color="black"))
a
ggsave("20210702_DPHLv2_transitions_4lib_compair.pdf", plot=a,device = NULL, width = 8, height = 8)



#DPHL1、PHL####
DPHLv1<-read.csv("20181226DPHL_v6.1.csv_library_0421.csv",stringsAsFactors = F,header = T,row.names = 1)
PHL<-read.csv("phl004_canonical_s64_osw_decoys.csv_library_0421.csv",stringsAsFactors = F,header = T,row.names = 1)
DPHLv1$transition_group_id<-gsub("^\\d+_","",DPHLv1$transition_group_id)
DPHLv1$ProteinName<-gsub("1/sp|","",DPHLv1$ProteinName,fixed = T)
DPHLv1$ProteinName<-gsub("\\|[a-zA-Z0-9]+_HUMAN","",DPHLv1$ProteinName,perl  = T)
DPHLv1<-DPHLv1[!grepl("iRT",DPHLv1$ProteinName),]
PHL$transition_group_id<-gsub("^\\d+_","",PHL$transition_group_id)
PHL$ProteinName<-gsub("1/","",PHL$ProteinName,fixed = T)
write.csv(DPHLv1,"20210421_DPHL1_use_matrix.csv")
write.csv(PHL,"20210421_PHL_use_matrix.csv")
total<-total<-data.frame(matrix(NA,nrow=3,ncol=2))
row.names(total)<-c("proteins","peptides","precursors")
colnames(total)<-c("DPHL1","Phl")
total[1,1]<-length(unique(DPHLv1$ProteinName))
total[2,1]<-length(unique(DPHLv1$FullUniModPeptideName))
total[3,1]<-length(unique(DPHLv1$transition_group_id))
total[1,2]<-length(unique(PHL$ProteinName))
total[2,2]<-length(unique(PHL$FullUniModPeptideName))
total[3,2]<-length(unique(PHL$transition_group_id))
write.csv(total,"20210526_DPHL2_sum_of_dphl1_phl.csv")


#The subset of DPHL2, DPHL1 and PHL without decory  and group###
ReFull<-read.table("reviewedfull_library_QC.tsv",sep="\t",header = T)
ReSemi<-read.table("reviewedsemi_library_QC.tsv",sep="\t",header = T)
IsoFull<-read.table("isoformfull_library_QC.tsv",sep="\t",header = T)
IsoSemi<-read.table("isoformsemi_library_QC.tsv",sep="\t",header = T)
DPHL1<-read.csv("20210421_DPHL1_use_matrix.csv",header = T,row.names = 1)
PHL<-read.csv("20210421_PHL_use_matrix.csv",header = T,row.names = 1)
# ReFull1<-ReFull[ReFull$Decoy==0,]
# ReSemi1<-ReSemi[ReSemi$Decoy==0,]
# IsoFull1<-IsoFull[IsoFull$Decoy==0,]
# IsoSemi1<-IsoSemi[IsoSemi$Decoy==0,]

#fig 2B
listinput <- list(reviewedFull=unique(ReFull$ProteinId),
                  reviwedSemi=unique(ReSemi$ProteinId),
                  isoformFull=unique(IsoFull$ProteinId),
                  isoformSemi=unique(IsoSemi$ProteinId),
                  DPHL1=unique(DPHL1$ProteinName),
                  PHL=unique(PHL$ProteinName))

pdf('20210705_DPHLv2_proteins_overlap_6_lib.pdf',height = 8,width = 8)
upset(fromList(listinput),nsets = 6, order.by = "freq",text.scale =1.5)
dev.off()

#fig 2D
listinput <- list(reviewedFull=unique(ReFull$ModifiedPeptideSequence),
                  reviwedSemi=unique(ReSemi$ModifiedPeptideSequence),
                  isoformFull=unique(IsoFull$ModifiedPeptideSequence),
                  isoformSemi=unique(IsoSemi$ModifiedPeptideSequence),
                  DPHL1=unique(DPHL1$FullUniModPeptideName),
                  PHL=unique(PHL$FullUniModPeptideName))


pdf('20210705_DPHLv2_peptides_overlap_6_lib.pdf',height = 8,width = 8)####
upset(fromList(listinput),nsets = 6, order.by = "freq",text.scale =1.5)
dev.off()

ReFull$TransitionGroupId<-paste(ReFull$ModifiedPeptideSequence,ReFull$PrecursorCharge,sep="_")
ReSemi$TransitionGroupId<-paste(ReSemi$ModifiedPeptideSequence,ReSemi$PrecursorCharge,sep="_")
IsoFull$TransitionGroupId<-paste(IsoFull$ModifiedPeptideSequence,IsoFull$PrecursorCharge,sep="_")
IsoSemi$TransitionGroupId<-paste(IsoSemi$ModifiedPeptideSequence,IsoSemi$PrecursorCharge,sep="_")


#fig 2F
listinput <- list(reviewedFull=unique(ReFull$TransitionGroupId),
                  reviwedSemi=unique(ReSemi$TransitionGroupId),
                  isoformFull=unique(IsoFull$TransitionGroupId),
                  isoformSemi=unique(IsoSemi$TransitionGroupId),
                  DPHL1=unique(DPHL1$transition_group_id),
                  PHL=unique(PHL$transition_group_id))


pdf('20210705_DPHLv2_precursors_overlap_6_lib.pdf',height = 8,width = 10)
upset(fromList(listinput),nsets = 6, order.by = "freq",text.scale =1.5)
dev.off()


#4
diff<-read.csv("20210422_DPHL2_6lib_diff.csv",header = T)
# difplot<-melt(diff,c("library"),na.rm=T)
a<-ggplot(diff[,c(1,2,5)],aes(x=library,y=proteins,fill=name))+geom_bar(stat="identity",width=0.7,position = position_dodge(width=0.9),
          colour="black")+geom_text(aes(label=proteins),size=4,position = position_dodge(0.9),vjust = -0.25)+guides(fill=guide_legend(reverse=F))+
  scale_fill_brewer(palette="RdBu")+annotate("text",x=-Inf,y=Inf,label="",hjust=-0.2,vjust=2)
a
ggsave("20210422_DPHLv2_4lib_diff_dphl1_phl_compair_proteins_1.pdf", plot=a,device = NULL, width = 8, height = 8)

a<-ggplot(diff[,c(1,3,5)],aes(x=library,y=peptides,fill=name))+geom_bar(stat="identity",width=0.7,position = position_dodge(width=0.9),
          colour="black")+geom_text(aes(label=peptides),size=4,position = position_dodge(0.9),vjust = -0.25)+guides(fill=guide_legend(reverse=F))+
  scale_fill_brewer(palette="RdBu")+annotate("text",x=-Inf,y=Inf,label="",hjust=-0.2,vjust=2)
a
ggsave("20210422_DPHLv2_4lib_diff_dphl1_phl_compair_peptides_1.pdf", plot=a,device = NULL, width = 8, height = 8)
a<-ggplot(diff[,c(1,4,5)],aes(x=library,y=precursors,fill=name))+geom_bar(stat="identity",width=0.7,position = position_dodge(width=0.9),
         colour="black")+geom_text(aes(label=precursors),size=4,position = position_dodge(0.9),vjust = -0.25)+guides(fill=guide_legend(reverse=F))+
  scale_fill_brewer(palette="RdBu")+annotate("text",x=-Inf,y=Inf,label="",hjust=-0.2,vjust=2)
a
ggsave("20210422_DPHLv2_4lib_diff_dphl1_phl_compair_precursor_1.pdf", plot=a,device = NULL, width = 8, height = 8)

######


library(circlize)
alllib<-read.csv("20210704_summary_of 6lib.csv",header = T)
# alllib<- alllib %>% mutate(library = fct_reorder(library,proteins ))
color<-rev(rainbow(length(alllib$library)))

#fig 2A
#protein
pdf("20210704_6lib_protein_number_compair.pdf",width =8,height = 8)
circos.par("start.degree" = 90, cell.padding = c(0, 0, 0, 0))
circos.initialize("a",xlim=c(0,16000))
circos.track(ylim = c(0.5, length(alllib$proteins)+0.5), track.height = 0.8, 
             bg.border = NA, panel.fun = function(x, y) {
               xlim = CELL_META$xlim
               circos.segments(rep(xlim[1], 6), 1:6,
                               rep(xlim[2], 6), 1:6,
                               col = "#CCCCCC")
               circos.rect(rep(0, 6), 1:6 - 0.45, alllib$proteins, 1:6 + 0.45,
                           col = color, border = "white")
               circos.text(rep(xlim[1], 6), 1:6, 
                           paste(alllib$library, " - ", alllib$proteins), 
                           facing = "downward", adj = c(1.05, 0.5), cex = 1.2) 
               breaks = seq(0, 16000, by = 1000)
               circos.axis(h = "top", major.at = breaks, labels = breaks, 
                           labels.cex = 1.2)
               
             })
circos.clear()
dev.off()

#fig 2C
#peptide
pdf("20210704_6lib_peptides_number_compair.pdf",width =8,height = 8)
circos.par("start.degree" = 90, cell.padding = c(0, 0, 0, 0))
circos.initialize("a",xlim=c(0,700000))
circos.track(ylim = c(0.5, length(alllib$peptides)+0.5), track.height = 0.8, 
             bg.border = NA, panel.fun = function(x, y) {
               xlim = CELL_META$xlim
               circos.segments(rep(xlim[1], 6), 1:6,
                               rep(xlim[2], 6), 1:6,
                               col = "#CCCCCC")
               circos.rect(rep(0, 6), 1:6 - 0.45, alllib$peptides, 1:6 + 0.45,
                           col = color, border = "white")
               circos.text(rep(xlim[1], 6), 1:6, 
                           paste(alllib$library, " - ", alllib$peptides), 
                           facing = "downward", adj = c(1.05, 0.5), cex = 1.2) 
               breaks = seq(0, 700000, by = 31000)
               circos.axis(h = "top", major.at = breaks, labels = breaks, 
                           labels.cex = 1.2)
               
             })
circos.clear()
dev.off()

#fig 2E
#precursor
pdf("20210704_6lib_precursors_number_compair.pdf",width =8,height = 8)
circos.par("start.degree" = 90, cell.padding = c(0, 0, 0, 0))
circos.initialize("a",xlim=c(0,900000))
circos.track(ylim = c(0.5, length(alllib$precursors)+0.5), track.height = 0.8, 
             bg.border = NA, panel.fun = function(x, y) {
               xlim = CELL_META$xlim
               circos.segments(rep(xlim[1], 6), 1:6,
                               rep(xlim[2], 6), 1:6,
                               col = "#CCCCCC")
               circos.rect(rep(0, 6), 1:6 - 0.45, alllib$precursors, 1:6 + 0.45,
                           col = color, border = "white")
               circos.text(rep(xlim[1], 6), 1:6, 
                           paste(alllib$library, " - ", alllib$precursors), 
                           facing = "downward", adj = c(1.05, 0.5), cex = 1.2) 
               breaks = seq(0, 900000, by = 50000)
               circos.axis(h = "top", major.at = breaks, labels = breaks, 
                           labels.cex = 1.2)
               
             })
circos.clear()
dev.off()


#fig 3A and 3B
#proteins and prptides overlap 
protsum<-read.csv("4lib_prot_overlap_numb.csv",header = T)
library(ggplot2)

a<-ggplot(protsum,aes(x=name,y=prot,fill=lable)) + geom_bar(stat="identity",position="stack")
ggsave("4lib_prot_overlap.pdf", plot=a,device = NULL, width = 12, height = 8)

pepsum<-read.csv("4lib_pep_overlap_numb.csv",header = T)
a<-ggplot(pepsum,aes(x=name,y=pep,fill=lable)) + geom_bar(stat="identity",position="stack")
ggsave("4lib_pep_overlap.pdf", plot=a,device = NULL, width = 12, height = 8)





#############################
#fig 3C
#spec list
# setwd("E:/Qsync/work/DPHL2/20210228_DPHL2_buildlib_QC")
hall<-read.csv("hallmark_gene2pro.csv",header = T)
ovay<-read.csv("ovary_spec_prot_list.csv",header = T)
esop<-read.csv("esophagus_spec_prot_list.csv",header = T)
FDA<-read.csv("FDA_gene_name.csv",header = T)
brain<-read.csv("brain_spec_prot_list.csv",header = T)
setwd("E:/Qsync/work/DPHL2/20210228_DPHL2_buildlib_QC/reviewedFull/20210607_reviewedfull_new/sub_lib")
lib<-read.csv("reviewedfull_lib_prot_gnen.csv",header = T)
hallP<-lib[lib$ProteinId%in%hall$x,]
ovayP<-lib[lib$ProteinId%in%ovay$x,]
esopP<-lib[lib$ProteinId%in%esop$x,]
FDAP<-lib[lib$GeneName%in%FDA$Gene,]
brainP<-lib[lib$ProteinId%in%brain$x,]
setwd("E:/Qsync/work/DPHL2/20210228_DPHL2_buildlib_QC/reviewedSemi/20210608_reviewedsemi_new/sub_lib")
lib<-read.csv("reviewedsemi_lib_prot_gnen.csv",header = T)
hallP<-lib[lib$ProteinId%in%hall$x,]
ovayP<-lib[lib$ProteinId%in%ovay$x,]
esopP<-lib[lib$ProteinId%in%esop$x,]
FDAP<-lib[lib$GeneName%in%FDA$Gene,]
brainP<-lib[lib$ProteinId%in%brain$x,]
setwd("E:/Qsync/work/DPHL2/20210228_DPHL2_buildlib_QC/isoformfull/20210608_isoformfull_new/sub_lib")
lib<-read.csv("isoformfull_lib_prot_gnen.csv",header = T)
hallP<-lib[lib$ProteinId%in%hall$x,]
ovayP<-lib[lib$ProteinId%in%ovay$x,]
esopP<-lib[lib$ProteinId%in%esop$x,]
FDAP<-lib[lib$GeneName%in%FDA$Gene,]
brainP<-lib[lib$ProteinId%in%brain$x,]
setwd("E:/Qsync/work/DPHL2/20210228_DPHL2_buildlib_QC/isoformsemi/20210609_isoformsemi_new/sub_lib")
lib<-read.csv("isoformsemi_lib_prot_gnen.csv",header = T)
hallP<-lib[lib$ProteinId%in%hall$x,]
ovayP<-lib[lib$ProteinId%in%ovay$x,]
esopP<-lib[lib$ProteinId%in%esop$x,]
FDAP<-lib[lib$GeneName%in%FDA$Gene,]
brainP<-lib[lib$ProteinId%in%brain$x,]
#DPHL1
setwd("E:/Qsync/work/DPHL2/DPHL1")
lib<-read.csv("dphl1_prot_gene_list.csv",header = T)
hallP<-lib[lib$uniprot%in%hall$x,]
ovayP<-lib[lib$uniprot%in%ovay$x,]
esopP<-lib[lib$uniprot%in%esop$x,]
FDAP<-lib[lib$gene%in%FDA$Gene,]
brainP<-lib[lib$uniprot%in%brain$x,]
####################

####################
#spec list coverage
lib<-read.csv("dphl2_speclist_overlap_summary.csv",header = T)

lib1<-melt(lib,value.name = "value")

library(ggplot2)

a<-ggplot(lib1, aes(x = variable, y = value,fill=name))+geom_bar(stat="identity",position=position_dodge2(0.75))
a
ggsave("dphl2_speclist_overlap_summary.pdf",plot=a,device = NULL, width = 8, height = 10)

###%
lib<-read.csv("dphl2_speclist_overlap_summary%.csv",header = T)
lib1<-melt(lib,value.name = "value")
library(ggplot2)
a<-ggplot(lib1, aes(x = variable, y = value,fill=name))+geom_bar(stat="identity",position=position_dodge2(0.75))+
  scale_y_continuous(expand = c(0, 0))+
  scale_fill_manual(values =c("#FC4E07","#00AFBB", "#E7B800","#4682B4" ,"#A193BF"))+theme_bw()+coord_flip(ylim = c(0,1))
a

ggsave("dphl2_speclist_overlap_summary1.pdf",plot=a,device = NULL, width = 8, height = 10)





