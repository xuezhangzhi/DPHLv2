type<-read.csv("rawdata_type.csv",header = T)
setwd("E:/Qsync/work/DPHL2/20210228_DPHL2_buildlib_QC/reviewedFull/20210607_reviewedfull_new/sub_lib")
revfullpsm<-read.csv("../psm.tsv",sep = "\t",header = T)

revfullpsm$name<-sapply(revfullpsm$Spectrum,function(e){unlist(strsplit(e,"\\."))[1]} ) 
revfull<-read.csv("../reviewedfull_library_QC.tsv",header = T,sep = "\t")
sel<-revfullpsm[,c(28,34)]
sel1<-merge(sel,type,by="name",all.x = T)
write.csv(sel1,"reviewedfull_psm_prot_file_all.csv",row.names = F)
a<-unique(sel1$type)
for (i in a) {
  indsel<-unique(sel1$Protein.ID[sel1$type==i])
  lib<-revfull[revfull$ProteinId%in%indsel,]
  write.table(lib,paste0(i,"_sub_library.tsv"),sep="\t",row.names=F,col.names=T,na="",quote=F)
}


setwd("E:/Qsync/work/DPHL2/20210228_DPHL2_buildlib_QC/isoformfull/20210608_isoformfull_new/sub_lib")

revfullpsm<-read.csv("../psm.tsv",sep = "\t",header = T)

revfullpsm$name<-sapply(revfullpsm$Spectrum,function(e){unlist(strsplit(e,"\\."))[1]} ) 
revfull<-read.csv("../isoformfull_library_QC.tsv",header = T,sep = "\t")
sel<-revfullpsm[,c(28,34)]
sel1<-merge(sel,type,by="name",all.x = T)
write.csv(sel1,"isoformfull_psm_prot_file_all.csv",row.names = F)
a<-unique(sel1$type)
for (i in a) {
  indsel<-unique(sel1$Protein.ID[sel1$type==i])
  lib<-revfull[revfull$ProteinId%in%indsel,]
  write.table(lib,paste0(i,"_isofull_sub_library.tsv"),sep="\t",row.names=F,col.names=T,na="",quote=F)
}


setwd("E:/Qsync/work/DPHL2/20210228_DPHL2_buildlib_QC/reviewedSemi/20210608_reviewedsemi_new/sub_lib")


revfullpsm<-read.csv("../psm.tsv",sep = "\t",header = T)

revfullpsm$name<-sapply(revfullpsm$Spectrum,function(e){unlist(strsplit(e,"\\."))[1]} ) 
revfull<-read.csv("../reviewedsemi_library_QC.tsv",header = T,sep = "\t")
sel<-revfullpsm[,c(28,34)]
sel1<-merge(sel,type,by="name",all.x = T)
write.csv(sel1,"reviewedsemi_psm_prot_file_all.csv",row.names = F)
a<-unique(sel1$type)
for (i in a) {
  indsel<-unique(sel1$Protein.ID[sel1$type==i])
  lib<-revfull[revfull$ProteinId%in%indsel,]
  write.table(lib,paste0(i,"_revsemi_sub_library.tsv"),sep="\t",row.names=F,col.names=T,na="",quote=F)
}

setwd("E:/Qsync/work/DPHL2/20210228_DPHL2_buildlib_QC/isoformsemi/20210609_isoformsemi_new/sub_lib")

revfullpsm<-read.csv("../psm.tsv",sep = "\t",header = T)

revfullpsm$name<-sapply(revfullpsm$Spectrum,function(e){unlist(strsplit(e,"\\."))[1]} ) 
revfull<-read.csv("../isoformsemi_library_QC.tsv",header = T,sep = "\t")
sel<-revfullpsm[,c(28,34)]
sel1<-merge(sel,type,by="name",all.x = T)
write.csv(sel1,"isoformsemi_psm_prot_file_all.csv",row.names = F)
a<-unique(sel1$type)
for (i in a) {
  indsel<-unique(sel1$Protein.ID[sel1$type==i])
  lib<-revfull[revfull$ProteinId%in%indsel,]
  write.table(lib,paste0(i,"_isosemi_sub_library.tsv"),sep="\t",row.names=F,col.names=T,na="",quote=F)
}



################################
#number of samples
readin<-read.csv("dphl2_raw_file_number.csv",header = T)
readin$type<-paste0(readin$type,":",readin$numb)
#
library(ggplot2)
#barplot
p<-ggplot(readin,aes(x=as.factor(type),y=numb))+geom_bar(stat="identity",fill=alpha("#f88421",0.7))#目前还是不太清楚stat参数的作用
p<-p+coord_polar()
ggsave("dphl2_sample_sub_number.pdf",p,height = 8,width = 8)


df1<-data.frame(individual=paste("Mister",seq(1,60),sep=""),
                value=rep(c(sample(60:100,9,replace=T),NA),6))
df1$id<-seq(1,nrow(df1))
df1
df1$angle<-df$angle1
df1$hjust<-df$hjust
df1
df1$fill<-c(rep("A",10),rep("B",10),rep("C",10),rep("D",10),rep("E",10),rep("F",10))
ggplot(df1,aes(x=as.factor(id),y=value))+
  geom_bar(stat="identity",aes(fill=fill))+
  coord_polar()+ylim(-100,120)+
  geom_text(aes(x=id,y=value+20,label=individual,
                angle=angle,hjust=hjust),size=3)+
  theme_minimal()+ylab("")+xlab("")+
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_blank(),
        panel.grid = element_blank(),
        legend.position="none")+
  scale_fill_manual(values=c("red","yellow","blue","green","orange","skyblue"))




library(ggplot2)
p <- ggplot(readin, aes(x=as.factor(id), y=prot_numb, fill=type)) +       # Note that id is a factor. If x is numeric, there is some space between the first bar
  geom_bar(stat="identity", alpha=0.5) +
  ylim(-10000,12000) +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    plot.margin = unit(rep(-1,4), "cm")
  ) +
  coord_polar() +
  geom_text(data=label_data, aes(x=id, y=prot_numb, label=num, hjust=hjust), color="black", fontface="bold",alpha=0.6, size=2.5, angle= label_data$angle, inherit.aes = FALSE )
p
ggsave("dphl2_4lib_sub_lib_prot_summary.pdf",p,height = 10,width = 10)


p <- ggplot(readin, aes(x=as.factor(id), y=spec_numb, fill=type)) +       # Note that id is a factor. If x is numeric, there is some space between the first bar
  geom_bar(stat="identity", alpha=0.5) +
  ylim(-200,200) +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    plot.margin = unit(rep(-1,4), "cm")
  ) +
  coord_polar() +
  geom_text(data=label_data, aes(x=id, y=spec_numb, label=num, hjust=hjust), color="black", fontface="bold",alpha=0.6, size=2.5, angle= label_data$angle, inherit.aes = FALSE )
p
ggsave("dphl2_4lib_sub_lib__spec_prot_summary.pdf",p,height = 8,width = 8)



#586 new samples
lib<-read.csv("4lib_586_add_prot_summary.csv",header = T)
View(lib)
df <- tibble(
  gene = factor(paste0("gene_", rep(1:16, 2)), levels = paste0("gene_", 16:1)),
  stat = c(seq(-10, -100, -10), seq(-90, -40, 10), seq(10, 100, 10), seq(90, 40, -10)),
  direct = rep(c("down", "up"), each=16)
)
ggplot(df, aes(gene, stat, fill = direct)) +
  geom_col() +
  coord_flip() +
  scale_y_continuous(breaks = seq(-100, 100, 20),
                     labels = c(seq(100, 0, -20), seq(20, 100, 20)))
rm(list = ls())
options(stringsAsFactors = F)
library("ggplot2")

#the protein number of 586 samples
lib<-read.csv("4lib_586_add_prot_summary.csv",header = T)
library("ggplot2")
library(reshape2)
lib1<-melt(lib,value.name = "lib")
a<-ggplot(lib,aes(x=number,y=name))+geom_segment(aes(yend=name),xend=0,colour="grey50",size=0.5)+
  geom_point(size=4,colour="red")+theme_bw()

a<-ggplot(lib,aes(x=lib,y=name))+geom_segment(aes(yend=name),xend=0,colour="grey50",size=0.5)+
  geom_point(size=4,colour="red")+theme_bw()
a
a<-ggplot(lib1,aes(x=lib,y=name))+geom_segment(aes(yend=name),xend=0,colour="grey50",size=0.5)+
  geom_point(size=4,colour="red")+theme_bw()
a
ggsave("dphl2_4lib_586file_add_prot_summary.pdf", plot=a,device = NULL, width = 8, height = 10)
a<-ggplot(lib1,aes(x=lib,y=name))+geom_segment(aes(yend=name),xend=0,colour="grey50",size=1.2)+
  geom_point(size=4,colour="red")+theme_bw()
a
ggsave("dphl2_4lib_586file_add_prot_summary.pdf", plot=a,device = NULL, width = 8, height = 10)
#
lib<-read.csv("4lib_586_add_prot_summary.csv",header = T)
lib1<-melt(lib,value.name = "lib")
a<-ggplot(lib1,aes(x=lib,y=name))+geom_segment(aes(yend=name),xend=0,colour="grey50",size=1.2)+
  geom_point(size=4,colour="red")+theme_bw()
a

lib1$name<-as.factor(lib1$name)
a<-ggplot(lib1,aes(x=lib,y=name))+geom_segment(aes(yend=name),xend=0,colour="grey50",size=1.2)+
  geom_point(size=4,colour="red")+theme_bw()
a
############################
#
setwd("E:/Qsync/work/DPHL2/20210228_DPHL2_buildlib_QC")
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
#
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
  scale_y_continuous(expand = c(0, 0))+#消除x轴与绘图区的间隙
  scale_fill_manual(values =c("#FC4E07","#00AFBB", "#E7B800","#4682B4" ,"#A193BF"))+theme_bw()+coord_flip(ylim = c(0,1))
a

ggsave("dphl2_speclist_overlap_summary1.pdf",plot=a,device = NULL, width = 8, height = 10)
################################
#
readin<-read.csv("dphl2_raw_file_number.csv",header = T)
readin$type<-paste0(readin$type,":",readin$numb)

#
library(ggplot2)

a<-ggplot(readin, aes(x = numb, y = as.factor(type)))+geom_bar(stat="identity",position=position_dodge2(0.75))
a

readin<-readin[order(readin$numb,decreasing = T),]

a<-ggplot(readin, aes(x =numb, y =  reorder(type,- numb))+geom_bar(stat="identity",position=position_dodge2(0.75)))
################################
#
readin<-read.csv("dphl2_raw_file_number.csv",header = T)
readin$type<-paste0(readin$type,":",readin$numb)
readin<-readin[order(readin$numb,decreasing = T),]


#
p<-ggplot(readin,aes(x=as.factor(type),y=numb))+geom_bar(stat="identity",fill=alpha("#f88421",0.7))#目前还是不太清楚stat参数的作用
p<-p+coord_polar()
p
readin$type<-as.factor(readin$type)
a<-ggplot(readin, aes(x = numb, y =type))+geom_bar(stat="identity",position=position_dodge2(0.75))+coord_flip()
a
a<-ggplot(readin, aes(x =type, y =numb))+geom_bar(stat="identity",position=position_dodge2(0.75))+coord_flip()
a

################################
#
readin<-read.csv("dphl2_raw_file_number.csv",header = T)
readin$type<-paste0(readin$type,":",readin$numb)

a<-ggplot(readin, aes(x = reorder(type, numb), y = numb)) +geom_col(fill ="#00AFBB")+
  geom_bar(stat="identity",position=position_dodge2(0.75))+coord_flip()+theme_bw()
a
ggsave("dphl2_sample_sub_number_2.pdf",plot=a,device = NULL, width = 8, height = 10)

############################################
#
library(tidyverse)
readin<-read.csv("dphl2_4lib_sub_prot_summary.csv",header = T)

################
#reviewedfull
readin<-read.csv("reviewedfull_sub_prot_summary.csv",header = T)

#protein identify
a<-ggplot(readin, aes(x = reorder(X, prot_numb), y = prot_numb)) +geom_col(fill ="#00AFBB")+
  geom_bar(stat="identity",position=position_dodge2(0.75))+coord_flip()+theme_bw()+ geom_text(aes(label=prot_numb),size=4,vjust=0.5)
a
ggsave("dphl2_reviewedfull_sub_lib_prot_number.pdf",plot=a,device = NULL, width = 8, height = 10)
a<-ggplot(readin, aes(x = reorder(X, spec_numb), y = spec_numb)) +geom_col(fill ="#00AFBB")+
  geom_bar(stat="identity",position=position_dodge2(0.75))+coord_flip()+theme_bw()+ geom_text(aes(label=spec_numb),size=4,vjust=0.5)
a
ggsave("dphl2_reviewedfull_sub_lib_specprot_number.pdf",plot=a,device = NULL, width = 8, height = 10)

#reviewedsemi
readin<-read.csv("reviewedsemi_sub_prot_summary.csv",header = T)
#protein identify
a<-ggplot(readin, aes(x = reorder(X, prot_numb), y = prot_numb)) +geom_col(fill ="#00AFBB")+
  geom_bar(stat="identity",position=position_dodge2(0.75))+coord_flip()+theme_bw()+ geom_text(aes(label=prot_numb),size=4,vjust=0.5)
a
ggsave("dphl2_reviewedsemi_sub_lib_prot_number.pdf",plot=a,device = NULL, width = 8, height = 10)

a<-ggplot(readin, aes(x = reorder(X, spec_numb), y = spec_numb)) +geom_col(fill ="#00AFBB")+
  geom_bar(stat="identity",position=position_dodge2(0.75))+coord_flip()+theme_bw()+ geom_text(aes(label=spec_numb),size=4,vjust=0.5)
a
ggsave("dphl2_reviewedsemi_sub_lib_specprot_number.pdf",plot=a,device = NULL, width = 8, height = 10)


#isoformfull
readin<-read.csv("isoformfull_sub_prot_summary.csv",header = T)
#protein identify
a<-ggplot(readin, aes(x = reorder(X, prot_numb), y = prot_numb)) +geom_col(fill ="#00AFBB")+
  geom_bar(stat="identity",position=position_dodge2(0.75))+coord_flip()+theme_bw()+ geom_text(aes(label=prot_numb),size=4,vjust=0.5)
a
ggsave("dphl2_isoformfull_sub_lib_prot_number.pdf",plot=a,device = NULL, width = 8, height = 10)
a<-ggplot(readin, aes(x = reorder(X, spec_numb), y = spec_numb)) +geom_col(fill ="#00AFBB")+
  geom_bar(stat="identity",position=position_dodge2(0.75))+coord_flip()+theme_bw()+ geom_text(aes(label=spec_numb),size=4,vjust=0.5)
a
ggsave("dphl2_isoformdfull_sub_lib_specprot_number.pdf",plot=a,device = NULL, width = 8, height = 10)

#isoformsemi
readin<-read.csv("isoformsemi_sub_prot_summary.csv",header = T)
#protein identify
a<-ggplot(readin, aes(x = reorder(X, prot_numb), y = prot_numb)) +geom_col(fill ="#00AFBB")+
  geom_bar(stat="identity",position=position_dodge2(0.75))+coord_flip()+theme_bw()+ geom_text(aes(label=prot_numb),size=4,vjust=0.5)
a
ggsave("dphl2_isoformsemi_sub_lib_prot_number.pdf",plot=a,device = NULL, width = 8, height = 10)
a<-ggplot(readin, aes(x = reorder(X, spec_numb), y = spec_numb)) +geom_col(fill ="#00AFBB")+
  geom_bar(stat="identity",position=position_dodge2(0.75))+coord_flip()+theme_bw()+ geom_text(aes(label=spec_numb),size=4,vjust=0.5)
a
ggsave("dphl2_isoformsemi_sub_lib_specprot_number.pdf",plot=a,device = NULL, width = 8, height = 10)




























