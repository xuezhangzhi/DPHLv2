rm(list = ls())
options(stringsAsFactors = F)
# first_category_name = list.files() 
crcdinf<-read.csv("20190522CRCD_sampleinf.csv",stringsAsFactors = F,header = T)
genename<-read.csv("uniprot-human-gene_name.csv",stringsAsFactors = F,header = T)
first_category_name<-c("reviewedfull","reviewedsemi","isoformfull","isoformsemi","dphl1","library_free")
filename<-paste0(first_category_name,"_crcd_diann_prot_90.csv")
for(i in filename ){
  protmat<-read.csv(i,header = T)
  sum(is.na(protmat))/(ncol(protmat)*nrow(protmat))
protmat[is.na(protmat)]<-min(protmat[,2:ncol(protmat)],na.rm = T)
  # protmat<-data.frame(t(readin))
  # protmat<-readin[,apply(!is.na(readin),2,sum)>nrow(readin)*0.2]
  for (j in 2:ncol(protmat)) {
    indsel<-which(genename$Entry==colnames(protmat)[j])
    if(length(indsel)>0){
      colnames(protmat)[j]<-paste(colnames(protmat)[j],genename$Gene[indsel],sep = "_")
    }
  }
  
  protmat2<-merge(protmat,crcdinf,by.x = "X",by.y = "batch_number",all.x = T)
  protmat2<-protmat2[,2:(ncol(protmat2)-7)]
  delmat<-protmat2[duplicated(protmat2$PatientId),]
  #mean of proteins
  MeanNA<-data.frame(matrix(NA,nrow=length(unique(protmat2$PatientId)),ncol=(ncol(protmat2)-1)))
  colnames(MeanNA)<-colnames(protmat2)[1:(ncol(protmat2)-1)]
  row.names(MeanNA)<-unique(protmat2$PatientId)
  for (j in 1:nrow(MeanNA)) {
    indsel<-which(protmat2$PatientId==row.names(MeanNA)[j])
    if(length(indsel)>1){
      MeanNA[j,]<-log2(apply(2^protmat2[indsel,1:(ncol(protmat2)-1)],2,mean,na.rm=T)) 
      
    }
    else
      MeanNA[j,]<-protmat2[indsel,1:(ncol(protmat2)-1)]
  }
  MeanNA$patID<-row.names(MeanNA)
  crcdinf1<-crcdinf[!duplicated(crcdinf$PatientId),]
  unimat<-merge(MeanNA,crcdinf1,by.x = "patID",by.y = "PatientId",all.x = T)
  unimat$pat<-sapply(strsplit(unimat$patID,""),function(e){e[1]})
  protmat1<-unimat[,c(2:(ncol(unimat)-9),ncol(unimat))]
  row.names(protmat1)<-unimat$patID
  readCont<-protmat1[grepl("N",protmat1$pat),]
  readCont<-readCont[,-ncol(readCont)]
  readPat<-protmat1[grepl("P",protmat1$pat),]
  readPat<-readPat[,-ncol(readPat)]
  
  FC<-vector();
  pvalue<-vector();
  for(x in 1:ncol(readCont)){
    if(sum(!is.na(readCont[,x]))>1&sum(!is.na(readPat[,x]))>1){
      
      FC[x]<-mean(2^readPat[,x],na.rm=T)/mean(2^readCont[,x],na.rm = T)
      pvalue[x]<-t.test(readPat[,x],readCont[,x], paxred = F,  var.equal = F)$p.value
    }
  }
  pvalueAd<-p.adjust(pvalue,method = "BH",n=sum(!is.na(pvalue)))
  FC_PAd<-cbind(FC,pvalueAd)
  row.names(FC_PAd)<-colnames(readCont)
  write.csv(FC_PAd,paste0("FC_PAd_",i))
  hist(pvalueAd)
  hist(log2(FC),breaks = 1000)
  dev.off()
  library(ggplot2)
  pdf(paste0(i,"_vocano.pdf"),width =8,height = 6)
  plot(log2(FC),-log10(pvalueAd),
       pch=19, col="gray",
       xlim=c(-10,10),
       main="patient vs control"
       ,xlab="log2 FC",ylab="-log10 pvalueAd",
       cex.axis = 1.5,cex.lab=1.5)
  up <- which(pvalueAd < 0.01 & FC > 2)
  down <- which(pvalueAd < 0.01 & FC <  0.5)
  points(log2(FC[up]), -log10(pvalueAd[up]), col=1, bg = "red", pch=21, cex=1)
  points(log2(FC[down]), -log10(pvalueAd[down]), col = "blue", pch=19,cex=1)
  abline(h=-log10(0.01),v=c(-log2(2),log2(2)),lty=2,lwd=1)
  dev.off()
  protUP<-FC_PAd[up,]
  protdown<-FC_PAd[down,]
  write.csv(protUP,paste0("protUP_",i))
  write.csv(protdown,paste0("protDown_",i))
  upprot<-cbind(protmat1[,colnames(protmat1)%in%row.names(protUP)],protmat1$pat)
  downprot<-cbind(protmat1[,colnames(protmat1)%in%row.names(protdown)],protmat1$pat)
  write.csv(upprot,paste0("protUP_matrix_",i))
  write.csv(downprot,paste0("protDown_matrix_",i))
  }

#4lib dia overlap differential proteins
#down
first_category_name<-c("reviewedfull","reviewedsemi","isoformfull","isoformsemi","dphl1","library_free")
filename<-paste0("protDown_",first_category_name,"_crcd_diann_prot_90.csv")
c<-list()
m<-1
for (i in filename) {
  indsel<-read.csv(i,header = T)
  c[[m]]<-indsel$X
  m<-m+1
}
uniP<-unique(c(unlist(c[[1]]),unlist(c[[2]]),unlist(c[[3]]),unlist(c[[4]]),unlist(c[[5]]),unlist(c[[6]])))
unimat<-data.frame(matrix(0,nrow = length(uniP),ncol = 6))
row.names(unimat)<-uniP
colnames(unimat)<-first_category_name
for (j in 1:ncol(unimat)) {
  for (i in 1:nrow(unimat)) {
    pr_this <- row.names(unimat)[i]
    if(pr_this%in%c[[j]]){
      unimat[i,j]<-1
    } 
}

}

###heatmap

#cycle heatmap
library(circlize)
library(grid)
library(ComplexHeatmap)
col_fun1 = colorRamp2(c(0,1), c( "#8F82BC","#F19E9E"))

pdf("20220103_DPHL2_diann_crcd_375Down_prot_circheatmap.pdf",width=8,height=8)
circos.heatmap(unimat,col =col_fun1,cluster =T,rownames.side = "outside",dend.side = "inside",bg.border = T)
circos.clear()
# lgd = Legend(title = "heatmap", col_fun = col_fun1)
# grid.draw(lgd)
# circos.clear()
dev.off()


#upset
library(UpSetR)

listinput <- list(reviewedFull=c[[1]],
                  reviwedSemi=c[[2]],
                  isoformFull=c[[3]],
                  isoformSemi=c[[4]],
                  DPHL1=c[[5]],
                library_free=c[[6]])
pdf('20220103_DPHL2_diann_crcd_down_prot.pdf',height = 8,width = 10)
upset(fromList(listinput),nsets = 6, order.by = "freq",text.scale =1.5)
dev.off()

###########################
# 6 lib compair
temp<-read.csv("6lib_summary_diann_result.csv",header = T)

library(RColorBrewer)
library(randomcoloR)
library(reshape2)
library(ggplot2)
library(plyr)
protsum<-melt(temp)

a<-ggplot(protsum,aes(x=X,y=value,fill=variable))+geom_bar(stat="identity",position = 'stack',colour="black")+
  geom_text(aes(label=value),size=4,position = position_stack())+guides(fill=guide_legend(reverse=F))+
  scale_fill_manual(values=c("burlywood1","gold2","gray74"))+theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(),
                                                  panel.background = element_blank(),axis.line = element_line(colour = "black"))

ggsave("46lib_prot_up_down.pdf", plot=a,device = NULL, width = 12, height = 8)




#venn
library(gplots)
library(VennDiagram)
venn1<-venn.diagram(listinput,resolution = 300, imagetype = "png",  cat.fontface=1,fontfamily=1,cat.cex=2,cex=2,
                    main="Protein down overlap of 6 librarys",main.cex = 2, main.fontface = 2, main.fontfamily = 3,filename = NULL)

# pdf(file="20220103_dphl_protein_down_overlap_lib.pdf")
grid.draw(venn1)#
dev.off()



unimat$type<-paste(unimat$reviewedfull,unimat$reviewedsemi,unimat$isoformfull,unimat$isoformsemi,unimat$dphl1,sep = "_")
unimat1<-unimat[unimat$type=="1_1_1_1_1_1",]
write.csv(unimat1,"20220103_DPHL2_6lib_overlap_prot_down.csv",row.names = T)



#UP
filename<-paste0("protUP_",first_category_name,"_crcd_diann_prot_90.csv")
c<-list()
m<-1
for (i in filename) {
  indsel<-read.csv(i,header = T)
  c[[m]]<-indsel$X
  m<-m+1
}
uniP<-unique(c(unlist(c[[1]]),unlist(c[[2]]),unlist(c[[3]]),unlist(c[[4]]),unlist(c[[5]]),unlist(c[[6]])))
unimat<-data.frame(matrix(0,nrow = length(uniP),ncol = 6))
row.names(unimat)<-uniP
colnames(unimat)<-first_category_name
for (j in 1:ncol(unimat)) {
  for (i in 1:nrow(unimat)) {
    pr_this <- row.names(unimat)[i]
    if(pr_this%in%c[[j]]){
      unimat[i,j]<-1
    } 
  }
  
}

###heatmap

#heatmap
library(circlize)
library(grid)
library(ComplexHeatmap)
col_fun1 = colorRamp2(c(0,1), c(  "#8F82BC","#F19E9E"))

pdf("20220103_DPHL2_diann_crcd_763UP_prot_circheatmap.pdf",width=8,height=8)
circos.heatmap(unimat,col =col_fun1,cluster =T,rownames.side = "outside",dend.side = "inside",bg.border = T)
circos.clear()
# lgd = Legend(title = "heatmap", col_fun = col_fun1)
# grid.draw(lgd)
# circos.clear()
dev.off()



#upset
library(UpSetR)

listinput <- list(reviewedFull=c[[1]],
                  reviwedSemi=c[[2]],
                  isoformFull=c[[3]],
                  isoformSemi=c[[4]],
                  DPHL1=c[[5]],
                  library_free=c[[6]])
pdf('20220103_DPHL2_diann_crcd_up_prot.pdf',height = 8,width = 10)
upset(fromList(listinput),nsets = 6, order.by = "freq",text.scale =1.5)
dev.off()






#venn
library(gplots)
library(VennDiagram)
venn1<-venn.diagram(listinput,resolution = 300, imagetype = "png",  cat.fontface=1,fontfamily=1,cat.cex=2,cex=2,
                    main="Protein up overlap of 6 librarys",main.cex = 2, main.fontface = 2, main.fontfamily = 3,filename =NULL)

grid.draw(venn1)#
dev.off()

unimat$type<-paste(unimat$reviewedfull,unimat$reviewedsemi,unimat$isoformfull,unimat$isoformsemi,unimat$dphl1,sep = "_")
unimat1<-unimat[unimat$type=="1_1_1_1_1_1",]
write.csv(unimat1,"20220103_DPHL2_6lib_overlap_prot_up.csv",row.names = T)



sum<-read.csv("20220103_6lib_analysis_result_summary.csv",header = T)

library(ggplot2)
library(reshape2)
library(dplyr)
data_m <- melt(sum[,1:3], id.vars=c("X"))
p <- ggplot(data_m, aes(x=X, y=value))+ geom_bar(stat="identity", position="dodge", aes(fill=variable))+geom_text(aes(label=value))
p

ggsave("20220103_DPHL2_6lib_result_summary.pdf", plot=p,device = NULL, width = 12, height = 8)




library(Biostrings)
fasta<-readBStringSet("2021-05-06-nodecoys-UP000005640_reviewed_human_20200401_no_irt_fragpipe.fasta",
                      format="fasta", nrec=-1L, skip=0L, seek.first.rec=FALSE, use.names=TRUE)

a<-data.frame(fasta@ranges@NAMES)

b<-data.frame(sapply(a$fasta.ranges.NAMES,function(e){unlist(strsplit(e," "))[1]}))#

write.csv(b,"reviewed_fasta_prot_list.csv",row.names = F)


fasta<-readBStringSet("2021-05-06-nodecoys-UP000005640_isoform_reviewed_human_20200401_no_irt_fragpipe.fasta",
                      format="fasta", nrec=-1L, skip=0L, seek.first.rec=FALSE, use.names=TRUE)

a<-data.frame(fasta@ranges@NAMES)

b<-data.frame(sapply(a$fasta.ranges.NAMES,function(e){unlist(strsplit(e," "))[1]}))#f

write.csv(b,"isoform_fasta_prot_list.csv",row.names = F)





########################################################################
RF<-read.csv("RF_CRCD_diff_prot.csv",header = T)
RS<-read.csv("RS_CRCD_diff_prot.csv",header = T)
IF<-read.csv("IF_CRCD_diff_prot.csv",header = T)
IS<-read.csv("IS_CRCD_diff_prot.csv",header = T)
DPHL1<-read.csv("DPHL1_CRCD_diff_prot.csv",header = T)
Libfree<-read.csv("Libfree_CRCD_diff_prot.csv",header = T)

#upset
library(UpSetR)

listinput <- list(reviewedFull=RF$RF,
                  reviwedSemi=RS$RS,
                  isoformFull=IF$IF,
                  isoformSemi=IS$IS,
                  DPHL1=DPHL1$DPHL1,
                  library_free=Libfree$lib_free)
pdf('20220727_DPHL2_diann_crcd_diff_prot.pdf',height = 8,width = 10)
upset(fromList(listinput),nsets = 6, order.by = "freq",text.scale =1.5)
dev.off()



#######################################################################
#fig 4A
#missing value
library(ggplot2)
missing<-read.csv("dphlv2_crcd_6lib_missing_value.csv",header = T)
missing$X<-factor(missing$X,levels=c("Libfree","DPHL1","IS","IF","RS","RF"))
fig<-ggplot(missing,aes(x=X,y=missing))+geom_bar(aes(y=missing,fill=X),stat = "identity",width=0.6)+coord_flip()+
  ylim(0,0.5)+theme_minimal()+theme(legend.position = "none")
ggsave("dphl2_crcd_missing_value.pdf",fig,width = 8,height = 8)











