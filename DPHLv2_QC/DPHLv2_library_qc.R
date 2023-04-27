options(stringsAsFactors=F)
d<-read.table("easypqp_lib_openswath.tsv",header=T,sep="\t")

charge<-unique(d$PrecursorCharge)
charge2<-d[d$PrecursorCharge==2,]
charge2<-charge2[!duplicated(charge2$ModifiedPeptideSequence),]
charge3<-d[d$PrecursorCharge==3,]
charge3<-charge3[!duplicated(charge3$ModifiedPeptideSequence),]

charge2cor<-charge2[charge2$ModifiedPeptideSequence%in%charge3$ModifiedPeptideSequence,]

charge3cor<-charge3[charge3$ModifiedPeptideSequence%in%charge2$ModifiedPeptideSequence,]

charge3cor<-charge3cor[match(charge2cor$ModifiedPeptideSequence,charge3cor$ModifiedPeptideSequence),]

corraction<-cor(charge2cor$NormalizedRetentionTime, charge3cor$NormalizedRetentionTime, use = "na.or.complete",method = c("pearson"))
z<-lm(charge2cor$NormalizedRetentionTime~charge3cor$NormalizedRetentionTime)
pdf("DPHL2_correlation of +2_+3_before_qc_1.pdf",width =8,height = 8)
plot(charge2cor$NormalizedRetentionTime, charge3cor$NormalizedRetentionTime,main = "correlation of +2/+3" ,type = "p" ,col = "#F08080", lwd = 3)
abline(z,lwd=3,col="#4F94CD")
legend(-200,350,paste0("R2=",corraction^0.5) , cex=1.5, col="black", pch=21, lty=3)
dev.off()

d2<-d[d$NormalizedRetentionTime>= -60 & d$NormalizedRetentionTime<=200,]
d2<-d2[d2$LibraryIntensity>10,]
a<-unique(d2$PrecursorCharge)
d3<-d2[!d2$PrecursorCharge==1,]
d3$match<-paste(d3$ModifiedPeptideSequence,d3$PrecursorCharge,sep="_")

d4<-d3[duplicated(d3$match),]
unimodprec<-unique(d4$match)
d5<-d3[d3$match%in%unimodprec,]
unisemi<-d5[!duplicated(d5$match),]
unipep<-unique(unisemi$ModifiedPeptideSequence)
onepep<-twopep<-thrpep<-data.frame(matrix(NA,nrow = 1,ncol=ncol(unisemi)))
colnames(onepep)<-colnames(twopep)<-colnames(thrpep)<-colnames(unisemi)
for (i in 1:length(unipep)) {
  indsel<-unisemi[which(unisemi$ModifiedPeptideSequence==unipep[i]),]
  if(nrow(indsel)<2){
    onepep<-rbind(onepep,indsel) 
  }else{
    if(nrow(indsel)==2){
      twopep<-rbind(twopep,indsel) 
    }else{
      thrpep<-rbind(thrpep,indsel)
    }
  }
}
onepep1<-onepep[-1,]
twopep1<-twopep[-1,]  
thrpep1<-thrpep[-1,]  
unitwo<-unique(twopep1$ModifiedPeptideSequence)
twogod1<-twogod2<-data.frame(matrix(NA,nrow = 1,ncol = ncol(twopep1)))
colnames(twogod1)<-colnames(twogod2)<-colnames(twopep1)
for (j in 1:length(unitwo)) {
  indsel<-twopep1[which(twopep1$ModifiedPeptideSequence==unitwo[j]),]
  diff<-indsel$NormalizedRetentionTime[1]-indsel$NormalizedRetentionTime[2]
  if(diff>-5&diff<5){
    twogod1<-rbind(twogod1,indsel)
  }else{
    pep1<-d5[d5$match%in%indsel$match[1],]
    pep2<-d5[d5$match%in%indsel$match[2],]
    indsel1<-indsel[which.max(c(sum(pep1$LibraryIntensity),sum(pep2$LibraryIntensity))),]
    twogod2<-rbind(twogod2,indsel1)
  }
  
}
twogod3<-twogod1[-1,]
twogod4<-twogod2[-1,]
twogod<-rbind(twogod3,twogod4)
unithr<-unique(thrpep1$ModifiedPeptideSequence)
thrgod<-thrgod1<-data.frame(matrix(NA,nrow = 1,ncol = ncol(thrpep1)))
colnames(thrgod)<-colnames(thrgod1)<-colnames(thrpep1)
for (k in 1:length(unithr)) {
  indsel<-thrpep1[which(thrpep1$ModifiedPeptideSequence==unithr[k]),]
  diff<-mean(indsel$NormalizedRetentionTime)
  indsel1<-indsel[which(indsel$NormalizedRetentionTime-diff>-5&indsel$NormalizedRetentionTime-diff<5),]
  thrgod<-rbind(thrgod,indsel1)
  if(nrow(indsel1)<=0){
    indsel2<-indsel[which(indsel$NormalizedRetentionTime-median(indsel$NormalizedRetentionTime)>
                            -5&indsel$NormalizedRetentionTime-median(indsel$NormalizedRetentionTime)<5),]
    
    thrgod1<-rbind(thrgod1,indsel2)
  }
}
thrgod2<-thrgod[-1,]
thrgod3<-thrgod1[-1,]

usepep<-rbind(onepep1,twogod,thrgod2,thrgod3)
reviewedpep<-d5[d5$match%in%usepep$match,]
d6<-reviewedpep[,c(-12,-13)]
a<-unique(d6$PrecursorCharge)
prec2<-d6[d6$PrecursorCharge==2,]
prec2<-prec2[!duplicated(prec2$ModifiedPeptideSequence),]
prec3<-d6[d6$PrecursorCharge==3,]
prec3<-prec3[!duplicated(prec3$ModifiedPeptideSequence),]

prec2cor<-prec2[prec2$ModifiedPeptideSequence%in%prec3$ModifiedPeptideSequence,]

prec3cor<-prec3[prec3$ModifiedPeptideSequence%in%prec2$ModifiedPeptideSequence,]

prec3cor<-prec3cor[match(prec2cor$ModifiedPeptideSequence,prec3cor$ModifiedPeptideSequence),]

corraction<-cor(prec2cor$NormalizedRetentionTime, prec3cor$NormalizedRetentionTime, use = "na.or.complete",method = c("pearson"))
z<-lm(prec2cor$NormalizedRetentionTime~prec3cor$NormalizedRetentionTime)
pdf("DPHL2_correlation of +2_+3.pdf",width =8,height = 8)
plot(prec2cor$NormalizedRetentionTime, prec3cor$NormalizedRetentionTime,main = "correlation of +2/+3" ,type = "p" ,col = "#F08080", lwd = 3)
abline(z,lwd=2,col="#4F94CD")
legend(-55,150,paste0("R2=",corraction^0.5) , cex=1.2, col="black", pch=21, lty=3)
dev.off()
write.table(d6,"easypqp_lib_openswathQC.tsv",sep="\t",row.names=F,col.names=T,na="",quote=F)

#open swath correlation of +2/+3####
options(stringsAsFactors=F)
revfull<-read.table("reviewedFull_easypqp_lib_openswath_opt_decoys.tsv",header=T,sep="\t")

d<-revfull[revfull$Decoy==0,]
charge<-unique(d$PrecursorCharge)
charge2<-d[d$PrecursorCharge==2,]
charge2<-charge2[!duplicated(charge2$ModifiedPeptideSequence),]
charge3<-d[d$PrecursorCharge==3,]
charge3<-charge3[!duplicated(charge3$ModifiedPeptideSequence),]

charge2cor<-charge2[charge2$ModifiedPeptideSequence%in%charge3$ModifiedPeptideSequence,]

charge3cor<-charge3[charge3$ModifiedPeptideSequence%in%charge2$ModifiedPeptideSequence,]

charge3cor<-charge3cor[match(charge2cor$ModifiedPeptideSequence,charge3cor$ModifiedPeptideSequence),]

corraction<-cor(charge2cor$NormalizedRetentionTime, charge3cor$NormalizedRetentionTime, use = "na.or.complete",method = c("pearson"))
z<-lm(charge2cor$NormalizedRetentionTime~charge3cor$NormalizedRetentionTime)
pdf("DPHL2_correlation of +2_+3_finalllib.pdf",width =8,height = 8)
plot(charge2cor$NormalizedRetentionTime, charge3cor$NormalizedRetentionTime,main = "correlation of +2/+3" ,type = "p" ,col = "#F08080", lwd = 3)
abline(z,lwd=2,col="#4F94CD")
legend(-50,150,paste0("R2=",corraction^0.5) , cex=1.2, col="black", pch=21, lty=3)
dev.off()



