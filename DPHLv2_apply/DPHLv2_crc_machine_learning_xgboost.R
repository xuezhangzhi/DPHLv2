#####################################################################
#Machine learning using 4lib overlap differential proteins

file_list_up<-list.files(path = "E:/Qsync/work/DPHL2/0415_DPHL2_apply_new/crcd/analysis4_90/analysis_90_no_gene_name/4lib_diff_prot/",pattern = "protUP")
prot_name_up<-c()
for (i in file_list_up) {
  indsel<-read.csv(i,header = T)
  prot_name_up<-c(prot_name_up,indsel$X)
  
}       
prot_num_up<-data.frame(table(prot_name_up))

prot_ovlap_up<-prot_num_up[prot_num_up$Freq==4,]


file_list_down<-list.files(path = "E:/Qsync/work/DPHL2/0415_DPHL2_apply_new/crcd/analysis4_90/analysis_90_no_gene_name/4lib_diff_prot/",pattern = "protDown_")
prot_name_down<-c()
for (i in file_list_down) {
  indsel<-read.csv(i,header = T)
  prot_name_down<-c(prot_name_down,indsel$X)
  
}       
prot_num_down<-data.frame(table(prot_name_down))

prot_ovlap_down<-prot_num_down[prot_num_down$Freq==4,]

###########################################################################
#reviewedfull matrix as input
refull<-read.csv("reviewedfull_uni_matrix.csv",header = T)

protmat_up<-refull[,colnames(refull)%in%prot_ovlap_up$prot_name_up]
protmat_down<-refull[,colnames(refull)%in%prot_ovlap_down$prot_name_down]
protmat<-cbind(protmat_up,protmat_down)
row.names(protmat)<-refull$X
protmat$label<-refull$pat
protmat$label<-gsub("N",0,protmat$label)
protmat$label<-gsub("P",1,protmat$label)
# protmat$label<-as.factor(protmat$label)
set.seed(123)
nn<-sample(241, 200, replace = FALSE, prob = NULL)
protmat<-protmat[,c(1427,1:1426)]
train_mat<-protmat[nn,]
test_mat<-protmat[-nn,]



library("xgboost")
library("Matrix")


#build model
accu<-c()
for (eta_s in seq(0.2,0.3,0.05)) {
  for (ss in seq(0.8,1,0.05)) {
    for (gam in seq(0.05,0.2,0.05)) {
      for (seeds in round(runif(10,1,1000))) {
        set.seed(seeds)
        nnn<-sample(200, 160, replace = FALSE, prob = NULL)
        protmat<-protmat[,c(1427,1:1426)]
        trains<-train_mat[nnn,]
        tests<-train_mat[-nnn,]
        trains_matrix <- sparse.model.matrix(label ~ ., data = trains)
        trains_label <- as.numeric(trains$label)
        trains_fin <- list(data=trains_matrix,label=trains_label) 
        dtrains <- xgb.DMatrix(data = trains_fin$data, label = trains_fin$label) 
       
      
          xgb <- xgboost(data = dtrains, eta=eta_s,objective='binary:logistic',nround=50,subsample=ss,
                         gamma = gam,scale_pos_weight=1,set.seed(1234),eval_metric='logloss')
          
          
          #Ranking of importance
          importance <- xgb.importance(trains_matrix@Dimnames[[2]], model = xgb)  
          # head(importance)
          # xgb.plot.importance(importance)
          #############################################################
          if(nrow(importance)>5 ){
            
          
          #feature select
          for (fea_selt in seq(5,20,1)) {
            feature_choose<-importance$Feature[1:fea_selt]
            # feature_choose<-importance$Feature[importance$Gain>fea_selt]
            trains2<-trains[,colnames(trains)%in%feature_choose]
            trains2$label<-trains$label
            tests2<-tests[,colnames(tests)%in%feature_choose]
            tests2$label<-tests$label
            
            
            train_matrix2 <- sparse.model.matrix(label ~ ., data = trains2)
            test_matrix2 <- sparse.model.matrix(label ~ ., data = tests2)
            
            train_label2 <- as.numeric(trains2$label)
            test_label2 <-  as.numeric(tests2$label)
            
            train_fin2 <- list(data=train_matrix2,label=train_label2) 
            test_fin2 <- list(data=test_matrix2,label=test_label2) 
           
            dtrain2 <- xgb.DMatrix(data = train_fin2$data, label = train_fin2$label) 
            dtest2 <- xgb.DMatrix(data = test_fin2$data, label = test_fin2$label)
            
            
            
            
            
            xgb2 <- xgboost(data = dtrain2, eta=eta_s,objective='binary:logistic', nround=50,subsample=ss,
                            gamma = gam,scale_pos_weight=1,set.seed(1234),eval_metric='logloss')
            
            
            
            
            ########################################
            #test
            
            
            pre_xgb = round(predict(xgb2,newdata = dtest2))
            pre_xgb_tes<-predict(xgb2,newdata = dtest2)
            aaa<-data.frame(table(test_label2,pre_xgb,dnn=c("true","pre")))
            acc_test<-sum(pre_xgb==test_label2)/length(test_label2)
            
            #ROC
            library(pROC)
            xgboost_roc <- roc(test_label2,as.numeric(pre_xgb_tes))
            
           
            auc_test<-xgboost_roc[["auc"]]
            
            
            
            accu <- rbind(accu,c(eta_s,ss,gam,seeds,fea_selt,acc_test,auc_test))
            
            }
          }
           
              
          
        }
        
      }
      
    }
    
  }
  

accu<-data.frame(accu)
colnames(accu)<-c("eta_s","ss","gam","seeds","fea_selt","acc_test","auc_test")

write.csv(accu,"20230407_dphl2_4lib_overlap_diffprot_RFmat_xgboost_search_cosvali.csv",row.names = F)

accu3<-accu[accu$acc_test==1&accu$auc_test==1,]

train_matrix <- sparse.model.matrix(label ~ ., data = train_mat)

train_label <- as.numeric(train_mat$label)

train_fin <- list(data=train_matrix,label=train_label)

dtrain <- xgb.DMatrix(data = train_fin$data, label = train_fin$label)

accu4<-feature_total<-c()
for (i in 1:nrow(accu3)) {
  
  xgb <- xgboost(data = dtrain, eta=accu3$eta_s[i],objective='binary:logistic',nround=50,subsample=accu3$ss[i],
                 gamma = accu3$gam[i],scale_pos_weight=1,set.seed(1234),eval_metric='logloss')
  
  
  #Ranking of importance 
  importance <- xgb.importance(train_matrix@Dimnames[[2]], model = xgb)  
  # head(importance)
  # xgb.plot.importance(importance)
  #############################################################
      feature_choose<-importance$Feature[1:accu3$fea_selt[i]]
      feature_total<-qpcR:::cbind.na(feature_total,feature_choose)
        
      # feature_choose<-importance$Feature[importance$Gain>fea_selt]
      trains2<-train_mat[,colnames(train_mat)%in%feature_choose]
      trains2$label<-train_mat$label
      tests2<-test_mat[,colnames(test_mat)%in%feature_choose]
      tests2$label<-test_mat$label
      
      
      train_matrix2 <- sparse.model.matrix(label ~ ., data = trains2)
      test_matrix2 <- sparse.model.matrix(label ~ ., data = tests2)
      
      train_label2 <- as.numeric(trains2$label)
      test_label2 <-  as.numeric(tests2$label)
      
      train_fin2 <- list(data=train_matrix2,label=train_label2) 
      test_fin2 <- list(data=test_matrix2,label=test_label2) 
      
      dtrain2 <- xgb.DMatrix(data = train_fin2$data, label = train_fin2$label) 
      dtest2 <- xgb.DMatrix(data = test_fin2$data, label = test_fin2$label)
      
      
      
      
      
      xgb2 <- xgboost(data = dtrain2, eta=accu3$eta_s[i],objective='binary:logistic',nround=50,subsample=accu3$ss[i],
                      gamma = accu3$gam[i],scale_pos_weight=1,set.seed(1234),eval_metric='logloss')
      
      
      ########################################
      #train
      
      
      pre_xgb_train = round(predict(xgb2,newdata = dtrain2))
      pre_xgb_tra<-predict(xgb2,newdata = dtrain2)
      aaa_train<-data.frame(table(train_label2,pre_xgb_train,dnn=c("true","pre")))
      acc_train<-sum(pre_xgb_train==train_label2)/length(train_label2)
      
      #ROC
      library(pROC)
      xgboost_roc_train <- roc(train_label2,as.numeric(pre_xgb_tra))
      
      
      auc_train<-xgboost_roc_train[["auc"]]
      
      
      #test
      
      
      pre_xgb = round(predict(xgb2,newdata = dtest2))
      pre_xgb_tes=predict(xgb2,newdata = dtest2)
      aaa<-data.frame(table(test_label2,pre_xgb,dnn=c("true","pre")))
      acc_test<-sum(pre_xgb==test_label2)/length(test_label2)
      
      #ROC
      library(pROC)
      xgboost_roc <- roc(test_label2,as.numeric(pre_xgb_tes))
      
    
      auc_test<-xgboost_roc[["auc"]]
      
      
      
      accu4 <- rbind(accu4,c(accu3$eta_s[i],accu3$ss[i],accu3$gam[i],accu3$fea_selt[i],acc_train,auc_train,acc_test,auc_test))
      
     
      
}
accu4<-data.frame(accu4)
names(accu4)<-c("eta_s","ss","gam","fea_selt","acc_train","auc_train","acc_test","auc_test")
write.csv(accu4,"20230408_dphl2_RF_overlap_diffprot_RFmat_xgboost_model_result.csv",row.names = T)
sum(accu4$acc_test>0.9)
sum(accu4$auc_test>0.9)
# accu5<-accu4(accu4$acc_train==1&accu4$auc_test==1,)
feature_total1<-data.frame(feature_total)
feature_total1<-data.frame(unlist(feature_total1))
feature_total1<-na.omit(feature_total1)
feature_tab<-data.frame(table(feature_total1))
write.csv(feature_tab,"20230408_dphl2_RF_overlap_diffprot_RFmat_xgboost_model_featuretable.csv",row.names = F)
write.csv(feature_total,"20230408_dphl2_RF_overlap_diffprot_RFmat_xgboost_model_featuretotal.csv",row.names = F)
#########################
#the select model
#32
i=609
xgb <- xgboost(data = dtrain, eta=accu3$eta_s[i],objective='binary:logistic',nround=50,subsample=accu3$ss[i],
               gamma = accu3$gam[i],scale_pos_weight=1,set.seed(1234),eval_metric='logloss')


#重要重要性排序 
importance <- xgb.importance(train_matrix@Dimnames[[2]], model = xgb)  
# head(importance)
# xgb.plot.importance(importance)
#############################################################
feature_choose<-importance$Feature[1:accu3$fea_selt[i]]
write.csv(feature_choose,"20230407_DPHL2_RF_xgboost_bast_model609_feature.csv",row.names = F)
 # feature_choose<-importance$Feature[importance$Gain>fea_selt]
fea<-importance[1:accu3$fea_selt[i],]
pdf("20230407_dphl2_RF_xgboost_overlap_bast_model609_feature.pdf",width = 8,height = 6)
plot(fea$Gain,xgap.axis = 1)
dev.off()

trains2<-train_mat[,feature_choose]
trains2$label<-train_mat$label
tests2<-test_mat[,feature_choose]
tests2$label<-test_mat$label


train_matrix2 <- sparse.model.matrix(label ~ ., data = trains2)
test_matrix2 <- sparse.model.matrix(label ~ ., data = tests2)

train_label2 <- as.numeric(trains2$label)
test_label2 <-  as.numeric(tests2$label)

train_fin2 <- list(data=train_matrix2,label=train_label2) 
test_fin2 <- list(data=test_matrix2,label=test_label2) 

dtrain2 <- xgb.DMatrix(data = train_fin2$data, label = train_fin2$label) 
dtest2 <- xgb.DMatrix(data = test_fin2$data, label = test_fin2$label)





xgb2 <- xgboost(data = dtrain2, eta=accu3$eta_s[i],objective='binary:logistic',nround=50,subsample=accu3$ss[i],
                gamma = accu3$gam[i],scale_pos_weight=1,set.seed(1234),eval_metric='logloss')



########################################
#train


pre_xgb_train = round(predict(xgb2,newdata = dtrain2))
pre_xgb_tra<-predict(xgb2,newdata = dtrain2)
aaa_train<-data.frame(table(train_label2,pre_xgb_train,dnn=c("true","pre")))
acc_train<-sum(pre_xgb_train==train_label2)/length(train_label2)

#ROC
library(pROC)
xgboost_roc_train <- roc(train_label2,as.numeric(pre_xgb_tra))
pdf("20230407_DPHL2_RF_xgboost_bastmosel609_train_ROC.pdf",width = 8,height = 8)
plot.roc(xgboost_roc_train,col = "blue3",ylim=c(0,1), print.auc=TRUE, print.thres="best",
         main="Train ROC",legacy.axes = TRUE,print.auc.cex=1.2)
dev.off()
auc_train<-xgboost_roc_train[["auc"]]

####################################################################################
#train predict 

predicted1<-cbind(trains2,pre_xgb_tra)
predicted1$mean<-apply(predicted1[,1:14],1,mean)
# row.names(predicted1)<-predicted1$batch
pdf("20221219_dphl2_RF_xgboost_bastmosel609_train_acc.pdf",width = 5.5,height = 5)
plot(predicted1$pre_xgb_tra,predicted1$mean,xlim = c(0,1))
points(predicted1$pre_xgb_tra[predicted1$label==0],predicted1$mean[predicted1$label==0], col="red", pch=19, cex=1)
points(predicted1$pre_xgb_tra[predicted1$label==1],predicted1$mean[predicted1$label==1], col = "blue", pch=19,cex=1)
abline(v=0.5,lty=2,lwd=1)
dev.off()


#test


pre_xgb = round(predict(xgb2,newdata = dtest2))
pre_xgb_tes = predict(xgb2,newdata = dtest2)
aaa<-data.frame(table(test_label2,pre_xgb,dnn=c("true","pre")))
acc_test<-sum(pre_xgb==test_label2)/length(test_label2)

#ROC
library(pROC)
xgboost_roc <- roc(test_label2,as.numeric(pre_xgb_tes))

pdf("20230407_DPHL2_RF_xgboost_bastmosel609_test_ROC.pdf",width = 8,height = 8)
plot.roc(xgboost_roc,col = "blue3",ylim=c(0,1), print.auc=TRUE, print.thres="best",
         main="Train ROC",legacy.axes = TRUE,print.auc.cex=1.2)
dev.off()
auc_test<-xgboost_roc[["auc"]]




####################################################################################
#test predict 

predicted1<-cbind(tests2,pre_xgb_tes)
predicted1$mean<-apply(predicted1[,1:14],1,mean)
# row.names(predicted1)<-predicted1$batch
pdf("20230407_DPHL2_RF_xgboost_bastmosel609_test_acc.pdf",width = 5.5,height = 5)
plot(predicted1$pre_xgb_tes,predicted1$mean,xlim = c(0,1))
points(predicted1$pre_xgb_tes[predicted1$label==0],predicted1$mean[predicted1$label==0], col="red", pch=19, cex=1)
points(predicted1$pre_xgb_tes[predicted1$label==1],predicted1$mean[predicted1$label==1], col = "blue", pch=19,cex=1)
abline(v=0.5,lty=2,lwd=1)
dev.off()

##############################################
################################################################
#predict using RS IF IS sets
#resemi
resemi<-read.csv("reviewedsemi_uni_matrix.csv",header = T)

# protmat_rs<-resemi[,colnames(resemi)%in%feature_choose]
protmat_rs<-resemi[,feature_choose]
row.names(protmat_rs)<-resemi$X
protmat_rs$label<-resemi$pat
protmat_rs$label<-gsub("N",0,protmat_rs$label)
protmat_rs$label<-gsub("P",1,protmat_rs$label)
protmat_rs$label<-as.factor(protmat_rs$label)
RS_matrix2 <- sparse.model.matrix(label ~ ., data = protmat_rs)
RS_label2 <- protmat_rs$label
RS_fin2 <- list(data=RS_matrix2,label=RS_label2) 
dRS2 <- xgb.DMatrix(data = RS_fin2$data, label = RS_fin2$label) 



pre_xgb_rs = round(predict(xgb2,newdata = dRS2))
pre_xgb_rs_1<-predict(xgb2,newdata = dRS2)
aaa_rs<-data.frame(table(RS_label2,pre_xgb_rs,dnn=c("true","pre")))
acc_rs<-sum(pre_xgb_rs==RS_label2)/length(RS_label2)

#ROC
library(pROC)
xgboost_roc_rs <- roc(RS_label2,as.numeric(pre_xgb_rs_1))
pdf("20230407_DPHL2_RF_xgboost_bastmosel609_RS_ROC.pdf",width = 8,height = 8)
plot.roc(xgboost_roc_rs,col = "blue3",ylim=c(0,1), print.auc=TRUE, print.thres="best",
         main="Train ROC",legacy.axes = TRUE,print.auc.cex=1.2)
dev.off()
auc_RS<-xgboost_roc_rs[["auc"]]

####################################################################################
#train predict 

predicted1<-cbind(protmat_rs,pre_xgb_rs_1)
predicted1$mean<-apply(predicted1[,1:14],1,mean)
# row.names(predicted1)<-predicted1$batch
pdf("20230407_dphl2_RF_xgboost_bastmosel609_RS_acc.pdf",width = 5.5,height = 5)
plot(predicted1$pre_xgb_rs_1,predicted1$mean,xlim = c(0,1))
points(predicted1$pre_xgb_rs_1[predicted1$label==0],predicted1$mean[predicted1$label==0], col="red", pch=19, cex=1)
points(predicted1$pre_xgb_rs_1[predicted1$label==1],predicted1$mean[predicted1$label==1], col = "blue", pch=19,cex=1)
abline(v=0.5,lty=2,lwd=1)
dev.off()



################################################################
#issemi
issemi<-read.csv("isoformsemi_uni_matrix.csv",header = T)

protmat_is<-issemi[,feature_choose]
row.names(protmat_is)<-issemi$X
protmat_is$label<-issemi$pat
protmat_is$label<-gsub("N",0,protmat_is$label)
protmat_is$label<-gsub("P",1,protmat_is$label)
protmat_is$label<-as.factor(protmat_is$label)

IS_matrix2 <- sparse.model.matrix(label ~ ., data = protmat_is)
IS_label2 <- protmat_is$label
IS_fin2 <- list(data=IS_matrix2,label=IS_label2) 
dIS2 <- xgb.DMatrix(data = IS_fin2$data, label = IS_fin2$label) 


pre_xgb_is = round(predict(xgb2,newdata = dIS2))
pre_xgb_is_1<-predict(xgb2,newdata = dIS2)
aaa_is<-data.frame(table(IS_label2,pre_xgb_is,dnn=c("true","pre")))
acc_is<-sum(pre_xgb_is==IS_label2)/length(IS_label2)

#ROC
library(pROC)
xgboost_roc_is <- roc(IS_label2,as.numeric(pre_xgb_is_1))
pdf("20230407_DPHL2_RF_xgboost_bastmosel609_IS_ROC.pdf",width = 8,height = 8)
plot.roc(xgboost_roc_is,col = "blue3",ylim=c(0,1), print.auc=TRUE, print.thres="best",
         main="Train ROC",legacy.axes = TRUE,print.auc.cex=1.2)
dev.off()
auc_IS<-xgboost_roc_is[["auc"]]

####################################################################################
#train predict 

predicted1<-cbind(protmat_is,pre_xgb_is_1)
predicted1$mean<-apply(predicted1[,1:14],1,mean)
# row.names(predicted1)<-predicted1$batch
pdf("20230407_dphl2_RF_xgboost_bastmosel609_IS_acc.pdf",width = 5.5,height = 5)
plot(predicted1$pre_xgb_is_1,predicted1$mean,xlim = c(0,1))
points(predicted1$pre_xgb_is_1[predicted1$label==0],predicted1$mean[predicted1$label==0], col="red", pch=19, cex=1)
points(predicted1$pre_xgb_is_1[predicted1$label==1],predicted1$mean[predicted1$label==1], col = "blue", pch=19,cex=1)
abline(v=0.5,lty=2,lwd=1)
dev.off()



################################################################
#isfull
isfull<-read.csv("isoformfull_uni_matrix.csv",header = T)

protmat_if<-isfull[,feature_choose]
row.names(protmat_if)<-isfull$X
protmat_if$label<-isfull$pat
protmat_if$label<-gsub("N",0,protmat_if$label)
protmat_if$label<-gsub("P",1,protmat_if$label)
protmat_if$label<-as.factor(protmat_if$label)

IF_matrix2 <- sparse.model.matrix(label ~ ., data = protmat_if)
IF_label2 <- protmat_if$label
IF_fin2 <- list(data=IF_matrix2,label=IF_label2) 
dIF2 <- xgb.DMatrix(data = IF_fin2$data, label = IS_fin2$label) 


pre_xgb_if = round(predict(xgb2,newdata = dIF2))
pre_xgb_if_1<-predict(xgb2,newdata = dIF2)
aaa_if<-data.frame(table(IF_label2,pre_xgb_if,dnn=c("true","pre")))
acc_if<-sum(pre_xgb_if==IF_label2)/length(IF_label2)

#ROC
library(pROC)
xgboost_roc_if <- roc(IF_label2,as.numeric(pre_xgb_if_1))
pdf("20230407_DPHL2_RF_xgboost_bastmosel609_IF_ROC.pdf",width = 8,height = 8)
plot.roc(xgboost_roc_if,col = "blue3",ylim=c(0,1), print.auc=TRUE, print.thres="best",
         main="IF ROC",legacy.axes = TRUE,print.auc.cex=1.2)
dev.off()
auc_IF<-xgboost_roc_if[["auc"]]

####################################################################################
#train predict 

predicted1<-cbind(protmat_if,pre_xgb_if_1)
predicted1$mean<-apply(predicted1[,1:14],1,mean)
# row.names(predicted1)<-predicted1$batch
pdf("20230407_dphl2_RF_xgboost_bastmosel609_IF_acc.pdf",width = 5.5,height = 5)
plot(predicted1$pre_xgb_if_1,predicted1$mean,xlim = c(0,1))
points(predicted1$pre_xgb_if_1[predicted1$label==0],predicted1$mean[predicted1$label==0], col="red", pch=19, cex=1)
points(predicted1$pre_xgb_if_1[predicted1$label==1],predicted1$mean[predicted1$label==1], col = "blue", pch=19,cex=1)
abline(v=0.5,lty=2,lwd=1)
dev.off()


##################################
#shaoyingkuan 20230406_DPHL2_IF_apply_CRCD_unimat.csv

Sh_prot<-read.csv("20230406_DPHL2_IF_apply_CRCB_unimat.csv",header = T,row.names = 1,check.names = F)


Sh_pat_inf<-read.csv("cac212240-sup-0009-tables4_2.csv",header = T)
Sh_prot$label<-Sh_pat_inf$label[match(row.names(Sh_prot),Sh_pat_inf$Patient.ID)]
Sh_prot1<-Sh_prot[!is.na(Sh_prot$label),]
a<-data.frame(apply(is.na(Sh_prot1), 2,sum))
b<-row.names(a)[a$apply.is.na.Sh_prot1...2..sum.<45]
Sh_prot1<-Sh_prot1[,b]
Sh_prot1[is.na(Sh_prot1)]<-min(Sh_prot1[1:6760],na.rm = T)
Sh_prot1[1:6760]<-log2(Sh_prot1[1:6760])

protmat_Sh<-Sh_prot1[,feature_choose]
protmat_Sh$label<-Sh_prot1$label
row.names(protmat_Sh)<-row.names(Sh_prot1)

Sh_matrix2 <- sparse.model.matrix(label ~ ., data = protmat_Sh)
Sh_label2 <- protmat_Sh$label
Sh_fin2 <- list(data=Sh_matrix2,label=Sh_label2) 
dSh2 <- xgb.DMatrix(data = Sh_fin2$data, label =Sh_fin2$label) 



pre_xgb_sh = round(predict(xgb2,newdata = dSh2))
pre_xgb_sh_1<-predict(xgb2,newdata = dSh2)
aaa_sh<-data.frame(table(Sh_label2,pre_xgb_sh,dnn=c("true","pre")))
acc_sh<-sum(pre_xgb_sh==Sh_label2)/length(Sh_label2)

#ROC
library(pROC)
xgboost_roc_sh <- roc(Sh_label2,as.numeric(pre_xgb_sh_1))
pdf("20230407_DPHL2_RF_xgboost_bastmosel609_Sh_ROC.pdf",width = 8,height = 8)
plot.roc(xgboost_roc_sh,col = "blue3",ylim=c(0,1), print.auc=TRUE, print.thres="best",
         main="Train ROC",legacy.axes = TRUE,print.auc.cex=1.2)
dev.off()
auc_sh<-xgboost_roc_sh[["auc"]]

####################################################################################
#train predict 散点图

predicted1<-cbind(protmat_Sh,pre_xgb_sh_1)
predicted1$mean<-apply(predicted1[,1:14],1,mean)
# row.names(predicted1)<-predicted1$batch
pdf("20221219_dphl2_RF_xgboost_bastmosel609_Sh_acc.pdf",width = 5.5,height = 5)
plot(predicted1$pre_xgb_sh_1,predicted1$mean,xlim = c(0,1))
points(predicted1$pre_xgb_sh_1[predicted1$label==0],predicted1$mean[predicted1$label==0], col="red", pch=19, cex=1)
points(predicted1$pre_xgb_sh_1[predicted1$label==1],predicted1$mean[predicted1$label==1], col = "blue", pch=19,cex=1)
abline(v=0.5,lty=2,lwd=1)
dev.off()



############################
#ROC
pdf("20230411_DPHLv2_ROC_6mat.pdf", width = 8, height = 6)
plot(xgboost_roc_train, print.auc = TRUE, print.thres = F, xlab = "1-Specificity", ylab = "Sensitivity",
     main = "ROC Curve of 5-fold", col = "red", 
     print.auc.x = 0.5, print.auc.y = 0.5,legacy.axes = TRUE)

plot(xgboost_roc, add = T, col = 'blue', print.auc = T, print.auc.x = 0.5, print.auc.y = 0.45,legacy.axes = TRUE)
plot(xgboost_roc_rs, add = T, col = 'pink', print.auc = T, print.auc.x = 0.5, print.auc.y = 0.45,legacy.axes = TRUE)
plot(xgboost_roc_if, add = T, col = 'green', print.auc = T, print.auc.x = 0.5, print.auc.y = 0.45,legacy.axes = TRUE)
plot(xgboost_roc_is, add = T, col = 'orange', print.auc = T, print.auc.x = 0.5, print.auc.y = 0.45,legacy.axes = TRUE)
plot(xgboost_roc_sh, add = T, col = 'purple', print.auc = T, print.auc.x = 0.5, print.auc.y = 0.45,legacy.axes = TRUE)


legend("bottomright", legend = c(paste0("Train:",auc_train),paste0("Test:",auc_test),paste0("RS:",auc_RS),paste0("IF:",auc_IF),
                                 paste0("IS:",auc_IS), paste0("Shao:",auc_sh)),
       col=c("red", "blue","pink",'green','orange',"black"), lwd=2)
dev.off()



















