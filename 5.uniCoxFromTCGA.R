library(glmnet)
library(clusterProfiler)
library(org.Hs.eg.db)
library(limma)
library(survival)
setwd("...\\findTCGA")

DEG1<-read.table("DEG_R_D_10X.txt",sep = "\t",header = T)
cyy<-function(x){
  up<-rownames(x)[which(x[,"p_val_adj"]<0.05 & x[,"avg_log2FC"]>0.5)]
  down<-rownames(x)[which(x[,"p_val_adj"]<0.05 & x[,"avg_log2FC"]<(-0.5))]
  return(list(up,down))
}

DEG11<-cyy(DEG1)

DEG2<-read.table("DEG_R_D_SS2.txt",sep = "\t",header = T)

DEG22<-cyy(DEG2)

inter_up<-intersect(DEG11[[1]],DEG22[[1]])#183
inter_down<-intersect(DEG11[[2]],DEG22[[2]])#84
gene<-c(inter_up,inter_down)

ov<-read.table("TCGA-OV.htseq_fpkm.tsv",sep = "\t",header = T,check.names = F)
ov_clin<-read.table("ov_clin.txt",sep = "\t",header = T)
ov[,1]<-substr(ov[,1],1,15)

s2e <- bitr(ov[,1], 
            fromType = "ENSEMBL",
            toType = "SYMBOL",
            OrgDb = org.Hs.eg.db)
colnames(ov)[1]<-"ENSEMBL"
ov1<-merge(s2e,ov,by="ENSEMBL")
ov1<-ov1[,-1]
ov_exp<-avereps(ov1[,-1],ID=ov1$SYMBOL)
ov_exp<-as.data.frame(ov_exp)
#write.table(ov_exp,"...\\ov_exp_all.txt",sep = "\t",quote = F)
#write.table(ov_exp1,"ov_exp.txt",sep = "\t",quote = F)
ov_exp1<-ov_exp[which(rownames(ov_exp)%in%gene),]
colnames(ov_exp1)<-substr(colnames(ov_exp1),1,15)

rm(list = ls())

getPrintFile<-function(coxph.fit,gene,fileName){
  result<-summary(coxph.fit);txt<-paste(fileName,"_cox.txt",sep="")  # print(result)
  header=t(c("geneName","coef","exp(coef)[HR]","se(coef)","Z[coef/se]","Pr(>|Z|)[pval]","lower.95","upper.95","Likelihood ratio test","Wald test","Log rank test"))
  count<-rownames(result[[7]])
  if(length(count)==1){
    ##单因素输�?
    aa<-result$coefficients
    bb<-result$conf.int
    temp<-t(c(gene,aa,bb[c(3,4)],result[[9]][3],result[[12]][3],result[[10]][3]))
  }else{
    ##多因素输�?
    aa<-result$coefficients
    bb<-result$conf.int
    n<-dim(aa)[1]
    temp<-data.frame(gene[1:n],aa,bb[,c(3,4)],rep(result[[9]][3],times=n),rep(result[[12]][3],times=n),rep(result[[10]][3],times=n))
  }
  if(file.exists(txt)==FALSE){
    write.table(header,txt,col.names=F,row.names=F,sep="\t",quote=F)
  }
  write.table(temp,txt,col.names=F,row.names=F,sep="\t",quote=F,append=TRUE)
}

setwd("F:\\survival\\findTCGA")
###246
ov_clin<-read.table("ov_clin.txt",sep = "\t",header = T,row.names = 1)
rownames(ov_clin)<-paste(rownames(ov_clin),"-01",sep = "")

ov_clin<-na.omit(ov_clin)
ov_exp1<-t(read.table("ov_exp.txt",sep = "\t",header = T,row.names = 1,check.names = F))

aa<-intersect(rownames(ov_exp1),rownames(ov_clin))

ov_clin<-ov_clin[which(rownames(ov_clin)%in%aa),];
ov_clin<-ov_clin[order(ov_clin[,2],decreasing = T),]

library(survival)
library(survminer)
library("ggpubr")
m=1

for (i in 1:10000) {
  ov_clin_test<-ov_clin[c(sample(c(0:370),185)),]
  ov_clin_test<-ov_clin_test[,c(2,1)]
  colnames(ov_clin_test)<-c("Survival","Events")
  ov_exp11<-ov_exp1[rownames(ov_clin_test),]
  #jsurv_ss <- survival::Surv( ov_clin[,"OS.Time"], ov_clin[,"OS"])
  
  ov<-cbind(ov_clin_test,ov_exp11)
  #ov[,3:101]
  BaSurv<-Surv(ov[,1],ov[,2])
  
  result_unicox<-matrix(0,246,11)
  colnames(result_unicox)<-c("geneName","coef","exp(coef)[HR]","se(coef)","Z[coef/se]","Pr(>|Z|)[pval]","lower.95","upper.95","Likelihood ratio test","Wald test","Log rank test")
  
  for(jj in 1:246){
    
    coxph.fit<-coxph(BaSurv~ov[,2+jj],data = ov)
    result=summary(coxph.fit)
    aa<-result$coefficients
    bb<-result$conf.int
    temp<-t(c(colnames(ov)[2+jj],aa,bb[c(3,4)],result[[9]][3],result[[12]][3],result[[10]][3]))
    result_unicox[jj,]<-temp
  }
  
  if(length(which(result_unicox[,"Pr(>|Z|)[pval]"]<0.05))<15 & length(which(result_unicox[,"Pr(>|Z|)[pval]"]<0.05))>6){
    
    a<-result_unicox[which(result_unicox[,"Pr(>|Z|)[pval]"]<0.05),1]
    
    fml<-as.formula(paste0('BaSurv~',paste0(a,collapse = "+")))
    coxph.fit<-coxph(fml,data=ov)
    
    write.table(ov_clin_test,paste(paste("clin",m,sep = "_"),".txt",sep=""),sep="\t",quote=F)
    
    getPrintFile(coxph.fit,a,paste("zzz",m,sep = "_"))
    m=m+1
  }
}

#rm(list = ls())
filezz<-list.files()

ov_clin<-read.table("ov_clin.txt",sep = "\t",header = T,row.names = 1)
rownames(ov_clin)<-paste(rownames(ov_clin),"-01",sep = "")

ov_clin<-na.omit(ov_clin)
ov_exp1<-t(read.table("ov_exp.txt",sep = "\t",header = T,row.names = 1,check.names = F))

aa<-intersect(rownames(ov_exp1),rownames(ov_clin))

ov_clin<-ov_clin[which(rownames(ov_clin)%in%aa),];
ov_clin<-ov_clin[order(ov_clin[,2],decreasing = T),]

phenotype<-read.table('...\\TCGA-OV.GDC_phenotype.tsv',sep = "\t",header = T)

phenotype[,1]<-substr(phenotype[,1],1,15)

phenotype<-phenotype[which(!duplicated(phenotype[,1])),]

rownames(phenotype)<-phenotype[,1]

for (j in 1:m) {
  samplez<-rownames(read.table(filezz[grep(paste(paste("_",j,sep = ""),".txt",sep = ""),filezz)],sep = "\t"))
  
  ov_clin_valid<-ov_clin[-which(rownames(ov_clin)%in%samplez),]
  
  ov_exp111<-t(ov_exp1)
  
  ov_exp111<-ov_exp111[,-which(colnames(ov_exp111)%in%samplez)]
  
  cox_cc<-read.table(filezz[grep(paste(paste("_",j,sep = ""),"_cox.txt",sep = ""),filezz)],sep = "\t",header = T,check.names = F)
  #colnames(cox_cc)
  aa<-cox_cc[,1]
  bb<-cox_cc
  exp_ov<-ov_exp111[which(rownames(ov_exp111)%in%aa),]
  aab<-cox_cc[,2]*exp_ov
  coln_z<-apply(aab,2,sum)
    ov_clin_valid$score<-coln_z[match(rownames(ov_clin_valid),names(coln_z))]

      aa<-intersect(rownames(ov_clin_valid),rownames(phenotype))
      ov_clin_valid<-ov_clin_valid[aa,]
      phenotype1<-phenotype[aa,]
      
      ov_clin_valid<-cbind(ov_clin_valid,phenotype1)
      ov_clin_valid1<-ov_clin_valid[,c("OS","OS.Time","score","neoplasm_histologic_grade","clinical_stage","age_at_initial_pathologic_diagnosis" )]
      
      ov_clin_valid1[grep("IIIA | IIIB | IIIC",ov_clin_valid1[,"clinical_stage"]),"clinical_stage"]<-3
      ov_clin_valid1[grep("IIA|IIB|IIC",ov_clin_valid1[,"clinical_stage"]),"clinical_stage"]<-2
      ov_clin_valid1[grep("IA|IB|IC",ov_clin_valid1[,"clinical_stage"]),"clinical_stage"]<-1
      ov_clin_valid1[grep("IV",ov_clin_valid1[,"clinical_stage"]),"clinical_stage"]<-4
      
      ov_clin_valid1[grep("G4",ov_clin_valid1[,"neoplasm_histologic_grade"]),"neoplasm_histologic_grade"]<-4
      ov_clin_valid1[grep("G1",ov_clin_valid1[,"neoplasm_histologic_grade"]),"neoplasm_histologic_grade"]<-1
      ov_clin_valid1[grep("G3",ov_clin_valid1[,"neoplasm_histologic_grade"]),"neoplasm_histologic_grade"]<-3
      ov_clin_valid1[grep("G2",ov_clin_valid1[,"neoplasm_histologic_grade"]),"neoplasm_histologic_grade"]<-2
      ov_clin_valid1[grep("GX|GB",ov_clin_valid1[,"neoplasm_histologic_grade"]),"neoplasm_histologic_grade"]<-0
      
      ov_clin_valid11<-apply(ov_clin_valid1,2,as.numeric)
      rownames(ov_clin_valid11)<-rownames(ov_clin_valid1)
      ov_clin_valid11<-as.data.frame(ov_clin_valid11)
      
      BaSurv<-Surv(ov_clin_valid11[,"OS.Time"],ov_clin_valid11[,"OS"])
      a=colnames(ov_clin_valid11)[3:6]
      
      fml<-as.formula(paste0('BaSurv~',paste0(a,collapse = "+")))
      coxph.fit<-coxph(fml,data=ov_clin_valid11)
      coxph.fit_uni<-coxph(Surv(ov_clin_valid11[,"OS.Time"],ov_clin_valid11[,"OS"])~ov_clin_valid11[,"score"],data=ov_clin_valid11)
      
      if(summary(coxph.fit_uni)[["coefficients"]][5]<0.05 && summary(coxph.fit)[["coefficients"]][1,5]<0.05){
        
          high<-names(which(coln_z>median(coln_z)))
          low<-names(which(coln_z<median(coln_z)))
          ov_clin_valid$type<-NA
          ov_clin_valid[which(rownames(ov_clin_valid)%in%high),"type"]<-"high"
          ov_clin_valid[which(rownames(ov_clin_valid)%in%low),"type"]<-"low"
          ov_clin_valid<-na.omit(ov_clin_valid)

          km<-survdiff(Surv(OS.Time, OS) ~ type, data=ov_clin_valid)
          p.val <- 1 - pchisq(km$chisq, length(km$n) - 1)
          if(p.val<0.05){
            
            filezz1<-list.files("D:\\xxx\\OV\\zzz")
            a<-filezz1[grep("GSE",filezz1)]
            b<-a[grep("txt",a)]
            for(ii in 1:length(b)){
              exp_zz<-read.table(paste("D:\\xxx\\OV\\zzz\\",b[ii],sep=""),sep = "\t",header = T,row.names = 1)
              inter_gene<-intersect(rownames(exp_zz),cox_cc[,1])
              exp_zz<-exp_zz[which(rownames(exp_zz)%in%inter_gene),]
              
              asdasd<-cox_cc[which(cox_cc[,1]%in%inter_gene),]
              expzz1<-exp_zz[match(asdasd[,1],rownames(exp_zz)),]
              expzz1<-t(na.omit(t(expzz1)))
              expzz11<-asdasd[,"coef"]*expzz1
              coln_z1<-apply(expzz11,2,sum)
              high<-names(which(coln_z1>median(coln_z1)))
              low<-names(which(coln_z1<median(coln_z1)))
              # "Survival","Events"
              
              clin<-read.csv(paste("D:\\xxx\\OV\\zzz\\", paste(substr(b[ii],1,nchar(b[ii])-8),"_survival.csv",sep = ""),sep=""),header = T)
              clin$score<-coln_z1[match(clin[,1],names(coln_z1))]
              
              coxpph_GEO<-coxph(Surv(month, status) ~ score, data=clin)
              if(summary(coxpph_GEO)[["coefficients"]][5]<0.05){
                clin$type<-NA
                clin[which(clin[,1]%in%high),"type"]<-"high"
                clin[which(clin[,1]%in%low),"type"]<-"low"
                clin<-na.omit(clin)
                km1<-survdiff(Surv(month, status) ~ type, data=clin)
                p.val1 <- 1 - pchisq(km1$chisq, length(km1$n) - 1)
                if(p.val1<0.05){
                  print(j)
              }
            }
          }
    }
  }
}
