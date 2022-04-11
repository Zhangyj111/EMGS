PD1 PDCD1
PDL1 CD274

rm(list = ls())

setwd("...\\Figure6")

exp<-read.table("ov_exp_all.txt",sep = "\t",header = T,row.names = 1,check.names = F)

exp1<-as.data.frame(t(exp[which(rownames(exp)%in%c("PDCD1","CD274")),]))

cox_cc<-read.table("...\\zzz_2501_cox.txt",sep = "\t",header = T,check.names = F)

exp_zz<-read.table("ov_exp.txt",sep = "\t",header = T,row.names = 1,check.names = F)

inter_gene<-intersect(rownames(exp_zz),cox_cc[,1])
exp_zz<-exp_zz[which(rownames(exp_zz)%in%inter_gene),]

asdasd<-cox_cc[which(cox_cc[,1]%in%inter_gene),]
expzz1<-exp_zz[match(asdasd[,1],rownames(exp_zz)),]
expzz1<-t(na.omit(t(expzz1)))
expzz11<-asdasd[,"coef"]*expzz1##EMGS score
coln_z1<-apply(expzz11,2,sum)
# "Survival","Events"

exp1$score<-coln_z1[match(rownames(exp1),names(coln_z1))]

library(ggplot2)

p1=ggplot(exp1,aes(x=score,y=PDCD1))+
  geom_point()+
geom_smooth(method = "lm")+
  theme_bw()

cor(exp1[,"score"],exp1[,"PDCD1"],method = "spearman")

cor(exp1[,"score"],exp1[,"CD274"],method = "spearman")

p2=ggplot(exp1,aes(x=score,y=CD274))+
  geom_point()+
  geom_smooth(method = "lm")+
  theme_bw()

wilcox.test(y, x, exact = TRUE)

a<-exp1[which(exp1[,"PDCD1"]>median(exp1[,"PDCD1"])),"score"]
b<-exp1[which(exp1[,"PDCD1"]<median(exp1[,"PDCD1"])),"score"]
wilcox.test(a, b, exact = TRUE)

boxplot(a,b,col = c("#DF7068","#23ABAF"),names=c("High risk","Low risk"))
c<-exp1[which(exp1[,"CD274"]>median(exp1[,"CD274"])),"score"]
d<-exp1[which(exp1[,"CD274"]<median(exp1[,"CD274"])),"score"]
wilcox.test(c,d, exact = TRUE)
boxplot(c,d,col = c("#DF7068","#23ABAF"),names=c("High risk","Low risk"))
