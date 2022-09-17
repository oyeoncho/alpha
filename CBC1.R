library(readxl)
library(tidyverse)
Cx <- read_excel("cervical cancer_new.xlsx", sheet = "Exosome") 
Cx <- Cx %>% mutate(stage_m=ifelse(stage=="1B"|stage=="2A"|stage=="2B", 1, 
                                   ifelse(stage=="3C1"|stage=="3A"|stage=="3B", 2, 3)))
Cx <- Cx %>% filter(EQD2>55) %>% filter(css==1 | fu_date>12) %>% filter(CTx=="cisplatin") 
Cx %>% filter(stage=="4B") #2


CBC <- read_excel("CBC.xlsx", sheet="Exosome") 
CBC <- CBC %>% as_tibble() %>% mutate(wbc=as.numeric(wbc), hb=as.numeric(hb), neu=as.numeric(neu), lym=as.numeric(lym),
                                      anc=wbc*neu*10, alc=wbc*lym*10, nlr=neu/lym) %>% arrange(date) 
CBC %>% ggplot(aes(x=date)) + geom_histogram(bins=40)
CBC <- CBC %>% filter( date<40)
CBC %>% ggplot(aes(x=date, y=alc)) +geom_point()

z1 <- c()
for(i in 1:41){
  z <- CBC %>% filter(No==i) %>% nrow()
  z1 <- rbind(z1, cbind(i, z))
}
z1 %>% as_tibble() %>% filter(z>0) %>% summary()


z1 <- c()
for(i in 1:325){
  z <- CBC %>% filter(No==i) %>% nrow()
  z1 <- rbind(z1, cbind(i, z))
}

z1 %>% as_tibble() %>% filter(z>0) %>% summary()

######### individual ALC
m1 <- c()

for (i in 1:41) {
  CBC1 <- CBC %>% filter(No==i) %>% filter(is.na(alc)!=T)
  alc <- c(); lm_alc <- c(); m <- c()  
  
  for (a in 1:20) {
    lm_alc <- lm(alc ~ exp(-date*(a/100)), data=CBC1)
    r_alc <- summary(lm_alc)$r.squared
    df <- c(r_alc, a); alc <- rbind(alc, df)
  }
  
  colnames(alc) <- c("R", "a")
  sel <- alc %>% as_tibble() %>% arrange(desc(R)) %>%   as.data.frame()  
  m2 <- c()
  for (j in 1:10){
    a <- sel[j,2]
    model <- lm(alc ~ exp(-date*(a/100)) , data=CBC1)
    a1 <- model$coefficients[2]
    e <-  model$coefficients[1]
    m <- c(i, sel[j,1], a, a1, e); m2 <- rbind(m2, m)
  }
  colnames(m2) <- c("No","R","a","a1","e")
  sel1 <- m2 %>% as_tibble() %>% filter(e>0 & a1>0) %>% arrange(a) %>% slice(1) %>% mutate(min_alc=min(CBC1$alc)) %>%  as.data.frame() 
  m1 <- rbind(m1,sel1)
  # garph
  a <- sel1[,3]
  a1 <- sel1[,4]
  e <- sel1[,5]
  time <- c(1:40)
  p <- data.frame(time, eval(expression(a1*exp(-time*(a/100))+e)))
  names(p) <- c("time", "value")
  plot(alc~date, data=CBC1,xlab="",ylab="", col="blue",  pch=18)
  mtext("Abolute lymphocyte counts (cells/uL)", side=2, line=2, cex=1.5)
  mtext("Time from RT", side=1, line=2, cex=1.5)
  lines(p$time, p$value, col="red", lwd=5)
}

m1 <- m1 %>% mutate(a2=a1+e, a2=round(a2,0), e=round(e,0), min_alc=round(min_alc,0), a=a/100)
Cx <-  Cx%>% inner_join(m1, by="No") 

summary(Cx)

Cx1 <- Cx %>% mutate(rtf=ifelse(`RT field`=="L5",0,1), recur_m=ifelse(recur1>0,1,0),
                     dm_m=ifelse(recur1==2,1,0), lp_m=ifelse(recur1==1,1,0), lpdm_m=ifelse(recur1==3,1,0),
                     a2_m=ifelse(a2<1393,1,0), e_m=ifelse(e<101,1,0), alpha=ifelse(a<0.07,1,0),
                     ALC0_m=ifelse(ALC0<1746 ,1,0), min_alc_m=ifelse(min_alc<203,1,0),
                     age_m=ifelse(Age<50,1,0), path_m=ifelse(pathology=="sqcc", 0, 1),
                     nlr_m=ifelse(NLR<2.72,0,1), risk=ifelse(nlr_m==1 | alpha==1,"high","low"),
                     risk1=ifelse(nlr_m==1 & alpha==1,"high",ifelse(risk=="low", "low", "intermediate")))

library(moonBook)
library(survival)
library(survminer)
Cx1 %>% filter(css==1) %>% select(fu_date)
summary(Cx1)
#######
fit=survfit(Surv(fu_date, css==1)~1, data=Cx1)
ggsurvplot(fit)
summary(fit, times=30) # 89%

Cx1$recur_date
fit=survfit(Surv(recur_date, recur_m==1)~1, data=Cx1)
summary(fit, times=30) # 76.7%%


Cx1$risk <- factor(Cx1$risk, levels=c("low", "high"))
fit=survfit(Surv(fu_date, css==1)~risk, data=Cx1)
ggsurvplot(fit, xlab="Months", ylab="Disease specific survival rate", pval=T, pval.size=8, fun="pct", conf.int=F, risk.table = T, size=1, linetype = "solid",
           palette = c("blue","red"), legend= "bottom", legend.title=expression(paste("NLR & ", alpha)),
           legend.lab=c("Both <2.72 and \u22650.07","Both \u22652.72 or <0.07"), 
           break.time.by=6, xlim=c(0,36), risk.table.height=0.3, ylim=c(0,100))
summary(fit, times=30) # 100, 84.9%

fit=survfit(Surv(fu_date, css==1)~risk1, data=Cx1)
ggsurvplot(fit, xlab="Months", ylab="Disease specific survival rate", pval=T, pval.size=8, fun="pct", conf.int=F, risk.table = T, size=1, linetype = "solid")

fit=survfit(Surv(fu_date, css==1)~alpha, data=Cx1)
ggsurvplot(fit, xlab="Months", ylab="Disease specific survival rate", pval=T, pval.size=8, fun="pct", conf.int=F, risk.table = T, size=1, linetype = "solid")

fit=survfit(Surv(fu_date, css==1)~nlr_m, data=Cx1)
ggsurvplot(fit, xlab="Months", ylab="Disease specific survival rate", pval=T, pval.size=8, fun="pct", conf.int=F, risk.table = T, size=1, linetype = "solid")

############3
library(edgeR)
library(lubridate)
library(Hmisc) # rcorr
library(leaps) # rugsubset
load(file='rt_fc1.Rdata')
miR <- rt_fc[["miR"]]; piR <- rt_fc[["piR"]]; snoR <- rt_fc[["snoR"]]
mR <- rt_fc[["mR"]]; lncR <- rt_fc[["lncR"]]; snR <- rt_fc[["snR"]]
tR <- rt_fc[["tR"]]; yR <- rt_fc[["yR"]]
cx <- Cx1 %>% mutate(ID=substr(SID, 1,8)) %>% select(ID, css, recur_m, a, a2, e, NLR)

##
mR <- mR %>% select(-contains("pvalue")) %>% 
  column_to_rownames(var="Gene_Symbol") %>% t %>% 
  as.data.frame() %>% rownames_to_column(var="ID") %>% 
  as_tibble() %>% mutate(ID=substr(ID, 1,8))  # key ID
mR <- mR %>% 
  inner_join(cx, by="ID") 
mRa <- mR %>% column_to_rownames(var="ID")
x<-select_if (mRa, is.numeric)
rx <- rcorr (as.matrix (x))$r
result <-ifelse (rcorr (as.matrix (x))$P<0.05, rx, 0)
mR_cor <- result %>% as_tibble() %>% mutate(Gene_Symbol = colnames(result)) %>% 
  select(Gene_Symbol, css:NLR) %>% filter(Gene_Symbol!="a"& Gene_Symbol!="a2"& Gene_Symbol!="e"& Gene_Symbol!="NLR"& Gene_Symbol!="css"& Gene_Symbol!="recur_m")

## selection according to correlation of alpha , NLR, css 
s <- mR_cor  %>% filter(abs(NLR) > 0 | abs(a) > 0) %>% 
  mutate(r=ifelse(a>0 & NLR>0, 1, ifelse(a<0 & NLR<0,1,0)), r1=ifelse(NLR>0 & css>0,1, ifelse(NLR<0 &css <0,1,
                                                                                              ifelse(a<0 & css>0, 1,
                                                                                                     ifelse(a>0 & css<0,1,0))))) %>% filter(r==0 & r1==1) 

s1 <- s %>% filter(abs(recur_m)>0 | abs(css)>0.4 ) %>% arrange(desc(abs(css)))
c <- mR %>% select(s1$Gene_Symbol,css)

c[is.na(c)] <- 0
g <-  regsubsets(css ~ .,
                 data = c,
                 nbest = 1,       # 1 best model for each number of predictors
                 nvmax = 15,    # NULL for no limit on number of variables
                 force.in = NULL, force.out = NULL,
                 method = "exhaustive")
out <- summary(g)

plot(g, scale = "adjr2")
polygon(c(2.5,3.5,3.5,2.5),c(1.5,1.5,18,18),col=adjustcolor("red",alpha=0.5), border=NA)
polygon(c(5.5,6.5,6.5,5.5),c(1.5,1.5,18,18),col=adjustcolor("red",alpha=0.5), border=NA)
polygon(c(12.5,13.5,13.5,12.5),c(2.5,2.5,18,18),col=adjustcolor("red",alpha=0.5), border=NA)
polygon(c(20.5,21.5,21.5,20.5),c(3.5,3.5,18,18),col=adjustcolor("red",alpha=0.5), border=NA)

plot(out$adjr2, type="o", xlab="", ylab="",xaxt="n",yaxt="n")
axis(side=1, at=seq(0,15,1),cex.axis=1.5)
axis(side=2, at=seq(0,1,0.05),cex.axis=1.5)
mtext(expression(paste("Adujsted ", R^2)), side=2, line=2.5, cex=1.5)
mtext("Number of selected mRNAs", side=1, line=2.5, cex=1.5)
abline(v=4, col="red", lwd=2, lty=1)
which.max(out$adjr2)
out$which[4,]

s %>% filter(Gene_Symbol=="CCDC113"|Gene_Symbol=="STX6"|Gene_Symbol=="E2F8"|Gene_Symbol=="ACOT9") # CCDC114 (ciliated cells): NLR, E2F8 (plasma cell): alpha, STX6 (Machrophage): alpha, ACOT9: monocyte (NLR))

mR1 <- mR %>% select(ID, E2F8, STX6, CCDC113, ACOT9)

library(ggpubr)
library(survminer)

test <- Cx1 %>% mutate(ID=substr(SID, 1,8))  %>% inner_join(mR1, by="ID") %>% 
  mutate(t1=E2F8+STX6, t1_m=ifelse(t1<0.2429,0,1), t2=CCDC113-ACOT9, t2_m=ifelse(t2 < -0.3184,0,1),t=CCDC113+E2F8+STX6-ACOT9, t_m=ifelse(t < -1.179, 0,1), 
         t3=STX6-ACOT9, t3_m = ifelse(t3<0.3657,0,1), rtf=`RT field`)

summary(test)

###########
out=mytable(risk~age_m+pathology+stage_m+rtf+duration+EQD2+ALC0_m+min_alc_m+a2_m+e_m+recur1+css+E2F8+STX6+CCDC113+ACOT9+t1+t2+t+t3, data=test, digits=2, method=3)
mycsv(out, file="table2.csv")
out

## zero number
test %>% filter(E2F8==0) #11
test %>% filter(STX6==0) #0
test %>% filter(CCDC113==0) #7
test %>% filter(ACOT9==0) #0
test %>% filter(t2==0) #0
test %>% filter(t3==0) #0

#############
test %>% mutate(`Disease specific death`=ifelse(css==1,"Yes","No"), `E2F8+STX6(log2FC)`=t1) %>% 
  ggboxplot(x="Disease specific death",y="E2F8+STX6(log2FC)", palette = "jco", add="jitter", color="Disease specific death")+stat_compare_means(label="p.signif", label.x=1.5, label.y=10, cex=10)+
  theme(axis.title =element_text(face="bold", size=20), axis.text = element_text(size=20), legend.position ="NA")+ coord_cartesian(ylim=c(-7,11))

fit=survfit(Surv(fu_date, css==1)~t1_m, data=test)
ggsurvplot(fit, xlab="Months", ylab="Disease specific survival rate", pval=T,pval.size=10, fun="pct", conf.int=F, risk.table = T, size=1, linetype = "solid",
           palette = c("blue","red"), legend= "bottom", legend.title="E2F8+STX6 (log2FC)", legend.lab=c("<0.2429", "\u22650.2429"),  break.time.by=6, xlim=c(0,36), ylim=c(0,100))
summary(fit, times=30) # 100% , 77.2%

test %>% mutate(`Disease specific death`=ifelse(css==1,"Yes","No"), `E2F8+STX6+CCDC113-ACOT9 (log2FC)`=t) %>% 
  ggboxplot(x="Disease specific death",y="E2F8+STX6+CCDC113-ACOT9 (log2FC)", palette = "jco", add="jitter", color="Disease specific death")+stat_compare_means(label="p.signif", label.x=1.5, label.y=10, cex=10)+
  theme(axis.title =element_text(face="bold", size=17), axis.text = element_text(size=20), legend.position = "NA")

fit=survfit(Surv(fu_date, css==1)~t_m, data=test)
ggsurvplot(fit, xlab="Months", ylab="Disease specific survival rate", pval=T, pval.size=10, fun="pct", conf.int=F, risk.table = T, size=1, linetype = "solid",
           palette = c("blue","red"), legend= "bottom", legend.title="E2F8+STX6+CCDC113-ACOT9 (log2FC)", legend.lab=c("< -1.179", "\u2265 -1.179"), break.time.by=6, xlim=c(0,36), ylim=c(0,100))
summary(fit, times=30) # 100%, 78.46%

test %>% mutate(`Disease specific death`=ifelse(css==1,"Yes","No"), `CCDC113-ACOT9 (log2FC)`=t2) %>% 
  ggboxplot(x="Disease specific death",y="CCDC113-ACOT9 (log2FC)", palette = "jco", add="jitter", color="Disease specific death")+stat_compare_means(label="p.signif", label.x=1.5, label.y=10, cex=10)+
  theme(axis.title =element_text(face="bold", size=20), axis.text = element_text(size=20), legend.position = "NA")+ coord_cartesian(ylim=c(-8,11))

fit=survfit(Surv(fu_date, css==1)~t2_m, data=test)
ggsurvplot(fit, xlab="Months", ylab="Disease specific survival rate", pval=T, pval.size=10, fun="pct", conf.int=F, risk.table = T, size=1, linetype = "solid",
           palette = c("blue","red"), legend= "bottom", legend.title="CCDC113-ACOT9 (log2FC)", legend.lab=c("< -0.3184", "\u2265 -0.3184"), break.time.by=6, xlim=c(0,36), ylim=c(0,100))
summary(fit, times=30) # 100%, 78.2%

test %>% mutate(`Disease specific death`=ifelse(css==1,"Yes","No"), `STX6-ACOT9 (log2FC)`=t3) %>% 
  ggboxplot(x="Disease specific death",y="STX6-ACOT9 (log2FC)", palette = "jco", add="jitter", color="Disease specific death")+stat_compare_means(label="p.signif", label.x=1.5, label.y=10, cex=10)+
  theme(axis.title =element_text(face="bold", size=20), axis.text = element_text(size=20), legend.position = "NA")+ coord_cartesian(ylim=c(-8,11))

fit=survfit(Surv(fu_date, css==1)~t3_m, data=test)
ggsurvplot(fit, xlab="Months", ylab="Disease specific survival rate", pval=T, pval.size=10, fun="pct", conf.int=F, risk.table = T, size=1, linetype = "solid",
           palette = c("blue","red"), legend= "bottom", legend.title="STX-ACOT9 (log2FC)", legend.lab=c("< 0.3657", "\u2265 0.3657"), break.time.by=6, xlim=c(0,36), ylim=c(0,100))
summary(fit, times=30) # 100%, 77.7%

mytable(t1_m~age_m+pathology+stage_m+rtf+duration+EQD2+ALC0_m+min_alc_m+a2_m+e_m+recur1+css+E2F8+STX6+CCDC113+ACOT9+t1+t2+t+alpha+a+NLR+nlr_m, data=test, digits=2, method=3)
mytable(risk~age_m+pathology+stage_m+rtf+duration+EQD2+ALC0_m+min_alc_m+a2_m+e_m+recur1+css+E2F8+STX6+CCDC113+ACOT9+t1+t2+t, data=test, digits=2, method=3)

summary(lm(a~t1, data=test))
summary(lm(a~STX6, data=test))
summary(lm(a~E2F8, data=test))

test %>% filter(risk=="high") %>% select(Name, Age, E2F8, STX6, CCDC113, ACOT9, t1, t2, t,a, NLR,css, recur1, fu_date) %>% arrange(desc(t)) 
test %>% filter(risk=="high") %>% select(Name, Age, E2F8, STX6, CCDC113, ACOT9, t1, t2, t,a, NLR,css, recur1, fu_date) %>% arrange(desc(t1)) 
test %>% filter(risk=="high") %>% select(Name, Age, E2F8, STX6, CCDC113, ACOT9, t1, t2, t,a, NLR,css, recur1, fu_date) %>% arrange(desc(t2))

