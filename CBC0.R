library(readxl)
library(tidyverse)

Cx <- read_excel("cervical cancer_new.xlsx", sheet = "CCRT") 
Cx <- Cx %>% mutate(stage_m=ifelse(stage=="1B"|stage=="2A"|stage=="2B", "IB-IIB", 
                                   ifelse(stage=="3C1"|stage=="3A"|stage=="3B", "IIIA-IIIC1", "IIIC2-IVB")))
## selection 
Cx <- Cx %>% filter(css==1 | fu_date>12) %>% filter(CTx=="Cisplatin") 
Cx %>% filter(stage=="4B") #5

# replace two missing values with median values for NLR and ALC 
Cx <- Cx %>% mutate(NLR=ifelse(is.na(NLR)==T, summary(Cx$NLR)[3], NLR), ALC0=ifelse(is.na(ALC0)==T, summary(Cx$ALC0)[3], ALC0) ) 

## ALC 1927 from 323
CBC <- read_excel("CBC.xlsx", sheet="CCRT") 
CBC <- CBC %>% as_tibble() %>% mutate(wbc=as.numeric(wbc), hb=as.numeric(hb), neu=as.numeric(neu), lym=as.numeric(lym),
                                      anc=wbc*neu*10, alc=wbc*lym*10, nlr=neu/lym) %>% arrange(date) 
CBC %>% ggplot(aes(x=date)) + geom_histogram(bins=30)+theme(axis.title =element_text(face="bold", size=20), axis.text = element_text(size=20))+
  xlab("Time from radiation therapy (days)")+ylab("Number of ALC measurements")+scale_x_continuous(breaks=seq(0,50,7))

CBC <- CBC %>% filter(date<40)  %>% filter(No !=133 & No!=271)

z1 <- c()
for(i in 1:325){
  z <- CBC %>% filter(No==i) %>% nrow()
  z1 <- rbind(z1, cbind(i, z))
}

z1 %>% as_tibble() %>% filter(z>0) %>% summary()

######### ALL ALC
alc <- c()
### fitting R2 value (alpha 0.01~0.2)
for (a in 1:20) {
  lm_alc <- lm(alc ~ exp(-date*(a/100)), data=CBC)
  r_alc <- summary(lm_alc)$r.squared
  df <- c(r_alc, a); alc <- rbind(alc, df)
}
colnames(alc) <- c("R", "a")
sel <- alc %>% as_tibble() %>%  arrange(desc(R)) %>%   as.data.frame()  
# top 10 of 20 R2 selection-> calculate a1, e for 10 alpha equation
m2 <- c()
for (j in 1:10){
  a <- sel[j,2]
  model <- lm(alc ~ exp(-date*(a/100)) , data=CBC)
  a1 <- model$coefficients[2]
  e <-  model$coefficients[1]
  m <- c(sel[j,1], a, a1, e); m2 <- rbind(m2, m)
}
colnames(m2) <- c("R","a","a1","e")
m2 %>%as_tibble() %>%  ggplot(aes(x=a,y=a1))+geom_point()
m2 %>%as_tibble() %>%  ggplot(aes(x=a,y=e))+geom_point()
# selection e (lower limit)>0 arrange to alpha
sel1 <- m2 %>% as_tibble() %>% filter(e>0) %>% arrange(a) %>% as.data.frame()  
sel1
CBC

## ALL ALC biggest a (A) & smallest a(B) of equations with top 10 R2 
plot(alc~date, data=CBC,xlab="",ylab="", col="gray",  pch=18, xaxt="n", yaxt="n", xlim=c(0,70))
axis(side=1, at=seq(0,70,7),cex.axis=1.5)
axis(side=2, at=seq(0,3600,200),cex.axis=1.5)
mtext(expression(paste("Abolute lymphocyte counts (cells/",mu*L,")")), side=2, line=2.5, cex=1.5)
mtext("Time from radiation therapy (days)", side=1, line=2.5, cex=1.5)
text(35,3200,expression(italic(ALC==a1*e^(-time*alpha)+e1)), cex=2)
time <- c(0:70)

p <- data.frame(time, eval(expression(1888*exp(-time*(0.15))+359)))
names(p) <- c("time", "value")
lines(p$time, p$value, col="blue", lwd=5, lty=1)
abline(h=359, col="blue", lwd=2, lty=2)
text(35,2800,expression(paste(bold("(A):"), italic(ALC==1888*e^(-time*0.15)+359))), cex=1.5, col="blue")
text(4.5,1600,expression(paste(bold("(A)"))),cex=1.5, col="blue")


p <- data.frame(time, eval(expression(1412*exp(-time*(0.06))+104)))
names(p) <- c("time", "value")
lines(p$time, p$value, col="red", lwd=5, lty=1)
abline(h=104, col="red", lwd=2, lty=2)
text(35,2600,expression(paste(bold("(B):"), italic(ALC==1412*e^(-time*0.06)+104))), cex=1.5, col="red")
text(12,900,expression(paste(bold("(B)"))),cex=1.5, col="red")

p <- data.frame(time, eval(expression(1888*exp(-time*(0.06))+359)))
names(p) <- c("time", "value")
lines(p$time, p$value, col="black", lwd=2, lty=1)
text(35,2400,expression(paste(bold("(C):"), italic(ALC==1888*e^(-time*0.06)+359))), cex=1.5, col="black")
text(15.5,1200,expression(paste(bold("(C)"))), cex=1.5, col="black")

#a2=a1+e1 (Estimated pre ALC)
arrows(0,1550,0,2150, col="black", lwd=3, angle=15)
text(3.5,2300,expression(italic(a2==a1+e1)), cex=1.5, col="black")

#e1 (Estimated ALC nadir)
arrows(70,105,70,360, col="black", lwd=4, angle=15)
text(70,450,expression(italic(e1)), cex=1.5, col="black")

text(35,2200,expression(paste(italic(alpha),bold("(slope)"),": (A) > (B)=(C)")), cex=1.5, col="black")
text(35,2000,"a2(estimated pre ALC): (B) < (A)=(C)", cex=1.5, col="black")
text(35,1800,"e1(estimated ALC nadir): (B) < (A)=(C)", cex=1.5, col="black")


######### individual ALC
CBC <- read_excel("CBC.xlsx", sheet="CCRT") 
CBC <- CBC %>% as_tibble() %>% mutate(wbc=as.numeric(wbc), hb=as.numeric(hb), neu=as.numeric(neu), lym=as.numeric(lym),
                                      anc=wbc*neu*10, alc=wbc*lym*10, nlr=neu/lym) %>% arrange(date) 
CBC %>% ggplot(aes(x=date)) + geom_histogram(bins=30)+theme(axis.title =element_text(face="bold", size=20), axis.text = element_text(size=20))+
  xlab("Time from radiation therapy (days)")+ylab("Number of ALC measurements")+scale_x_continuous(breaks=seq(0,50,7))
CBC <- CBC %>% filter(date<40)

m1 <- c()

for (i in 1:325) {
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
  sel1 <- m2 %>% as_tibble() %>% filter(e>0 & a1>0)  %>% arrange(a) %>% slice(1) %>% mutate(min_alc=min(CBC1$alc)) %>%  as.data.frame() 
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

m1 <- m1 %>%  mutate(a2=a1+e, a2=round(a2,0), e=round(e,0), min_alc=round(min_alc,0), a=a/100)
Cx <-  Cx %>% inner_join(m1, by="No") %>% mutate(recur_m=ifelse(recur1>0,1,0))

library(ggpubr)

b1 <- Cx %>% ggplot(aes(R))+geom_density(alpha=.2, fill="#999999")+ 
  ylab("Density")+ xlab(expression(R^2))+  scale_fill_grey() + theme_classic()+
  theme(axis.title =element_text(face="bold", size=20), axis.text = element_text(size=20))

b2 <- Cx %>%  ggplot(aes(x=a))+geom_density(alpha=.2, fill="#999999")+ 
  ylab("Density")+ xlab(expression(alpha))+  scale_fill_grey() + theme_classic()+
  theme(axis.title =element_text(face="bold", size=20), axis.text = element_text(size=20))+
  geom_vline(aes(xintercept=median(a)), color="red", linetype="dashed", size=1)+
  annotate("text", 0.08, 15, label="Median", cex=8, color="black")

a <- Cx %>% select(ALC0) %>% mutate(ALC=ALC0, type=rep("pre ALC", nrow(Cx))) %>% select(-ALC0)
b <- Cx %>% select(a2) %>% mutate(ALC=a2, type=rep("Estimated pre ALC (a2)", nrow(Cx))) %>% select(-a2)
b3 <- a %>% bind_rows(b) %>% ggplot(aes(ALC, fill=type))+ geom_density(alpha=.2)+ 
  ylab("Density")+ xlab(expression(paste("Abolute lymphocyte counts (cells/",mu*L,")")))+  theme_classic()+
  theme(axis.title =element_text(face="bold", size=20), axis.text = element_text(size=20), 
        legend.position = "top", legend.title =element_text(size=15), legend.text =element_text(size=15) ) +
  geom_vline(aes(xintercept=1629), color="red", linetype="dashed", size=1)+
  geom_vline(aes(xintercept=1884), color="blue", linetype="dashed", size=1)+
  annotate("text", 1400, 0.0004, label="Median (a2)", cex=8, color="black")+
    annotate("text", 2200, 0.0006, label="Median (pre ALC)", cex=8, color="black")

a <- Cx %>% select(min_alc) %>% mutate(ALC=min_alc, type=rep("Min ALC", nrow(Cx))) %>% select(-min_alc)
b <- Cx %>% select(e) %>% mutate(ALC=e, type=rep("Estimated ALC nadir (e1)", nrow(Cx))) %>% select(-e)
b4 <- a %>% bind_rows(b) %>% ggplot(aes(ALC, fill=type))+ geom_density(alpha=.2)+ 
  ylab("Density")+ xlab(expression(paste("Abolute lymphocyte counts (cells/",mu*L,")")))+  theme_classic()+
  theme(axis.title =element_text(face="bold", size=20), axis.text = element_text(size=20), 
        legend.position = "top", legend.title =element_text(size=15), legend.text =element_text(size=15) ) +
geom_vline(aes(xintercept=144), color="red", linetype="dashed", size=1)+
geom_vline(aes(xintercept=276), color="blue", linetype="dashed", size=1)+
annotate("text", 100, 0.003, label="Median (e1)", cex=8, color="black")+
annotate("text", 350, 0.002, label="Median (Min ALC)", cex=8, color="black")

b5 <- Cx %>% ggplot(aes(x=ALC0, y=a2))+geom_point()+geom_smooth(method = "lm", se=F)+
  ylab(expression(paste("Estimated pre ALC (cells/",mu*L,")")))+ xlab(expression(paste("pre ALC (cells/",mu*L,")")))+  theme_classic()+
  annotate("text", 2500, 4000, label=expression(paste(italic(R^2==0.2648))), cex=8, color="black")+
  theme(axis.title =element_text(face="bold", size=20), axis.text = element_text(size=20))

b6 <- Cx %>% ggplot(aes(x=min_alc, y=e))+geom_point()+geom_smooth(method = "lm", se=F)+
  ylab(expression(paste("Estimated ALC nadir (cells/",mu*L,")")))+ xlab(expression(paste("MIn ALC (cells/",mu*L,")")))+  theme_classic()+
  annotate("text", 400, 600, label=expression(paste(italic(R^2==0.2913))), cex=8, color="black")+
  theme(axis.title =element_text(face="bold", size=20), axis.text = element_text(size=20))
  
ggarrange(b1,b2,b3,b5,b4,b6, labels=c("A","B","C","D","E","F"), font.label = list(size = 30, color = "black"),  nrow=3, ncol=2)


# mean 30
a <- Cx %>% arrange(recur_date) %>% filter(css==1) %>% select(No, recur_date, recur1, fu_date, css) %>% 
  mutate(progression=ifelse(recur1==1, "LP", ifelse(recur1==2, "DM", "LP+DM")), progression=factor(progression, levels = c("LP", "DM", "LP+DM"))) %>% 
  ggplot(aes(x=recur_date, y=fu_date, color=progression)) + geom_point(size=3)+ scale_color_manual(values=c("gray", "blue","red"))+theme_bw()+
  ylab("Time taken for disease specific deaths (Months)")+ xlab("Time taken for cancer progressions (Months)")+
  theme(axis.title =element_text(face="bold", size=20), axis.text = element_text(size=20), legend.position = "top")+
  scale_x_continuous(breaks = seq(0, 150,5))+scale_y_continuous(breaks = seq(0, 120,5))+
  geom_hline(aes(yintercept=mean(fu_date)), color="red", linetype="dashed", size=1)+
  annotate("text", 45, 32, label="Average of time taken for disease specific deaths", cex=5, color="red")+
  geom_segment(aes(x=80,y=30,xend=80,yend=50), arrow=arrow(), color="red", size=2)+
  geom_segment(aes(x=78,y=30,xend=78,yend=10), arrow=arrow(), color="red", size=2)+
  geom_text(aes(x=80,53), label="Non-aggressive group", color="black", size=5)+
  geom_text(aes(x=78,7), label="Aggressive group", color="black", size=5)

b <- Cx %>% filter(css==1) %>% ggplot(aes(fu_date))+geom_histogram(aes(y=..density..), color="black", fill="white", bins=50)+geom_density(alpha=.2, fill="#999999")+ 
  ylab("Density")+ xlab("Time taken for disease specific deaths (Months)")+  scale_fill_grey() + theme_classic()+ 
  theme(axis.title =element_text(face="bold", size=20), axis.text = element_text(size=20))+
  geom_vline(aes(xintercept=mean(fu_date)), color="red", linetype="dashed", size=1)+
  annotate("text", 38, 0.03, label="Average", cex=8, color="red")+
  scale_x_continuous(expand=expansion(mult=c(0,0)), breaks = seq(0, 150,5))+scale_y_continuous(expand=expansion(mult=c(0,0)))+
  geom_segment(aes(x=30,y=0.04,xend=50,yend=0.04), arrow=arrow(), color="red", size=2)+
  geom_segment(aes(x=30,y=0.039,xend=10,yend=0.039), arrow=arrow(), color="red", size=2)+
  geom_text(aes(x=48,0.042), label="Non-aggressive group", color="black", size=5)+
  geom_text(aes(x=12,0.041), label="Aggressive group", color="black", size=5)

# aggressive css median 15, recur 5.5, non-aggressive 50, 12
c <- Cx %>% filter(css==1) %>% mutate(Group=ifelse(fu_date<30 & css==1, "aggressive", ifelse(fu_date>=30 &css==1, "non-aggressive", "NA"))) %>%
  ggplot(aes(fu_date, fill=Group))+geom_density(alpha=0.2)+ ylab("Density of disease specific deaths")+ xlab("Time taken for disease specific deaths (Months)")+   
  theme_bw()+theme(axis.title =element_text(face="bold", size=20), axis.text = element_text(size=20), legend.position = "top")+
  scale_x_continuous(expand=expansion(mult=c(0,0)), breaks=seq(0,120,5))+scale_y_continuous(expand=expansion(mult=c(0,0)))+
  geom_vline(aes(xintercept=15), color="red", linetype="dashed", size=1)+ geom_vline(aes(xintercept=50), color="red", linetype="dashed", size=1)+
  annotate("text", 15, 0.03, label="Median", cex=8, color="black")+ annotate("text", 50, 0.02, label="Median", cex=8, color="black")

d <- Cx %>% filter(css==1) %>% mutate(Group=ifelse(fu_date<30 & css==1, "aggressive", ifelse(fu_date>=30 &css==1, "non-aggressive", "NA"))) %>% 
  ggplot(aes(recur_date, fill=Group))+geom_density()+ ylab("Density of cancer progressions")+ xlab("Time taken for cancer progressions (Months)")+   
  theme_bw()+theme(axis.title =element_text(face="bold", size=20), axis.text = element_text(size=20), legend.position = "top")+
  scale_x_continuous(expand=expansion(mult=c(0,0)), breaks=seq(0,120,5))+scale_y_continuous(expand=expansion(mult=c(0,0)))+
  geom_vline(aes(xintercept=5.5), color="red", linetype="dashed", size=1)+ geom_vline(aes(xintercept=12), color="red", linetype="dashed", size=1)+
  annotate("text", 5, 0.06, label="Median", cex=8, color="black")+ annotate("text", 12, 0.03, label="Median", cex=8, color="black")

ggarrange(a,b,c,d,  labels=c("A","B","C","D"), font.label = list(size = 30, color = "black"), nrow=2, ncol=2)

########
Cx1 <- Cx %>%  mutate(Group1=ifelse(fu_date<30 & css==1, 1,0), Group2=ifelse(fu_date>=30 &css==1, 1, 0), 
                     rtf=`RT field`, dm_m=ifelse(recur1==2,1,0), lp_m=ifelse(recur1==1,1,0), lpdm_m=ifelse(recur1==3,1,0),
                     a2_m=ifelse(a2<1629,1,0),e_m=ifelse(e<144,1,0), alpha=ifelse(a<0.08,1,0),
                     ALC0_m=ifelse(ALC0<1884,1,0), min_alc_m=ifelse(min_alc<276,1,0),
                     age_m=ifelse(Age<50,1,0), path_m=ifelse(pathology=="sqcc", 0, 1), 
                     nlr_m=ifelse(NLR<2.43,0,1), risk=ifelse(nlr_m==1 & alpha==1,"high",
                                                           ifelse(nlr_m==0 & alpha==0,"low","intermediate")))
summary(Cx1)
#########################################################
library(moonBook)
library(survival)
library(survminer)

## 69 patients with DS deaths
t <- Cx1 %>%  filter(Group1==1 | Group2==1) %>% mutate(all=ifelse(Group1==1, "Aggressive", "Non-aggressive"))
out=mytable(all~age_m+recur1+dm_m+lp_m+lpdm_m+path_m+stage_m+NLR+ALC0+min_alc+a+a2+e, data=t, digits=2)
mycsv(out, file="reference.csv")
out


###### 8OS 73.7% 
fit=survfit(Surv(fu_date, survival==1)~1, data=Cx1)
ggsurvplot(fit, xlab="Months", ylab="Overall survival rate", fun="pct", conf.int=T, risk.table = T, size=1, break.time.by=12, xlim=c(0,130))
summary(fit, times=96)

###### 8DSS 75.2% 5DSS 80.1%
fit=survfit(Surv(fu_date, css==1)~1, data=Cx1)
ggsurvplot(fit, xlab="Months", ylab="Diease specific survival rate", fun="pct", conf.int=T, risk.table = T, size=1, break.time.by=12, xlim=c(0,130))
summary(fit, times=96)
summary(fit, times=60)

##  8DSS stage 1B~IIB (86) : 82.8%, 3A-3C1 (190): 75.3%, 3C2-4B (47): 61.7%
fit=survfit(Surv(fu_date, css==1)~stage_m, data=Cx1)
summary(fit, times=96)
ggsurvplot(fit, xlab="Months", ylab="Disease specific survival rate", pval=T, fun="pct", conf.int=F, risk.table = T, size=1, linetype = "strata",
           palette = c("red","black","blue"), legend= "bottom", legend.title="Stage", legend.lab=c("IB-IIB", "IIIA-IIIC1", "IIIC2-IVB"), break.time.by=12, xlim=c(0,120))

### 8PFS 65.7 % 
fit=survfit(Surv(recur_date, recur_m==1)~1, data=Cx1)
ggsurvplot(fit, xlab="Months", ylab="Progression free survival rate", fun="pct", conf.int=T, risk.table = T, size=1, break.time.by=12, xlim=c(0,120))
summary(fit, times=96)

### 8 PFS stage 1B~IIB (86) : 5 yrs PFS 80.2%, 3A-3C1 (190): 65.3%, 3C2-4B (47): 41.3% 
fit=survfit(Surv(recur_date, recur_m==1)~stage_m, data=Cx1)
summary(fit, times=96)
ggsurvplot(fit, xlab="Months", ylab="Progression free survival rate", fun="pct", conf.int=F, risk.table = T, size=1, break.time.by=12, xlim=c(0,96))


### 2.5DSS-A 86.2%
fit=survfit(Surv(fu_date, Group1==1)~1, data=Cx1)
summary(fit, times=30)
ggsurvplot(fit, xlab="Months", ylab="Diease specific survival rate (aggressive)", fun="pct", conf.int=T, risk.table = T, size=1, break.time.by=6,xlim=c(0,30))

# 8DSS-NA 87.2% 
fit=survfit(Surv(fu_date, Group2==1)~1, data=Cx1)
summary(fit, times=96)
ggsurvplot(fit, xlab="Months", ylab="Diease specific survival rate (Non-aggressive)", fun="pct", conf.int=T, risk.table = T, size=1, break.time.by=12, xlim=c(30, 120))

## Cox

#DSS
Cx_s <- Cx1 %>% mutate(NLR=ifelse(nlr_m==1,"\u22652.43","< 2.43"), `FIGO stage`=stage_m, Pathology=ifelse(path_m==1, "Non-sqcc", "Sqcc"), 
                       Pathology=factor(Pathology, levels=c("Sqcc", "Non-sqcc")), `pre ALC(a2)`=ifelse(a2_m==1, "<1629", "\u22651629"),
                       `pre ALC(a2)`=factor(`pre ALC(a2)`, levels = c("\u22651629", "<1629")), `ALC nadir(e1)`=ifelse(e_m==1, "<144", "\u2265144"),
                       `ALC nadir(e1)`=factor(`ALC nadir(e1)`, levels = c("\u2265144", "<144")),
                       alpha=ifelse(a<0.08, "<0.08", "\u22650.08"), alpha= factor(alpha, levels=c("\u22650.08","<0.08")),
                       `pre ALC`=ifelse(ALC0_m==1, "<1884", "\u22651884"), `pre ALC`=factor(`pre ALC`, levels = c("\u22651884", "<1884")),
                       `min ALC`=ifelse(min_alc_m==1, "<276", "\u2265276"), `min ALC`=factor(`min ALC`, levels = c("\u2265276", "<276")),
                       Age=ifelse(age_m==1,"<50","\u226550"), Age=factor(Age,levels = c("\u226550","<50")))



Cx_s <- Cx_s %>% mutate(TS = Surv(fu_date, css==1)) 
out=mycph(TS~ Age+Pathology+NLR+`FIGO stage`+`pre ALC`+`min ALC`+ alpha+ `pre ALC(a2)`+ `ALC nadir(e1)`, data=Cx_s)
HRplot(out, type=2, show.CI=TRUE, main="Hazard ratios of all individual variables for DSS")
result<-coxph(Surv(fu_date,css)
              ~Age+Pathology+NLR+`FIGO stage`+`pre ALC`+`min ALC`+ alpha+ `pre ALC(a2)`+ `ALC nadir(e1)`,data=Cx_s)    
finalmodel<-step(result, direction = 'backward')
summary(finalmodel)
ggforest(finalmodel, data=Cx_s, fontsize = 1.5)

# PFS
Cx_s <- Cx_s %>% mutate(TS = Surv(recur_date,recur_m==1)) 
out=mycph(TS~ Age+Pathology+NLR+`FIGO stage`+`pre ALC`+`min ALC`+ alpha+ `pre ALC(a2)`+ `ALC nadir(e1)`, data=Cx_s)
HRplot(out, type=2, show.CI=TRUE, main="Hazard ratios of all individual variables for PFS")
result<-coxph(Surv(recur_date,recur_m)
              ~Age+Pathology+NLR+`FIGO stage`+`pre ALC`+`min ALC`+ alpha+ `pre ALC(a2)`+ `ALC nadir(e1)`,
              data=Cx_s)    
finalmodel<-step(result, direction = 'backward')
summary(finalmodel)
ggforest(finalmodel, data=Cx_s, fontsize = 1.5)

####### aggressive group
Cx_s <- Cx_s %>% mutate(TS = Surv(fu_date,Group1==1))
out=mycph(TS~ Age+Pathology+NLR+`FIGO stage`+`pre ALC`+`min ALC`+ alpha+ `pre ALC(a2)`+ `ALC nadir(e1)`, data=Cx_s)
HRplot(out, type=2, show.CI=TRUE, main="Hazard ratios of all individual variables for DSS (aggressive)")
#result<-coxph(Surv(fu_date,Group1)
#              ~Age+Pathology+NLR+`FIGO stage`+`min ALC`+`pre ALC`+ alpha+ `pre ALC(a2)`+ `ALC nadir(e1)`,
#              data=Cx_s)    
# min ALC - different form expected or UV analysis in MV - noisy factor  
result<-coxph(Surv(fu_date,Group1)
              ~Age+Pathology+NLR+`FIGO stage`+`pre ALC`+ alpha+ `pre ALC(a2)`+ `ALC nadir(e1)`,
              data=Cx_s)    
finalmodel<-step(result, direction = 'backward')
summary(finalmodel)
ggforest(finalmodel, data=Cx_s, fontsize = 1.5)

##non-aggressive group
Cx_sa <- Cx_s %>% mutate(alpha=ifelse(a>=0.08,"\u22650.08","<0.08"), TS = Surv(fu_date,Group2==1))
out=mycph(TS~ Age+Pathology+NLR+`FIGO stage`+`pre ALC`+`min ALC`+ alpha+ `pre ALC(a2)`+ `ALC nadir(e1)`, data=Cx_sa)
HRplot(out, type=2, show.CI=TRUE, main="Hazard ratios of all individual variables for DSS (non-aggressive)")
result<-coxph(Surv(fu_date,Group2)
              ~Age+Pathology+NLR+`FIGO stage`+`pre ALC`+`min ALC`+ alpha+ `pre ALC(a2)`+ `ALC nadir(e1)`,  data=Cx_sa)    
finalmodel<-step(result, direction = 'backward')
summary(finalmodel)
ggforest(finalmodel, data=Cx_sa,  fontsize = 1.5)

###
## stage 3C2-4B
Cx_s1 <-  Cx_s %>% filter(`FIGO stage`=="IIIC2-IVB")
Cx_s1 <- Cx_s1 %>% mutate(TS = Surv(fu_date,Group1==1))
out=mycph(TS~ Age+Pathology+NLR+`pre ALC`+`min ALC`+ alpha+ `pre ALC(a2)`+ `ALC nadir(e1)`, data=Cx_s1)
HRplot(out, type=2, show.CI=TRUE, main="Hazard ratios of all individual variables for DSS (aggressive) in stage IIIC2-IVB")
result<-coxph(Surv(fu_date,Group1)
              ~Age+Pathology+NLR+`pre ALC`+ alpha+ `pre ALC(a2)`+ `ALC nadir(e1)`,
              data=Cx_s1)    
finalmodel<-step(result, direction = 'backward')
summary(finalmodel)
ggforest(finalmodel, data=Cx_s1, fontsize = 1.5)

## stage 3A-3C1
Cx_s1 <-  Cx_s %>% filter(stage_m=="IIIA-IIIC1")
Cx_s1 <- Cx_s1 %>% mutate(TS = Surv(fu_date,Group1==1))
out=mycph(TS~ Age+Pathology+NLR+`pre ALC`+`min ALC`+ alpha+ `pre ALC(a2)`+ `ALC nadir(e1)`, data=Cx_s1)
HRplot(out, type=2, show.CI=TRUE, main="Hazard ratios of all individual variables for DSS (aggressive) in stage IIIA-IIIC1")
result<-coxph(Surv(fu_date,Group1)
              ~Age+Pathology+NLR+`pre ALC`+ alpha+ `pre ALC(a2)`+ `ALC nadir(e1)`,
              data=Cx_s1)    
finalmodel<-step(result, direction = 'backward')
summary(finalmodel)
ggforest(finalmodel, data=Cx_s1, fontsize = 1.5)

## KM
fit=survfit(Surv(fu_date, Group1==1)~alpha, data=Cx1)
a1 <- ggsurvplot(fit, xlab="Months", ylab="Disease specific survival rate (aggressive)", pval=T, pval.size=8, fun="pct", conf.int=F, risk.table = T, size=1, linetype = "solid",
           palette = c("blue","red"), legend= "top", legend.title=expression(alpha), legend.lab=c("\u22650.08", "<0.08"), 
           break.x.by=6, xlim=c(0,36), ylim=c(0,100) )
summary(fit, times=30) # 90.3 81.8 (2.5DSS-A, alpha)

fit=survfit(Surv(fu_date, Group2==1)~alpha, data=Cx1)
a2 <- ggsurvplot(fit, xlab="Months", ylab="Disease specific survival rate (non-aggressive)", pval=T, pval.size=8,fun="pct", conf.int=F, risk.table = T, size=1, linetype = "solid",
           palette = c("blue","red"), legend= "top", legend.title=expression(alpha), legend.lab=c("\u22650.08", "<0.08"), 
           break.time.by=12, xlim=c(0,100),ylim=c(0,100))
summary(fit, times=96) #82.4 92.8 (8DSS-NA, alpha)

fit=survfit(Surv(fu_date, css==1)~alpha, data=Cx1)
a3 <- ggsurvplot(fit, xlab="Months", ylab="Disease specific survival rate", pval=T, pval.size=8,fun="pct", conf.int=F, risk.table = T, size=1, linetype = "solid",
           palette = c("blue","red"), legend= "top", legend.title=expression(alpha), legend.lab=c("\u22650.08", "<0.08"), 
           break.time.by=12, xlim=c(0,100),ylim=c(0,100))
summary(fit, times=96)# 74.4 75.9 (8DSS, alpha)

##################
fit=survfit(Surv(fu_date, Group1==1)~nlr_m, data=Cx1)
n1 <- ggsurvplot(fit, xlab="Months", ylab="Disease specific survival rate (aggressive)", pval=T, pval.size=8, fun="pct", conf.int=F, risk.table = T, size=1, linetype = "solid",
                palette = c("blue","red"), legend= "top", legend.title="NLR", legend.lab=c("<2.43", "\u22652.43"), 
                break.time.by=6, xlim=c(0,36), ylim=c(0,100))
summary(fit, times=30) # 92.4 80.1 (2.5DSS-A, nlr)

fit=survfit(Surv(fu_date, Group2==1)~nlr_m, data=Cx1)
n2 <- ggsurvplot(fit, xlab="Months", ylab="Disease specific survival rate (non-aggressive)", pval=T, pval.size=8, fun="pct", conf.int=F, risk.table = T, size=1, linetype = "solid",
           palette = c("blue","red"), legend= "top", legend.title="NLR", legend.lab=c("<2.43", "\u22652.43"), 
           break.time.by=12, xlim=c(0,100), ylim=c(0,100))
summary(fit, times=96) #88.7 85.5 (8DSS-NA , nlr)

fit=survfit(Surv(fu_date, css==1)~nlr_m, data=Cx1)
n3 <- ggsurvplot(fit, xlab="Months", ylab="Disease specific survival rate", pval=T, pval.size=8, fun="pct", conf.int=F, risk.table = T, size=1, linetype = "solid",
           palette = c("blue","red"), legend= "top", legend.title="NLR", legend.lab=c("<2.43", "\u22652.43"), 
           break.time.by=12, xlim=c(0,100), ylim=c(0,100))
summary(fit, times=96) # 82 68.47 (8DSS)


##########
Cx1 <- Cx1 %>% mutate(risk1=ifelse(nlr_m==0 & alpha==0, 0, ifelse(nlr_m==0 & alpha==1, 1, ifelse(nlr_m==1 & alpha==0, 3, 4))))
fit=survfit(Surv(fu_date, Group1==1)~risk1, data=Cx1)
r1 <- ggsurvplot(fit, xlab="Months", ylab="Disease specific survival rate (aggressive)", pval=T, pval.size=8,fun="pct", conf.int=F, risk.table = T, size=1, linetype = c("dashed","solid","dashed","solid"),
           palette = c("dark blue","blue","dark red","red"), legend= "bottom", legend.title=expression(paste("NLR & ", alpha)),
           legend.lab=c("<2.43 and \u22650.08", "<2.43 and <0.08", "\u22652.43 and \u22650.08",  "\u22652.43 and <0.08" ), 
           break.time.by=6, xlim=c(0,36), risk.table.height=0.3, ylim=c(0,100))
summary(fit, times=30)#93.9, 90.8, p=0.45 / 86.7, 73 p=0.028


fit=survfit(Surv(fu_date, Group2==1)~risk1, data=Cx1)
r2 <- ggsurvplot(fit, xlab="Months", ylab="Disease specific survival rate (non-aggressive)", pval=T, pval.size=8,fun="pct", conf.int=F, risk.table = T, size=1, linetype = c("dashed","solid","dashed","solid"),
           palette = c("dark blue","blue","dark red","red"), legend= "bottom", legend.title=expression(paste("NLR & ", alpha)),
           legend.lab=c("<2.43 and \u22650.08", "<2.43 and <0.08", "\u22652.43 and \u22650.08",  "\u22652.43 and <0.08" ), 
           break.time.by=12, xlim=c(0,100), risk.table.height=0.3, ylim=c(0,100))
summary(fit, times=96) # 86.2, 91.3, p=0.21 /78.5, 94.8 p=0.046 (8DSS-NA)
r2

fit=survfit(Surv(fu_date, css==1)~risk1, data=Cx1)
r3 <- ggsurvplot(fit, xlab="Months", ylab="Disease specific survival rate", pval=T, pval.size=8,fun="pct", conf.int=F, risk.table = T, size=1, linetype = c("dashed","solid","dashed","solid"),
           palette = c("dark blue","blue","dark red","red"), legend= "bottom", legend.title=expression(paste("NLR & ", alpha)),
           legend.lab=c("<2.43 and \u22650.08", "<2.43 and <0.08", "\u22652.43 and \u22650.08",  "\u22652.43 and <0.08" ), 
           break.time.by=12, xlim=c(0,100), risk.table.height=0.3, ylim=c(0,100))
summary(fit, times=96) # 81, 83 p=0.78 //68.1, 69.2 p=0.5 (8DSS)
r3

## detailed p vlaues
test <- Cx1 %>% filter(risk1>1)
fit=survfit(Surv(fu_date, Group1==1)~risk1, data=test)
ggsurvplot(fit, xlab="Months", ylab="Disease specific survival rate (aggressive)", pval=T, pval.size=8,fun="pct", conf.int=F, risk.table = T, size=1)
summary(fit, times=30)
fit=survfit(Surv(fu_date, Group2==1)~risk1, data=test)
ggsurvplot(fit, xlab="Months", ylab="Disease specific survival rate (aggressive)", pval=T, pval.size=8,fun="pct", conf.int=F, risk.table = T, size=1)
summary(fit, times=60)
fit=survfit(Surv(fu_date, css==1)~risk1, data=test)
ggsurvplot(fit, xlab="Months", ylab="Disease specific survival rate (aggressive)", pval=T, pval.size=8,fun="pct", conf.int=F, risk.table = T, size=1)
summary(fit, times=96)

test <- Cx1 %>% filter(risk1<2)
fit=survfit(Surv(fu_date, Group1==1)~risk1, data=test)
ggsurvplot(fit, xlab="Months", ylab="Disease specific survival rate (aggressive)", pval=T, pval.size=8,fun="pct", conf.int=F, risk.table = T, size=1)
summary(fit, times=30)
fit=survfit(Surv(fu_date, Group2==1)~risk1, data=test)
ggsurvplot(fit, xlab="Months", ylab="Disease specific survival rate (aggressive)", pval=T, pval.size=8,fun="pct", conf.int=F, risk.table = T, size=1)
summary(fit, times=60)
fit=survfit(Surv(fu_date, css==1)~risk1, data=test)
ggsurvplot(fit, xlab="Months", ylab="Disease specific survival rate (aggressive)", pval=T, pval.size=8,fun="pct", conf.int=F, risk.table = T, size=1)
summary(fit, times=96)


### table
summary(Cx1)
Cx1$rtf <- ifelse(Cx1$rtf=="L5", "pelvis", "EFRT")
out=mytable(alpha~age_m+pathology+stage_m+duration+rtf+EQD2+ALC0_m+min_alc_m+nlr_m+a2_m+e_m+recur1+Group1+Group2+css, data=Cx1, digits=2, method=3)
mycsv(out, file="table_alpha.csv")
mytable(~age_m+pathology+stage_m+duration+rtf+EQD2+ALC0_m+min_alc_m+nlr_m+a2_m+e_m+recur1+Group1+Group2+css, data=Cx1, digits=2, method=3)
out=mytable(nlr_m~age_m+pathology+stage_m+ALC0_m+min_alc_m+a2_m+alpha+e_m, data=Cx1, digits=2, method=3)
mycsv(out, file="table_nlr.csv")
mytable(risk1~age_m+pathology+stage_m+duration+rtf+EQD2+ALC0_m+min_alc_m+nlr_m+a2_m+e_m+recur1+Group1+Group2+css, data=Cx1, digits=2, method=3)

