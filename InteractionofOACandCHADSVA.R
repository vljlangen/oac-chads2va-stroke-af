# =============================================================================
# Interaction of OAC and CHA2DS2-VA score — analysis code
# For use with data: IncidentCohort.sav (place in working directory)
# Requires: R with packages haven, dplyr, survival, knitr, gridExtra, Epi,
#           cohorttools, ggplot2, tidyr; optional boot for NNT CIs
# =============================================================================

library(haven)
library(dplyr)
library(survival)
library(knitr)
library(gridExtra)
library(Epi)
library(cohorttools)
library(ggplot2)
library(tidyr)

### All OACs ###
Ranalyysit.dt <- data.frame(read_sav("IncidentCohort.sav"))
# Original data variable names (from IncidentCohort.sav) — English meanings:
#   DateISorLoppuOrDeath = date of IS or end of follow-up or death
#   LastAKdateplus120days = last anticoagulation (AK) date + 120 days
#   DateISorLoppuOrDeath_taiOACloppu120pv = earlier of IS/end/death vs OAC end + 120 days
#   VKA_aloitus = VKA start; ostodate = purchase date
#   kuolpvmSPSSdate = date of death; SukupuoliBin = sex
#   Incometertiles = income tertiles; Koulutus1perus2toinen3korkeinUusi = education level
head(Ranalyysit.dt)

# Ensure both columns are in Date format (if not already)
Ranalyysit.dt$DateISorLoppuOrDeath <- as.Date(Ranalyysit.dt$DateISorLoppuOrDeath)
Ranalyysit.dt$LastAKdateplus120days <- as.Date(Ranalyysit.dt$LastAKdateplus120days)

# Create the new variable: the earlier of the two dates
Ranalyysit.dt$DateISorLoppuOrDeath_taiOACloppu120pv <- pmin(
  Ranalyysit.dt$DateISorLoppuOrDeath, 
  Ranalyysit.dt$LastAKdateplus120days, 
  na.rm = TRUE
)

# Ensure the dates are properly formatted
Ranalyysit.dt$ISfirstdateAfterCohortfirsteverandrecurrent <- as.Date(Ranalyysit.dt$ISfirstdateAfterCohortfirsteverandrecurrent)
Ranalyysit.dt$DateISorLoppuOrDeath_taiOACloppu120pv <- as.Date(Ranalyysit.dt$DateISorLoppuOrDeath_taiOACloppu120pv)

# Create the new variable with the given condition
Ranalyysit.dt$ISaftercohortallBUT_BEFORE_OAC_END120pv <- ifelse(
  Ranalyysit.dt$ISaftercohortall == 1 &
    Ranalyysit.dt$ISfirstdateAfterCohortfirsteverandrecurrent <= Ranalyysit.dt$DateISorLoppuOrDeath_taiOACloppu120pv,
  1,
  0
)

# Ensure both variables are Date objects
Ranalyysit.dt$LastAKdateplus120days <- as.Date(Ranalyysit.dt$LastAKdateplus120days)
Ranalyysit.dt$CohortEntryDate <- as.Date(Ranalyysit.dt$CohortEntryDate)



Ranalyysit.dt$LastAKdateplus120daysNEWnodatesbeforecohortentry <- if_else(
  Ranalyysit.dt$LastAKdateplus120days >= Ranalyysit.dt$CohortEntryDate,
  Ranalyysit.dt$LastAKdateplus120days,
  as.Date(NA)
)


# Show count of 1s and 0s
table(Ranalyysit.dt$ISaftercohortallBUT_BEFORE_OAC_END120pv)


# Create a new Date variable and fill it only where VKA_aloitus == 1 (VKA start)
Ranalyysit.dt$ostodateFirstVKA <- as.Date(NA)  # Initialize with NA
Ranalyysit.dt$ostodateFirstVKA[Ranalyysit.dt$VKA_aloitus == 1] <- Ranalyysit.dt$ostodate[Ranalyysit.dt$VKA_aloitus == 1]  # ostodate = purchase date


with(Ranalyysit.dt,range(cal.yr(DateISorLoppuOrDeath_taiOACloppu120pv,format = "%Y-%m-%d") ))
with(Ranalyysit.dt,range(cal.yr(kuolpvmSPSSdate,format = "%Y-%m-%d"),na.rm = TRUE))  # kuolpvmSPSSdate = date of death

apu1<-with(Ranalyysit.dt,factor(rep("SOF",nrow(Ranalyysit.dt)),levels=c("SOF","EOF","IS")))
apu2<-with(Ranalyysit.dt,factor(ifelse(ISaftercohortallBUT_BEFORE_OAC_END120pv==1,"IS","EOF"),levels=c("SOF","EOF","IS")))

table(apu1,useNA = "always")
table(apu2,useNA = "always")

tmp.dt<-Lexis(entry=list(age=Age,
                         fu=0,
                         per=cal.yr(CohortEntryDate,format = "%Y-%m-%d")),
              duration = cal.yr(DateISorLoppuOrDeath_taiOACloppu120pv,format = "%Y-%m-%d")-
                cal.yr(CohortEntryDate,format = "%Y-%m-%d"),
              entry.status = apu1,
              exit.status = apu2,
              id=SID,
              data=Ranalyysit.dt)
# NOTE: Dropping  249  rows with duration of follow up < tol

summary(tmp.dt)
timeScales(tmp.dt)

tmp.dt1<-cutLexis(data = tmp.dt,
                  cut=cal.yr(tmp.dt$ostodate,format = "%Y-%m-%d"),
                  timescale = "per",
                  new.state = "AK",
                  new.scale = "AK.Start"
)
summary(tmp.dt1)


tmp.dt2<-cutLexis(data = tmp.dt1,
                  cut=cal.yr(tmp.dt1$LastAKdateplus120daysNEWnodatesbeforecohortentry,format = "%Y-%m-%d"),
                  timescale = "per",
                  new.state = "AK.quitted",
                  new.scale = "AK.quit"
)

summary(tmp.dt2)


boxesLx(tmp.dt2,show.persons=FALSE)






# Split by calendar year (per)
range(tmp.dt2$per)
tmp.dt3<-splitLexis(lex=tmp.dt2,time.scale = "per",breaks = c(2009,2011,2013,2015,2017))


# Year as factor
tmp.dt3$per.c<-timeBand(lex = tmp.dt3,time.scale = "per",type="factor")
# Year numeric
tmp.dt3$per.num<-with(tmp.dt3,per+lex.dur/2)

# Age as factor, age in middle of time slice
apu<-with(tmp.dt3,cut(age+lex.dur/2,c(0,40,45,50,55,60,65,70,75,80,85,90,95,100,Inf)))
tmp.dt3$age.c<-apu
table(apu,tmp.dt3$lex.Xst)

# Variable names
names(tmp.dt3)

# Sex (original data column SukupuoliBin)
apu<-factor(unclass(tmp.dt3$SukupuoliBin),levels=0:1,labels =c("female","male"))
tmp.dt3$sex<-apu

tmp.dt3$sex <- relevel(tmp.dt3$sex, ref = "male")

# Other variables
apu<-unclass(tmp.dt3$HyperlipidemiaBOAC);table(apu)
tmp.dt3$Hyperlipidemia<-factor(apu,levels=0:1,labels =c("no","yes"))

apu<-unclass(tmp.dt3$HypertensionBOAC);table(apu)
tmp.dt3$Hypertension<-factor(apu,levels=0:1,labels =c("no","yes"))

apu<-unclass(tmp.dt3$DiabetesBOAC);table(apu)
tmp.dt3$Diabetes<-factor(apu,levels=0:1,labels =c("no","yes"))

apu<-unclass(tmp.dt3$CongestiveHeartFailureBOAC);table(apu)
tmp.dt3$HF<-factor(apu,levels=0:1,labels =c("no","yes"))

apu<-unclass(tmp.dt3$AbnormalLiverFunctionBOAC);table(apu)
tmp.dt3$LiverVT<-factor(apu,levels=0:1,labels =c("no","yes"))

apu<-unclass(tmp.dt3$AbnormalRenalFunctionBOAC);table(apu)
tmp.dt3$RenalVT<-factor(apu,levels=0:1,labels =c("no","yes"))

apu<-unclass(tmp.dt3$BleedingsBOAC);table(apu)
tmp.dt3$Bleeding<-factor(apu,levels=0:1,labels =c("no","yes"))



apu<-unclass(tmp.dt3$AlcoholBOAC);table(apu)  # AlcoholBOAC = alcohol use (from data)
tmp.dt3$alcohol<-factor(apu,levels=0:1,labels =c("no","yes"))

apu<-unclass(tmp.dt3$AnyVascularDiseaseBOAC);table(apu)
tmp.dt3$anyvasc<-factor(apu,levels=0:1,labels =c("no","yes"))

apu<-unclass(tmp.dt3$CancerBeforeOrAtCohort);table(apu)  # Cancer before or at cohort entry
tmp.dt3$cancer<-factor(apu,levels=0:1,labels =c("no","yes"))

apu<-unclass(tmp.dt3$DementiaBOAC);table(apu)
tmp.dt3$dementia<-factor(apu,levels=0:1,labels =c("no","yes"))

apu<-unclass(tmp.dt3$PsychiatricDiseaseBOAC);table(apu)  # Psychiatric disease (from data)
tmp.dt3$psychiatric<-factor(apu,levels=0:1,labels =c("no","yes"))

apu<-unclass(tmp.dt3$Incometertiles);table(apu)  # Incometertiles = income tertiles (from data)
tmp.dt3$income_tertile<-factor(apu,levels=1:3,labels =c("low","mid", "high"))

apu<-unclass(tmp.dt3$IschemicStrokeBOAC);table(apu)
tmp.dt3$stroke<-factor(apu,levels=0:1,labels =c("no","yes"))

apu<-unclass(tmp.dt3$Koulutus1perus2toinen3korkeinUusi);table(apu)  # Education level (from data)
tmp.dt3$education<-factor(apu,levels=1:3,labels =c("low","mid", "high"))

apu<-unclass(tmp.dt3$low2moderate3high);table(apu)
tmp.dt3$low2moderate3high<-factor(apu,levels=1:3,labels =c("low","mid", "high"))

apu<-unclass(tmp.dt3$CoronaryHeartDiseaseBOAC);table(apu)
tmp.dt3$MCC<-factor(apu,levels=0:1,labels =c("no","yes"))

apu<-unclass(tmp.dt3$MyocardialInfarctionBOAC);table(apu)
tmp.dt3$MI<-factor(apu,levels=0:1,labels =c("no","yes"))

tmp.dt3 <- subset(tmp.dt3, CHADSVAbaseline > 0)

apu<-unclass(tmp.dt3$CHADSVAbaseline);table(apu)
tmp.dt3$CHADSVAbaseline<-factor(apu,levels=1:8)



#------------------------------
#------------------------------
# Rate table
#------------------------------

tmp.rt<-mkratetable(Surv(lex.dur,lex.Xst=="IS")~Relevel(lex.Cst,list(1,2:3,4))+CHADSVAbaseline,
                    data=tmp.dt3,scale=100,add.RR = TRUE)

sink(file="RateTableCHADSVA.html")
knitr::kable(tmp.rt,caption="Rate table (1/1000 person years)",format="html",digits=3)
sink()

#------------------------------
# Rate table: time WITH AK
#------------------------------
tmp.dt.AK <- subset(tmp.dt3, lex.Cst %in% c("AK"))  # keep only time in AK state

rate_AK <- mkratetable(
  Surv(lex.dur, lex.Xst == "IS") ~  CHADSVAbaseline,
  data = tmp.dt.AK,
  scale = 100,       # rates per 1000 person-years
  add.RR = TRUE
)

#------------------------------
# Rate table: time WITHOUT AK
#------------------------------
tmp.dt.noAK <- subset(tmp.dt3, lex.Cst %in% c("SOF"))  # keep only time off AK

rate_noAK <- mkratetable(
  Surv(lex.dur, lex.Xst == "IS") ~  CHADSVAbaseline,
  data = tmp.dt.noAK,
  scale = 100,       # rates per 1000 person-years
  add.RR = TRUE
)

sink(file = "RateTableIS_yesANDnoAK.html")
knitr::kable(rate_noAK, caption = "Rate table for time OFF anticoagulation (per 100 PY)", format = "html", digits = 3)
knitr::kable(rate_AK, caption = "Rate table for time ON anticoagulation (per 100 PY)", format = "html", digits = 3)
sink()


####regressions###

tmp.m1 <- glm(cbind(lex.Xst=="IS",lex.dur) ~ Relevel(lex.Cst,list(1,2:3,4))+CHADSVAbaseline,
              data=tmp.dt3,family="poisreg")

sink(file="AdjustedAnalysisCHADSVA.html")

knitr::kable(round(ci.exp(tmp.m1),2),caption="IRR from Poisson regression model",format="html")

sink()


tmp.dt3$lex.Cst.uusi<-droplevels(tmp.dt3$lex.Cst)

mkratetable(Surv(lex.dur,lex.Xst=="IS")~ lex.Cst.uusi+sex+CHADSVAbaseline,data=tmp.dt3,scale=100)

###with all variables######
tmp.m1 <- glm(cbind(lex.Xst=="IS",lex.dur) ~ anyvasc+income_tertile+stroke+Hypertension+RenalVT+LiverVT+HF+alcohol+Bleeding+Diabetes+cancer+dementia+psychiatric+Hyperlipidemia+age.c+lex.Cst.uusi+sex+CHADSVAbaseline,
              data=tmp.dt3,
              family="poisreg",maxit = 200)


tmp.m2 <- update(tmp.m1, ~ .+CHADSVAbaseline:lex.Cst.uusi)

anova(tmp.m1,tmp.m2,test="Chisq")

round(ci.exp(tmp.m2),3)


tmp.df<-with(tmp.dt3,
             expand.grid( age.c=sort(unique(age.c))[1],
                          sex=sort(unique(sex))[1],
                          income_tertile=sort(unique(income_tertile))[1],
                          anyvasc=sort(unique(anyvasc))[1],
                          Hyperlipidemia=sort(unique(Hyperlipidemia))[1],
                          psychiatric=sort(unique(psychiatric))[1],
                          dementia=sort(unique(dementia))[1],
                          cancer=sort(unique(cancer))[1],
                          Diabetes=sort(unique(Diabetes))[1],
                          alcohol=sort(unique(alcohol))[1],
                          Bleeding=sort(unique(Bleeding))[1],
                          HF=sort(unique(HF))[1],
                          LiverVT=sort(unique(LiverVT))[1],
                          RenalVT=sort(unique(RenalVT))[1],
                          Hypertension=sort(unique(Hypertension))[1],
                          stroke=sort(unique(stroke))[1],
                          CHADSVAbaseline=sort(unique(CHADSVAbaseline)),
                          
                          
                          lex.Cst.uusi=sort(unique(lex.Cst.uusi))))

tmp.mtr<-model.matrix(~anyvasc+stroke+income_tertile+Hypertension+RenalVT+LiverVT+HF+alcohol+Bleeding+Diabetes+cancer+dementia+psychiatric+Hyperlipidemia+age.c+lex.Cst.uusi+sex + CHADSVAbaseline + CHADSVAbaseline:lex.Cst.uusi,data=tmp.df)


apu.all<-as.data.frame(ci.exp(tmp.m2,ctr.mat=(tmp.mtr[-(1:8),]-tmp.mtr[(1:8),])))

names(apu.all)<-c("MRR","MRR.lo","MRR.hi");

apu.all$CHADSVAbaseline<-c("1","2", "3", "4","5","6","7","8")


fig_fullyadjusted <- ggplot(apu.all, aes(x = CHADSVAbaseline, y = MRR, group = 1)) +  # Specify 'group = 1' to connect all points with a single line
  geom_pointrange(aes(ymin = MRR.lo, ymax = MRR.hi), fatten = 2.5, lwd = 1.5, color = "#0D2C54") +
  geom_line(color = "#0D2C54") +  # Add a line connecting the points
  geom_hline(yintercept = 1, linetype = 3, color = "red") +
labs(y = "Adjusted incidence rate ratio", x = "CHA2DS2-VA score", title = " ") +
  ggtitle("C") + 
  theme_bw() +
  theme(
    panel.grid.major = element_line(color = "#D3D3D3", linetype = "dotted"),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5),
    axis.title = element_text(size = 12),
    plot.title = element_text(size = 14, face = "bold")
  )

fig_fullyadjusted




###unadjusted######
tmp.m1unadj <- glm(cbind(lex.Xst=="IS",lex.dur) ~ lex.Cst.uusi+CHADSVAbaseline,
              data=tmp.dt3,
              family="poisreg",maxit = 200)


tmp.m2unadj <- update(tmp.m1unadj, ~ .+CHADSVAbaseline:lex.Cst.uusi)

anova(tmp.m1unadj,tmp.m2unadj,test="Chisq")

round(ci.exp(tmp.m2unadj),3)


tmp.dfunadj<-with(tmp.dt3,
             expand.grid( 
                          CHADSVAbaseline=sort(unique(CHADSVAbaseline)),
                          
                          
                          lex.Cst.uusi=sort(unique(lex.Cst.uusi))))

tmp.mtrunadj<-model.matrix(~lex.Cst.uusi+ CHADSVAbaseline + CHADSVAbaseline:lex.Cst.uusi,data=tmp.dfunadj)


apu.allunadj<-as.data.frame(ci.exp(tmp.m2unadj,ctr.mat=(tmp.mtrunadj[-(1:8),]-tmp.mtrunadj[(1:8),])))

names(apu.allunadj)<-c("MRR","MRR.lo","MRR.hi");

apu.allunadj$CHADSVAbaseline<-c("1","2", "3", "4","5","6","7","8")


figunadjusted <- ggplot(apu.allunadj, aes(x = CHADSVAbaseline, y = MRR, group = 1)) +  # Specify 'group = 1' to connect all points with a single line
  geom_pointrange(aes(ymin = MRR.lo, ymax = MRR.hi), fatten = 2.5, lwd = 1.5, color = "#0D2C54") +
  geom_line(color = "#0D2C54") +  # Add a line connecting the points
  geom_hline(yintercept = 1, linetype = 3, color = "red") +
  labs(y = "Adjusted incidence rate ratio", x = "CHA2DS2-VA score", title = " ") +
  ggtitle("C") + 
  theme_bw() +
  theme(
    panel.grid.major = element_line(color = "#D3D3D3", linetype = "dotted"),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5),
    axis.title = element_text(size = 12),
    plot.title = element_text(size = 14, face = "bold")
  )

figunadjusted



# Adjusted plot
fig_fullyadjusted <- ggplot(apu.all, aes(x = CHADSVAbaseline, y = MRR, group = 1)) +
  geom_pointrange(aes(ymin = MRR.lo, ymax = MRR.hi),
                  fatten = 2.5, lwd = 1.2, color = "#0D2C54") +
  geom_line(color = "#0D2C54", linewidth = 0.8) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  labs(
    x = "CHA₂DS₂-VA score",
    y = "Incidence rate ratio",
    title = "Adjusted analysis"
  ) +
  scale_y_continuous(
    limits = c(0.2, 1.2),
    breaks = seq(0.2, 1.2, 0.2)
  ) +
  theme_bw(base_size = 12, base_family = "sans") +
  theme(
    legend.position = "none",
    panel.grid.major = element_line(color = "grey80", linetype = "dotted"),
    panel.grid.minor = element_blank(),
    axis.title = element_text(size = 13, face = "bold"),
    axis.text = element_text(size = 11),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    panel.border = element_rect(color = "black", size = 0.6)
  )

# Unadjusted plot
figunadjusted <- ggplot(apu.allunadj, aes(x = CHADSVAbaseline, y = MRR, group = 1)) +
  geom_pointrange(aes(ymin = MRR.lo, ymax = MRR.hi),
                  fatten = 2.5, lwd = 1.2, color = "#0D2C54") +
  geom_line(color = "#0D2C54", linewidth = 0.8) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  labs(
    x = "CHA₂DS₂-VA score",
    y = "Incidence rate ratio",
    title = "Unadjusted analysis"
  ) +
  scale_y_continuous(
    limits = c(0.2, 1.2),
    breaks = seq(0.2, 1.2, 0.2)
  ) +
  theme_bw(base_size = 12, base_family = "sans") +
  theme(
    legend.position = "none",
    panel.grid.major = element_line(color = "grey80", linetype = "dotted"),
    panel.grid.minor = element_blank(),
    axis.title = element_text(size = 13, face = "bold"),
    axis.text = element_text(size = 11),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    panel.border = element_rect(color = "black", size = 0.6)
  )

# Combine into one figure with 2 panels
grid.arrange(figunadjusted, fig_fullyadjusted, ncol = 2)

# Save as JPEG
jpeg("Combined_MRR_Figures.jpg", width = 12, height = 6, units = "in", res = 300)  # 300 dpi for publication quality
grid.arrange(figunadjusted, fig_fullyadjusted, ncol = 2)
dev.off()

# Save as PDF
pdf("Combined_MRR_Figures.pdf", width = 12, height = 6)
grid.arrange(figunadjusted, fig_fullyadjusted, ncol = 2)
dev.off()

######## NUMBER NEEDED TO TREAT #########

# 1. Merge datasets as before
dat <- apu.all %>%
  mutate(CHADSVAbaseline = as.numeric(as.character(CHADSVAbaseline))) %>%
  inner_join(
    as.data.frame(rate_noAK) %>%
      transmute(
        CHADSVAbaseline = as.numeric(as.character(CHADSVAbaseline)),
        rate = rate,   # events per 100 patient-years
        low  = low,
        high = high
      ),
    by = "CHADSVAbaseline"
  ) %>%
  mutate(
    logIRR    = log(MRR),
    SE_logIRR = (log(MRR.hi) - log(MRR.lo)) / (2*1.96),
    SE_rate   = (high - low) / (2*1.96)
  )

# 2. Bootstrap function with 10-year scaling
boot_one <- function(row, R = 5000) {
  
  IRR_star <- rlnorm(R, meanlog = row$logIRR, sdlog = row$SE_logIRR)
  rate_star <- rnorm(R, mean = row$rate, sd = row$SE_rate)
  rate_star[rate_star < 0] <- 0
  
  # ARR scaled to 10-year treatment per patient
  ARR_star <- (rate_star * 10 / 100) * (1 - IRR_star)
  NNT_star <- 1 / ARR_star
  
  tibble(
    ARR_lo  = quantile(ARR_star, 0.025, na.rm = TRUE),
    ARR_hi  = quantile(ARR_star, 0.975, na.rm = TRUE),
    NNT_lo  = quantile(NNT_star, 0.025, na.rm = TRUE),
    NNT_hi  = quantile(NNT_star, 0.975, na.rm = TRUE)
  )
}

# 3. Apply bootstrap per row
boot_results <- dat %>%
  rowwise() %>%
  mutate(
    ARR  = (rate * 10 / 100) * (1 - MRR),  # 10-year ARR per patient
    NNT  = 1 / ARR,
    boot = list(boot_one(cur_data()))
  ) %>%
  unnest(boot) %>%
  ungroup()

# 4. Final table
NNT_table <- boot_results %>%
  select(
    CHADSVAbaseline,
    rate, low, high,
    MRR, MRR.lo, MRR.hi,
    ARR, ARR_lo, ARR_hi,
    NNT, NNT_lo, NNT_hi
  )


# 1. Merge datasets as before
dat <- apu.allunadj %>%
  mutate(CHADSVAbaseline = as.numeric(as.character(CHADSVAbaseline))) %>%
  inner_join(
    as.data.frame(rate_noAK) %>%
      transmute(
        CHADSVAbaseline = as.numeric(as.character(CHADSVAbaseline)),
        rate = rate,   # events per 100 patient-years
        low  = low,
        high = high
      ),
    by = "CHADSVAbaseline"
  ) %>%
  mutate(
    logIRR    = log(MRR),
    SE_logIRR = (log(MRR.hi) - log(MRR.lo)) / (2*1.96),
    SE_rate   = (high - low) / (2*1.96)
  )

# 2. Bootstrap function with 10-year scaling
boot_one <- function(row, R = 5000) {
  
  IRR_star <- rlnorm(R, meanlog = row$logIRR, sdlog = row$SE_logIRR)
  rate_star <- rnorm(R, mean = row$rate, sd = row$SE_rate)
  rate_star[rate_star < 0] <- 0
  
  # ARR scaled to 10-year treatment per patient
  ARR_star <- (rate_star * 10 / 100) * (1 - IRR_star)
  NNT_star <- 1 / ARR_star
  
  tibble(
    ARR_lo  = quantile(ARR_star, 0.025, na.rm = TRUE),
    ARR_hi  = quantile(ARR_star, 0.975, na.rm = TRUE),
    NNT_lo  = quantile(NNT_star, 0.025, na.rm = TRUE),
    NNT_hi  = quantile(NNT_star, 0.975, na.rm = TRUE)
  )
}

# 3. Apply bootstrap per row
boot_results <- dat %>%
  rowwise() %>%
  mutate(
    ARR  = (rate * 10 / 100) * (1 - MRR),  # 10-year ARR per patient
    NNT  = 1 / ARR,
    boot = list(boot_one(cur_data()))
  ) %>%
  unnest(boot) %>%
  ungroup()

# 4. Final table
NNT_tableunadj <- boot_results %>%
  select(
    CHADSVAbaseline,
    rate, low, high,
    MRR, MRR.lo, MRR.hi,
    ARR, ARR_lo, ARR_hi,
    NNT, NNT_lo, NNT_hi
  )


print(NNT_tableunadj)

print(NNT_table)


##### ABSOLUTE RATE DIFFERENCE #####

tmp.dt4<-tmp.dt3

# 1) Fit Poisson model with the exact same covariates as your IRR model
mod <- glm(
  I(lex.Xst == "IS") ~ CHADSVAbaseline * lex.Cst.uusi 
  
  #+per.c+ RenalVT + income_tertile + LiverVT + alcohol + Bleeding +
  #cancer + dementia + psychiatric + Hyperlipidemia + sex 
  
  ,
  family = poisson(),
  offset = log(lex.dur),
  data = tmp.dt4
)

# 2) Build newdata for predictions
cov_ref <- tmp.dt4 %>%
  summarise(across(where(is.factor), ~ levels(.x)[1]),
            across(where(is.numeric), ~ mean(.x, na.rm = TRUE)))

newdata <- expand.grid(
  CHADSVAbaseline = levels(tmp.dt4$CHADSVAbaseline),
  lex.Cst.uusi = levels(tmp.dt4$lex.Cst.uusi)
)

# Add constant covariates from cov_ref
for (v in names(cov_ref)) {
  if (!v %in% names(newdata)) {
    newdata[[v]] <- cov_ref[[v]]
  }
}

# Match factor levels
for (v in names(newdata)) {
  if (is.factor(tmp.dt4[[v]])) {
    newdata[[v]] <- factor(newdata[[v]], levels = levels(tmp.dt4[[v]]))
  }
}

# 3) Predict adjusted rates (per 100 PY)
pred <- predict(mod, newdata, type = "link", se.fit = TRUE)

newdata <- newdata %>%
  mutate(
    fit_rate = exp(pred$fit) * 100,
    lower = exp(pred$fit - 1.96 * pred$se.fit) * 100,
    upper = exp(pred$fit + 1.96 * pred$se.fit) * 100
  )

# 4) Calculate ARD vs reference treatment
ref_lev <- levels(tmp.dt4$lex.Cst.uusi)[1]

ref_rates <- newdata %>%
  filter(lex.Cst.uusi == ref_lev) %>%
  select(CHADSVAbaseline, fit_rate_ref = fit_rate)

ard_df <- newdata %>%
  left_join(ref_rates, by = "CHADSVAbaseline") %>%
  mutate(
    ARD = fit_rate - fit_rate_ref,
    ARD_lower = lower - fit_rate_ref,
    ARD_upper = upper - fit_rate_ref
  )

# 5) Plot ARDs
# 5) Plot ARDs with same style as IRR plots
# ARD plot styled like the unadjusted MRR plot
ARDtype1 <- ggplot(
  ard_df %>% filter(lex.Cst.uusi != ref_lev),
  aes(x = CHADSVAbaseline, y = ARD, group = lex.Cst.uusi)
) +
  geom_pointrange(
    aes(ymin = ARD_lower, ymax = ARD_upper),
    fatten = 2.5, lwd = 1.2, color = "#0D2C54", 
    position = position_dodge(width = 0.4)
  ) +
  geom_line(color = "#0D2C54", linewidth = 0.8, position = position_dodge(width = 0.4)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  labs(
    x = "CHA₂DS₂-VA score",
    y = "Rate difference\n(events per 100 PY)",
    title = "Adjusted absolute rate difference"
  ) +
  scale_y_continuous(
    limits = c(-3, 0),  # adjust as needed
    breaks = seq(-3, 0, 0.5)
  ) +
  theme_bw(base_size = 12, base_family = "sans") +
  theme(
    legend.position = "none",
    panel.grid.major = element_line(color = "grey80", linetype = "dotted"),
    panel.grid.minor = element_blank(),
    axis.title = element_text(size = 13, face = "bold"),
    axis.text = element_text(size = 11),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    panel.border = element_rect(color = "black", size = 0.6)
  )


ARDtype1


# Save ARDtype1 as JPEG
jpeg("ARDtype1_Figure.jpg", width = 6, height = 6, units = "in", res = 300)  # 6x6 inches, 300 dpi
print(ARDtype1)
dev.off()

###VKA and NOAC###



with(Ranalyysit.dt,range(cal.yr(DateISorLoppuOrDeath_taiOACloppu120pv,format = "%Y-%m-%d") ))
with(Ranalyysit.dt,range(cal.yr(kuolpvmSPSSdate,format = "%Y-%m-%d"),na.rm = TRUE))  # kuolpvmSPSSdate = date of death

apu1<-with(Ranalyysit.dt,factor(rep("SOF",nrow(Ranalyysit.dt)),levels=c("SOF","EOF","IS")))
apu2<-with(Ranalyysit.dt,factor(ifelse(ISaftercohortallBUT_BEFORE_OAC_END120pv==1,"IS","EOF"),levels=c("SOF","EOF","IS")))

table(apu1,useNA = "always")
table(apu2,useNA = "always")

tmp.dt<-Lexis(entry=list(age=Age,
                         fu=0,
                         per=cal.yr(CohortEntryDate,format = "%Y-%m-%d")),
              duration = cal.yr(DateISorLoppuOrDeath_taiOACloppu120pv,format = "%Y-%m-%d")-
                cal.yr(CohortEntryDate,format = "%Y-%m-%d"),
              entry.status = apu1,
              exit.status = apu2,
              id=SID,
              data=Ranalyysit.dt)
# NOTE: Dropping  249  rows with duration of follow up < tol

summary(tmp.dt)
timeScales(tmp.dt)

tmp.dt1<-cutLexis(data = tmp.dt,
                  cut=cal.yr(tmp.dt$ostodateFirstNoac,format = "%Y-%m-%d"),
                  timescale = "per",
                  new.state = "NOAC",
                  new.scale = "NOAC.Start"
)
summary(tmp.dt1)

tmp.dt2<-cutLexis(data = tmp.dt1,
                  cut=cal.yr(tmp.dt1$ostodateFirstVKA,format = "%Y-%m-%d"),
                  timescale = "per",
                  new.state = "VKA",
                  new.scale = "VKA.Start"
)
summary(tmp.dt2)

tmp.dt3<-cutLexis(data = tmp.dt2,
                  cut=cal.yr(tmp.dt2$LastAKdateplus120daysNEWnodatesbeforecohortentry,format = "%Y-%m-%d"),
                  timescale = "per",
                  new.state = "AK.quitted",
                  new.scale = "AK.quit"
)

summary(tmp.dt3)

boxesLx(tmp.dt2,show.persons=FALSE)

tmp.dt3<-tmp.dt2

# Split by calendar year (per)
range(tmp.dt3$per)
tmp.dt4<-splitLexis(lex=tmp.dt3,time.scale = "per",breaks = c(2009,2011,2013,2015,2017))


# Year as factor
tmp.dt4$per.c<-timeBand(lex = tmp.dt4,time.scale = "per",type="factor")
# Year numeric
tmp.dt4$per.num<-with(tmp.dt4,per+lex.dur/2)

# Age as factor, age in middle of time slice
apu<-with(tmp.dt4,cut(age+lex.dur/2,c(0,40,45,50,55,60,65,70,75,80,85,90,95,100,Inf)))
tmp.dt4$age.c<-apu
table(apu,tmp.dt4$lex.Xst)

# Variable names
names(tmp.dt4)

# Sex (original data column SukupuoliBin)
apu<-factor(unclass(tmp.dt4$SukupuoliBin),levels=0:1,labels =c("female","male"))
tmp.dt4$sex<-apu

tmp.dt4$sex <- relevel(tmp.dt4$sex, ref = "male")

# Other variables
apu<-unclass(tmp.dt4$HyperlipidemiaBOAC);table(apu)
tmp.dt4$Hyperlipidemia<-factor(apu,levels=0:1,labels =c("no","yes"))

apu<-unclass(tmp.dt4$HypertensionBOAC);table(apu)
tmp.dt4$Hypertension<-factor(apu,levels=0:1,labels =c("no","yes"))

apu<-unclass(tmp.dt4$DiabetesBOAC);table(apu)
tmp.dt4$Diabetes<-factor(apu,levels=0:1,labels =c("no","yes"))

apu<-unclass(tmp.dt4$CongestiveHeartFailureBOAC);table(apu)
tmp.dt4$HF<-factor(apu,levels=0:1,labels =c("no","yes"))

apu<-unclass(tmp.dt4$AbnormalLiverFunctionBOAC);table(apu)
tmp.dt4$LiverVT<-factor(apu,levels=0:1,labels =c("no","yes"))

apu<-unclass(tmp.dt4$AbnormalRenalFunctionBOAC);table(apu)
tmp.dt4$RenalVT<-factor(apu,levels=0:1,labels =c("no","yes"))

apu<-unclass(tmp.dt4$BleedingsBOAC);table(apu)
tmp.dt4$Bleeding<-factor(apu,levels=0:1,labels =c("no","yes"))



apu<-unclass(tmp.dt4$AlcoholBOAC);table(apu)  # AlcoholBOAC = alcohol use (from data)
tmp.dt4$alcohol<-factor(apu,levels=0:1,labels =c("no","yes"))

apu<-unclass(tmp.dt4$AnyVascularDiseaseBOAC);table(apu)
tmp.dt4$anyvasc<-factor(apu,levels=0:1,labels =c("no","yes"))

apu<-unclass(tmp.dt4$CancerBeforeOrAtCohort);table(apu)  # Cancer before or at cohort entry
tmp.dt4$cancer<-factor(apu,levels=0:1,labels =c("no","yes"))

apu<-unclass(tmp.dt4$DementiaBOAC);table(apu)
tmp.dt4$dementia<-factor(apu,levels=0:1,labels =c("no","yes"))

apu<-unclass(tmp.dt4$PsychiatricDiseaseBOAC);table(apu)  # Psychiatric disease (from data)
tmp.dt4$psychiatric<-factor(apu,levels=0:1,labels =c("no","yes"))

apu<-unclass(tmp.dt4$Incometertiles);table(apu)  # Incometertiles = income tertiles (from data)
tmp.dt4$income_tertile<-factor(apu,levels=1:3,labels =c("low","mid", "high"))

apu<-unclass(tmp.dt4$IschemicStrokeBOAC);table(apu)
tmp.dt4$stroke<-factor(apu,levels=0:1,labels =c("no","yes"))

apu<-unclass(tmp.dt4$Koulutus1perus2toinen3korkeinUusi);table(apu)  # Education level (from data)
tmp.dt4$education<-factor(apu,levels=1:3,labels =c("low","mid", "high"))

apu<-unclass(tmp.dt4$low2moderate3high);table(apu)
tmp.dt4$low2moderate3high<-factor(apu,levels=1:3,labels =c("low","mid", "high"))

apu<-unclass(tmp.dt4$CoronaryHeartDiseaseBOAC);table(apu)
tmp.dt4$MCC<-factor(apu,levels=0:1,labels =c("no","yes"))

apu<-unclass(tmp.dt4$MyocardialInfarctionBOAC);table(apu)
tmp.dt4$MI<-factor(apu,levels=0:1,labels =c("no","yes"))

tmp.dt4 <- subset(tmp.dt4, CHADSVAbaseline > 0)

apu<-unclass(tmp.dt4$CHADSVAbaseline);table(apu)
tmp.dt4$CHADSVAbaseline<-factor(apu,levels=1:8)



#------------------------------
#------------------------------
# Rate table
#------------------------------

tmp.rt<-mkratetable(Surv(lex.dur,lex.Xst=="IS")~Relevel(lex.Cst,list(1,2,3,4))+CHADSVAbaseline,
                    data=tmp.dt4,scale=100,add.RR = TRUE)

sink(file="RateTableCHADSVA.html")
knitr::kable(tmp.rt,caption="Rate table (1/1000 person years)",format="html",digits=3)
sink()

tmp.m1 <- glm(cbind(lex.Xst=="IS",lex.dur) ~ Relevel(lex.Cst,list(1,2,3,4))+CHADSVAbaseline,
              data=tmp.dt4,family="poisreg")

sink(file="AdjustedAnalysisCHADSVANOACVKA.html")

knitr::kable(round(ci.exp(tmp.m1),2),caption="IRR from Poisson regression model",format="html")

sink()


tmp.dt4$lex.Cst.uusi<-droplevels(tmp.dt4$lex.Cst)

mkratetable(Surv(lex.dur,lex.Xst=="IS")~ lex.Cst.uusi+sex+CHADSVAbaseline,data=tmp.dt4,scale=100)

###with all variables######
tmp.m1 <- glm(cbind(lex.Xst=="IS",lex.dur) ~ anyvasc+income_tertile+stroke+Hypertension+RenalVT+LiverVT+HF+alcohol+Bleeding+Diabetes+cancer+dementia+psychiatric+Hyperlipidemia+age.c+lex.Cst.uusi+sex+CHADSVAbaseline,
              data=tmp.dt4,
              family="poisreg",maxit = 200)


tmp.m2 <- update(tmp.m1, ~ .+CHADSVAbaseline:lex.Cst.uusi)

anova(tmp.m1,tmp.m2,test="Chisq")

round(ci.exp(tmp.m2),3)


tmp.df<-with(tmp.dt4,
             expand.grid( age.c=sort(unique(age.c))[1],
                          sex=sort(unique(sex))[1],
                          income_tertile=sort(unique(income_tertile))[1],
                          anyvasc=sort(unique(anyvasc))[1],
                          Hyperlipidemia=sort(unique(Hyperlipidemia))[1],
                          psychiatric=sort(unique(psychiatric))[1],
                          dementia=sort(unique(dementia))[1],
                          cancer=sort(unique(cancer))[1],
                          Diabetes=sort(unique(Diabetes))[1],
                          alcohol=sort(unique(alcohol))[1],
                          Bleeding=sort(unique(Bleeding))[1],
                          HF=sort(unique(HF))[1],
                          LiverVT=sort(unique(LiverVT))[1],
                          RenalVT=sort(unique(RenalVT))[1],
                          Hypertension=sort(unique(Hypertension))[1],
                          stroke=sort(unique(stroke))[1],
                          CHADSVAbaseline=sort(unique(CHADSVAbaseline)),
                          
                          
                          lex.Cst.uusi=sort(unique(lex.Cst.uusi))))

tmp.mtr<-model.matrix(~anyvasc+stroke+income_tertile+Hypertension+RenalVT+LiverVT+HF+alcohol+Bleeding+Diabetes+cancer+dementia+psychiatric+Hyperlipidemia+age.c+lex.Cst.uusi+sex + CHADSVAbaseline + CHADSVAbaseline:lex.Cst.uusi,data=tmp.df)


########
# Get indices for each treatment group
idx_noac <- tmp.df$lex.Cst.uusi == "NOAC"
idx_sof  <- tmp.df$lex.Cst.uusi == "SOF"
idx_VKA  <- tmp.df$lex.Cst.uusi == "VKA"


# Make sure they align by CHADSVAbaseline
stopifnot(all(tmp.df$CHADSVAbaseline[idx_noac] == tmp.df$CHADSVAbaseline[idx_sof]))

# Create contrast matrix: NOAC vs SOF
ctr.mat <- tmp.mtr[idx_noac, ] - tmp.mtr[idx_sof, ]

# Compute IRRs with CIs
apu.all <- as.data.frame(ci.exp(tmp.m2, ctr.mat = ctr.mat))


apu.all$CHADSVAbaseline <- tmp.df$CHADSVAbaseline[idx_noac]
print(apu.all)





NOAC<-ggplot(apu.all, aes(x = factor(CHADSVAbaseline), y = `exp(Est.)`)) +
  geom_point() +
  geom_errorbar(aes(ymin = `2.5%`, ymax = `97.5%`), width = 0.2) +
  labs(
    x = "CHADSVAbaseline Score",
    y = "IRR (NOAC vs SOF)",
    title = "Adjusted IRRs of NOAC vs SOF across CHADSVAbaseline"
  ) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  theme_minimal()
# Ensure VKA and SOF CHADSVAbaseline levels match
stopifnot(all(tmp.df$CHADSVAbaseline[idx_VKA] == tmp.df$CHADSVAbaseline[idx_sof]))

# Create contrast matrix: VKA vs SOF
ctr.mat_vka <- tmp.mtr[idx_VKA, ] - tmp.mtr[idx_sof, ]

# Compute IRRs with CIs
apu.vka <- as.data.frame(ci.exp(tmp.m2, ctr.mat = ctr.mat_vka))

# Add CHADSVAbaseline levels to the result
apu.vka$CHADSVAbaseline <- tmp.df$CHADSVAbaseline[idx_VKA]
print(apu.vka)

VKA<- ggplot(apu.vka, aes(x = factor(CHADSVAbaseline), y = `exp(Est.)`)) +
  geom_point() +
  geom_errorbar(aes(ymin = `2.5%`, ymax = `97.5%`), width = 0.2) +
  labs(
    x = "CHADSVAbaseline Score",
    y = "IRR (VKA vs SOF)",
    title = "Adjusted IRRs of VKA vs SOF across CHADSVAbaseline"
  ) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  theme_minimal()

########

# Drop unused factor levels just in case
tmp.dt4$lex.Cst.uusi <- droplevels(tmp.dt4$lex.Cst)

# Fit model with interaction and full adjustment
tmp.m1 <- glm(cbind(lex.Xst == "IS", lex.dur) ~
                anyvasc+income_tertile+stroke
              +Hypertension+RenalVT+LiverVT+HF+alcohol+Bleeding
              +Diabetes+cancer+dementia+psychiatric+Hyperlipidemia
              +age.c+lex.Cst.uusi+sex+CHADSVAbaseline 
                ,
              data = tmp.dt4,
              family = "poisreg", maxit = 200)

tmp.m2 <- update(tmp.m1, ~ . + CHADSVAbaseline:lex.Cst.uusi)
anova(tmp.m1, tmp.m2, test = "Chisq")

tmp.df <- with(tmp.dt4,
               expand.grid(
                 age.c=sort(unique(age.c))[1],
                 sex=sort(unique(sex))[1],
                 income_tertile=sort(unique(income_tertile))[1],
                 anyvasc=sort(unique(anyvasc))[1],
                 Hyperlipidemia=sort(unique(Hyperlipidemia))[1],
                 psychiatric=sort(unique(psychiatric))[1],
                 dementia=sort(unique(dementia))[1],
                 cancer=sort(unique(cancer))[1],
                 Diabetes=sort(unique(Diabetes))[1],
                 alcohol=sort(unique(alcohol))[1],
                 Bleeding=sort(unique(Bleeding))[1],
                 HF=sort(unique(HF))[1],
                 LiverVT=sort(unique(LiverVT))[1],
                 RenalVT=sort(unique(RenalVT))[1],
                 Hypertension=sort(unique(Hypertension))[1],
                 stroke=sort(unique(stroke))[1],
                 CHADSVAbaseline=sort(unique(CHADSVAbaseline)),
                 
                 
                 lex.Cst.uusi=sort(unique(lex.Cst.uusi))))


tmp.mtr <- model.matrix(~ 
                          anyvasc+stroke+income_tertile
                        +Hypertension+RenalVT+LiverVT+HF
                        +alcohol+Bleeding+Diabetes+cancer+dementia
                        +psychiatric+Hyperlipidemia+age.c
                        +lex.Cst.uusi+sex + CHADSVAbaseline + CHADSVAbaseline:lex.Cst.uusi,
                        data = tmp.df)


# Identify rows by treatment group
idx_noac <- tmp.df$lex.Cst.uusi == "NOAC"
idx_sof  <- tmp.df$lex.Cst.uusi == "SOF"
idx_vka  <- tmp.df$lex.Cst.uusi == "VKA"

# Ensure same CHADSVAbaseline alignment
stopifnot(all(tmp.df$CHADSVAbaseline[idx_noac] == tmp.df$CHADSVAbaseline[idx_sof]))
stopifnot(all(tmp.df$CHADSVAbaseline[idx_vka]  == tmp.df$CHADSVAbaseline[idx_sof]))

# Create contrast matrices
ctr.noac <- tmp.mtr[idx_noac, ] - tmp.mtr[idx_sof, ]
ctr.vka  <- tmp.mtr[idx_vka, ]  - tmp.mtr[idx_sof, ]

# Compute IRRs
irr.noac <- as.data.frame(ci.exp(tmp.m2, ctr.mat = ctr.noac))
irr.vka  <- as.data.frame(ci.exp(tmp.m2, ctr.mat = ctr.vka))

# Label CHADSVAbaseline and comparison
irr.noac$CHADSVAbaseline <- tmp.df$CHADSVAbaseline[idx_noac]
irr.vka$CHADSVAbaseline  <- tmp.df$CHADSVAbaseline[idx_vka]

irr.noac$Comparison <- "NOAC vs. no anticoagulation"
irr.vka$Comparison  <- "VKA vs. no anticoagulation"

# Standardize column names
names(irr.noac)[1:3] <- c("IRR", "IRR.lo", "IRR.hi")
names(irr.vka)[1:3]  <- c("IRR", "IRR.lo", "IRR.hi")


# Combine into one dataframe
irr.all <- rbind(irr.noac, irr.vka)

# Plot
figVKAandNOACadjusted <- ggplot(
  irr.all, 
  aes(x = factor(CHADSVAbaseline), y = IRR, color = Comparison, group = Comparison)
) +
  geom_point(position = position_dodge(width = 0.4), size = 3) +
  geom_errorbar(
    aes(ymin = IRR.lo, ymax = IRR.hi),
    position = position_dodge(width = 0.4),
    width = 0.15
  ) +
  geom_line(position = position_dodge(width = 0.4), linewidth = 0.7) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "black") +
  labs(
    x = "CHA₂DS₂-VA Score",
    y = "Incidence rate ratio",
    title = "Adjusted analysis"
  ) +
  scale_y_continuous(
    limits = c(0.2, 1.4),
    breaks = seq(0.2, 1.4, 0.2)
  ) +
  theme_bw(base_size = 12, base_family = "sans") +
  theme(
    legend.title = element_blank(),
    legend.position = "top",
    legend.key = element_blank(),
    panel.grid.major = element_line(color = "grey85", linetype = "dotted"),
    panel.grid.minor = element_blank(),
    axis.title = element_text(size = 13, face = "bold"),
    axis.text = element_text(size = 11),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    panel.border = element_rect(color = "black", size = 0.6)
  )

print(figVKAandNOACadjusted)


# Fit model with interaction and full adjustment
tmp.m1 <- glm(cbind(lex.Xst == "IS", lex.dur) ~
                lex.Cst.uusi+CHADSVAbaseline 
              ,
              data = tmp.dt4,
              family = "poisreg", maxit = 200)

tmp.m2 <- update(tmp.m1, ~ . + CHADSVAbaseline:lex.Cst.uusi)
anova(tmp.m1, tmp.m2, test = "Chisq")

tmp.df <- with(tmp.dt4,
               expand.grid(
                 
                 CHADSVAbaseline=sort(unique(CHADSVAbaseline)),
                 
                 
                 lex.Cst.uusi=sort(unique(lex.Cst.uusi))))


tmp.mtr <- model.matrix(~ 
                          
                        lex.Cst.uusi + CHADSVAbaseline + CHADSVAbaseline:lex.Cst.uusi,
                        data = tmp.df)


# Identify rows by treatment group
idx_noac <- tmp.df$lex.Cst.uusi == "NOAC"
idx_sof  <- tmp.df$lex.Cst.uusi == "SOF"
idx_vka  <- tmp.df$lex.Cst.uusi == "VKA"

# Ensure same CHADSVAbaseline alignment
stopifnot(all(tmp.df$CHADSVAbaseline[idx_noac] == tmp.df$CHADSVAbaseline[idx_sof]))
stopifnot(all(tmp.df$CHADSVAbaseline[idx_vka]  == tmp.df$CHADSVAbaseline[idx_sof]))

# Create contrast matrices
ctr.noac <- tmp.mtr[idx_noac, ] - tmp.mtr[idx_sof, ]
ctr.vka  <- tmp.mtr[idx_vka, ]  - tmp.mtr[idx_sof, ]

# Compute IRRs
irr.noac <- as.data.frame(ci.exp(tmp.m2, ctr.mat = ctr.noac))
irr.vka  <- as.data.frame(ci.exp(tmp.m2, ctr.mat = ctr.vka))

# Label CHADSVAbaseline and comparison
irr.noac$CHADSVAbaseline <- tmp.df$CHADSVAbaseline[idx_noac]
irr.vka$CHADSVAbaseline  <- tmp.df$CHADSVAbaseline[idx_vka]

irr.noac$Comparison <- "NOAC vs. no anticoagulation"
irr.vka$Comparison  <- "VKA vs. no anticoagulation"

# Standardize column names
names(irr.noac)[1:3] <- c("IRR", "IRR.lo", "IRR.hi")
names(irr.vka)[1:3]  <- c("IRR", "IRR.lo", "IRR.hi")


# Combine into one dataframe
irr.all <- rbind(irr.noac, irr.vka)

# Plot

figVKAandNOACunadjusted <- ggplot(
  irr.all, 
  aes(x = factor(CHADSVAbaseline), y = IRR, color = Comparison, group = Comparison)
) +
  geom_point(position = position_dodge(width = 0.4), size = 3) +
  geom_errorbar(
    aes(ymin = IRR.lo, ymax = IRR.hi),
    position = position_dodge(width = 0.4),
    width = 0.15
  ) +
  geom_line(position = position_dodge(width = 0.4), linewidth = 0.7) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "black") +
  labs(
    x = "CHA₂DS₂-VA Score",
    y = "Incidence rate ratio",
    title = "Unadjusted analysis"
  ) +
  scale_y_continuous(
    limits = c(0.2, 1.4),
    breaks = seq(0.2, 1.4, 0.2)
  ) +
  theme_bw(base_size = 12, base_family = "sans") +
  theme(
    legend.title = element_blank(),
    legend.position = "top",
    legend.key = element_blank(),
    panel.grid.major = element_line(color = "grey85", linetype = "dotted"),
    panel.grid.minor = element_blank(),
    axis.title = element_text(size = 13, face = "bold"),
    axis.text = element_text(size = 11),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    panel.border = element_rect(color = "black", size = 0.6)
  )

print(figVKAandNOACunadjusted)



grid.arrange(
  figVKAandNOACunadjusted,
  figVKAandNOACadjusted, 
  ncol = 2
)

# Save VKA and NOAC figure as JPEG
jpeg("VKA_NOAC_IRR_Figures.jpg", width = 12, height = 6, units = "in", res = 300)  # 12x6 inches, 300 dpi
grid.arrange(
  figVKAandNOACunadjusted,
  figVKAandNOACadjusted, 
  ncol = 2
)
dev.off()

# Save VKA and NOAC figure as PDF
pdf("VKA_NOAC_IRR_Figures.pdf", width = 12, height = 6)  # size in inches

grid.arrange(
  figVKAandNOACunadjusted,
  figVKAandNOACadjusted,
  ncol = 2
)

dev.off()

######## NUMBER NEEDED TO TREAT #########

# 1. Merge datasets as before
dat <- apu.all %>%
  mutate(CHADSVAbaseline = as.numeric(as.character(CHADSVAbaseline))) %>%
  inner_join(
    as.data.frame(rate_noAK) %>%
      transmute(
        CHADSVAbaseline = as.numeric(as.character(CHADSVAbaseline)),
        rate = rate,   # events per 100 patient-years
        low  = low,
        high = high
      ),
    by = "CHADSVAbaseline"
  ) %>%
  mutate(
    logIRR    = log(MRR),
    SE_logIRR = (log(MRR.hi) - log(MRR.lo)) / (2*1.96),
    SE_rate   = (high - low) / (2*1.96)
  )

# 2. Bootstrap function with 10-year scaling
boot_one <- function(row, R = 5000) {
  
  IRR_star <- rlnorm(R, meanlog = row$logIRR, sdlog = row$SE_logIRR)
  rate_star <- rnorm(R, mean = row$rate, sd = row$SE_rate)
  rate_star[rate_star < 0] <- 0
  
  # ARR scaled to 10-year treatment per patient
  ARR_star <- (rate_star * 10 / 100) * (1 - IRR_star)
  NNT_star <- 1 / ARR_star
  
  tibble(
    ARR_lo  = quantile(ARR_star, 0.025, na.rm = TRUE),
    ARR_hi  = quantile(ARR_star, 0.975, na.rm = TRUE),
    NNT_lo  = quantile(NNT_star, 0.025, na.rm = TRUE),
    NNT_hi  = quantile(NNT_star, 0.975, na.rm = TRUE)
  )
}

# 3. Apply bootstrap per row
boot_results <- dat %>%
  rowwise() %>%
  mutate(
    ARR  = (rate * 10 / 100) * (1 - MRR),  # 10-year ARR per patient
    NNT  = 1 / ARR,
    boot = list(boot_one(cur_data()))
  ) %>%
  unnest(boot) %>%
  ungroup()

# 4. Final table
NNT_table <- boot_results %>%
  select(
    CHADSVAbaseline,
    rate, low, high,
    MRR, MRR.lo, MRR.hi,
    ARR, ARR_lo, ARR_hi,
    NNT, NNT_lo, NNT_hi
  )


# 1. Merge datasets as before
dat <- apu.allunadj %>%
  mutate(CHADSVAbaseline = as.numeric(as.character(CHADSVAbaseline))) %>%
  inner_join(
    as.data.frame(rate_noAK) %>%
      transmute(
        CHADSVAbaseline = as.numeric(as.character(CHADSVAbaseline)),
        rate = rate,   # events per 100 patient-years
        low  = low,
        high = high
      ),
    by = "CHADSVAbaseline"
  ) %>%
  mutate(
    logIRR    = log(MRR),
    SE_logIRR = (log(MRR.hi) - log(MRR.lo)) / (2*1.96),
    SE_rate   = (high - low) / (2*1.96)
  )

# 2. Bootstrap function with 10-year scaling
boot_one <- function(row, R = 5000) {
  
  IRR_star <- rlnorm(R, meanlog = row$logIRR, sdlog = row$SE_logIRR)
  rate_star <- rnorm(R, mean = row$rate, sd = row$SE_rate)
  rate_star[rate_star < 0] <- 0
  
  # ARR scaled to 10-year treatment per patient
  ARR_star <- (rate_star * 10 / 100) * (1 - IRR_star)
  NNT_star <- 1 / ARR_star
  
  tibble(
    ARR_lo  = quantile(ARR_star, 0.025, na.rm = TRUE),
    ARR_hi  = quantile(ARR_star, 0.975, na.rm = TRUE),
    NNT_lo  = quantile(NNT_star, 0.025, na.rm = TRUE),
    NNT_hi  = quantile(NNT_star, 0.975, na.rm = TRUE)
  )
}

# 3. Apply bootstrap per row
boot_results <- dat %>%
  rowwise() %>%
  mutate(
    ARR  = (rate * 10 / 100) * (1 - MRR),  # 10-year ARR per patient
    NNT  = 1 / ARR,
    boot = list(boot_one(cur_data()))
  ) %>%
  unnest(boot) %>%
  ungroup()

# 4. Final table
NNT_tableunadj <- boot_results %>%
  select(
    CHADSVAbaseline,
    rate, low, high,
    MRR, MRR.lo, MRR.hi,
    ARR, ARR_lo, ARR_hi,
    NNT, NNT_lo, NNT_hi
  )


print(NNT_tableunadj)

print(NNT_table)
