library(tidyverse)
library(mgcv)
library(gratia)
library(itsadug)
library(readxl)
library(bbmle)

d <- read_excel("HR_combined.xlsx", sheet = "nodead")

## turn variables into factors
## center (shift) continuous predictors
## create interaction treatment by egg sequence
d <- d %>% mutate(NestID=factor(NestID),
                  egg_ID=factor(egg_ID),
                  treat = factor(treat),
                  egg_seq = factor(egg_seq),
                  year=factor(year),
                  sex = factor(sex),
                  x = Inc_day-6,
                  egg_temp_c = egg_temp - mean(egg_temp),
                  TbyES = interaction(treat,egg_seq))

k_value <- 5

m0 <- gamm(HR ~ year + sex + s(egg_temp_c, k=k_value) + TbyES + s(x, by = TbyES, k=k_value) + 
             s(x, egg_ID, bs = "fs", k=k_value) + s(NestID, bs = "re"),
           correlation = corAR1(form = ~ x | egg_ID), 
           data = d,
           method = "ML")

draw(m0$gam)
summary(m0$gam)
## Error, not sure why

summary(m0$lme)
## Phy1 = 0.12, small

## Model without autocorrelation
m1 <- gamm(HR ~ year + sex + s(egg_temp_c, k=k_value) + TbyES + s(x, by = TbyES, k=k_value) + 
             s(x, egg_ID, bs = "fs", k=k_value) + s(NestID, bs = "re"),
           data = d,
           method = "ML")

summary(m1$gam)

AICctab(m0$lme,m1$lme)
## No need for autocorrelation

## This model doesn't work for k_value > 5 :
m2 <- gam(HR ~ year + sex + s(egg_temp_c, k=k_value) + TbyES + s(x, by = TbyES, k=k_value) + 
            s(x, egg_ID, bs = "fs", k=k_value) + s(NestID, bs = "re"),
          data = d,
          method = "ML")

summary(m2)
draw(m2)
anova(m2)

## Build more models ...

m3 <- gam(HR ~ year + sex + s(egg_temp_c, k=k_value) + TbyES + s(x, k=k_value) + 
            s(x, egg_ID, bs = "fs", k=k_value) + s(NestID, bs = "re"),
          data = d,
          method = "ML")

AICctab(m2,m3)
## 

m4 <- gam(HR ~ year + sex + s(egg_temp_c, k=k_value) + s(x, k=k_value) + 
            s(x, egg_ID, bs = "fs", k=k_value) + s(NestID, bs = "re"),
          data = d,
          method = "ML")

AICctab(m2,m3,m4)

m5 <- gam(HR ~ year + sex + s(egg_temp_c, k=k_value) + treat + s(x, k=k_value) + 
            s(x, egg_ID, bs = "fs", k=k_value) + s(NestID, bs = "re"),
          data = d,
          method = "ML")


diff <- difference_smooths(m2,smooth="s(x)")
draw(diff)



#############################################################################################
# Original code
#############################################################################################

#correct autocorelation
data <- start_event(data, event = "NestID")

Null <- gam(HR ~ s(Inc_day)+year+
              s(Inc_day, NestID, bs = "fs", k=5, m=1),
            rho = 0.912, data = data, AR.start = data$start.event,  method = "ML")
Null.temp <- gam(HR ~ egg_temp + year+sex+
                   s(Inc_day)+
                   s(Inc_day, NestID, bs = "fs", k=5, m=1),
                 rho = 0.912, data = data, AR.start = data$start.event,  method = "ML")
M <- gam(HR ~ egg_temp +year+
           s(Inc_day, by = TreatEgg_seq)+
           s(Inc_day, NestID, bs = "fs", k=5, m=1),
         rho = 0.912, data = data, AR.start = data$start.event,  method = "ML")


M.T <- gam(HR ~ egg_temp +year+
             s(Inc_day, by = treat)+
             s(Inc_day, NestID, bs = "fs", k=5, m=1),
           rho = 0.912, data = data, AR.start = data$start.event,  method = "ML")


M.S <- gam(HR ~ egg_temp +year+
             s(Inc_day, by = egg_seq)+
             s(Inc_day, NestID, bs = "fs", k=5, m=1),
           rho = 0.912, data = data, AR.start = data$start.event,  method = "ML")


M.TS <- gam(HR ~ egg_temp +year+
              s(Inc_day, by = egg_seq)+
              s(Inc_day, by = treat)+
              s(Inc_day, NestID, bs = "fs", k=5, m=1),
            rho = 0.912, data = data, AR.start = data$start.event,  method = "ML")

M0 <- gam(HR ~ treat*egg_seq + egg_temp +year+sex+
            s(Inc_day, by = TreatEgg_seq)+
            s(Inc_day, NestID, bs = "fs", k=5, m=1),
          rho = 0.912, data = data, AR.start = data$start.event,  method = "ML")


M0.T <- gam(HR ~ treat*egg_seq + egg_temp +year+sex+
              s(Inc_day, by = treat)+
              s(Inc_day, NestID, bs = "fs", k=5, m=1),
            rho = 0.912, data = data, AR.start = data$start.event,  method = "ML")


M0.S <- gam(HR ~ treat*egg_seq + egg_temp +year+sex+
              s(Inc_day, by = egg_seq)+
              s(Inc_day, NestID, bs = "fs", k=5, m=1),
            rho = 0.912, data = data, AR.start = data$start.event,  method = "ML")


M0.TS <- gam(HR ~ treat*egg_seq  + egg_temp +year+sex+
               s(Inc_day, by = egg_seq)+
               s(Inc_day, by = treat)+
               s(Inc_day, NestID, bs = "fs", k=5, m=1),
             rho = 0.912, data = data, AR.start = data$start.event,  method = "ML")


M1 <- gam(HR ~ treat + egg_seq + egg_temp +year+
            s(Inc_day, by = TreatEgg_seq)+
            s(Inc_day, NestID, bs = "fs", k=5, m=1),
          rho = 0.912, data = data, AR.start = data$start.event,  method = "ML")


M1.T <- gam(HR ~ treat + egg_seq + egg_temp +year+sex+
              s(Inc_day, by = treat)+
              s(Inc_day, NestID, bs = "fs", k=5, m=1),
            rho = 0.912, data = data, AR.start = data$start.event,  method = "ML")


M1.S <- gam(HR ~ treat + egg_seq + egg_temp +year+sex+
              s(Inc_day, by = egg_seq)+
              s(Inc_day, NestID, bs = "fs", k=5, m=1),
            rho = 0.912, data = data, AR.start = data$start.event,  method = "ML")


M1.TS <- gam(HR ~ treat + egg_seq + egg_temp +year+sex+
               s(Inc_day, by = egg_seq)+
               s(Inc_day, by = treat)+
               s(Inc_day, NestID, bs = "fs", k=5, m=1),
             rho = 0.912, data = data, AR.start = data$start.event,  method = "ML")



M2 <- gam(HR ~ treat + egg_temp +year+
            s(Inc_day, by = TreatEgg_seq)+
            s(Inc_day, NestID, bs = "fs", k=5, m=1),
          rho = 0.912, data = data, AR.start = data$start.event,  method = "ML")


M2.T <- gam(HR ~ treat + egg_temp +year+sex+
              s(Inc_day, by = treat)+
              s(Inc_day, NestID, bs = "fs", k=5, m=1),
            data = data, AR.start = data$start.event,  method = "ML")


M2.S <- gam(HR ~ treat + egg_temp +year+
              s(Inc_day, by = egg_seq)+
              s(Inc_day, NestID, bs = "fs", k=5, m=1),
            rho = 0.912, data = data, AR.start = data$start.event,  method = "ML")


M2.TS <- gam(HR ~ treat + egg_temp +year+
               s(Inc_day, by = egg_seq)+
               s(Inc_day, by = treat)+
               s(Inc_day, NestID, bs = "fs", k=5, m=1),
             rho = 0.912, data = data, AR.start = data$start.event,  method = "ML")



M3 <- gam(HR ~ egg_seq + egg_temp +year+
            s(Inc_day, by = TreatEgg_seq)+
            s(Inc_day, NestID, bs = "fs", k=5, m=1),
          rho = 0.912, data = data, AR.start = data$start.event,  method = "ML")


M3.T <- gam(HR ~ egg_seq + egg_temp +year+
              s(Inc_day, by = treat)+
              s(Inc_day, NestID, bs = "fs", k=5, m=1),
            rho = 0.912, data = data, AR.start = data$start.event,  method = "ML")


M3.S <- gam(HR ~ egg_seq + egg_temp +year+sex+
              s(Inc_day, by = egg_seq)+
              s(Inc_day, NestID, bs = "fs", k=5, m=1),
            rho = 0.912, data = data, AR.start = data$start.event,  method = "ML")


M3.TS <- gam(HR ~ egg_seq + egg_temp +year+
               s(Inc_day, by = egg_seq)+
               s(Inc_day, by = treat)+
               s(Inc_day, NestID, bs = "fs", k=5, m=1),
             rho = 0.912, data = data, AR.start = data$start.event,  method = "ML")

M4 <- gam(HR ~ TreatEgg_seq + egg_temp +year+sex+
            s(Inc_day, by = TreatEgg_seq)+
            s(Inc_day, NestID, bs = "fs", k=5, m=1),
          rho = 0.912, data = data, AR.start = data$start.event,  method = "ML")


M4.T <- gam(HR ~ TreatEgg_seq + egg_temp +year+sex+
              s(Inc_day, by = treat)+
              s(Inc_day, NestID, bs = "fs", k=5, m=1),
            rho = 0.912, data = data, AR.start = data$start.event,  method = "ML")


M4.S <- gam(HR ~ TreatEgg_seq + egg_temp +year+sex+
              s(Inc_day, by = egg_seq)+
              s(Inc_day, NestID, bs = "fs", k=5, m=1),
            rho = 0.912, data = data, AR.start = data$start.event,  method = "ML")


M4.TS <- gam(HR ~ TreatEgg_seq  + egg_temp +year+sex+
               s(Inc_day, by = egg_seq)+
               s(Inc_day, by = treat)+
               s(Inc_day, NestID, bs = "fs", k=5, m=1),
             rho = 0.912, data = data, AR.start = data$start.event,  method = "ML")




compareML(M4, M2.T)
compareML(M4, M3.S)
compareML(M0.TS, M4)
compareML(M0.T, M4)
compareML(M0.S, M4)


aic <- AIC(Null.temp, Null, M, M.T, M.S, M.TS, M0, M0.T, M0.S, M0.TS, M1, M1.T, M1.S, M1.TS, M2, M2.T, M2.S, M2.TS, M3, M3.T, M3.S, M3.TS, M4, M4.T, M4.S, M4.TS)


m0.TS <- gam(HR ~ treat*egg_seq  + egg_temp +year+sex+
               s(Inc_day, by = egg_seq)+
               s(Inc_day, by = treat)+
               s(Inc_day, NestID, bs = "fs", k=5, m=1),
             rho = 0.912, data = data, AR.start = data$start.event)
summary(m0.TS)

