# This script fits the models that test the impact of helpers on fitness
# Jeremy Summers
# July 2025

library(tidyverse)
library(glmmTMB)
library(usdm)
library(DHARMa)

# load dataset
load("data/FitnessDataSet_20250714.rdata")

# JUVENILE SURVIVAL

# convert NAs to 0
Offspring.data <- Offspring.input %>% 
  # remove years with >10% territories with unsexed helpers
  filter(!Year %in% c(1994, 1997, 1999)) %>%
  # set NAs to 0 for records without helpers
  replace_na(list(mean.relate = 0, sex.ratio = 0))

# create models
# Base model
O.null.model <- glmmTMB(surv ~ num.helpers + 
                          pedF + HatchDate + ha + pairYear +
                          (1|Terr) + (1|Year) + (1|NatalNest),
                        family = "binomial",
                        data = Offspring.data)

# Non-linear model
O.null.2.model <- glmmTMB(surv ~ num.helpers + I(num.helpers^2) +
                          pedF + HatchDate + ha + pairYear +
                          (1|Terr) + (1|Year) + (1|NatalNest),
                        family = "binomial",
                        data = Offspring.data)

# Relatedness model
O.relate.model <- glmmTMB(surv ~ num.helpers +
                            num.helpers:mean.relate +
                     pedF + HatchDate + ha + pairYear +
                     (1|Terr) + (1|Year) + (1|NatalNest),
                   family = "binomial",
                   data = Offspring.data)

# Sex model
O.sex.model <- glmmTMB(surv ~ num.helpers + 
                            num.helpers:sex.ratio +
                            pedF + HatchDate + ha + pairYear +
                         (1|Terr) + (1|Year) + (1|NatalNest),
                          family = "binomial",
                          data = Offspring.data)

# compare models
anova(O.null.model, O.null.2.model)

anova(O.null.model, O.relate.model)

anova(O.null.model, O.sex.model)

anova(O.null.model, O.null.2.model,
      O.relate.model, O.sex.model)

# test model fit
O.null.sim <- simulateResiduals(O.null.model)
testResiduals(O.null.sim)
# no significant issues

O.null.2.sim <- simulateResiduals(O.null.2.model)
testResiduals(O.null.2.sim)
# no significant issues

O.relate.sim <- simulateResiduals(O.relate.model)
testResiduals(O.relate.sim)
# no significant issues

O.sex.sim <- simulateResiduals(O.sex.model)
testResiduals(O.sex.sim)
# no significant issues

### POST-HOC TERRITORY AREA INTERACTION

# Sex-only model with territory interaction
O.ha.sex.model <- glmmTMB(surv ~ num.helpers*ha + 
                                  num.helpers:sex.ratio*ha +
                          pedF + HatchDate + pairYear + 
                          (1|Terr) + (1|Year) + (1|NatalNest),
                        family = "binomial",
                        data = Offspring.data)

# Relatedness-only model with territory interaction
O.ha.relate.model <- glmmTMB(surv ~ num.helpers*ha +
                                  num.helpers:mean.relate*ha +
                          pedF + HatchDate + pairYear + 
                          (1|Terr) + (1|Year) + (1|NatalNest),
                        family = "binomial",
                        data = Offspring.data)

# Base-only model with territory interaction
O.null.ha.model <- glmmTMB(surv ~ num.helpers*ha +
                          pedF + HatchDate + pairYear + 
                          (1|Terr) + (1|Year) + (1|NatalNest),
                        family = "binomial",
                        data = Offspring.data)

# Non-linear model with territory interaction
O.null.2.ha.model <- glmmTMB(surv ~ num.helpers*ha + I(num.helpers^2)*ha +
                             pedF + HatchDate + pairYear + 
                             (1|Terr) + (1|Year) + (1|NatalNest),
                           family = "binomial",
                           data = Offspring.data)

# compare models
anova(O.null.model, O.null.ha.model, O.ha.relate.model,
      O.ha.sex.model, O.null.2.ha.model)

# test model fit
O.null.ha.sim <- simulateResiduals(O.null.ha.model)
testResiduals(O.null.ha.sim)
# no significant issues

O.ha.sex.sim <- simulateResiduals(O.ha.sex.model)
testResiduals(O.ha.sex.sim)
# no significant issues

O.ha.relate.sim <- simulateResiduals(O.ha.relate.model)
testResiduals(O.ha.sex.sim)
# no significant issues

O.null.2.ha.sim <- simulateResiduals(O.null.2.ha.model)
testResiduals(O.null.2.ha.sim)
# no significant issues

# BREEDER SURVIVAL

# separate dataset by sex
Breeders.Census.F <- filter(Breeders.input, Sex == "F") %>%
  # remove years with >10% territories with unsexed helpers
  filter(!Year %in% c(1994, 1997, 1999)) %>%
  # set NAs to 0 for records without helpers
  replace_na(list(mean.relate = 0, sex.ratio = 0))

Breeders.Census.M <- filter(Breeders.input, Sex == "M") %>%
  # remove years with >10% territories with unsexed helpers
  filter(!Year %in% c(1994, 1997, 1999)) %>%
  # set NAs to 0 for records without helpers
  replace_na(list(mean.relate = 0, sex.ratio = 0))

# create half-normal random effect prior
Bprior <- data.frame(prior = "normal(0, 1)",
                     class = "ranef",
                     coef = "")

# create models
# Base models
B.F.null.model <- glmmTMB(surv ~ num.helpers + 
                            pedF + age + I(age^2) + ha +
                            (1|Terr) + (1|Year) + (1|USFWS),
                          family = "binomial",
                          data = Breeders.Census.F,
                          prior = Bprior)

B.M.null.model <- glmmTMB(surv ~ num.helpers + 
                            pedF + age + I(age^2) + ha +
                            (1|Terr) + (1|Year) + (1|USFWS),
                          family = "binomial",
                          data = Breeders.Census.M,
                          prior = Bprior)

# Non-linear models
B.F.null.2.model <- glmmTMB(surv ~ num.helpers + I(num.helpers^2) +
                            pedF + age + I(age^2) + ha +
                              (1|Terr) + (1|Year) + (1|USFWS),
                          family = "binomial",
                          data = Breeders.Census.F,
                          prior = Bprior)

B.M.null.2.model <- glmmTMB(surv ~ num.helpers + I(num.helpers^2) +
                            pedF + age + I(age^2) + ha +
                            (1|Terr) + (1|Year) + (1|USFWS),
                          family = "binomial",
                          data = Breeders.Census.M,
                          prior = Bprior)

# Relatedness models
B.F.relate.model <- glmmTMB(surv ~ num.helpers + 
                       num.helpers:mean.relate +
                       pedF + age + I(age^2) + ha +
                       (1|Terr) + (1|Year) + (1|USFWS),
                     family = "binomial",
                     data = Breeders.Census.F,
                     prior = Bprior)

B.M.relate.model <- glmmTMB(surv ~ num.helpers + 
                       num.helpers:mean.relate +
                       pedF + age + I(age^2) + ha +
                       (1|Terr) + (1|Year) + (1|USFWS),
                     family = "binomial",
                     data = Breeders.Census.M,
                     prior = Bprior)

# Sex models
B.F.sex.model <- glmmTMB(surv ~ num.helpers + 
                       num.helpers:sex.ratio +
                       pedF + age + I(age^2) + ha +
                       (1|Terr) + (1|Year) + (1|USFWS),
                     family = "binomial",
                     data = Breeders.Census.F,
                     prior = Bprior)

B.M.sex.model <- glmmTMB(surv ~ num.helpers + 
                       num.helpers:sex.ratio +
                       pedF + age + I(age^2) + ha +
                       (1|Terr) + (1|Year) + (1|USFWS),
                     family = "binomial",
                     data = Breeders.Census.M,
                     prior = Bprior)

# compare models
anova(B.F.null.model, B.F.null.2.model,
      B.F.relate.model, B.F.sex.model)

anova(B.M.null.model, B.M.null.2.model,
      B.M.relate.model, B.M.sex.model)

# test model fit
B.F.null.sim <- simulateResiduals(B.F.null.model)
testResiduals(B.F.null.sim)
# DHARMa has issues with binomial for testing outliers,
# it is recommended to re-run with type = "bootstrap"
testOutliers(B.F.null.sim, type = "bootstrap")
# no significant issues

B.F.null.2.sim <- simulateResiduals(B.F.null.2.model)
testResiduals(B.F.null.2.sim)
# DHARMa has issues with binomial for testing outliers,
# it is recommended to re-run with type = "bootstrap"
testOutliers(B.F.null.2.sim, type = "bootstrap")
# no significant issues

B.F.relate.sim <- simulateResiduals(B.F.relate.model)
testResiduals(B.F.relate.sim)
# DHARMa has issues with binomial for testing outliers,
# it is recommended to re-run with type = "bootstrap"
testOutliers(B.F.relate.sim, type = "bootstrap")
# no significant issues

B.F.sex.sim <- simulateResiduals(B.F.sex.model)
testResiduals(B.F.sex.sim)
# DHARMa has issues with binomial for testing outliers,
# it is recommended to re-run with type = "bootstrap"
testOutliers(B.F.relate.sim, type = "bootstrap")
# no significant issues

B.M.null.sim <- simulateResiduals(B.M.null.model)
testResiduals(B.M.null.sim)
# DHARMa has issues with binomial for testing outliers,
# it is recommended to re-run with type = "bootstrap"
testOutliers(B.M.null.sim, type = "bootstrap")
# no significant issues

B.M.null.2.sim <- simulateResiduals(B.M.null.2.model)
testResiduals(B.M.null.2.sim)
# DHARMa has issues with binomial for testing outliers,
# it is recommended to re-run with type = "bootstrap"
testOutliers(B.M.null.2.sim, type = "bootstrap")
# no significant issues

B.M.relate.sim <- simulateResiduals(B.M.relate.model)
testResiduals(B.M.relate.sim)
# DHARMa has issues with binomial for testing outliers,
# it is recommended to re-run with type = "bootstrap"
testOutliers(B.M.null.sim, type = "bootstrap")
# no significant issues

B.M.sex.sim <- simulateResiduals(B.M.sex.model)
testResiduals(B.M.sex.sim)
# DHARMa has issues with binomial for testing outliers,
# it is recommended to re-run with type = "bootstrap"
testOutliers(B.M.null.sim, type = "bootstrap")
# no significant issues

# NESTLING PRODUCTION
NestPairs.data <- replace_na(NestPairs.input, 
                             # set NAs to 0 for records without helpers
                             list(mean.relate = 0, sex.ratio = 0)) %>%
  # remove years with >10% territories with unsexed helpers
  filter(!Year %in% c(1994, 1997, 1999))

# test zero-inflation
ARS.null.model.poisson <- glmmTMB(ARS ~ num.helpers + 
                                    pairYear + ha + 
                                    pair.pedF + Month +
                                    (1|Year) + (1|Terr) + (1|pairID),
                                  family = "poisson",
                                  data = NestPairs.data)

testZeroInflation(ARS.null.model.poisson)
# ratioObsSim = 1.71, p-value < 2.2e-16

ARS.null.model.poisson.zi <- glmmTMB(ARS ~ num.helpers + 
                                       pairYear + ha + 
                                       pair.pedF + Month +
                                       (1|Year) + (1|Terr) + (1|pairID),
                                     zi = ~ 1,
                                     family = "poisson",
                                     data = NestPairs.data)

testZeroInflation(ARS.null.model.poisson.zi)
# ratioObsSim = 1.0217, p-value = 0.584

# test dispersion
testDispersion(ARS.null.model.poisson.zi)
# dispersion = 0.6832, p-value < 2.2e-16

# Base model
ARS.null.model <- glmmTMB(ARS ~ num.helpers + 
                       pairYear + ha + 
                       pair.pedF + Month + (1|Year) + (1|Terr) + (1|pairID),
                     zi = ~ 1,
                     family = "genpois",
                     data = NestPairs.data)

testZeroInflation(ARS.null.model)
# ratioObsSim = 1.0041, p-value = 0.912

testDispersion(ARS.null.model)
# dispersion = 1.0041, p-value = 0.904

# create remaining models

# Non-linear model
ARS.null.2.model <- glmmTMB(ARS ~ num.helpers + I(num.helpers^2) +
                            pairYear + ha + 
                            pair.pedF + Month +
                            (1|Year) + (1|Terr) + (1|pairID),
                          zi = ~ 1,
                          family = "genpois",
                          data = NestPairs.data)

# Relatedness model
ARS.relate.model <- glmmTMB(ARS ~ num.helpers + 
                       num.helpers:mean.relate +
                       pairYear + ha + 
                       pair.pedF + Month +
                       (1|Year) + (1|Terr) + (1|pairID),
                     zi = ~ 1,
                     family = "genpois",
                     data = NestPairs.data)

# Sex model
ARS.sex.model <- glmmTMB(ARS ~ num.helpers +
                       num.helpers:sex.ratio +
                       pairYear + ha + 
                       pair.pedF + Month +
                       (1|Year) + (1|Terr) + (1|pairID),
                     zi = ~ 1,
                     family = "genpois",
                     data = NestPairs.data)


# compare models
anova(ARS.null.model, ARS.null.2.model, ARS.relate.model,
      ARS.sex.model)

# test model fit
ARS.null.sim <- simulateResiduals(ARS.null.model)
testResiduals(ARS.null.sim)
# no significant issues

ARS.null.2.sim <- simulateResiduals(ARS.null.2.model)
testResiduals(ARS.null.2.sim)
# no significant issues

ARS.relate.sim <- simulateResiduals(ARS.relate.model)
testResiduals(ARS.relate.sim)
# no significant issues

ARS.sex.sim <- simulateResiduals(ARS.sex.model)
testResiduals(ARS.sex.sim)
# no significant issues

# save results
save(O.null.model, O.null.2.model, O.relate.model, O.sex.model,
     O.null.sim, O.null.2.sim, O.relate.sim, O.sex.sim,
     O.null.ha.model, O.null.ha.sim,
     O.ha.relate.model, O.ha.relate.sim,
     O.ha.sex.model, O.ha.sex.sim,
     O.null.2.ha.model, O.null.2.ha.sim,
     B.M.null.model, B.M.null.2.model, B.M.relate.model, B.M.sex.model,
     B.F.null.model, B.F.null.2.model, B.F.relate.model, B.F.sex.model,
     B.M.null.sim, B.M.null.2.sim, B.M.relate.sim, B.M.sex.sim,
     B.F.null.sim, B.F.null.2.sim, B.F.relate.sim, B.F.sex.sim,
     ARS.null.model, ARS.relate.model, ARS.sex.model, ARS.null.2.model,
     ARS.null.sim, ARS.null.2.sim, ARS.relate.sim, ARS.sex.sim,
     Offspring.data,
     Breeders.Census.F, Breeders.Census.M, NestPairs.data,
     file = "data/HelperModels_20250714.rdata")
