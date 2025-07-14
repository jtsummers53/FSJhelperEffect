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
  replace_na(list(mean.relate = 0, sex.ratio = 0)) %>%
  # scale non-focal fixed effects
  mutate_at(c("ha", "HatchDate", "pairYear", "pedF"), scale)

# create models
# Base model
O.null.model <- glmmTMB(surv ~ num.helpers + 
                          pedF + HatchDate + ha + pairYear +
                          (1|Terr) + (1|Year) + (1|NatalNest),
                        family = "binomial",
                        data = Offspring.data)

# Base model with non-linear effects
O.null.2.model <- glmmTMB(surv ~ num.helpers + I(num.helpers^2) +
                          pedF + HatchDate + ha + pairYear +
                          (1|Terr) + (1|Year) + (1|NatalNest),
                        family = "binomial",
                        data = Offspring.data)

# Full model
O.model <- glmmTMB(surv ~ num.helpers + 
                     num.helpers:mean.relate +
                     num.helpers:sex.ratio +
                     pedF + HatchDate + ha + pairYear +
                     (1|Terr) + (1|Year) + (1|NatalNest),
                   family = "binomial",
                   data = Offspring.data)

# Full model with non-linear effects
O.2.model <- glmmTMB(surv ~ num.helpers + I(num.helpers^2) +
                       num.helpers:mean.relate +
                       I(num.helpers^2):mean.relate +
                       num.helpers:sex.ratio +
                       I(num.helpers^2):sex.ratio +
                       pedF + HatchDate + ha + pairYear + 
                       (1|Terr) + (1|Year) + (1|NatalNest),
                     family = "binomial",
                     data = Offspring.data)

# compare models
anova(O.null.model, O.null.2.model, O.model, O.2.model)

# test model fit
O.null.sim <- simulateResiduals(O.null.model)
testResiduals(O.null.sim)
# no significant issues

O.null.2.sim <- simulateResiduals(O.null.2.model)
testResiduals(O.null.2.sim)
# no significant issues

O.sim <- simulateResiduals(O.model)
testResiduals(O.sim)
# no significant issues

O.2.sim <- simulateResiduals(O.2.model)
testResiduals(O.2.sim)
# no significant issues

### POST-HOC TERRITORY AREA INTERACTION

# Full model with territory interaction
O.ha.model <- glmmTMB(surv ~ num.helpers*ha + 
                          num.helpers:mean.relate*ha +
                          num.helpers:sex.ratio*ha +
                          pedF + HatchDate + pairYear + 
                          (1|Terr) + (1|Year) + (1|NatalNest),
                        family = "binomial",
                        data = Offspring.data)

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

# compare models
anova(O.null.model, O.null.ha.model, O.ha.relate.model,
      O.ha.sex.model,
      O.ha.model)

# test model fit
O.null.ha.sim <- simulateResiduals(O.null.ha.model)
testResiduals(O.null.ha.sim)
# no significant issues

O.ha.sim <- simulateResiduals(O.ha.model)
testResiduals(O.ha.sim)
# no significant issues

O.ha.sex.sim <- simulateResiduals(O.ha.sex.model)
testResiduals(O.ha.sex.sim)
# no significant issues

O.ha.relate.sim <- simulateResiduals(O.ha.relate.model)
testResiduals(O.ha.sex.sim)
# no significant issues

# BREEDER SURVIVAL

# separate dataset by sex
Breeders.Census.F <- filter(Breeders.input, Sex == "F") %>%
  # remove years with >10% territories with unsexed helpers
  filter(!Year %in% c(1994, 1997, 1999)) %>%
  # set NAs to 0 for records without helpers
  replace_na(list(mean.relate = 0, sex.ratio = 0)) %>%
  # scale non-focal fixed effects
  mutate_at(c("ha", "pedF", "age"), scale)

Breeders.Census.M <- filter(Breeders.input, Sex == "M") %>%
  # remove years with >10% territories with unsexed helpers
  filter(!Year %in% c(1994, 1997, 1999)) %>%
  # set NAs to 0 for records without helpers
  replace_na(list(mean.relate = 0, sex.ratio = 0)) %>%
  # scale non-focal fixed effects
  mutate_at(c("ha", "pedF", "age"), scale)

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

# Base models with non-linear effects
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

# Full models
B.F.model <- glmmTMB(surv ~ num.helpers + 
                       num.helpers:mean.relate +
                       num.helpers:sex.ratio +
                       pedF + age + I(age^2) + ha +
                       (1|Terr) + (1|Year) + (1|USFWS),
                     family = "binomial",
                     data = Breeders.Census.F,
                     prior = Bprior)

B.M.model <- glmmTMB(surv ~ num.helpers + 
                       num.helpers:mean.relate +
                       num.helpers:sex.ratio +
                       pedF + age + I(age^2) + ha +
                       (1|Terr) + (1|Year) + (1|USFWS),
                     family = "binomial",
                     data = Breeders.Census.M,
                     prior = Bprior)

# Full models with non-linear effects
B.F.2.model <- glmmTMB(surv ~ num.helpers + I(num.helpers^2) +
                         num.helpers:mean.relate +
                         I(num.helpers^2):mean.relate +
                         num.helpers:sex.ratio +
                         I(num.helpers^2):sex.ratio +
                         pedF + age + I(age^2) + ha +
                         (1|Terr) + (1|Year) + (1|USFWS),
                       family = "binomial",
                       data = Breeders.Census.F,
                       prior = Bprior)

B.M.2.model <- glmmTMB(surv ~ num.helpers + I(num.helpers^2) +
                         num.helpers:mean.relate +
                         I(num.helpers^2):mean.relate +
                         num.helpers:sex.ratio +
                         I(num.helpers^2):sex.ratio +
                         pedF + age + I(age^2) + ha +
                          (1|Terr) + (1|Year) + (1|USFWS),
                       family = "binomial",
                       data = Breeders.Census.M,
                       prior = Bprior)

# compare models
anova(B.F.null.model, B.F.model, B.F.2.model, B.F.null.2.model)
anova(B.M.null.model, B.M.model, B.M.2.model, B.M.null.2.model)

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

B.F.sim <- simulateResiduals(B.F.model)
testResiduals(B.F.sim)
# DHARMa has issues with binomial for testing outliers,
# it is recommended to re-run with type = "bootstrap"
testOutliers(B.F.sim, type = "bootstrap")
# no significant issues

B.F.2.sim <- simulateResiduals(B.F.2.model)
testResiduals(B.F.2.sim)
# DHARMa has issues with binomial for testing outliers,
# it is recommended to re-run with type = "bootstrap"
testOutliers(B.F.2.sim, type = "bootstrap")
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

B.M.sim <- simulateResiduals(B.M.model)
testResiduals(B.M.sim)
# DHARMa has issues with binomial for testing outliers,
# it is recommended to re-run with type = "bootstrap"
testOutliers(B.M.sim, type = "bootstrap")
# no significant issues

B.M.2.sim <- simulateResiduals(B.M.2.model)
testResiduals(B.M.2.sim)
# DHARMa has issues with binomial for testing outliers,
# it is recommended to re-run with type = "bootstrap"
testOutliers(B.M.2.sim, type = "bootstrap")
# no significant issues

# NESTLING PRODUCTION
NestPairs.data <- replace_na(NestPairs.input, 
                             # set NAs to 0 for records without helpers
                             list(mean.relate = 0, sex.ratio = 0)) %>%
  # remove years with >10% territories with unsexed helpers
  filter(!Year %in% c(1994, 1997, 1999)) %>%
  # scale non-focal fixed effects
  mutate_at(c("pairYear", "ha", "pair.pedF", "Month"), scale)

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

# Base model with non-linear effects
ARS.null.2.model <- glmmTMB(ARS ~ num.helpers + I(num.helpers^2) +
                            pairYear + ha + 
                            pair.pedF + Month +
                            (1|Year) + (1|Terr) + (1|pairID),
                          zi = ~ 1,
                          family = "genpois",
                          data = NestPairs.data)

# Full model
ARS.model <- glmmTMB(ARS ~ num.helpers + 
                       num.helpers:mean.relate +
                       num.helpers:sex.ratio +
                       pairYear + ha + 
                       pair.pedF + Month +
                       (1|Year) + (1|Terr) + (1|pairID),
                     zi = ~ 1,
                     family = "genpois",
                     data = NestPairs.data)

# Full model with non-linear effects
ARS.2.model <- glmmTMB(ARS ~ num.helpers + I(num.helpers^2) +
                       num.helpers:mean.relate +
                       I(num.helpers^2):mean.relate +
                       num.helpers:sex.ratio +
                       I(num.helpers^2):sex.ratio +
                       pairYear + ha + 
                       pair.pedF + Month +
                       (1|Year) + (1|Terr) + (1|pairID),
                     zi = ~ 1,
                     family = "genpois",
                     data = NestPairs.data)

# compare models
anova(ARS.null.model, ARS.model, ARS.2.model, ARS.null.2.model)

# test model fit
ARS.null.sim <- simulateResiduals(ARS.null.model)
testResiduals(ARS.null.sim)
# no significant issues

ARS.null.2.sim <- simulateResiduals(ARS.null.2.model)
testResiduals(ARS.null.2.sim)
# no significant issues

ARS.sim <- simulateResiduals(ARS.model)
testResiduals(ARS.sim)
# no significant issues

ARS.2.sim <- simulateResiduals(ARS.2.model)
testResiduals(ARS.2.sim)
# no significant issues

# save results
save(O.null.model, O.model, O.2.model, O.null.2.model,
     O.null.sim, O.sim, O.2.sim, O.null.2.sim,
     O.null.ha.model, O.null.ha.sim,
     O.ha.model, O.ha.sim,
     O.ha.relate.model, O.ha.relate.sim,
     O.ha.sex.model, O.ha.sex.sim,
     B.M.null.model, B.M.model, B.M.2.model, B.M.null.2.model,
     B.F.null.model, B.F.model, B.F.2.model, B.F.null.2.model,
     B.M.null.sim, B.M.sim, B.M.2.sim, B.M.null.2.sim,
     B.F.null.sim, B.F.sim, B.F.2.sim, B.F.null.2.sim,
     ARS.null.model, ARS.model, ARS.2.model, 
     ARS.null.sim, ARS.sim, ARS.2.sim, ARS.null.2.model, ARS.null.2.sim,
     Offspring.data,
     Breeders.Census.F, Breeders.Census.M, NestPairs.data,
     file = "data/HelperModels_20250714.rdata")
