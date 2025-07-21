# Jeremy Summers
# This script generates tables and figures for the sex-specific effects of helpers paper
# July 2025

library(glmmTMB)
library(flextable)
library(usdm)
library(modelsummary)
library(jstable)
library(ggeffects)
library(kableExtra)
library(patchwork)
library(tidyverse)

options(knitr.kable.NA = '')

# load models and dataset
load("data/HelperModels_20250714.rdata")
load("data/FitnessDataSet_20250714.rdata")

#### SURVIVAL

# predict survival for male and female breeders
# predict for sex based female model
Surv.pred <- ggpredict(B.F.sex.model,
                         c("num.helpers [all]",
                           # predict between scenarios with all male (0)
                           # and all female (1) helpers
                           "sex.ratio [0, 1]"),
                       # bias_correction adjusts values by random effects
                         bias_correction = TRUE) %>%
  data.frame() %>%
  mutate(Bsex = "F", var = "sex") %>%
  # repeat for sex based male model
  bind_rows(ggpredict(B.M.sex.model,
                      c("num.helpers [all]",
                        "sex.ratio [0, 1]"),
                      bias_correction = TRUE) %>%
              data.frame() %>%
              mutate(Bsex = "M", var = "sex")) %>%
  # repeat for relatedness based female model
  bind_rows(ggpredict(B.F.relate.model,
                      c("num.helpers [all]",
                        # predict between scenarios with all related (0)
                        # and all unrelated (1) helpers
                        "mean.relate [0, 1]"),
                      bias_correction = TRUE) %>%
              data.frame() %>%
              mutate(Bsex = "F", var = "relate")) %>%
  # repeat for relatedness based male model
  bind_rows(ggpredict(B.M.relate.model,
                      c("num.helpers [all]",
                        # predict between scenarios with all related (0)
                        # and all unrelated (1) helpers
                        "mean.relate [0, 1]"),
                      bias_correction = TRUE) %>%
              data.frame() %>%
              mutate(Bsex = "M", var = "relate")) %>%
  # repeat for base male model
  bind_rows(ggpredict(B.M.null.model,
                      c("num.helpers [all]"),
                      bias_correction = TRUE) %>%
              data.frame() %>%
              mutate(Bsex = "M", group = "B", var = "null")) %>%
  # repeat for base female model
  bind_rows(ggpredict(B.F.null.model,
                      c("num.helpers [all]"),
                      bias_correction = TRUE) %>%
              data.frame() %>%
              mutate(Bsex = "F", group = "B", var = "null")) %>%
  # adjust labels for plotting
  mutate(num.helpers = x,
         group = if_else(var == "relate",
                         if_else(group == "0", "U", "R"),
                         if_else(var == "sex",
                                 if_else(group == "0", "M", "F"),
                                 "B")),
         var = factor(var, levels = c("null", "sex", "relate")))

# plot breeder survival prediction curves
Breed_Surv.plot <- ggplot(Surv.pred %>%
                            # limit number of helpers to maximum
                            # number of single-sex helpers observed (4)
                      filter(num.helpers <= 5),
                    aes(x = num.helpers, y = predicted,
                        col = group, fill = group)) +
  geom_line(linewidth = 0.5) +
  # add 95% confidence intervals
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.2, col = NA) +
  scale_color_manual(values = c(F = "#FF2B0F", M = "#009AC5", B = "#754D98",
                                R = "#2DA93C", U = "#8D5802"),
                     labels = c(F = "Female", M = "Male", B = "All",
                                R = "Related", U = "Unrelated")) +
  scale_fill_manual(values = c(F = "#FF2B0F", M = "#009AC5", B = "#754D98",
                               R = "#2DA93C", U = "#8D5802"),
                    labels = c(F = "Female", M = "Male", B = "All",
                               R = "Related", U = "Unrelated")) +
  # add horizontal line showing predicted survival with 0 helpers
  geom_hline(data = group_by(Surv.pred, Bsex) %>%
               summarize(zero.helpers = predicted[1]),
             aes(yintercept = zero.helpers), col = "black", lty = "dashed") +
  labs(x = "# of Helpers", y = "Predicted Breeder Survival",
       fill = "Helpers",
       col = "Helpers") +
  facet_grid(Bsex~var,
             labeller = labeller(Bsex = c(F = "Female Breeder Survival",
                                          M = "Male Breeder Survival"),
                                 var = c(null = "Total Helper #",
                                        sex = "Helper Sex",
                                      relate = "Helper\nRelatedness"))) +
  theme_minimal() +
  theme(text = element_text(size = 12),
        panel.border = element_rect(fill = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

ggsave(filename = "figures/Figure1_20250714.pdf",
       Breed_Surv.plot,
       width = 5,
       height = 5,
       units = "in")

ggsave(filename = "figures/Figure1_20250714.png",
       Breed_Surv.plot,
       width = 5,
       height = 4,
       units = "in",
       dpi = 500)

# Follow-up offspring territory area analyses

# predict offspring survival
# predict for sex based model
O.pred <- ggpredict(O.sex.model,
                       c("num.helpers [all]",
                         # predict between scenarios with all male (0)
                         # and all female (1) helpers
                         "sex.ratio [0, 1]"),
                       # bias_correction adjusts values by random effects
                       bias_correction = TRUE) %>%
  data.frame() %>%
  mutate(var = "sex") %>%
  # repeat for relatedness based model
  bind_rows(ggpredict(O.relate.model,
                      c("num.helpers [all]",
                        # predict between scenarios with all related (0)
                        # and all unrelated (1) helpers
                        "mean.relate [0, 1]"),
                      bias_correction = TRUE) %>%
              data.frame() %>%
              mutate(var = "relate")) %>%
  # repeat for base model
  bind_rows(ggpredict(O.null.model,
                      c("num.helpers [all]"),
                      bias_correction = TRUE) %>%
              data.frame() %>%
              mutate(group = "B", var = "null")) %>%
  # adjust labels for plotting
  mutate(num.helpers = x,
         group = if_else(var == "relate",
                         if_else(group == "0", "U", "R"),
                         if_else(var == "sex",
                                 if_else(group == "0", "M", "F"),
                                 "B")),
         var = factor(var, levels = c("null", "sex", "relate")))

# predict offspring survival using best performing model
O.surv.ha <- ggpredict(O.ha.sex.model, 
                       # predict over standardized 1st quartile,
                       # mean, and 3rd quartile values
                       c("num.helpers [all]", "sex.ratio [0, 1]",
                                           "ha [9.8, 14.1, 17.7]"),
                       bias_correction = TRUE) %>%
  data.frame() %>%
  mutate(num.helpers = x,
         group = if_else(group == "0", "M",
                        if_else(group == "1", "F", group)),
         ha = facet)

# create offspring survival plot
Off_Surv.plot <- ggplot(O.pred %>%
                            filter(num.helpers <= 5),
                          aes(x = num.helpers, y = predicted,
                              col = group, fill = group)) +
  geom_line(linewidth = 0.5) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.2, col = NA) +
  scale_color_manual(values = c(F = "#FF2B0F", M = "#009AC5", B = "#754D98",
                                R = "#2DA93C", U = "#8D5802"),
                     labels = c(F = "Female", M = "Male", B = "All",
                                R = "Related", U = "Unrelated")) +
  scale_fill_manual(values = c(F = "#FF2B0F", M = "#009AC5", B = "#754D98",
                               R = "#2DA93C", U = "#8D5802"),
                    labels = c(F = "Female", M = "Male", B = "All",
                               R = "Related", U = "Unrelated")) +
  geom_hline(data = O.pred %>%
               summarize(zero.helpers = predicted[1]),
             aes(yintercept = zero.helpers), col = "black", lty = "dashed") +
  labs(x = "# of Helpers", y = "Predicted\nOffspring Survival",
       fill = "Helpers",
       col = "Helpers") +
  facet_grid(.~var,
             labeller = labeller(var = c(null = "Total Helper #",
                                         sex = "Helper Sex",
                                         relate = "Helper\nRelatedness"))) +
  theme_minimal() +
  theme(text = element_text(size = 12),
        panel.border = element_rect(fill = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

# create offspring survival plot predicted over territory area
Off_Surv_ha.plot <- ggplot(O.surv.ha %>%
                           filter(num.helpers <= 5),
                         aes(x = num.helpers, y = predicted,
                             col = group, fill = group)) +
  geom_line(linewidth = 0.5) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.2, col = NA) +
  scale_color_manual(values = c(F = "#FF2B0F", M = "#009AC5", B = "#754D98",
                                R = "#2DA93C", U = "#8D5802"),
                     labels = c(F = "Female", M = "Male", B = "All",
                                R = "Related", U = "Unrelated")) +
  scale_fill_manual(values = c(F = "#FF2B0F", M = "#009AC5", B = "#754D98",
                               R = "#2DA93C", U = "#8D5802"),
                    labels = c(F = "Female", M = "Male", B = "All",
                               R = "Related", U = "Unrelated")) +
  geom_hline(data = group_by(O.surv.ha, ha) %>%
               summarize(zero.helpers = predicted[1]),
             aes(yintercept = zero.helpers), col = "black", lty = "dashed") +
  labs(x = "# of Helpers", y = "Predicted\nOffspring Survival",
       fill = "Helpers",
       col = "Helpers") +
  facet_grid(.~ha,
             labeller = labeller(ha = c("9.8" = "Small Territory",
                                          "14.1" = "Average Territory",
                                          "17.7" = "Large Territory"))) +
  theme_minimal() +
  theme(text = element_text(size = 12),
        panel.border = element_rect(fill = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none")

# combine offspring survival plots
Off.plot <- Off_Surv.plot/Off_Surv_ha.plot +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "A")

ggsave(filename = "figures/Figure2_20250714.pdf",
       Off.plot,
       width = 6.96,
       height = 5,
       units = "in")

ggsave(filename = "figures/Figure2_20250714.png",
       Off.plot,
       width = 6.96,
       height = 5,
       units = "in",
       dpi = 500)

# ARS plot

# predict # of offspring produced (annual reproductive success, ARS)
# using best performing base model
ARS.pred <- ggpredict(ARS.null.2.model,
                       c("num.helpers [all]"),
                       bias_correction = TRUE,
                      # prediction type incorporates zero-inflated model
                      # to predict mean # of offspring produced including
                      # the abundance of zeros
                      type = "zero_inflated") %>%
  data.frame() %>%
  mutate(var = "null", group = "B") %>%
  # repeat with sex based model
  bind_rows(ggpredict(ARS.sex.model,
                      c("num.helpers [all]",
                        "sex.ratio [0, 1]"),
                      bias_correction = TRUE,
                      type = "zero_inflated") %>%
              data.frame() %>%
              mutate(var = "sex")) %>%
  # repeat with relatedness based model
  bind_rows(ggpredict(ARS.relate.model,
                      c("num.helpers [all]",
                        "mean.relate [0, 1]"),
                      bias_correction = TRUE,
                      type = "zero_inflated") %>%
              data.frame() %>%
              mutate(var = "relate")) %>%
  mutate(num.helpers = x,
         group = if_else(var == "relate",
                         if_else(group == "0", "U", "R"),
                         if_else(var == "sex",
                                 if_else(group == "0", "M", "F"),
                                 "B")),
         var = factor(var, levels = c("null", "sex", "relate")))

# plot prediction of # of offspring produced
ARS.plot <- ggplot(ARS.pred %>%
                            filter(num.helpers <= 5),
                          aes(x = num.helpers, y = predicted,
                              col = group, fill = group)) +
  geom_line(linewidth = 0.5) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.2, col = NA) +
  scale_color_manual(values = c(F = "#FF2B0F", M = "#009AC5", B = "#754D98",
                                R = "#2DA93C", U = "#8D5802"),
                     labels = c(F = "Female", M = "Male", B = "All",
                                R = "Related", U = "Unrelated")) +
  scale_fill_manual(values = c(F = "#FF2B0F", M = "#009AC5", B = "#754D98",
                               R = "#2DA93C", U = "#8D5802"),
                    labels = c(F = "Female", M = "Male", B = "All",
                               R = "Related", U = "Unrelated")) +
  geom_hline(data = group_by(ARS.pred, var) %>%
               summarize(zero.helpers = predicted[1]),
             aes(yintercept = zero.helpers), col = "black", lty = "dashed") +
  labs(x = "# of Helpers", y = "Predicted # of Nestlings Produced",
       fill = "Helpers",
       col = "Helpers") +
  facet_grid(.~var,
             labeller = labeller(var = c(null = "Total Helper #",
                                         sex = "Helper Sex",
                                         relate = "Helper\nRelatedness"))) +
  theme_minimal() +
  theme(text = element_text(size = 12),
        panel.border = element_rect(fill = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

ggsave(filename = "figures/FigureS3_20250714.pdf",
       ARS.plot,
       width = 5,
       height = 3,
       units = "in")

ggsave(filename = "figures/FigureS3_20250714.png",
       ARS.plot,
       width = 5,
       height = 3,
       units = "in",
       dpi = 500)

### CREATE MODEL TABLES
# create key to label variables
BSurvModelKey <- c("(Intercept)" = "Intercept",
                  "age" = "Breeder Age",
                  "I(age^2)" = "Breeder Age^2",
                  "pairYear" = "Pair Experience",
                  "pedF" = "Inbreeding Coef.",
                  "ha" = "Territory Area",
                  "num.helpers" = "# of Helpers",
                  "I(num.helpers^2)" = "# of Helpers^2",
                  "num.helpers:sex.ratio" = "Helper Sex Ratio",
                  "num.helpers:mean.relate" = "Related Helper Ratio",
                  "SD (Intercept Year)" = "Year",
                  "SD (Intercept Terr)" = "Territory ID",
                  "SD (Intercept USFWS)" = "Individual ID")

# breeder survival models
modelsummary(list("Base" = B.M.null.model,
                  "Non-Linear" = B.M.null.2.model,
                  "Helper Sex" = B.M.sex.model,
                  "Helper Relatedness" = B.M.relate.model),
             stars = c("*" = 0.05, "**" = 0.01, "***" = 0.001),
             coef_map = BSurvModelKey,
             gof_omit = "ICC",
             output = "flextable") %>%
  hline(c(18, 21)) %>%
  autofit() %>%
  save_as_docx(path = "figures/BS_M_ModelTable_20250714.docx")

modelsummary(list("Base" = B.F.null.model,
                  "Non-Linear" = B.F.null.2.model,
                  "Helper Sex" = B.F.sex.model,
                  "Helper Relatedness" = B.F.relate.model),
             stars = c("*" = 0.05, "**" = 0.01, "***" = 0.001),
             coef_map = BSurvModelKey,
             gof_omit = "ICC",
             output = "flextable") %>%
  hline(c(18, 21)) %>%
  autofit() %>%
  save_as_docx(path = "figures/BS_F_ModelTable_20250714.docx")

# create key for offspring survival model
OSurvModelKey <- c("(Intercept)" = "Intercept",
                   "HatchDate" = "Hatch Day",
                   "pairYear" = "Pair Experience",
                   "pedF" = "Inbreeding Coef.",
                   "ha" = "Territory Area",
                   "num.helpers" = "# of Helpers",
                   "I(num.helpers^2)" = "# of Helpers^2",
                   "num.helpers:sex.ratio" = "Helper Sex Ratio",
                   "num.helpers:mean.relate" = "Related Helper Ratio",
                   "SD (Intercept Year)" = "Year",
                   "SD (Intercept Terr)" = "Territory ID",
                   "SD (Intercept NatalNest)" = "Nest ID")

# juvenile survival models
modelsummary(list("Base" = O.null.model,
                  "Non-Linear" = O.null.2.model,
                  "Helper Sex" = O.sex.model,
                  "Helper Relatedness" = O.relate.model),
             stars = c("*" = 0.05, "**" = 0.01, "***" = 0.001),
             coef_map = OSurvModelKey,
             gof_omit = "ICC",
             output = "flextable") %>%
  hline(c(18, 21)) %>%
  autofit() %>%
  save_as_docx(path = "figures/OSModelTable_20250714.docx")

# create key for offspring production (ARS) model
ARSModelKey <- c("(Intercept)" = "Intercept",
                 "Month" = "Breeding Start Month",
                 "pairYear" = "Pair Experience",
                 "pair.pedF" = "Breeder Relatedness",
                 "ha" = "Territory Area",
                 "num.helpers" = "# of Helpers",
                 "I(num.helpers^2)" = "# of Helpers^2",
                 "num.helpers:sex.ratio" = "Helper Sex Ratio",
                 "num.helpers:mean.relate" = "Related Helper Ratio",
                 "SD (Intercept Year)" = "Year",
                 "SD (Intercept Terr)" = "Territory ID",
                 "SD (Intercept pairID)" = "Pair ID")

modelsummary(list("Base" = ARS.null.model,
                  "Non-Linear" = ARS.null.2.model,
                  "Helper Sex" = ARS.sex.model,
                  "Helper Relatedness" = ARS.relate.model),
             stars = c("*" = 0.05, "**" = 0.01, "***" = 0.001),
             coef_map = ARSModelKey,
             gof_omit = "ICC",
             shape = component + term + statistic ~ model,
             group_map = c("conditional" = "Conditional",
                           "zero_inflated" = "Zero-Inflated",
                           "dispersion" = "Dispersion"),
             output = "flextable") %>%
  hline(c(18, 21, 24)) %>%
  autofit() %>%
  save_as_docx(path = "figures/ARSModelTable_20250714.docx")

# create model summary for post-hoc offspring survival models
OSurvModelKey.ha <- c("(Intercept)" = "Intercept",
                      "HatchDate" = "Hatch Day",
                   "pairYear" = "Pair Experience",
                   "pedF" = "Inbreeding Coef.",
                   "ha" = "Territory Area",
                   "num.helpers" = "# of Helpers",
                   "I(num.helpers^2)" = "# of Helpers^2",
                   "num.helpers:sex.ratio" = "Helper Sex Ratio",
                   "num.helpers:mean.relate" = "Related Helper Ratio",
                  "num.helpers:ha" = "# of Helpers\nTerritory Area Interaction",
            "ha:I(num.helpers^2)" = "# of Helper^2\nTerritory Area Interaction",
  "num.helpers:ha:sex.ratio" = "Helper Sex Ratio\nTerritory Area Interaction",
"num.helpers:ha:mean.relate" = "Related Helper Ratio\nTerritory Area Interaction",
                   "SD (Intercept Year)" = "Year",
                   "SD (Intercept Terr)" = "Territory ID",
                   "SD (Intercept NatalNest)" = "Nest ID")

# juvenile survival models
modelsummary(list("Base" = O.null.ha.model,
                  "Non-Linear" = O.null.2.ha.model,
                  "Helper Sex" = O.ha.sex.model,
                  "Helper Relatedness" = O.ha.relate.model),
             stars = c("*" = 0.05, "**" = 0.01, "***" = 0.001),
             coef_map = OSurvModelKey.ha,
             gof_omit = "ICC",
             output = "flextable") %>%
  hline(c(26, 29)) %>%
  autofit() %>%
  save_as_docx(path = "figures/OS_ha_ModelTable_20250714.docx")


# Make overall AIC table
## make funciton to quickly extract AIC
AIC_summary <- function(..., response){
  AIC(...) %>%
    mutate(response = response,
           model = c("Base", "Non-Linear", "Helper Sex", 
                     "Helper Relatedness"))
}

# create table
AIC_summary(B.M.null.model, B.M.null.2.model, B.M.sex.model, B.M.relate.model,
            response = "Male Breeder Survival") %>%
  bind_rows(AIC_summary(B.F.null.model, B.F.null.2.model, 
                        B.F.sex.model, B.F.relate.model,
                        response = "Female Breeder Survival")) %>%
  bind_rows(AIC_summary(O.null.model, O.null.2.model, 
                        O.sex.model, O.relate.model,
                        response = "Offspring Survival")) %>%
  bind_rows(AIC_summary(O.null.ha.model, O.null.2.ha.model, 
                        O.ha.sex.model, O.ha.relate.model,
                        response = "Offspring Survival") %>%
              mutate(model = paste("Area X", model))) %>%
  bind_rows(AIC_summary(ARS.null.model, ARS.null.2.model, 
                        ARS.sex.model, ARS.relate.model,
                        response = "Nestling Production")) %>%
  group_by(response) %>%
  # calculate delta AIC for each set of models
  mutate(dAIC = round(AIC - min(AIC), 2),
         AIC = round(AIC, 2),
         response = factor(response, levels = c("Male Breeder Survival", 
                                                "Female Breeder Survival", 
                                                "Offspring Survival", 
                                                "Nestling Production"))) %>%
  arrange(response, dAIC) %>%
  select("Response" = response, 
         "Model" = model, 
         "Degrees of Freedom" = df, AIC, dAIC) %>%
  as_grouped_data("Response") %>%
  flextable() %>%
  hline(c(5, 10, 19)) %>%
  autofit() %>%
  save_as_docx(path = "figures/AICTable_20250714.docx")


#### Helper summary statistics

# number of territories with helpers
filter(Breeders.input,
       !Year %in% c(1994, 1997, 1999),
       num.helpers > 0) %>%
  {n_distinct(.$TerrYr)}
# 787

filter(Breeders.input,
       !Year %in% c(1994, 1997, 1999)) %>%
  {n_distinct(.$TerrYr)}
# 1505

# 52.3% of territories have helpers

# average number of helpers
filter(Breeders.input,
       !Year %in% c(1994, 1997, 1999),
       num.helpers > 0) %>%
  distinct(TerrYr, .keep_all = TRUE) %>%
  {mean(.$num.helpers)}

# average number = 1.84

# proportion of related helpers overall
filter(Breeders.input,
       !Year %in% c(1994, 1997, 1999),
       num.helpers > 0) %>%
  mutate(related.helpers = mean.relate*num.helpers) %>%
  {sum(.$related.helpers)/sum(.$num.helpers)}

# 68% of helpers are related to their breeders

filter(Offspring.input %>% group_by(NatalNest) %>%
         summarize(Year = Year[1], num.helpers = num.helpers[1],
                   mean.relate = max(mean.relate)),
       !Year %in% c(1994, 1997, 1999),
       num.helpers > 0) %>%
  mutate(related.helpers = mean.relate*num.helpers) %>%
  {sum(.$related.helpers)/sum(.$num.helpers)}

# 84% of helpers are related to their breeders


### Correlations with helper presence

# check for differences in relatedness between helpers and male and females
wilcox.test(mean.relate ~ Sex, 
            Breeders.input %>% 
              filter(!Year %in% c(1994, 1997, 1999)) %>% 
              filter(!is.na(mean.relate)))
# W = 273546, p-value = 0.0005

# calculate mean relatedness between breeders and helpers
group_by(Breeders.input %>% 
           filter(!Year %in% c(1994, 1997, 1999)) %>% 
           filter(!is.na(mean.relate)), Sex) %>%
  summarize(relation.mean = mean(mean.relate), relation.sd = sd(mean.relate))
# F mean = 0.620, sd = 0.452
# M mean = 0.699, sd = 0.421

# check correlation between breeder age and number of helpers
cor.test(filter(Breeders.input %>% 
                  filter(!Year %in% c(1994, 1997, 1999)))$age, 
         filter(Breeders.input %>% 
                  filter(!Year %in% c(1994, 1997, 1999)))$num.helpers,
    method = "spearman")
# rho = 0.2352, p-value < 2.2e-16

# check correlation between breeder age and proportion of related helpers
cor.test(filter(Breeders.input %>% 
                  filter(!Year %in% c(1994, 1997, 1999)))$age, 
         filter(Breeders.input %>% 
                  filter(!Year %in% c(1994, 1997, 1999)))$mean.relate,
         method = "spearman")
# rho = 0.4434, p-value < 2.2e-16

# check correlation between breeder age and sex ratio of helpers
cor.test(filter(Breeders.input %>% 
                  filter(!Year %in% c(1994, 1997, 1999)))$age, 
         filter(Breeders.input %>% 
                  filter(!Year %in% c(1994, 1997, 1999)))$sex.ratio,
         method = "spearman")
# rho = -0.0382, p-value = 0.1329

# check correlation between territory size and number of helpers
cor.test(filter(Breeders.input %>% 
                  filter(!Year %in% c(1994, 1997, 1999)))$ha, 
         filter(Breeders.input %>% 
                  filter(!Year %in% c(1994, 1997, 1999)))$num.helpers,
         method = "spearman")
# rho = 0.1935, p-value < 2.2e-16

# check correlation between territory size and proportion of related helpers
cor.test(filter(Breeders.input %>% 
                  filter(!Year %in% c(1994, 1997, 1999)))$ha, 
         filter(Breeders.input %>% 
                  filter(!Year %in% c(1994, 1997, 1999)))$mean.relate,
         method = "spearman")
# rho = 0.0381, p-value = 0.134

# check correlation between territory size and sex ratio of helpers
cor.test(filter(Breeders.input %>% 
                  filter(!Year %in% c(1994, 1997, 1999)))$ha, 
         filter(Breeders.input %>% 
                  filter(!Year %in% c(1994, 1997, 1999)))$sex.ratio,
         method = "spearman")
# rho = 0.0656, p-value = 0.0098

# check correlation between number of helpers and
# proportion of related helpers to nestlings
group_by(Offspring.input %>% 
           filter(!Year %in% c(1994, 1997, 1999)), NatalNest) %>%
  summarize(num.helpers = num.helpers[1], mean.relate = max(mean.relate)) %>%
  {cor.test(.$num.helpers, .$mean.relate, method = "spearman")}
# rho = -0.0568, p-value = 0.2199

# check correlation between number of helpers and hatch date of nestlings
group_by(Offspring.input %>% 
           filter(!Year %in% c(1994, 1997, 1999)), NatalNest) %>%
  summarize(num.helpers = num.helpers[1], HatchDate = HatchDate[1]) %>%
  {cor.test(.$num.helpers, .$HatchDate, method = "spearman")}
# rho = -0.1317, p-value = 2.152e-06

# check correlation between proportion of related helpers and hatch date of
# nestlings
group_by(Offspring.input %>% 
           filter(!Year %in% c(1994, 1997, 1999)), NatalNest) %>%
  summarize(mean.relate = mean.relate[1], HatchDate = HatchDate[1]) %>%
  {cor.test(.$mean.relate, .$HatchDate, method = "spearman")}
# rho = -0.0101, p-value = 0.8275

# check correlation between sex ratio of helpers and hatch date of
# nestlings
group_by(Offspring.input %>% 
           filter(!Year %in% c(1994, 1997, 1999)), NatalNest) %>%
  summarize(sex.ratio = sex.ratio[1], HatchDate = HatchDate[1]) %>%
  {cor.test(.$sex.ratio, .$HatchDate, method = "spearman")}
# rho = 0.0575, p-value = 0.2143

### SUPPLEMENTAL FIGURES

## Create plot showing offspring survival with territory area and distribution
Off.surv.ha.plot <- ggplot(Offspring.input %>% 
                             filter(!Year %in% c(1994, 1997, 1999)) %>%
                             group_by(NatalNest) %>%
                             summarize(num.helpers = num.helpers[1],
                                       ha = ha[1],
                                       surv = mean(surv),
                                       num.juv = n_distinct(USFWS))) +
  geom_point(aes(x = ha, y = surv)) +
  geom_smooth(aes(x = ha, y = surv, weight = num.juv)) +
  # scale y axis of histogram to fit under data
  geom_histogram(aes(x = ha, y = 2*after_stat(density)), 
                 binwidth = 1, alpha = 0.6) +
  labs(x = "Territory Area (ha)", y = "Mean Offspring Survival\nPer Nest") +
  # add vertical lines indicating quartiles
  geom_vline(xintercept = c(9.8, 14.1, 17.7)) +
  theme_minimal() +
  theme(text = element_text(size = 12),
        panel.border = element_rect(fill = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

ggsave(filename = "figures/FigureS2_20250714.pdf",
       Off.surv.ha.plot,
       width = 5,
       height = 3,
       units = "in")

ggsave(filename = "figures/FigureS2_20250714.png",
       Off.surv.ha.plot,
       width = 5,
       height = 3,
       units = "in",
       dpi = 500)

# create figure that shows distribution of helpers
HelperCountData <- Breeders.input %>% 
  filter(!Year %in% c(1994, 1997, 1999)) %>%
  replace_na(list(sex.ratio = 0, mean.relate = 0)) %>%
  # calculate number of helpers of each sex and that are related or unrelated
  mutate(female.helpers = num.helpers*sex.ratio,
         male.helpers = num.helpers*(1-sex.ratio),
         related.helpers = num.helpers*mean.relate,
         unrelated.helpers = num.helpers*(1-mean.relate)) %>%
  # combine with information calculated from offspring survival data
  bind_rows(group_by(Offspring.data, 
                     Year, USFWS = NatalNest, TerrYr, Terr) %>%
              summarize(num.helpers = num.helpers[1],
                        sex.ratio = sex.ratio[1],
                        mean.relate = mean.relate[1],
                        female.helpers = num.helpers*sex.ratio,
                        male.helpers = num.helpers*(1-sex.ratio),
                        related.helpers = num.helpers*mean.relate,
                        unrelated.helpers = num.helpers*(1-mean.relate)) %>%
              mutate(Sex = "N")) %>%
  pivot_longer(c(female.helpers, male.helpers), names_to = "HSex",
               values_to = "sex.helpers") %>%
  pivot_longer(c(related.helpers, unrelated.helpers), names_to = "Hrelate",
               values_to = "relate.helpers")

# plot distribution of male and female helper counts
HelpCount.Sex <- ggplot(HelperCountData, aes(x = as.factor(sex.helpers), 
                                             fill = HSex)) + 
  geom_bar(position = position_dodge()) +
  labs(x = "# of Helpers", y = "# of Records", fill = "Helper Sex") +
  scale_fill_manual(values = c(female.helpers = "#FF2B0F", male.helpers = "#009AC5"),
                     labels = c(female.helpers = "Female", male.helpers = "Male")) +
  facet_wrap(.~Sex, labeller = labeller(Sex = c(F = "Female Breeders",
                                               M = "Male Breeders",
                                               N = "Nests with Nestlings"))) +
  theme_minimal() +
  theme(text = element_text(size = 12),
        panel.border = element_rect(fill = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

# plot distribution of related and unrelated helper counts
HelpCount.relate <- ggplot(HelperCountData, aes(x = as.factor(relate.helpers), 
                                                fill = Hrelate)) + 
  geom_bar(position = position_dodge()) +
  labs(x = "# of Helpers", y = "# of Records", fill = "Helper Relatedness") +
  scale_fill_manual(values = c(unrelated.helpers = "#8D5802", 
                               related.helpers = "#2DA93C"),
                    labels = c(unrelated.helpers = "Unrelated", 
                               related.helpers = "Related")) +
  facet_wrap(.~Sex, labeller = labeller(Sex = c(F = "Female Breeders",
                                                M = "Male Breeders",
                                                N = "Nests with Nestlings"))) +
  theme_minimal() +
  theme(text = element_text(size = 12),
        panel.border = element_rect(fill = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

# combine helper count plots
HelpCount.plot <- HelpCount.Sex/HelpCount.relate + plot_annotation(tag_levels = "A")


ggsave(filename = "figures/FigureS1_20250714.pdf",
       HelpCount.plot,
       width = 6.96,
       height = 5,
       units = "in")

ggsave(filename = "figures/FigureS1_20250714.png",
       HelpCount.plot,
       width = 6.5,
       height = 5,
       units = "in",
       dpi = 500)
