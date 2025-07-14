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
# predict for best performing full female model
Surv.pred <- ggpredict(B.F.model,
                         c("num.helpers [all]",
                           # predict between scenarios with all male (0)
                           # and all female (1) helpers
                           "sex.ratio [0, 1]",
                           # predict between scenarios with all unrelated (0)
                           # and all related (1) helpers
                           "mean.relate [0, 1]"),
                       # bias_correction adjusts values by random effects
                         bias_correction = TRUE) %>%
  data.frame() %>%
  mutate(Bsex = "F") %>%
  # repeat for best performing full male model
  bind_rows(ggpredict(B.M.2.model,
                      c("num.helpers [all]",
                        "sex.ratio [0, 1]",
                        "mean.relate [0, 1]"),
                      bias_correction = TRUE) %>%
              data.frame() %>%
              mutate(Bsex = "M")) %>%
  # predict for best performing full offspring model
  bind_rows(ggpredict(O.model,
                      c("num.helpers [all]",
                        "sex.ratio [0, 1]",
                        "mean.relate [0, 1]"),
                      bias_correction = TRUE) %>%
              data.frame() %>%
              mutate(Bsex = "A")) %>%
  # predict for base male breeder model
  bind_rows(ggpredict(B.M.null.model,
                      c("num.helpers [all]"),
                      bias_correction = TRUE) %>%
              data.frame() %>%
              mutate(Bsex = "M", group = "B")) %>%
  # predict for base female breeder model
  bind_rows(ggpredict(B.F.null.model,
                      c("num.helpers [all]"),
                      bias_correction = TRUE) %>%
              data.frame() %>%
              mutate(Bsex = "F", group = "B")) %>%
  # predict for base offspring model
  bind_rows(ggpredict(O.null.model,
                      c("num.helpers [all]"),
                      bias_correction = TRUE) %>%
              data.frame() %>%
              mutate(Bsex = "A", group = "B", )) %>%
  # adjust labels for plotting
  mutate(num.helpers = x,
         Hsex = if_else(group == "0", "M",
                        if_else(group == "1", "F", group)),
         mean.relate = if_else(is.na(facet), "A",
                               if_else(facet == "0", "U", "R")))

# plot breeder survival prediction curves
Breed_Surv.plot <- ggplot(Surv.pred %>%
                            # limit number of helpers to maximum
                            # number of single-sex helpers observed (4)
                      filter(num.helpers <= 4, Bsex != "A"),
                    aes(x = num.helpers, y = predicted,
                        col = Hsex, fill = Hsex)) +
  geom_line(linewidth = 0.5) +
  # add 95% confidence intervals
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.2, col = NA) +
  scale_color_manual(values = c(F = "#FF2B0F", M = "#009AC5", B = "#754D98"),
                     labels = c(F = "Female", M = "Male", B = "Both")) +
  scale_fill_manual(values = c(F = "#FF2B0F", M = "#009AC5", B = "#754D98"),
                    labels = c(F = "Female", M = "Male", B = "Both")) +
  # add horizontal line showing predicted survival with 0 helpers
  geom_hline(data = group_by(Surv.pred %>% filter(Bsex != "A"), Bsex) %>%
               summarize(zero.helpers = predicted[1]),
             aes(yintercept = zero.helpers), col = "black", lty = "dashed") +
  labs(x = "# of Helpers", y = "Predicted Breeder Survival",
       fill = "Helper Sex",
       col = "Helper Sex") +
  facet_grid(Bsex~mean.relate,
             labeller = labeller(Bsex = c(F = "Female Breeder Survival",
                                          M = "Male Breeder Survival",
                                          A = "Juvenile Survival"),
                                 mean.relate = c(U = "Unrelated Helpers",
                                                 R = "Related Helpers",
                                                 A = "All Helpers"))) +
  theme_minimal() +
  theme(text = element_text(size = 9),
        panel.border = element_rect(fill = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

ggsave(filename = "figures/Figure1_20250714.pdf",
       Breed_Surv.plot,
       width = 6.5,
       height = 5,
       units = "in")

ggsave(filename = "figures/Figure1_20250714.png",
       Breed_Surv.plot,
       width = 6.5,
       height = 5,
       units = "in",
       dpi = 500)

# Follow-up offspring territory area analyses

# predict offspring survival based on territory area
O.surv.ha <- ggpredict(O.ha.sex.model, 
                       # predict over standardized 1st quartile,
                       # mean, and 3rd quartile values
                       c("num.helpers [all]", "sex.ratio [0, 1]",
                                           "ha [-0.717, 0, 0.583]"),
                       bias_correction = TRUE) %>%
  data.frame() %>%
  # predict offspring survival based on territory area using base model
  bind_rows(ggpredict(O.null.ha.model,
                      c("num.helpers [all]", "ha [-0.717, 0, 0.583]"),
                      bias_correction = TRUE) %>%
              data.frame() %>%
              mutate(facet = group, group = "B")) %>%
  mutate(num.helpers = x,
         Hsex = if_else(group == "0", "M",
                        if_else(group == "1", "F", group)),
         ha = facet)

# create offspring survival plot
Off_Surv.plot <- ggplot(Surv.pred %>%
                            filter(num.helpers <= 4, Bsex == "A"),
                          aes(x = num.helpers, y = predicted,
                              col = Hsex, fill = Hsex)) +
  geom_line(linewidth = 0.5) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.2, col = NA) +
  scale_color_manual(values = c(F = "#FF2B0F", M = "#009AC5", B = "#754D98"),
                     labels = c(F = "Female", M = "Male", B = "Both")) +
  scale_fill_manual(values = c(F = "#FF2B0F", M = "#009AC5", B = "#754D98"),
                    labels = c(F = "Female", M = "Male", B = "Both")) +
  geom_hline(data = group_by(Surv.pred %>% filter(Bsex == "A"), Bsex) %>%
               summarize(zero.helpers = predicted[1]),
             aes(yintercept = zero.helpers), col = "black", lty = "dashed") +
  labs(x = "# of Helpers", y = "Predicted Offspring Survival",
       fill = "Helper Sex",
       col = "Helper Sex") +
  facet_grid(.~mean.relate,
             labeller = labeller(mean.relate = c(U = "Unrelated Helpers",
                                                 R = "Related Helpers",
                                                 A = "All Helpers"))) +
  theme_minimal() +
  theme(text = element_text(size = 9),
        panel.border = element_rect(fill = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

# create offspring survival plot predicted over territory area
Off_Surv_ha.plot <- ggplot(O.surv.ha %>%
                           filter(num.helpers <= 4, Hsex != "B"),
                         aes(x = num.helpers, y = predicted,
                             col = Hsex, fill = Hsex)) +
  geom_line(linewidth = 0.5) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.2, col = NA) +
  scale_color_manual(values = c(F = "#FF2B0F", M = "#009AC5", B = "#754D98"),
                     labels = c(F = "Female", M = "Male", B = "Both")) +
  scale_fill_manual(values = c(F = "#FF2B0F", M = "#009AC5", B = "#754D98"),
                    labels = c(F = "Female", M = "Male", B = "Both")) +
  geom_hline(data = group_by(O.surv.ha, ha) %>%
               summarize(zero.helpers = predicted[1]),
             aes(yintercept = zero.helpers), col = "black", lty = "dashed") +
  labs(x = "# of Helpers", y = "Predicted Offspring Survival",
       fill = "Helper Sex",
       col = "Helper Sex") +
  facet_grid(.~ha,
             labeller = labeller(ha = c("-0.717" = "Small Territory",
                                          "0" = "Average Territory",
                                          "0.583" = "Large Territory"))) +
  theme_minimal() +
  theme(text = element_text(size = 9),
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
       width = 6.5,
       height = 5,
       units = "in")

ggsave(filename = "figures/Figure2_20250714.png",
       Off.plot,
       width = 6.5,
       height = 5,
       units = "in",
       dpi = 500)

# ARS plot

# predict # of offspring produced (annual reproductive success, ARS)
# using best performing model
ARS.pred <- ggpredict(ARS.null.2.model,
                       c("num.helpers [all]"),
                       bias_correction = TRUE,
                      # prediction type incorporates zero-inflated model
                      # to predict mean # of offspring produced including
                      # the abundance of zeros
                      type = "zero_inflated") %>%
  data.frame() %>%
  mutate(model = "B", group = "B") %>%
  # repeat with best performing full model
  bind_rows(ggpredict(ARS.2.model,
                      c("num.helpers [all]",
                        "sex.ratio [0, 1]",
                        "mean.relate [0, 1]"),
                      bias_correction = TRUE,
                      type = "zero_inflated") %>%
              data.frame() %>%
              mutate(model = "F")) %>%
  mutate(num.helpers = x,
         Hsex = if_else(group == "0", "M",
                        if_else(group == "1", "F", group)),
         mean.relate = if_else(is.na(facet), "A",
                               if_else(facet == "0", "U", "R")))

# plot prediction of # of offspring produced
ARS.plot <- ggplot(ARS.pred %>%
                            filter(num.helpers <= 4),
                          aes(x = num.helpers, y = predicted,
                              col = Hsex, fill = Hsex)) +
  geom_line(linewidth = 0.5) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.2, col = NA) +
  scale_color_manual(values = c(F = "#FF2B0F", M = "#009AC5", B = "#754D98"),
                     labels = c(F = "Female", M = "Male", B = "Both")) +
  scale_fill_manual(values = c(F = "#FF2B0F", M = "#009AC5", B = "#754D98"),
                    labels = c(F = "Female", M = "Male", B = "Both")) +
  geom_hline(data = group_by(ARS.pred, model) %>%
               summarize(zero.helpers = predicted[1]),
             aes(yintercept = zero.helpers), col = "black", lty = "dashed") +
  labs(x = "# of Helpers", y = "Predicted # of Nestlings Produced",
       fill = "Helper Sex",
       col = "Helper Sex") +
  facet_grid(.~mean.relate,
             labeller = labeller(mean.relate = c(U = "Unrelated Helpers",
                                                 R = "Related Helpers",
                                                 A = "All Helpers"))) +
  theme_minimal() +
  theme(text = element_text(size = 9),
        panel.border = element_rect(fill = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

ggsave(filename = "figures/FigureS3_20250714.pdf",
       ARS.plot,
       width = 6.5,
       height = 3,
       units = "in")

ggsave(filename = "figures/FigureS3_20250714.png",
       ARS.plot,
       width = 6.5,
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
                  "num.helpers:mean.relate" = "Related Helper Ratio",
                  "I(num.helpers^2):mean.relate" = "Related Helper Ratio\nNon-Linear Term",
                  "num.helpers:sex.ratio" = "Helper Sex Ratio",
                  "I(num.helpers^2):sex.ratio" = "Helper Sex Ratio\n Non-Linear Term",
                  "SD (Intercept Year)" = "Year",
                  "SD (Intercept Terr)" = "Territory ID",
                  "SD (Intercept USFWS)" = "Individual ID")

# breeder survival models
modelsummary(list("Linear" = B.M.null.model,
                  "Non-Linear" = B.M.null.2.model,
                  "Linear" = B.M.model,
                  "Non-Linear" = B.M.2.model),
             stars = c("*" = 0.05, "**" = 0.01, "***" = 0.001),
             coef_map = BSurvModelKey,
             gof_omit = "ICC",
             output = "flextable") %>%
  hline(c(22, 25)) %>%
  add_header_row(top = TRUE, values = c("", "Base Model", "Full Model"), 
                 colwidths = c(1, 2, 2)) %>%
  autofit() %>%
  save_as_docx(path = "figures/BS_M_ModelTable_20250714.docx")

modelsummary(list("Linear" = B.F.null.model,
                  "Non-Linear" = B.F.null.2.model,
                  "Linear" = B.F.model,
                  "Non-Linear" = B.F.2.model),
             stars = c("*" = 0.05, "**" = 0.01, "***" = 0.001),
             coef_map = BSurvModelKey,
             gof_omit = "ICC",
             output = "flextable") %>%
  hline(c(22, 25)) %>%
  add_header_row(top = TRUE, values = c("", "Base Model", "Full Model"), 
                 colwidths = c(1, 2, 2)) %>%
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
                   "num.helpers:mean.relate" = "Related Helper Ratio",
                   "I(num.helpers^2):mean.relate" = "Related Helper Ratio\nNon-Linear Term",
                   "num.helpers:sex.ratio" = "Helper Sex Ratio",
                   "I(num.helpers^2):sex.ratio" = "Helper Sex Ratio\n Non-Linear Term",
                   "SD (Intercept Year)" = "Year",
                   "SD (Intercept Terr)" = "Territory ID",
                   "SD (Intercept NatalNest)" = "Nest ID")

# juvenile survival models
modelsummary(list("Linear" = O.null.model,
                  "Non-Linear" = O.null.2.model,
                  "Linear" = O.model,
                  "Non-Linear" = O.2.model),
             stars = c("*" = 0.05, "**" = 0.01, "***" = 0.001),
             coef_map = OSurvModelKey,
             gof_omit = "ICC",
             output = "flextable") %>%
  hline(c(22, 27)) %>%
  add_header_row(top = TRUE, values = c("", "Base Model", "Full Model"), 
                 colwidths = c(1, 2, 2)) %>%
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
                 "num.helpers:mean.relate" = "Related Helper Ratio",
                 "I(num.helpers^2):mean.relate" = "Related Helper Ratio\nNon-Linear Term",
                 "num.helpers:sex.ratio" = "Helper Sex Ratio",
                 "I(num.helpers^2):sex.ratio" = "Helper Sex Ratio\n Non-Linear Term",
                 "SD (Intercept Year)" = "Year",
                 "SD (Intercept Terr)" = "Territory ID",
                 "SD (Intercept pairID)" = "Pair ID")

modelsummary(list("Linear" = ARS.null.model,
                  "Non-Linear" = ARS.null.2.model,
                  "Linear" = ARS.model,
                  "Non-Linear" = ARS.2.model),
             stars = c("*" = 0.05, "**" = 0.01, "***" = 0.001),
             coef_map = ARSModelKey,
             gof_omit = "ICC",
             shape = component + term + statistic ~ model,
             group_map = c("conditional" = "Conditional",
                           "zero_inflated" = "Zero-Inflated",
                           "dispersion" = "Dispersion"),
             output = "flextable") %>%
  hline(c(22, 25, 28)) %>%
  add_header_row(top = TRUE, values = c("", "", "Base Model", "Full Model"), 
                 colwidths = c(1, 1, 2, 2)) %>%
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
                   "num.helpers:mean.relate" = "Related Helper Ratio",
                   "I(num.helpers^2):mean.relate" = "Related Helper Ratio\nNon-Linear Term",
                   "num.helpers:sex.ratio" = "Helper Sex Ratio",
                   "I(num.helpers^2):sex.ratio" = "Helper Sex Ratio\n Non-Linear Term",
                  "num.helpers:ha" = "# of Helpers\nTerritory Area Interaction",
"num.helpers:ha:mean.relate" = "Related Helper Ratio\nTerritory Area Interaction",
"num.helpers:ha:sex.ratio" = "Helper Sex Ratio\nTerritory Area Interaction",
                   "SD (Intercept Year)" = "Year",
                   "SD (Intercept Terr)" = "Territory ID",
                   "SD (Intercept NatalNest)" = "Nest ID")

# juvenile survival models
modelsummary(list("Null Model" = O.null.ha.model,
                  "Full Model" = O.ha.model,
                  "Relatedness\nModel" = O.ha.relate.model,
                  "Sex Model" = O.ha.sex.model),
             stars = c("*" = 0.05, "**" = 0.01, "***" = 0.001),
             coef_map = OSurvModelKey.ha,
             gof_omit = "ICC",
             output = "flextable") %>%
  hline(c(22, 25)) %>%
  autofit() %>%
  save_as_docx(path = "figures/OS_ha_ModelTable_20250714.docx")



#### Helper summary statistics

# number of territories with helpers
filter(Breeders.input,
       !Year %in% c(1994, 1997, 1999),
       num.helpers > 0) %>%
  {n_distinct(.$TerrYr)}
# 784

filter(Breeders.input,
       !Year %in% c(1994, 1997, 1999)) %>%
  {n_distinct(.$TerrYr)}
# 1488

# 52.6% of territories have helpers

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

### Correlations with helper presence

# check for differences in relatedness between helpers and male and females
wilcox.test(mean.relate ~ Sex, 
            Breeders.input %>% 
              filter(!Year %in% c(1994, 1997, 1999)) %>% 
              filter(!is.na(mean.relate)))
# W = 265401, p-value = 0.0009

# calculate mean relatedness between breeders and helpers
group_by(Breeders.input %>% 
           filter(!Year %in% c(1994, 1997, 1999)) %>% 
           filter(!is.na(mean.relate)), Sex) %>%
  summarize(relation.mean = mean(mean.relate), relation.sd = sd(mean.relate))
# F mean = 0.624, sd = 0.450
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
# rho = 0.4380, p-value < 2.2e-16

# check correlation between breeder age and sex ratio of helpers
cor.test(filter(Breeders.input %>% 
                  filter(!Year %in% c(1994, 1997, 1999)))$age, 
         filter(Breeders.input %>% 
                  filter(!Year %in% c(1994, 1997, 1999)))$sex.ratio,
         method = "spearman")
# rho = -0.0419, p-value = 0.1015

# check correlation between territory size and number of helpers
cor.test(filter(Breeders.input %>% 
                  filter(!Year %in% c(1994, 1997, 1999)))$ha, 
         filter(Breeders.input %>% 
                  filter(!Year %in% c(1994, 1997, 1999)))$num.helpers,
         method = "spearman")
# rho = 0.1882, p-value < 2.2e-16

# check correlation between territory size and proportion of related helpers
cor.test(filter(Breeders.input %>% 
                  filter(!Year %in% c(1994, 1997, 1999)))$ha, 
         filter(Breeders.input %>% 
                  filter(!Year %in% c(1994, 1997, 1999)))$mean.relate,
         method = "spearman")
# rho = 0.0359, p-value = 0.1502

# check correlation between territory size and sex ratio of helpers
cor.test(filter(Breeders.input %>% 
                  filter(!Year %in% c(1994, 1997, 1999)))$ha, 
         filter(Breeders.input %>% 
                  filter(!Year %in% c(1994, 1997, 1999)))$sex.ratio,
         method = "spearman")
# rho = 0.0647, p-value = 0.01146

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
  theme(text = element_text(size = 9),
        panel.border = element_rect(fill = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

ggsave(filename = "figures/FigureS2_20250714.pdf",
       Off.surv.ha.plot,
       width = 6.5,
       height = 3,
       units = "in")

ggsave(filename = "figures/FigureS2_20250714.png",
       Off.surv.ha.plot,
       width = 6.5,
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
  theme(text = element_text(size = 9),
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
  theme(text = element_text(size = 9),
        panel.border = element_rect(fill = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

# combine helper count plots
HelpCount.plot <- HelpCount.Sex/HelpCount.relate + plot_annotation(tag_levels = "A")


ggsave(filename = "figures/FigureS1_20250714.pdf",
       HelpCount.plot,
       width = 6.5,
       height = 5,
       units = "in")

ggsave(filename = "figures/FigureS1_20250714.png",
       HelpCount.plot,
       width = 6.5,
       height = 5,
       units = "in",
       dpi = 500)
