library(tidyverse)
library(ggpubr)
library(lmerTest)
library(emmeans)
cbbPalette <- c(
  "#E69F00", # orange
  "#56B4E9", # sky blue
  "#009E73", # bluish green
  "#F0E442", # yellow
  "#0072B2", # blue
  "#D55E00", # vermilion
  "#CC79A7",   # reddish purple
  "#000000"# black
)

# 2) pick two that you like – e.g. orange & sky‐blue:
myCols <- c(
  "Inoculated"     = cbbPalette[1],
  "Non-inoculated" = cbbPalette[2]
)
#Read
combined <- read.csv("Final_in_planta_validation_sheet.csv",
                     na.strings = c("", "NA", "#NUM!")) %>%
  rename(
    Kec12_qpcr = kec12_LOG.1.transformcorrected,
    Kec13_qpcr = kec13_LOG.1.transformcorrected,
    CFU = LOG.1.transform            # Rename CFU log column
  ) %>%
  mutate(
    Trial       = factor(Trial),                    # Convert Trial to factor for random effect
    Chamber     = factor(Chamber),
    Replication = factor(Replication),              # Replication as factor for grouping
    Time_days   = as.numeric(Time_days),            # Make sure Time is numeric
    Time_days_f = factor(Time_days)                 # Also prepare Time as factor for ANOVA
  )


str(combined[, c("Treatment_name","Trial","Chamber","Time_days","CFU","Kec12_qpcr","Kec13_qpcr")])
summary(combined$CFU)


# ==== Plot KEC12 ==== Final###
plot_kec12 <- combined %>%
  select(CFU, Kec12_qpcr, Treatment_name, Trial, Replication, Time_days) %>%
  pivot_longer(cols = c(CFU, Kec12_qpcr), names_to = "Metric", values_to = "value") %>%
  mutate(Method = recode(Metric,
                         "CFU" = "Plating (CFU)",
                         "Kec12_qpcr" = "qPCR KEC12")) %>%
  ggplot(aes(x = Time_days,
             y = value,
             color = Treatment_name,
             shape = Method,
             group = interaction(Treatment_name, Method))) +
  stat_summary(fun      = mean,      geom = "point",    size  = 2.5) +  # mean points
  stat_summary(fun.data = mean_se,   geom = "errorbar", width = 0.0) +  # SE bars
  stat_summary(fun      = mean,      geom = "line",     size  = 0.8, linetype="solid") +  # mean lines
  scale_shape_discrete(
    name   = "Method",
    labels = c("CFU", "qPCR"))+
  scale_color_manual(
    name   = "Treatment",
    values = myCols
  ) +
  labs(
    x        = "Days after inoculation",
    y        = expression(Log[10]~"DH44 per g root"),
    color    = "Treatment",
    shape    = "Method",
    linetype = "Method"
  ) +
  theme_classic(base_size = 10) +
  ggtitle("KEC12") +
  theme(
    plot.title = element_text(hjust = 0.5)  # <- center the title
  )
plot_kec12

# ==== Plot KEC13 ====
plot_kec13 <- combined %>%
  select(CFU, Kec13_qpcr, Treatment_name, Trial, Replication, Time_days) %>%
  pivot_longer(cols = c(CFU, Kec13_qpcr), names_to = "Metric", values_to = "value") %>%
  mutate(Method = recode(Metric,
                         "CFU" = "Plating (CFU)",
                         "Kec13_qpcr" = "qPCR KEC13")) %>%
  ggplot(aes(x = Time_days,
             y = value,
             color = Treatment_name,
             shape = Method,
             group = interaction(Treatment_name, Method))) +
  stat_summary(fun      = mean,      geom = "point",    size  = 2.5) +  # mean points
  stat_summary(fun.data = mean_se,   geom = "errorbar", width = 0.0) +  # SE bars
  stat_summary(fun      = mean,      geom = "line",     size  = 0.8, linetype="solid") +  # mean lines
  scale_shape_discrete(
    name   = "Method",
    labels = c("CFU", "qPCR"))+
  scale_color_manual(
    name   = "Treatment",
    values = myCols
  ) +
  labs(
    x        = "Days after inoculation",
    y        = expression(Log[10]~"DH44 per g root"),
    color    = "Treatment",
    shape    = "Method",
    linetype = "Method"
  ) +
  theme_classic(base_size = 10) +
  ggtitle("KEC13") +
  theme(
    plot.title = element_text(hjust = 0.5)  # <- center the title
  )

plot_kec13


# KEC13
# ==========================================================
kec13.lm <- combined %>%
  select(CFU, Kec13_qpcr, Treatment_name, Replication, Trial, Chamber, Time_days) %>%
  pivot_longer(cols = c(CFU, Kec13_qpcr),
               names_to = "method") %>%
  mutate(Time_days = as.factor(Time_days))   # Treat Time as factor

# Fit the model with Trial as random effect
mod_kec13 <- lmer(value ~ method * Time_days * Treatment_name + (1 | Trial), data = kec13.lm)
summary(mod_kec13)
anova(mod_kec13)
write.csv(as.data.frame(anova(mod_kec13)), "SuppTable1_anova_KEC13.csv")
VarCorr(mod_kec13)
confint(mod_kec13, method = "Wald")

# Estimated marginal means — pooled over time
emm_kec13 <- emmeans(mod_kec13, ~ method | Treatment_name)

emm_kec13
pairs(emm_kec13, adjust = "tukey")
confint(pairs(emm_kec13, adjust = "tukey"))
write.csv(as.data.frame(confint(pairs(emm_kec13, adjust = "tukey"))), "Table2_KEC13.csv")

# Per-timepoint contrasts (supplementary)
pairs(emmeans(mod_kec13, ~ method | Time_days | Treatment_name), adjust = "tukey")
write.csv(as.data.frame(pairs(emmeans(mod_kec13, ~ method | Time_days | Treatment_name), adjust = "tukey")), "TableS3_KEC13.csv")

# ==========================================================
# KEC12
# ==========================================================
kec12.lm <- combined %>%
  select(CFU, Kec12_qpcr, Treatment_name, Replication, Trial, Chamber, Time_days) %>%
  pivot_longer(cols = c(CFU, Kec12_qpcr),
               names_to = "method") %>%
  mutate(Time_days = as.factor(Time_days))

mod_kec12 <- lmer(value ~ method * Time_days * Treatment_name + (1 | Trial), data = kec12.lm)
summary(mod_kec12)
anova(mod_kec12)
write.csv(as.data.frame(anova(mod_kec12)), "SuppTable1_anova_KEC12.csv")
VarCorr(mod_kec12)
confint(mod_kec12, method = "Wald")

emm_kec12 <- emmeans(mod_kec12, ~ method | Treatment_name)
emm_kec12
pairs(emm_kec12, adjust = "tukey")
confint(pairs(emm_kec12, adjust = "tukey"))

write.csv(as.data.frame(confint(pairs(emm_kec12, adjust = "tukey"))), "Table2_KEC12.csv")
pairs(emmeans(mod_kec12, ~ method | Time_days | Treatment_name), adjust = "tukey")
write.csv(as.data.frame(pairs(emmeans(mod_kec12, ~ method | Time_days | Treatment_name), adjust = "tukey")), "TableS3_KEC12.csv")

# KEC12 vs KEC13 — do the two primer pairs differ?
# ==========================================================
all.lm <- combined %>%
  select(CFU, Kec12_qpcr, Kec13_qpcr, Treatment_name, Replication,
         Trial, Chamber, Time_days) %>%
  pivot_longer(cols = c(CFU, Kec12_qpcr, Kec13_qpcr), names_to = "assay") %>%
  mutate(Time_days = as.factor(Time_days),
         assay = factor(assay, levels = c("CFU", "Kec12_qpcr", "Kec13_qpcr")))

mod_all <- lmer(value ~ assay * Time_days * Treatment_name + (1 | Trial), data = all.lm)
anova(mod_all)

emm_all <- emmeans(mod_all, ~ assay | Treatment_name)
emm_all
pairs(emm_all, adjust = "tukey")          # KEC12 vs KEC13 is the third contrast
confint(pairs(emm_all, adjust = "tukey"))




# Chamber variance check — nested models, not reported.
# Fitted only to justify omitting chamber from the final models.
# ==========================================================
chk_kec13 <- lmer(value ~ method * Time_days * Treatment_name + (1 | Trial/Chamber),
                  data = kec13.lm)
summary(chk_kec13)
as.data.frame(VarCorr(chk_kec13))[, c("grp", "vcov")]
isSingular(chk_kec13)

chk_kec12 <- lmer(value ~ method * Time_days * Treatment_name + (1 | Trial/Chamber),
                  data = kec12.lm)
summary(chk_kec12)
as.data.frame(VarCorr(chk_kec12))[, c("grp", "vcov")]
isSingular(chk_kec12)


#model check

anova(chk_kec13, mod_kec13, refit = FALSE)
anova(chk_kec12, mod_kec12, refit = FALSE)

# Tests each random effect by dropping it
ranova(chk_kec13)
ranova(chk_kec12)

write.csv(as.data.frame(VarCorr(chk_kec12)), "S2_var_KEC12.csv")
write.csv(as.data.frame(VarCorr(chk_kec13)), "S2_var_KEC13.csv")
write.csv(as.data.frame(anova(chk_kec12, mod_kec12, refit = FALSE)), "S2_LRT_KEC12.csv")
write.csv(as.data.frame(anova(chk_kec13, mod_kec13, refit = FALSE)), "S2_LRT_KEC13.csv")

#note: significant variance of chamber so must include. also the trial itself have nonvaria so.
# ==========================================================
# KEC13 — final model
# ==========================================================
mod_kec13 <- lmer(value ~ method * Time_days * Treatment_name + (1 | Trial:Chamber),
                  data = kec13.lm)
summary(mod_kec13)
anova(mod_kec13)
VarCorr(mod_kec13)

emm_kec13 <- emmeans(mod_kec13, ~ method | Treatment_name)
emm_kec13
pairs(emm_kec13, adjust = "tukey")
confint(pairs(emm_kec13, adjust = "tukey"))

write.csv(as.data.frame(anova(mod_kec13)), "SuppTable1_anova_KEC13.csv")
write.csv(as.data.frame(confint(pairs(emm_kec13, adjust = "tukey"))), "Table2_KEC13.csv")
write.csv(as.data.frame(pairs(emm_kec13, adjust = "tukey")), "Table2_P_KEC13.csv")


# ==========================================================
# KEC12 — final model
# ==========================================================
mod_kec12 <- lmer(value ~ method * Time_days * Treatment_name + (1 | Trial:Chamber),
                  data = kec12.lm)
summary(mod_kec12)
anova(mod_kec12)
VarCorr(mod_kec12)

emm_kec12 <- emmeans(mod_kec12, ~ method | Treatment_name)
emm_kec12
pairs(emm_kec12, adjust = "tukey")
confint(pairs(emm_kec12, adjust = "tukey"))

write.csv(as.data.frame(anova(mod_kec12)), "SuppTable1_anova_KEC12.csv")
write.csv(as.data.frame(confint(pairs(emm_kec12, adjust = "tukey"))), "Table2_KEC12.csv")
write.csv(as.data.frame(pairs(emm_kec12, adjust = "tukey")), "Table2_P_KEC12.csv")


# ==========================================================
# KEC12 vs KEC13
# ==========================================================
mod_all <- lmer(value ~ assay * Time_days * Treatment_name + (1 | Trial:Chamber),
                data = all.lm)
anova(mod_all)

emm_all <- emmeans(mod_all, ~ assay | Treatment_name)
emm_all
pairs(emm_all, adjust = "tukey")
confint(pairs(emm_all, adjust = "tukey"))

with(combined, t.test(Kec12_qpcr, Kec13_qpcr, paired = TRUE))

write.csv(as.data.frame(anova(mod_all)), "SuppTable_anova_all_assays.csv")
write.csv(as.data.frame(confint(pairs(emm_all, adjust = "tukey"))), "Table_primer_comparison.csv")
write.csv(as.data.frame(pairs(emm_all, adjust = "tukey")), "Table_primer_comparison_P.csv")



# Per-timepoint contrasts — now required, since the 3-way is significant
pairs(emmeans(mod_kec12, ~ method | Time_days | Treatment_name), adjust = "tukey")
pairs(emmeans(mod_kec13, ~ method | Time_days | Treatment_name), adjust = "tukey")

write.csv(as.data.frame(pairs(emmeans(mod_kec12, ~ method | Time_days | Treatment_name),
                              adjust = "tukey")), "TableS3_by_time_KEC12.csv")
write.csv(as.data.frame(pairs(emmeans(mod_kec13, ~ method | Time_days | Treatment_name),
                              adjust = "tukey")), "TableS3_by_time_KEC13.csv")


# Population estimates
pop_est <- combined %>%
  pivot_longer(c(CFU, Kec12_qpcr, Kec13_qpcr), names_to = "Method", values_to = "Value") %>%
  group_by(Treatment_name, Method) %>%
  summarise(n = n(), mean = mean(Value), sd = sd(Value),
            min = min(Value), max = max(Value), .groups = "drop")

pop_est
write.csv(pop_est, "SuppTable2_population_estimates.csv", row.names = FALSE)

# Final plots (merged time points)


stat13 <- ggplot(kec13.lm, aes(x = Treatment_name, y = value, fill = method)) +  
  stat_summary(fun = mean, geom = "bar", position = "dodge", width = 0.5) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2, position = position_dodge(0.5)) +
  scale_fill_manual(name = "Method",
                    values = c(CFU = cbbPalette[1], Kec13_qpcr = cbbPalette[2]),
                    labels = c("CFU", "qPCR")) +
  labs(x = "Treatment", y = expression(Log[10]~"DH44 per g root") )+
  ggtitle("KEC13") +
  theme_classic(base_size = 10) +
  theme(plot.title = element_text(hjust = 0.5)) 
 # geom_pwc(aes(group = method), method = "t_test", label = "p.adj.signif") +
 # geom_pwc(aes(group = Treatment_name), method = "t_test", label = "p.adj.signif",y.position = c(9.5, 7),   # Adjust spacing if needed
       #    bracket.nudge.y = 0.1,
        #   tip.length = 0.01,
        #   position = position_dodge(0.5))
stat13

stat12 <- ggplot(kec12.lm, aes(x = Treatment_name, y = value, fill = method)) +  
  stat_summary(fun = mean, geom = "bar", position = "dodge", width = 0.5) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2, position = position_dodge(0.5)) +
  scale_fill_manual(name = "Method",
                    values = c(CFU = cbbPalette[1], Kec12_qpcr = cbbPalette[2]),
                    labels = c("CFU", "qPCR")) +
  labs(x = "Treatment", y = expression(Log[10]~"DH44 per g root")) +
  ggtitle("KEC12") +
  theme_classic(base_size = 10) +
  theme(plot.title = element_text(hjust = 0.5)) 
  #geom_pwc(aes(group = method), method = "t_test", label = "p.adj.signif") +
  #geom_pwc(aes(group = Treatment_name), method = "t_test", label = "p.adj.signif",y.position = c(9.5, 7),   # Adjust spacing if needed
         #  bracket.nudge.y = 0.1,
          # tip.length = 0.01,
          # position = position_dodge(0.5))
stat12

my_theme <- theme_classic(base_size = 14) +
  theme(
    axis.title     = element_text(size = 14),
    axis.text      = element_text(size = 14),
    legend.title   = element_text(size = 14, face = "bold"),
    legend.text    = element_text(size = 14),
    plot.title     = element_text(size = 14, face = "bold", hjust = 0.5)
  )


stat13     <- stat13     + my_theme
stat12     <- stat12     + my_theme
plot_kec12 <- plot_kec12 + my_theme
plot_kec13 <- plot_kec13 + my_theme

row_bottom <- ggarrange(stat12, stat13,
                     labels = c("C", "D"),
                     common.legend = TRUE, legend = "bottom")

row_top <- ggarrange(plot_kec12, plot_kec13,
                        labels = c("A", "B"),
                        common.legend = TRUE, legend = "top")

Treatment_plot <- ggarrange(row_top, row_bottom, nrow = 2, heights = c(1, 1))
Treatment_plot

ggsave(
  filename = "Validation_Treatment_Level_Plot.png",  # path/filename
  plot     = Treatment_plot,          # the plot object
  width    = 8,                # width in inches
  height   = 9,                # height in inches
  bg = "white",
  dpi      = 600               # dots per inch (300–600 is standard for publication)
)

ggsave(
  filename = "Validation_Treatment_Level_Plot.pdf",  # path/filename
  plot     = Treatment_plot,          # the plot object
  width    = 8,                # width in inches
  height   = 9,                # height in inches
  bg = "white",
  dpi      = 600               # dots per inch (300–600 is standard for publication)
)



summary_all <- combined %>%
  select(Treatment_name, Time_days, CFU, Kec12_qpcr, Kec13_qpcr) %>%
  pivot_longer(cols = c(CFU, Kec12_qpcr, Kec13_qpcr),
               names_to = "Method",
               values_to = "Value") %>%
  group_by(Treatment_name, Time_days, Method) %>%
  summarise(
    mean_log10 = mean(Value, na.rm = TRUE),
    sd_log10   = sd(Value, na.rm = TRUE),
    n          = dplyr::n(),
    .groups = "drop"
  )

summary_all

