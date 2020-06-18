################# HPV Project  #################

library(ggplot2)
library(foreign)
library(tidyverse)
library(sandwich)
library(miceadds)
library(boot)
library(grid)
set.seed(060320)

getwd()
setwd("D:/Google Drive/STOP-HPV_Cost-effectiveness-analysis")


## importing data
hpv_outcome <- read.csv("data/SiteLevel_20May2020.csv")
hpv_cost <- read.dta("data/final_cost_data.dta") 

# All variable names in small letters
colnames(hpv_outcome) <- tolower(colnames(hpv_outcome))

# Eliminating feedback period
hpv_outcome <- subset(hpv_outcome, period == "Communication" | period == "Baseline")

##### Preparing Data ######
# Generate variables for # of missed and # non-missed
hpv_outcome <- hpv_outcome %>% 
          mutate(well_init_missed_n = well_init*well_init_d,  
                 well_init_non_missed_n = well_init_d - well_init_missed_n,
                 well_sub_missed_n = well_sub*well_sub_d,
                 well_sub_non_missed_n = well_sub_d - well_sub_missed_n, 
                 sick_init_missed_n = sick_init*sick_init_d,
                 sick_init_non_missed_n = sick_init_d - sick_init_missed_n,
                 sick_sub_missed_n = sick_sub*sick_sub_d,
                 sick_sub_non_missed_n = sick_sub_d - sick_sub_missed_n)

# Creating a variable ID to merge with 
hpv_outcome <- hpv_outcome %>% mutate(merge_var = site)
hpv_cost <- hpv_cost %>% mutate(merge_var = practice_id)

# Merging, cleaning and adjusting for no cost during baseline
hpv_cea <- merge(hpv_outcome, hpv_cost, by = "merge_var")

hpv_cea <- select(hpv_cea, -c(clinicians, visits, patients, merge_var, practice_id))

hpv_cea <- hpv_cea %>% mutate(cost_physicians_hr = 
                                replace(cost_physicians_hr, period == "Baseline", 0))
hpv_cea <- hpv_cea %>% mutate(cost_np_hr = 
                                replace(cost_np_hr, period == "Baseline", 0))
hpv_cea <- hpv_cea %>% mutate(cost_pa_hr = 
                                replace(cost_pa_hr, period == "Baseline", 0))
hpv_cea <- hpv_cea %>% mutate(cost_aprn_hr = 
                                replace(cost_aprn_hr, period == "Baseline", 0))
hpv_cea <- hpv_cea %>% mutate(cost_rn_hr = 
                                replace(cost_rn_hr, period == "Baseline", 0))
hpv_cea <- hpv_cea %>% mutate(cost_ma_hr = 
                                replace(cost_ma_hr, period == "Baseline", 0))
hpv_cea <- hpv_cea %>% mutate(cost_lpn_hr = 
                                replace(cost_lpn_hr, period == "Baseline", 0))
hpv_cea <- hpv_cea %>% mutate(cost_manager_hr = 
                                replace(cost_manager_hr, period == "Baseline", 0))
hpv_cea <- hpv_cea %>% mutate(cost_oth_hr = 
                                replace(cost_oth_hr, period == "Baseline", 0))
hpv_cea <- hpv_cea %>% mutate(cost_non_dollars_1 = 
                                replace(cost_non_dollars_1, period == "Baseline", 0))
hpv_cea <- hpv_cea %>% mutate(total_cost = 
                                replace(total_cost, period == "Baseline", 0))

write_csv(hpv_cea, path = "data/hpv_cea_06-03-20.csv")

##### Descriptives #####
well_init_descriptive <- hpv_cea %>% group_by(group, period) %>%
      summarize(well_init_visit = mean(well_init_d),
                well_init_mo = mean(well_init_missed_n),
                well_init_mo_rate = mean(well_init))
well_init_descriptive
write_csv(well_init_descriptive, path = "output/well_init_summary.csv")


well_sub_descriptive <- hpv_cea %>% group_by(group, period) %>%
  summarize(well_sub_visit = mean(well_sub_d),
            well_sub_mo = mean(well_sub_missed_n),
            well_sub_mo_rate = mean(well_sub))
well_sub_descriptive
write_csv(well_sub_descriptive, path = "output/well_sub_summary.csv")


sick_init_descriptive <- hpv_cea %>% group_by(group, period) %>%
  summarize(sick_init_visit = mean(sick_init_d),
            sick_init_mo = mean(sick_init_missed_n),
            sick_init_mo_rate = mean(sick_init))
sick_init_descriptive
write_csv(sick_init_descriptive, path = "output/sick_init_summary.csv")


sick_sub_descriptive <- hpv_cea %>% group_by(group, period) %>%
  summarize(sick_sub_visit = mean(sick_sub_d),
            sick_sub_mo = mean(sick_sub_missed_n),
            sick_sub_mo_rate = mean(sick_sub))
sick_sub_descriptive
write_csv(sick_sub_descriptive, path = "output/sick_sub_summary.csv")

cost_descriptive <- hpv_cea %>% group_by(group, period) %>%
  summarize(average_cost <- mean(total_cost))
cost_descriptive


hpv_cea1 <-
a<- (hpv_cea$cost_physicians_hr*82.00 + 
  hpv_cea$cost_np_hr*51.46 + 
  hpv_cea$cost_pa_hr*52.22 + 
  hpv_cea$cost_aprn_hr*51.46 + 
  hpv_cea$cost_rn_hr*34.48 + 
  hpv_cea$cost_ma_hr*16.16 + 
  hpv_cea$cost_lpn_hr*22.23 + 
  hpv_cea$cost_manager_hr*47.95 + 
  hpv_cea$cost_oth_hr*16.76 + 
  hpv_cea$cost_non_dollars_1) / 0.686

##### Assessing Effectiveness   ##########

# Function for the delta method (diff-in-diff of predicted probabilities)
delta_method <-  function(missed, seized, data, i){
  data_sample <- data[i,]
  out_var <- c(missed, seized)
  grouped_outcome <- data_sample[out_var]
  grouped_outcome <- data.matrix(grouped_outcome)
  fit_logit <- glm.cluster(grouped_outcome ~ as.factor(group)*as.factor(period), 
                           family = binomial (link = logit), cluster = "site", 
                           data = data_sample)
  
  # Calculating odds and probabilites for the delta method
  ln_odds_post_treat <- as.numeric(coef(fit_logit)[1]) + as.numeric(coef(fit_logit)[2])*1 + 
    as.numeric(coef(fit_logit)[3])*1 + as.numeric(coef(fit_logit)[4])*1*1
  ln_odds_pre_treat <- as.numeric(coef(fit_logit)[1]) + as.numeric(coef(fit_logit)[2])*1 + 
    as.numeric(coef(fit_logit)[3])*0 + as.numeric(coef(fit_logit)[4])*1*0
  ln_odds_post_control <- as.numeric(coef(fit_logit)[1])+ as.numeric(coef(fit_logit)[2])*0 + 
    as.numeric(coef(fit_logit)[3])*1 + as.numeric(coef(fit_logit)[4])*0*1
  ln_odds_pre_control <- as.numeric(coef(fit_logit)[1])+ as.numeric(coef(fit_logit)[2])*0 + 
    as.numeric(coef(fit_logit)[3])*0 + as.numeric(coef(fit_logit)[4])*0*0
  
  prob_post_treat = exp(ln_odds_post_treat) / (1 + exp(ln_odds_post_treat))
  prob_pre_treat = exp(ln_odds_pre_treat) / (1 + exp(ln_odds_pre_treat))
  prob_post_control = exp(ln_odds_post_control) / (1 + exp(ln_odds_post_control))
  prob_pre_control = exp(ln_odds_pre_control) / (1 + exp(ln_odds_pre_control))
  
  delta <- (prob_post_treat - prob_pre_treat) - (prob_post_control - prob_pre_control)
  delta
  
  return(delta)
  
}

# Estimating ATE 
well_init_ate <-  delta_method(missed = "well_init_missed_n", seized = "well_init_non_missed_n",
                                     data = hpv_cea)

well_sub_ate <-  delta_method(missed = "well_sub_missed_n", seized = "well_sub_non_missed_n", 
                                    data = hpv_cea)

sick_init_ate <- delta_method(missed = "sick_init_missed_n", seized = "sick_init_non_missed_n",
                              data = hpv_cea)

sick_sub_ate <-  delta_method(missed = "sick_sub_missed_n", seized = "sick_sub_non_missed_n", 
                              data = hpv_cea)

table_ate <- matrix(NA, nrow = 2, ncol = 2)
rownames(table_ate) <- c("Well-child Care", "Acute/Chronic")
colnames(table_ate) <- c("Initial HPV Dose", "Subsequent HPV Dose")
table_ate[1,1] <- well_init_ate*100
table_ate[1,2] <- well_sub_ate*100
table_ate[2,1] <- sick_init_ate*100
table_ate[2,2] <- sick_sub_ate*100
table_ate

##### Bootstrap function for delta and costs ########
  
delta_method_boot <- function(missed, seized, visits, data, n_boot){
  
  out_var <- c(missed, seized)
  grouped_outcome <- data[out_var]
  
  wide_data <- reshape(data, 
                       direction='wide', 
                       timevar='period',
                       idvar= c('site', 'group'))
  
  # Bootstrapping
  bootstrap_output <- matrix(NA, ncol = 4, nrow = n_boot)
  colnames(bootstrap_output) <- c("delta", "cost", "averted_MO", "ce_ratio")
  
  for(i in 1:n_boot){
    boot_wide <- wide_data[c(sample(wide_data$site[wide_data$group == 1], sum(wide_data$group == 1), replace = TRUE),
                             sample(wide_data$site[wide_data$group == 2], sum(wide_data$group == 2), replace = TRUE)),]      
    # Wide to long
    boot_sample <- reshape(boot_wide, 
                           direction='long', 
                           varying= c("cost_aprn_hr.Baseline", "cost_lpn_hr.Baseline",
                                      "cost_ma_hr.Baseline", "cost_manager_hr.Baseline",
                                      "cost_non_dollars_1.Baseline", "cost_np_hr.Baseline",
                                      "cost_oth_hr.Baseline", "cost_pa_hr.Baseline",
                                      "cost_physicians_hr.Baseline", "cost_rn_hr.Baseline",
                                      "sick_init.Baseline", "sick_init_d.Baseline",
                                      "sick_init_missed_n.Baseline", "sick_init_non_missed_n.Baseline",
                                      "sick_sub.Baseline", "sick_sub_d.Baseline", 
                                      "sick_sub_missed_n.Baseline",  "sick_sub_non_missed_n.Baseline",
                                      "total_cost.Baseline", "well_init.Baseline",
                                      "well_init_d.Baseline","well_init_missed_n.Baseline", 
                                      "well_init_non_missed_n.Baseline","well_sub.Baseline",
                                      "well_sub_d.Baseline", "well_sub_missed_n.Baseline",
                                      "well_sub_non_missed_n.Baseline",
                                      "cost_aprn_hr.Communication", "cost_lpn_hr.Communication",
                                      "cost_ma_hr.Communication", "cost_manager_hr.Communication",
                                      "cost_non_dollars_1.Communication", "cost_np_hr.Communication",  
                                      "cost_oth_hr.Communication", "cost_pa_hr.Communication",
                                      "cost_physicians_hr.Communication","cost_rn_hr.Communication", 
                                      "sick_init.Communication", "sick_init_d.Communication",  
                                      "sick_init_missed_n.Communication", "sick_init_non_missed_n.Communication",
                                      "sick_sub.Communication", "sick_sub_d.Communication",
                                      "sick_sub_missed_n.Communication", "sick_sub_non_missed_n.Communication",
                                      "total_cost.Communication", "well_init.Communication",
                                      "well_init_d.Communication", "well_init_missed_n.Communication",
                                      "well_init_non_missed_n.Communication", "well_sub.Communication",
                                      "well_sub_d.Communication", "well_sub_missed_n.Communication",
                                      "well_sub_non_missed_n.Communication"), 
                           timevar = "period",
                           times = c("Baseline", "Communication"),
                           v.names = c("cost_aprn_hr", "cost_lpn_hr",
                                       "cost_ma_hr", "cost_manager_hr",
                                       "cost_non_dollars_1", "cost_np_hr", 
                                       "cost_oth_hr", "cost_pa_hr", 
                                       "cost_physicians_hr", "cost_rn_hr", 
                                       "sick_init", "sick_init_d", 
                                       "sick_init_missed_n", "sick_init_non_missed_n",
                                       "sick_sub", "sick_sub_d",
                                       "sick_sub_missed_n", "sick_sub_non_missed_n",
                                       "total_cost", "well_init",
                                       "well_init_d", "well_init_missed_n",
                                       "well_init_non_missed_n", "well_sub",
                                       "well_sub_d", "well_sub_missed_n",
                                       "well_sub_non_missed_n"),
                           idvar = 'site', 
                           new.row.names= c(1:96))
    
    # Applying delta method to obtain effectiveness measure
    boot_grouped_outcome <- boot_sample[out_var]
    boot_grouped_outcome <- data.matrix(boot_grouped_outcome)
    boot_delta <- delta_method(missed, seized, data = boot_sample)
    # Calculating cost
    boot_cost <- (boot_sample$cost_physicians_hr*82.00 + 
                    boot_sample$cost_np_hr*51.46 + 
                    boot_sample$cost_pa_hr*52.22 + 
                    boot_sample$cost_aprn_hr*51.46 + 
                    boot_sample$cost_rn_hr*34.48 + 
                    boot_sample$cost_ma_hr*16.16 + 
                    boot_sample$cost_lpn_hr*22.23 + 
                    boot_sample$cost_manager_hr*47.95 + 
                    boot_sample$cost_oth_hr*16.76) / 0.686 +
                    boot_sample$cost_non_dollars_1
    # CE ratio
    vis_n <- c(visits, "group", "period")
    mod_boot_sample <- boot_sample[which(boot_sample$group == 1 & 
                                           boot_sample$period == "Communication"), ]
    boot_visits <- mod_boot_sample[visits]
    boot_visits <- data.matrix(boot_visits)
    post_test_visits <- mean(boot_visits)
    ce_ratio <- -2*diff(tapply(boot_cost, boot_sample$group, mean))/(boot_delta*post_test_visits)
    # Creating matrix output
    bootstrap_output[i,1] <- boot_delta
    bootstrap_output[i,2] <- -2*diff(tapply(boot_cost, boot_sample$group, mean))
    bootstrap_output[i,3] <- boot_delta*post_test_visits
    bootstrap_output[i,4] <- ce_ratio
  }
  return(bootstrap_output)
}


# Well-child care & Initial Dose
  
  well_init_boot <- delta_method_boot(missed = "well_init_missed_n", seized = "well_init_non_missed_n",
                                      visits = "well_init_d", data = hpv_cea, n_boot = 10000)
  
    # Transforming output in data frame
  well_init_boot <- data.frame(well_init_boot)
  
  # Histogram of Bootstrap estimates
  hist(well_init_boot$delta*100, col ="lightblue",
       main = "Bootstrap estimates of ATE\nWell-child Care & Initial HPV Dose", 
       xlab = "Percentage Point Change in MOs")
  abline(v = 0, col ="red", lwd = 3)
  
  # Summary of estimates
  summary(well_init_boot$cost)
  
  # Obtaining CI bounds
  well_init_ci_lower_lim <- quantile(well_init_boot$delta, probs=0.025)
  well_init_ci_upper_lim <- quantile(well_init_boot$delta, probs=0.975)

  well_init_ci_lower_lim 
  well_init_ci_upper_lim 
  
# Well-child care & Subsequent HPV dose  
  
  well_sub_boot <- delta_method_boot(missed = "well_sub_missed_n", seized = "well_sub_non_missed_n",
                                      visits = "well_sub_d", data = hpv_cea, n_boot = 10000)
  
    # Transforming output in data frame
  well_sub_boot <- data.frame(well_sub_boot)
  
  # Histogram of Bootstrap estimates
  hist(well_sub_boot$delta*100, col ="lightyellow",
       main = "Bootstrap estimates of ATE\nWell-child Care & Subsequent HPV Dose", 
       xlab = "Percentage Point Change in MOs")
  abline(v = 0, col ="red", lwd = 3)
  
  # Summary of estimates
  summary(well_sub_boot$delta)
  
  # Obtaining CI bounds
  well_sub_ci_lower_lim <- quantile(well_sub_boot$delta, probs=0.025)
  well_sub_ci_upper_lim <- quantile(well_sub_boot$delta, probs=0.975)

  
# Acute/Chronic & Initial Dose
  
  sick_init_boot <- delta_method_boot(missed = "sick_init_missed_n", seized = "sick_init_non_missed_n",
                                      visits = "sick_init_d", data = hpv_cea, n_boot = 10000)
  
    # Transforming output in data frame
  sick_init_boot <- data.frame(sick_init_boot)
  
  # Histogram of Bootstrap estimates
  hist(sick_init_boot$delta*100, col ="lightgreen",
       main = "Bootstrap estimates of ATE\nAcute/Chronic Care & Initial HPV Dose", 
       xlab = "Percentage Point Change in MOs")
  abline(v = 0, col ="red", lwd = 3)
  
  # Summary of estimates
  summary(sick_init_boot$delta)
  
  # Obtaining CI bounds
  sick_init_ci_lower_lim <- quantile(sick_init_boot$delta, probs=0.025)
  sick_init_ci_upper_lim <- quantile(sick_init_boot$delta, probs=0.975)


    
# Acute/Chronic & Subsequent HPV Dose
  
  sick_sub_boot <- delta_method_boot(missed = "sick_sub_missed_n", seized = "sick_sub_non_missed_n",
                                      visits = "sick_sub_d", data = hpv_cea, n_boot = 10000)
  
  # Transforming output in data frame
  sick_sub_boot <- data.frame(sick_sub_boot)
  
  # Histogram of Bootstrap estimates
  hist(sick_sub_boot$delta*100, col ="peachpuff3",
       main = "Bootstrap estimates of ATE\nAcute/Chronic Care & Subsequent HPV Dose", 
       xlab = "Percentage Point Change in MOs")
  abline(v = 0, col ="red", lwd = 3)
  
  # Summary of estimates
  summary(sick_sub_boot$delta)
  
  # Obtaining CI bounds
  sick_sub_ci_lower_lim <- quantile(sick_sub_boot$delta, probs=0.025)
  sick_sub_ci_upper_lim <- quantile(sick_sub_boot$delta, probs=0.975)
  
  
  
# Obtaining figures of table 3
  table3 <- data.frame(variable = c("Well-child Care - Initial HPV Dose",
                                     "Well-child Care - Subsequent HPV Dose",
                                     "Acute/Chronic - Initial HPV Dose", 
                                     "Acute/Chronic - Subsequent HPV Dose"), 
                       point_estimates = c(well_init_ate*100,
                                           well_sub_ate*100, 
                                           sick_init_ate*100, 
                                           sick_sub_ate*100),
                       ci_lower_lim = c(well_init_ci_lower_lim*100, 
                                        well_sub_ci_lower_lim*100, 
                                        sick_init_ci_lower_lim*100, 
                                        sick_sub_ci_lower_lim*100),
                       ci_upper_lim = c(well_init_ci_upper_lim*100, 
                                        well_sub_ci_upper_lim*100, 
                                        sick_init_ci_upper_lim*100, 
                                        sick_sub_ci_upper_lim*100) 
                       )                       
  
  write_csv(table3, path = "output/table3.csv")
  
##### CEA  #############
  
  #### Table for MOs and Cost
  
  # Calculating CI bounds for MOs
  MO_well_init_ci_lower_lim <- quantile(well_init_boot$averted_MO, probs=0.025)
  MO_well_init_ci_upper_lim <- quantile(well_init_boot$averted_MO, probs=0.975)
  
  MO_well_sub_ci_lower_lim <- quantile(well_sub_boot$averted_MO, probs=0.025)
  MO_well_sub_ci_upper_lim <- quantile(well_sub_boot$averted_MO, probs=0.975)
  
  MO_sick_init_ci_lower_lim <- quantile(sick_init_boot$averted_MO, probs=0.025)
  MO_sick_init_ci_upper_lim <- quantile(sick_init_boot$averted_MO, probs=0.975)
  
  MO_sick_sub_ci_lower_lim <- quantile(sick_sub_boot$averted_MO, probs=0.025)
  MO_sick_sub_ci_upper_lim <- quantile(sick_sub_boot$averted_MO, probs=0.975)
  
  # Calculating bounds for costs
  cost_well_init_ci_lower_lim <- quantile(well_init_boot$cost, probs=0.025)
  cost_well_init_ci_upper_lim <- quantile(well_init_boot$cost, probs=0.975)
  
  cost_well_sub_ci_lower_lim <- quantile(well_sub_boot$cost, probs=0.025)
  cost_well_sub_ci_upper_lim <- quantile(well_sub_boot$cost, probs=0.975)
  
  cost_sick_init_ci_lower_lim <- quantile(sick_init_boot$cost, probs=0.025)
  cost_sick_init_ci_upper_lim <- quantile(sick_init_boot$cost, probs=0.975)
  
  cost_sick_sub_ci_lower_lim <- quantile(sick_sub_boot$cost, probs=0.025)
  cost_sick_sub_ci_upper_lim <- quantile(sick_sub_boot$cost, probs=0.975)
  
  
  table_cea <- data.frame(outcome = c("Well-child Care - Initial HPV Dose",
                                    "Well-child Care - Subsequent HPV Dose",
                                    "Acute/Chronic - Initial HPV Dose", 
                                    "Acute/Chronic - Subsequent HPV Dose"), 
                       averted_MOs = c(well_init_ate*well_init_descriptive$well_init_visit[2],
                                           well_sub_ate*well_sub_descriptive$well_sub_visit[2], 
                                           sick_init_ate*sick_init_descriptive$sick_init_visit[2], 
                                           sick_sub_ate*sick_sub_descriptive$sick_sub_visit[2]),
                       MO_ci_lower_lim = c(MO_well_init_ci_lower_lim, 
                                           MO_well_sub_ci_lower_lim, 
                                           MO_sick_init_ci_lower_lim, 
                                           MO_sick_sub_ci_lower_lim),
                       MO_ci_upper_lim = c(MO_well_init_ci_upper_lim, 
                                           MO_well_sub_ci_upper_lim, 
                                           MO_sick_init_ci_upper_lim, 
                                           MO_sick_sub_ci_upper_lim),
                       cost = cost_descriptive$`average_cost <- mean(total_cost)`[2], 
                       cost_ci_lower_lim = c(cost_well_init_ci_lower_lim, 
                                           cost_well_sub_ci_lower_lim, 
                                           cost_sick_init_ci_lower_lim, 
                                           cost_sick_sub_ci_lower_lim),
                       cost_ci_upper_lim = c(cost_well_init_ci_upper_lim, 
                                           cost_well_sub_ci_upper_lim, 
                                           cost_sick_init_ci_upper_lim, 
                                           cost_sick_sub_ci_upper_lim)
                       )    
  
  table_cea <- table_cea %>% mutate(ce_ratio = cost/averted_MOs)
  
  write_csv(table_cea, path = "output/table_cea.csv")
  
  ####  CE planes ####
  
  # Well-child care & Initial HPV Dose
  well_init_boot <- well_init_boot %>%
    mutate(ne_quadrant = if_else(cost > 0 & averted_MO > 0, 1, 0)) 
  
  percent_ne_quadrant_mo_well_1st <- sum(well_init_boot$ne_quadrant)/100
  percent_ne_quadrant_mo_well_1st
  
  CEplane_mo_well_1st <-ggplot(well_init_boot, 
                               aes(x = averted_MO, 
                                   y = cost)) +
    xlim(-10, 40) +
    ylim(-1000, 4000) +
    geom_vline(xintercept = 0) + 
    geom_hline(yintercept = 0) + 
    geom_point(size = 0.5, colour = "blue") +
    labs(x = "Absolute Change in MOs", y = "Cost difference") +
    ggtitle("Cost-Effectiveness Plane\n \n Well-child Care & Initial HPV Dose")
  
  CEplane_mo_well_1st
  ggsave(filename = "CEplane_mo_well_1st.png", plot = CEplane_mo_well_1st, path = "output")
  
  # Well-child care & Subsequent HPV Dose
  well_sub_boot <- well_sub_boot %>%
    mutate(ne_quadrant = if_else(cost > 0 & averted_MO > 0, 1, 0)) 
  
  percent_ne_quadrant_mo_well_sub <- sum(well_sub_boot$ne_quadrant)/100
  percent_ne_quadrant_mo_well_sub
  
  CEplane_mo_well_sub <-ggplot(well_sub_boot, 
                               aes(x = averted_MO, 
                                   y = cost)) +
    xlim(-10, 12) +
    ylim(-1000, 4000) +
    geom_vline(xintercept = 0) + 
    geom_hline(yintercept = 0) + 
    geom_point(size = 0.5, colour = "blue") +
    labs(x = "Absolute Change in MOs", y = "Cost difference") +
    ggtitle("Cost-Effectiveness Plan\n \n Well-child Care & Subsequent HPV Dose")
  
  CEplane_mo_well_sub
  ggsave(filename = "CEplane_mo_well_sub.png", plot = CEplane_mo_well_sub, path = "output")
  
  
  
  # Sick care & Initial HPV Dose
  sick_init_boot <- sick_init_boot %>%
    mutate(ne_quadrant = if_else(cost > 0 & averted_MO > 0, 1, 0)) 
  
  percent_ne_quadrant_mo_sick_1st <- sum(sick_init_boot$ne_quadrant)/100
  percent_ne_quadrant_mo_sick_1st
  
  CEplane_mo_sick_1st <-ggplot(sick_init_boot, 
                               aes(x = averted_MO, 
                                   y = cost)) +
    xlim(-5, 5) +
    ylim(-1000, 4000) +
    geom_vline(xintercept = 0) + 
    geom_hline(yintercept = 0) + 
    geom_point(size = 0.5, colour = "blue") +
    labs(x = "Absolute Change in MOs", y = "Cost difference") +
    ggtitle("Cost-Effectiveness Plan\n \n Acute/Chronic Care & Initial HPV Dose")
  
  CEplane_mo_sick_1st
  ggsave(filename = "CEplane_mo_sick_1st.png", plot = CEplane_mo_sick_1st, path = "output")
  
  
  # Sick care & Subsequent HPV Dose
  sick_sub_boot <- sick_sub_boot %>%
    mutate(ne_quadrant = if_else(cost > 0 & averted_MO > 0, 1, 0)) 
  
  percent_ne_quadrant_mo_sick_sub <- sum(sick_sub_boot$ne_quadrant)/100
  percent_ne_quadrant_mo_sick_sub
  
  CEplane_mo_sick_sub <-ggplot(sick_sub_boot, 
                               aes(x = averted_MO, 
                                   y = cost)) +
    xlim(-15, 15) +
    ylim(-1000, 4000) +
    geom_vline(xintercept = 0) + 
    geom_hline(yintercept = 0) + 
    geom_point(size = 0.5, colour = "blue") +
    labs(x = "Absolute Change in MOs", y = "Cost difference") +
    ggtitle("Cost-Effectiveness Plan\n \n Acute/Chronic Care & Subsequent HPV Dose")
  
  CEplane_mo_sick_sub
  ggsave(filename = "CEplane_mo_sick_sub.png", plot = CEplane_mo_sick_sub, path = "output")
  
  
  
  
