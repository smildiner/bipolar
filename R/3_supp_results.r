#==============================================================================#
#   Figures and results for paper2
#==============================================================================#

# .libPaths(c("X:\\My Documents\\R\\win-library\\4.1",.libPaths()))

# library(devtools)
# devtools::install_github("smildiner/mHMMbayes", ref = "develop")

# Preferred parallel library
library(tidyverse)
library(mHMMbayes)
library(viridis)
library(simHMM)

# Load utility functions
source("R/utils.r")
source("R/missingness/missingness_utils.r")

# load data
data <- haven::read_sav("ESM bipolar cleaned data_UU project.sav")
covariates <- readxl::read_xlsx("Patient demographics.xlsx")

#------------------------------------------------------------------------------#

# Elongate data with "missing" night occasions
extended_data <- data %>%
  group_by(patient_id) %>%
  summarise(maxdays = max(dayno)) %>%
  group_by(patient_id) %>%
  summarise(dayno = rep(1:maxdays,each=8), beepno = rep(1:8, maxdays))

train_df <- train_df_night <- left_join(extended_data, data) %>%
  dplyr::select(patient_id, time,
                bs_diary_5, bs_diary_13, bs_diary_22,
                bs_diary_15, bs_diary_9, bs_diary_10,
                bs_diary_7, bs_diary_11, bs_diary_17,
                bs_diary_8, bs_diary_14, bs_diary_16) %>%
  group_by(patient_id) %>%
  mutate(patient_id = cur_group_id(), time = row_number()) %>%
  ungroup() %>%
  arrange(patient_id, time) %>%
  dplyr::select(-time) %>%
  as.matrix()

# Elongate data with "missing" night occasions
train_df <- data %>%
  dplyr::select(patient_id, time,
                bs_diary_5, bs_diary_13, bs_diary_22,
                bs_diary_15, bs_diary_9, bs_diary_10,
                bs_diary_7, bs_diary_11, bs_diary_17,
                bs_diary_8, bs_diary_14, bs_diary_16) %>%
  group_by(patient_id) %>%
  mutate(patient_id = cur_group_id()) %>%
  ungroup() %>%
  arrange(patient_id, time) %>%
  dplyr::select(-time) %>%
  as.matrix()

# Load output
out <- readRDS(paste0("outputs/results/truncated/out_cont_gamma_prior_emiss_prior_m4_12dv_it4000_c",2,"_night.rds"))


#------------------------------------------------------------------------------#
# Check covariates by group:

bipolar_type <- read_csv(file = "bipolar_types.csv")

covariates_df <- inner_join(covariates,
                            bipolar_type %>%
                              rename("patient_ID" = "patient_id") %>%
                              mutate(diagnosis = ifelse(diagnosis == "bipolar type I",0,1))) %>%
  filter(patient_ID %in% unique(data$patient_id)) %>%
  group_by(patient_ID) %>%
  mutate(patient_ID = cur_group_id(),
         `sex (M=1)` = `sex (M=1)`-1) %>%
  ungroup() %>%
  rename("patient_id" = "patient_ID") %>%
  arrange(patient_id)



# Boxplot
covariates_df %>%
  mutate(cluster = case_when(patient_id %in% c(15,12,19,17,1,3) ~ "cluster_1",
                             patient_id %in% c(6,5,7,18,16,13) ~ "cluster_2",
                             patient_id %in% c(20,2,14,11,9,8,4,10) ~ "cluster_3")) %>%
  gather(covariate, value, -patient_id, -cluster) %>%
  filter(covariate %in% c("age","age of onset")) %>%
  ggplot(aes(x = cluster, y = value)) +
  geom_boxplot() +
  facet_wrap(covariate~.) +
  theme_minimal()

# Dotplot
covariates_df %>%
  mutate(cluster = case_when(patient_id %in% c(15,12,19,17,1,3) ~ "cluster_1",
                             patient_id %in% c(6,5,7,18,16,13) ~ "cluster_2",
                             patient_id %in% c(20,2,14,11,9,8,4,10) ~ "cluster_3")) %>%
  gather(covariate, value, -patient_id, -cluster) %>%
  filter(!covariate %in% c("age","age of onset")) %>%
  mutate(value = factor(value)) %>%
  ggplot(aes(x = cluster, y = value, fill = value)) +
  geom_dotplot(binaxis='y', stackdir='center',
               # position=position_dodge(0.8),
               dotsize = 2) +
  facet_wrap(covariate~.) +
  theme_minimal()

# Proportions
covariates_df %>%
  mutate(cluster = case_when(patient_id %in% c(15,12,19,17,1,3) ~ "pattern 1",
                             patient_id %in% c(6,5,7,18,16,13) ~ "pattern 2",
                             patient_id %in% c(20,2,14,11,9,8,4,10) ~ "pattern 3")) %>%
  gather(covariate, value, -patient_id, -cluster) %>%
  filter(!covariate %in% c("age","age of onset","sex (M=1)")) %>%
  mutate(value = factor(value, levels = 0:1, labels = c("No","Yes")),
         covariate = factor(covariate,
                            levels = c("anti-epileptic (1=yes)","antidepressant (1=yes)","antipsychotic (1=yes)",
                                       "diagnosis","litium (1=yes)", "personality_disorder (1=yes)"),
                            labels = c("anti-epileptic use", "antidepressant use","antipsychotic use",
                                       "bipolar type-I/II\ndiagnosis","litium use","personality disorder\ndiagnosis"))) %>%
  ggplot(aes(x = cluster, fill = value)) +
  geom_bar() +
  scale_fill_viridis_d(option = "viridis") +
  facet_wrap(covariate~.) +
  theme_minimal() +
  xlab("")

# Frequency
covariates_df %>%
  mutate(cluster = case_when(patient_id %in% c(15,12,19,17,1,3) ~ "cluster_1",
                             patient_id %in% c(6,5,7,18,16,13) ~ "cluster_2",
                             patient_id %in% c(20,2,14,11,9,8,4,10) ~ "cluster_3")) %>%
  gather(covariate, value, -patient_id, -cluster) %>%
  filter(!covariate %in% c("age","age of onset")) %>%
  group_by(covariate, cluster, value) %>%
  # summarise(mean_cov = mean(value)) %>%
  summarise(N_cov = n()) %>%
  spread(key = cluster, value = N_cov, fill = 0)
  # spread(key = cluster, value = mean_cov)

# Proportion
covariates_df %>%
  mutate(cluster = case_when(patient_id %in% c(15,12,19,17,1,3) ~ "cluster_1",
                             patient_id %in% c(6,5,7,18,16,13) ~ "cluster_2",
                             patient_id %in% c(20,2,14,11,9,8,4,10) ~ "cluster_3")) %>%
  gather(covariate, value, -patient_id, -cluster) %>%
  filter(!covariate %in% c("age","age of onset")) %>%
  group_by(covariate, cluster) %>%
  summarise(mean_cov = mean(value)) %>%
  spread(key = cluster, value = mean_cov)





#------------------------------------------------------------------------------#
# Plot durations:

# Get decoding
states <- vit_mHMM_cont_mar(object = out, s_data = train_df_night, burn_in = 1000) %>%
  as.data.frame() %>%
  gather(subject, state) %>%
  # drop_na() %>%
  mutate(subject = factor(subject, levels = paste0("Subj_",1:20))) %>%
  group_by(subject) %>%
  mutate(subject = cur_group_id()) %>%
  ungroup() %>%
  rename("patient_id" = "subject") %>%
  group_by(patient_id) %>%
  mutate(occasion = row_number()) %>%
  drop_na()

# State proportion by pattern
states %>%
  mutate(pattern = case_when(patient_id %in% c(15,12,19,17,1,3) ~ "pattern_1",
                           patient_id %in% c(6,5,7,18,16,13) ~ "pattern_2",
                           patient_id %in% c(20,2,14,11,9,8,4,10) ~ "pattern_3")) %>%
  group_by(pattern, state) %>%
  summarise(n_state = n()) %>%
  summarise(prop_state = round(n_state/sum(n_state),2)) %>%
  ungroup() %>%
  # mutate(area = prop_state/max(prop_state)*28)
  mutate(diameter = sqrt(prop_state/max(prop_state)*20**2*pi/(pi))+5)

25**2*pi

# # Expected duration given transitions
# get_durations(x = out)
# 
# m <- out$input$m
# n_dep <- out$input$n_dep
# 
# 1/(1-diag(matrix(as.numeric(apply(out$gamma_prob_bar[1001:2000,],2,median)), nrow=4, byrow = TRUE)))
# 1/(1-solve(t(diag(4) - matrix(as.numeric(apply(out$gamma_prob_bar[1001:2000,],2,median)), nrow = 4, byrow = TRUE) + 1), rep(1, 4)))
# 

subj_tpm <- do.call(rbind, lapply(1:length(out$PD_subj),function(s){
  out$PD_subj[[s]] %>%
    as.data.frame() %>%
    dplyr::select("S1toS1":"S4toS4") %>%
    mutate(patient_id = s, iter = row_number()) %>%
    filter(iter > 1000) %>%
    dplyr::select(-iter) %>%
    gather(state, value, -patient_id) %>%
    group_by(patient_id, state) %>%
    summarise(prob = median(value)) %>%
    spread(key = state, value = prob)
} ))

# Group-level: empirical duration on decoded data
group_state_duration <- states %>%
  drop_na() %>%
  group_by(patient_id) %>%
  summarise(decoding = rle(as.character(state))$values,
            duration = rle(as.character(state))$lengths,
            switches = length(rle(as.character(state)))) %>%
  group_by(decoding) %>% 
  summarise(mean_dur = mean(duration),
            sd_dur =  sd(duration),
            median_dur = median(duration),
            cci_lwr = quantile(duration, 0.025), cci_upr = quantile(duration, 0.975),
            min = min(duration),
            max = max(duration)) %>%
  mutate(decoding = factor(decoding, levels = 1:4, labels = c("euthymic","manic","mixed","depressive") ))

# Patient-level: empirical duration on decoded data
states %>%
  drop_na() %>%
  group_by(patient_id) %>%
  # mutate(state = factor(state, levels = 1:4)) %>%
  summarise(decoding = rle(as.character(state))$values,
            duration = rle(as.character(state))$lengths,
            switches = length(rle(as.character(state)))) %>%
  ungroup() %>%
  group_by(patient_id, decoding) %>% 
  summarise(mean_dur = mean(duration),
            sd_dur = sd(duration),
            median_dur = median(duration),
            cci_lwr = quantile(duration, 0.025), cci_upr = quantile(duration, 0.975),
            min = min(duration),
            max = max(duration)) %>%
  group_by(decoding) %>%
  summarise(mean_mean_dur = mean(mean_dur),
            mean_sd_dur = sd(mean_dur),
            mean_median_dur = median(mean_dur),
            mean_cci_lwr = quantile(mean_dur, 0.025), mean_cci_upr = quantile(mean_dur, 0.975),
            mean_min = min(mean_dur),
            mean_max = max(mean_dur))

# Patient-level: empirical duration on decoded data, by patient
patient_state_duration <- states %>%
  drop_na() %>%
  group_by(patient_id) %>%
  # mutate(state = factor(state, levels = 1:4)) %>%
  summarise(decoding = rle(as.character(state))$values,
            duration = rle(as.character(state))$lengths,
            switches = length(rle(as.character(state)))) %>%
  ungroup() %>%
  group_by(patient_id, decoding) %>% 
  summarise(mean_dur = mean(duration),
            sd_dur = sd(duration),
            median_dur = median(duration),
            cci_lwr = quantile(duration, 0.025), cci_upr = quantile(duration, 0.975),
            min = min(duration),
            max = max(duration)) %>%
  mutate(decoding = factor(decoding, levels = 1:4, labels = c("euthymic","manic","mixed","depressive") ))

# Compare durations by group
left_join(patient_state_duration[,1:3],
          covariates_df %>%
            mutate(cluster = case_when(patient_id %in% c(15,12,19,17,1,3) ~ "cluster_1",
                                       patient_id %in% c(6,5,7,18,16,13) ~ "cluster_2",
                                       patient_id %in% c(20,2,14,11,9,8,4,10) ~ "cluster_3")) %>%
            gather(covariate, cov_value, -patient_id, -cluster) %>%
            filter(!covariate %in% c("age","age of onset"))) %>%
  ggplot(aes(x = cluster, y = mean_dur)) +
  geom_boxplot() +
  geom_jitter() +
  facet_wrap(decoding~.) +
  coord_cartesian(ylim = c(0,50)) +
  theme_minimal()

left_join(patient_state_duration[,1:3],
          covariates_df %>%
            mutate(cluster = case_when(patient_id %in% c(15,12,19,17,1,3) ~ "cluster_1",
                                       patient_id %in% c(6,5,7,18,16,13) ~ "cluster_2",
                                       patient_id %in% c(20,2,14,11,9,8,4,10) ~ "cluster_3")) %>%
            gather(covariate, cov_value, -patient_id, -cluster) %>%
            filter(!covariate %in% c("age","age of onset"))) %>%
  group_by(cluster, decoding) %>%
  summarise(median_dur = median(mean_dur)) %>%
  spread(cluster, median_dur)

left_join(patient_state_duration[,1:3],
          covariates_df %>%
            mutate(cluster = case_when(patient_id %in% c(15,12,19,17,1,3) ~ "cluster_1",
                                       patient_id %in% c(6,5,7,18,16,13) ~ "cluster_2",
                                       patient_id %in% c(20,2,14,11,9,8,4,10) ~ "cluster_3")) %>%
            gather(covariate, cov_value, -patient_id, -cluster) %>%
            filter(!covariate %in% c("age","age of onset"))) %>%
  group_by(cluster, decoding) %>%
  summarise(mean_dur = mean(mean_dur)) %>%
  spread(cluster, mean_dur)

left_join(patient_state_duration[,1:3],
          covariates_df %>%
            mutate(cluster = case_when(patient_id %in% c(15,12,19,17,1,3) ~ "cluster_1",
                                       patient_id %in% c(6,5,7,18,16,13) ~ "cluster_2",
                                       patient_id %in% c(20,2,14,11,9,8,4,10) ~ "cluster_3")) %>%
            gather(covariate, cov_value, -patient_id, -cluster) %>%
            filter(!covariate %in% c("age","age of onset"))) %>%
  group_by(cluster, decoding) %>%
  summarise(sd_dur = sd(mean_dur)) %>%
  spread(cluster, sd_dur)



# Switches by cluster
# Number of transitions
# Switches on non-missing data
na_idx <- apply(is.na(train_df),1,any)
switches_table_omit <- drop_na(states)[!na_idx,] %>%
  group_by(patient_id) %>%
  summarise(n_occ = length(state),
            n_switch = length(rle(as.character(state))$values),
            rel_switch = length(rle(as.character(state))$values)/length(state))

left_join(switches_table_omit,
          covariates_df %>%
            mutate(cluster = case_when(patient_id %in% c(15,12,19,17,1,3) ~ "cluster_1",
                                       patient_id %in% c(6,5,7,18,16,13) ~ "cluster_2",
                                       patient_id %in% c(20,2,14,11,9,8,4,10) ~ "cluster_3")) %>%
            gather(covariate, cov_value, -patient_id, -cluster) %>%
            filter(!covariate %in% c("age","age of onset"))) %>%
  ggplot(aes(x = cluster, y = rel_switch)) +
  geom_boxplot() +
  # geom_dotplot(binaxis = "y") +
  geom_jitter() +
  theme_minimal()

left_join(switches_table_omit,
          covariates_df %>%
            mutate(cluster = case_when(patient_id %in% c(15,12,19,17,1,3) ~ "cluster_1",
                                       patient_id %in% c(6,5,7,18,16,13) ~ "cluster_2",
                                       patient_id %in% c(20,2,14,11,9,8,4,10) ~ "cluster_3")) %>%
            gather(covariate, cov_value, -patient_id, -cluster) %>%
            filter(!covariate %in% c("age","age of onset"))) %>%
  group_by(cluster) %>%
  summarise(mean_rel_switch = mean(rel_switch)) %>%
  spread(cluster, mean_rel_switch)

left_join(switches_table_omit,
          covariates_df %>%
            mutate(cluster = case_when(patient_id %in% c(15,12,19,17,1,3) ~ "cluster_1",
                                       patient_id %in% c(6,5,7,18,16,13) ~ "cluster_2",
                                       patient_id %in% c(20,2,14,11,9,8,4,10) ~ "cluster_3")) %>%
            gather(covariate, cov_value, -patient_id, -cluster) %>%
            filter(!covariate %in% c("age","age of onset"))) %>%
  group_by(cluster) %>%
  summarise(sd_rel_switch = sd(rel_switch)) %>%
  spread(cluster, sd_rel_switch)


# Geometric state durations:
geometric_patient_durations <- as.data.frame(cbind(1:20,t(apply(subj_tpm[,-1],1,function(s) 1/(1-diag(matrix(s, byrow = TRUE, nrow = 4))) ))))
names(geometric_patient_durations) <- c("patient_id","neutral","elevated","mixed","lowered")
geometric_patient_durations <- geometric_patient_durations %>%
  gather(decoding, median_dur, -patient_id)

# Geometric duration
geometric_group_state_duration <- get_durations(out, burnin = 1000) %>%
  rename("decoding"="parameter", "median_dur" = "median_val") %>%
  mutate(decoding = factor(decoding, levels = paste0("S",1:4,"toS",1:4), labels = c("neutral", "elevated","mixed","lowered")))

# Plot empirical duration
p <- ggplot(data = patient_state_duration %>%
              mutate(pattern = case_when(patient_id %in% c(15,12,19,17,1,3) ~ "pattern 1",
                                                                    patient_id %in% c(6,5,7,18,16,13) ~ "pattern 2",
                                                                    patient_id %in% c(20,2,14,11,9,8,4,10) ~ "pattern 3"))) +
  geom_line(aes(x = decoding, y = median_dur*3, group = patient_id), linetype = "dashed", alpha = 0.2) +
  geom_jitter(aes(x = decoding, y = median_dur*3, colour = decoding), width = 0.05, alpha = 0.75) +
  # geom_pointrange(data = group_state_duration,
  #                 aes(x = decoding, y = mean_dur*3,
  #                     ymin = cci_lwr*3, ymax = cci_upr*3)) +
  # geom_point(data = group_state_duration,
  #            aes(x = decoding, y = median_dur*3),
  #            # shape = 1,
  #            size = 3) +
  # geom_point(data = geometric_group_state_duration,
  #            aes(x = decoding, y = median_dur*3),
  #            # shape = 1,
  #            size = 3) +
  # geom_pointrange(data = group_state_duration,
  #                 aes(x = decoding, y = mean_dur*3,
  #                     ymin = (mean_dur-sd_dur)*3, ymax = (mean_dur+sd_dur)*3)) +
  scale_color_viridis(discrete = TRUE,
                      alpha=0.6,
                      direction = -1,
                      option = "plasma") +
  geom_hline(yintercept = 0) +
  facet_grid(.~pattern) +
  theme_minimal() +
  # scale_y_log10() +
  # ylab(label = "Empirical state duration (hours)") +
  ylab(label = "State duratio (hours)") +
  xlab(label = "Mood state") +
  # ggtitle("Empirical state duration") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(legend.position = "none") +
  theme(strip.text.y = element_text(angle = 0)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  coord_cartesian(ylim = c(0,45))

p

#------------------------------------------------------------------------------#
# Plot normalized off-diagonal transitions:

subj_tpm[,paste0("S",1:4,"toS",1:4)] <- 0
subj_tpm[,-1] <- apply(subj_tpm[,-1],1,function(s) matrix(s, nrow = 4, byrow = TRUE) %>%
        as.data.frame() %>%
        apply(.,1,function(r) r/sum(r))) %>%
  t()

p <- subj_tpm %>%
  gather(gamma, value, -patient_id) %>%
  filter(!gamma %in% paste0("S",1:4,"toS",1:4)) %>%
  separate(gamma, c("from","to"),"to") %>%
  mutate(from = factor(from, labels = paste0("from\n",c("neutral","elevated","mixed","lowered"))),
         to = factor(to, labels = paste0("to ",c("neutral","elevated","mixed","lowered"))),
         patient_id = factor(patient_id, levels = c(6,5,7,18,16,13,
                                                    15,12,19,17,1,3,
                                                    20,2,14,11,9,8,4,10)),
         pattern = case_when(patient_id %in% c(6,5,7,18,16,13) ~ "Neutral/Lowered",
                             patient_id %in% c(15,12,19,17,1,3) ~ "Neutral/Elevated",
                             patient_id %in% c(20,2,14,11,9,8,4,10) ~ "Mixed/Elevated/Lowered"),
         pattern = factor(pattern, levels = c("Neutral/Lowered", "Neutral/Elevated", "Mixed/Elevated/Lowered"))) %>%
  ggplot(aes(x=to, y=forcats::fct_rev(patient_id), fill = value)) +
  geom_tile() +
  geom_text(aes(label = format(round(value, digits = 2), nsmall = 1), color = value^3)) +
  scale_color_gradient(low = "white", high = "black", guide = NULL) +
  scale_fill_viridis_c(option = "rocket") +
  # facet_wrap(from~., ncol = 4, scales = "free") +
  facet_grid(pattern~from, scales = "free",space = "free_y") +
  theme_minimal() +
  # xlab("To state") +
  xlab(element_blank()) +
  ylab("Patient ID") +
  labs(fill = "Relative\nprob.") +
  theme(axis.text.x = element_text(angle = 45, hjust=1)) +
  theme(legend.key.height= unit(2.5, 'cm'),
        legend.key.width= unit(0.3, 'cm'))

ggsave(p, filename = "outputs/figures/off_diag_heatmap.pdf", dpi = 600)
ggsave(p, filename = "outputs/figures/off_diag_heatmap.jpeg", dpi = 600)
ggsave(p, filename = "outputs/figures/off_diag_heatmap.png", dpi = 600)


# Variability off diagonals:

# Mean
subj_tpm %>%
  gather(gamma, value, -patient_id) %>%
  filter(!gamma %in% paste0("S",1:4,"toS",1:4)) %>%
  mutate(patient_id = factor(patient_id, levels = c(15,12,19,17,1,3,
                                                    6,5,7,18,16,13,
                                                    20,2,14,11,9,8,4,10)),
         pattern = case_when(patient_id %in% c(15,12,19,17,1,3) ~ "pattern 1",
                             patient_id %in% c(6,5,7,18,16,13) ~ "pattern 2",
                             patient_id %in% c(20,2,14,11,9,8,4,10) ~ "pattern 3")) %>%
  group_by(pattern,gamma) %>%
  summarise(value = mean(value)) %>%
  spread(pattern, value)

# SD
subj_tpm %>%
  gather(gamma, value, -patient_id) %>%
  filter(!gamma %in% paste0("S",1:4,"toS",1:4)) %>%
  mutate(patient_id = factor(patient_id, levels = c(15,12,19,17,1,3,
                                                    6,5,7,18,16,13,
                                                    20,2,14,11,9,8,4,10)),
         pattern = case_when(patient_id %in% c(15,12,19,17,1,3) ~ "pattern 1",
                             patient_id %in% c(6,5,7,18,16,13) ~ "pattern 2",
                             patient_id %in% c(20,2,14,11,9,8,4,10) ~ "pattern 3")) %>%
  group_by(pattern,gamma) %>%
  summarise(value = sd(value)) %>%
  spread(pattern, value)



#------------------------------------------------------------------------------#
# Plot missing data by state:

# Train data ommitting the night
train_df_omit <- data %>%
  dplyr::select(patient_id, time,
                bs_diary_5, bs_diary_13, bs_diary_22,
                bs_diary_15, bs_diary_9, bs_diary_10,
                bs_diary_7, bs_diary_11, bs_diary_17,
                bs_diary_8, bs_diary_14, bs_diary_16) %>%
  group_by(patient_id) %>%
  mutate(patient_id = cur_group_id()) %>%
  ungroup() %>%
  arrange(patient_id, time) %>%
  dplyr::select(-time) %>%
  as.matrix()

# Get decoding
states <- vit_mHMM_cont_mar(object = out, s_data = train_df_omit, burn_in = 1000) %>%
  as.data.frame() %>%
  gather(subject, state) %>%
  # drop_na() %>%
  mutate(subject = factor(subject, levels = paste0("Subj_",1:20))) %>%
  group_by(subject) %>%
  mutate(subject = cur_group_id()) %>%
  ungroup() %>%
  rename("patient_id" = "subject") %>%
  group_by(patient_id) %>%
  mutate(occasion = row_number())

data_labelled <- left_join(as.data.frame(train_df_omit) %>% group_by(patient_id) %>%
                             mutate(occasion = row_number()), states) %>%
  group_by(patient_id) %>%
  mutate(occasion = row_number()) %>%
  as.data.frame()

# Bonus: states by missingness
decode_mis <- data_labelled %>%
  mutate(miss = is.na(bs_diary_5)) %>%
  group_by(miss, state) %>%
  summarise(n_state = n()) %>%
  summarise(prop_state = n_state/sum(n_state)) %>%
  ungroup() %>%
  mutate(state = rep(c("euthymic","manic","mixed","depressive"),2),
         state = factor(state, levels = c("euthymic","manic","mixed","depressive"))) %>%
  ungroup() %>%
  spread(miss, prop_state) %>%
  rename("prop_miss" = `FALSE`,
         "prop_obs" = `TRUE`)

# Calculate missing by state and plot it
group_miss <- data_labelled %>%
  mutate(miss = is.na(bs_diary_5),
         patient_id = "Overall") %>%
  group_by(patient_id, state) %>%
  summarise(prop_miss = mean(miss))

patient_miss <- data_labelled %>%
  mutate(miss = is.na(bs_diary_5),
         patient_id = as.character(patient_id)) %>%
  group_by(patient_id, state) %>%
  summarise(prop_miss = mean(miss))

p <- bind_rows(group_miss, patient_miss) %>%
  mutate(state = factor(state, levels = 1:4,
                        labels = c("neutral", "elevated", "mixed", "lowered")),
         patient_id = factor(patient_id, levels = c("Overall", as.character(1:20)))) %>%
  ggplot(aes(x=state, y=patient_id, fill = prop_miss)) +
  geom_tile() +
  geom_text(aes(label = format(round(prop_miss, digits = 2), nsmall = 2), color = prop_miss^10)) +
  scale_color_gradient(low = "white", high = "black", guide = NULL) +
  scale_fill_viridis_c(option = "rocket") +
  theme_minimal() +
  xlab("Mood state") +
  ylab("Patient ID") +
  labs(fill = "Prop.\nmissing") +
  theme(legend.key.height= unit(2.8, 'cm'),
        legend.key.width= unit(0.3, 'cm'))

ggsave(p, filename = "outputs/figures/prop_missing.pdf", dpi = 600)
ggsave(p, filename = "outputs/figures/prop_missing.jpeg", dpi = 600)


#------------------------------------------------------------------------------#
# Plot linked emissions:

# Figure R.1: group-level emissions

# Extract emission means
emiss_mu_bar <-  do.call(rbind, lapply(1:length(out$emiss_mu_bar), function(q){
  out$emiss_mu_bar[[q]] %>%
    as.data.frame() %>%
    mutate(iter = row_number(),
           dep = names(out$emiss_mu_bar)[q]) %>%
    gather(mu, value, -iter, -dep)
}))

# Add labels for the constructs

dep_vars <- c("bs_diary_5", "bs_diary_13", "bs_diary_22",
              "bs_diary_15", "bs_diary_9", "bs_diary_10")
man_vars <- c("bs_diary_7", "bs_diary_11", "bs_diary_17",
              "bs_diary_8", "bs_diary_14", "bs_diary_16")

# Plot

# Patient-specific values instead of group-level MAP
subj_emiss <- do.call(rbind, lapply(1:length(out$PD_subj),function(s){
  out$PD_subj[[s]] %>%
    as.data.frame() %>%
    dplyr::select("dep1_mu_S1":"dep12_mu_S4") %>%
    mutate(patient_id = s, iter = row_number()) %>%
    filter(iter > 1000) %>%
    dplyr::select(-iter) %>%
    gather(dep, value, -patient_id) %>%
    group_by(patient_id, dep) %>%
    summarise(value = median(value)) %>%
    separate(dep, into = c("dep","mu"), sep = "_mu_") %>%
    mutate(dep = factor(dep, levels = paste0("dep",1:12), labels = c("bs_diary_5", "bs_diary_13", "bs_diary_22",
                                                                     "bs_diary_15", "bs_diary_9", "bs_diary_10",
                                                                     "bs_diary_7", "bs_diary_11", "bs_diary_17",
                                                                     "bs_diary_8", "bs_diary_14", "bs_diary_16")))
  
} ))

# Boxplot
safe_colorblind_palette <- c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499", 
                             "#44AA99", "#999933", "#882255", "#661100", "#6699CC", "#888888")

# Subject and Group level parameters
subj_emiss_nice <- subj_emiss %>%
  mutate(construct = case_when(dep %in% dep_vars ~ "Depression items",
                               dep %in% man_vars ~ "Mania items"),
         construct = factor(construct, levels = c("Mania items", "Depression items")),
         dep = factor(dep,
                      levels = c("bs_diary_5", "bs_diary_13", "bs_diary_22",
                                 "bs_diary_15", "bs_diary_9", "bs_diary_10",
                                 "bs_diary_7", "bs_diary_11", "bs_diary_17",
                                 "bs_diary_8", "bs_diary_14", "bs_diary_16"),
                      labels = c("down","dread rest of day","worry",
                                 "inadequate", "tired", "content",
                                 "agitated", "irritated", "switch and focus",
                                 "extremely well", "full of ideas", "thoughts are racing")),
         mu = factor(mu, levels = paste0("S",1:4), labels = c("neutral","elevated","mixed","lowered"))) %>%
  rename("state" = "mu")

group_emiss_nice <- emiss_mu_bar %>%
  filter(iter > 1000) %>%
  group_by(dep, mu) %>%
  summarise(group_mean = mean(value),
            group_median = median(value),
            group_sd = sd(value),
            group_CCIlwr = quantile(value, 0.025),
            group_CCIupr = quantile(value, 0.975)) %>%
  mutate(construct = case_when(dep %in% dep_vars ~ "Depression items",
                               dep %in% man_vars ~ "Mania items"),
         construct = factor(construct, levels = c("Mania items", "Depression items")),
         dep = factor(dep,
                      levels = c("bs_diary_5", "bs_diary_13", "bs_diary_22",
                                 "bs_diary_15", "bs_diary_9", "bs_diary_10",
                                 "bs_diary_7", "bs_diary_11", "bs_diary_17",
                                 "bs_diary_8", "bs_diary_14", "bs_diary_16"),
                      labels = c("down","dread rest of day","worry",
                                 "inadequate", "tired", "content",
                                 "agitated", "irritated", "switch and focus",
                                 "extremely well", "full of ideas", "thoughts are racing")),
         mu = factor(mu, levels = paste0("mu_",1:4), labels = c("neutral","elevated","mixed","lowered"))) %>%
  rename("state" = "mu")

p <- ggplot(data = subj_emiss_nice) +
  geom_line(aes(x = state, y = value, group = patient_id), alpha = 0.3) +
  geom_jitter(aes(x = state, y = value, colour = state), width = 0.1, alpha = 0.6) +
  geom_pointrange(data = group_emiss_nice,
                  aes(x = state, y = group_mean,
                      ymin = group_CCIlwr, ymax = group_CCIupr)) +
  scale_color_viridis_d(option = "plasma", direction = -1) +
  geom_hline(yintercept = 0) +
  facet_wrap(construct+dep~., scales = "fixed", nrow = 2) +
  theme_minimal() +
  ylab(label = "EMA score") +
  xlab(label = "") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(legend.position = "none") +
  # ggtitle(label = "Composition of momentary mood states: group-level and\npatient-specific item emission scores by state") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(strip.text.y = element_text(angle = 0))  

p

ggsave(p, filename = "outputs/figures/linked_emiss.pdf", dpi = 600)
ggsave(p, filename = "outputs/figures/linked_emiss.jpeg", dpi = 600)


#------------------------------------------------------------------------------#
# Calculate number of transitions to state for diagrams:

# Switches by cluster
# Number of transitions
# Switches
switches_table <- states %>%
  group_by(patient_id) %>%
  summarise(n_occ = length(state),
            state = factor(rle(as.character(state))$values,
                           levels = 1:4, labels = c("neutral","elevated",'mixed',"lowered"))) %>%
  group_by(patient_id, state) %>%
  summarise(n_switch_to = n()) %>%
  spread(state, n_switch_to, fill = 0) %>%
  gather(state, value, -patient_id)

switches_table

# Relative numbers
switches_table %>%
  mutate(pattern = case_when(patient_id %in% c(15,12,19,17,1,3) ~ "pattern_1",
                             patient_id %in% c(6,5,7,18,16,13) ~ "pattern_2",
                             patient_id %in% c(20,2,14,11,9,8,4,10) ~ "pattern_3")) %>%
  group_by(pattern, state) %>%
  summarise(mean_switches_to_state = mean(value)) %>%
  mutate(state = factor(state, levels = c("neutral","elevated","mixed","lowered"))) %>%
  group_by(pattern) %>%
  mutate(mean_switches_to_state = mean_switches_to_state/sum(mean_switches_to_state)) %>%
  spread(pattern,mean_switches_to_state)

# Mean transitions
switches_table %>%
  mutate(pattern = case_when(patient_id %in% c(15,12,19,17,1,3) ~ "pattern_1",
                             patient_id %in% c(6,5,7,18,16,13) ~ "pattern_2",
                             patient_id %in% c(20,2,14,11,9,8,4,10) ~ "pattern_3")) %>%
  group_by(pattern, state) %>%
  summarise(mean_switches_to_state = mean(value)) %>%
  mutate(state = factor(state, levels = c("neutral","elevated","mixed","lowered"))) %>%
  spread(pattern, mean_switches_to_state)

# SD transitions
switches_table %>%
  mutate(pattern = case_when(patient_id %in% c(15,12,19,17,1,3) ~ "pattern_1",
                             patient_id %in% c(6,5,7,18,16,13) ~ "pattern_2",
                             patient_id %in% c(20,2,14,11,9,8,4,10) ~ "pattern_3")) %>%
  group_by(pattern, state) %>%
  summarise(sd_switches_to_state = sd(value)) %>%
  mutate(state = factor(state, levels = c("neutral","elevated","mixed","lowered"))) %>%
  spread(pattern, sd_switches_to_state)

# Plots
library(diagram)

subj_tpm_off <- subj_tpm
subj_tpm_off[,paste0("S",1:4,"toS",1:4)] <- 0
subj_tpm_off[,-1] <- apply(subj_tpm_off[,-1],1,function(s) matrix(s, nrow = 4, byrow = TRUE) %>%
                         as.data.frame() %>%
                         apply(.,1,function(r) r/sum(r))) %>%
  t()

# Variability off diagonals:

# Mean
subj_tpm_off <- subj_tpm_off %>%
  gather(gamma, value, -patient_id) %>%
  # filter(!gamma %in% paste0("S",1:4,"toS",1:4)) %>%
  mutate(patient_id = factor(patient_id, levels = c(15,12,19,17,1,3,
                                                    6,5,7,18,16,13,
                                                    20,2,14,11,9,8,4,10)),
         pattern = case_when(patient_id %in% c(15,12,19,17,1,3) ~ "pattern_1",
                             patient_id %in% c(6,5,7,18,16,13) ~ "pattern_2",
                             patient_id %in% c(20,2,14,11,9,8,4,10) ~ "pattern_3")) %>%
  group_by(pattern,gamma) %>%
  summarise(value = mean(value)) %>%
  spread(pattern, value)

gamma_1 <- t(matrix( round(pull(subj_tpm_off[,2]),2), nrow = 4, ncol = 4, byrow = TRUE,
                 dimnames = list(c("neutral","elevated","mixed","lowered"),c("neutral","elevated","mixed","lowered"))))
gamma_2 <- t(matrix( round(pull(subj_tpm_off[,3]),2), nrow = 4, ncol = 4, byrow = TRUE,
                     dimnames = list(c("neutral","elevated","mixed","lowered"),c("neutral","elevated","mixed","lowered"))))
gamma_3 <- t(matrix( round(pull(subj_tpm_off[,4]),2), nrow = 4, ncol = 4, byrow = TRUE,
                     dimnames = list(c("neutral","elevated","mixed","lowered"),c("neutral","elevated","mixed","lowered"))))


p1 <- plotmat(gamma_1, relsize = 1.1, box.col=c("yellow","red","purple","blue"), pos = c(2,2),
        arr.lwd = gamma_1*10, arr.pos = 0.65,
        arr.width = 0.2, shadow.size = 0)
p2 <- plotmat(gamma_2, relsize = 1.1, box.col=c("yellow","red","purple","blue"), pos = c(2,2),
        arr.lwd = gamma_2*10, arr.pos = 0.65,
        arr.width = 0.2, shadow.size = 0)
p3 <- plotmat(gamma_3, relsize = 1.1, box.col=c("yellow","red","purple","blue"), pos = c(2,2),
        arr.lwd = gamma_3*10, arr.pos = 0.65,
        arr.width = 0.2, shadow.size = 0)



