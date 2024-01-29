#==============================================================================#
#   Figures and results for paper2
#==============================================================================#

# library(devtools)
# devtools::install_github("smildiner/simHMM", ref = "dev")
# devtools::install_github("smildiner/mHMMbayes", ref = "develop")

# Preferred parallel library
library(tidyverse)
library(mHMMbayes)
library(viridis)
library(simHMM)
library(cowplot)
library(lme4)
library(lmerTest)
library(cowplot)

# Load utility functions
source("R/utils.R")
source("R/utils_missingness.r")

# load data
data <- haven::read_sav("data/ESM bipolar cleaned data_UU project.sav")

#------------------------------------------------------------------------------#

# Put data in right format and elongate data with "missing" night occasions
extended_data <- data %>%
  group_by(patient_id) %>%
  summarise(maxdays = max(dayno)) %>%
  group_by(patient_id) %>%
  summarise(dayno = rep(1:maxdays,each=8), beepno = rep(1:8, maxdays))

train_df <- left_join(extended_data, data) %>%
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

#------------------------------------------------------------------------------#

# Load output
out <- readRDS(paste0("results/out_cont_gamma_prior_emiss_prior_m4_12dv_it4000_c",2,"_night.rds"))


#==============================================================================#
# Main results (figures and tables)
#==============================================================================#

#------------------------------------------------------------------------------#
# Figure: Group- and patient-level emission means, linked by patient:

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

# Subject and Group level parameters
subj_emiss_nice <- subj_emiss %>%
  mutate(construct = case_when(dep %in% dep_vars ~ "Depression items",
                               dep %in% man_vars ~ "Mania items"),
         construct = factor(construct, levels = c("Depression items","Mania items")),
         dep = factor(dep,
                      levels = c("bs_diary_5", "bs_diary_13", "bs_diary_22",
                                 "bs_diary_15", "bs_diary_9", "bs_diary_10",
                                 "bs_diary_7", "bs_diary_11", "bs_diary_17",
                                 "bs_diary_8", "bs_diary_14", "bs_diary_16"),
                      labels = c("down","dread rest of day","worry",
                                 "inadequate", "tired", "content",
                                 "agitated", "irritated", "switch and focus",
                                 "extremely well", "full of ideas", "racing thoughts")),
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
         construct = factor(construct, levels = c("Depression items","Mania items")),
         dep = factor(dep,
                      levels = c("bs_diary_5", "bs_diary_13", "bs_diary_22",
                                 "bs_diary_15", "bs_diary_9", "bs_diary_10",
                                 "bs_diary_7", "bs_diary_11", "bs_diary_17",
                                 "bs_diary_8", "bs_diary_14", "bs_diary_16"),
                      labels = c("down","dread rest of day","worry",
                                 "inadequate", "tired", "content",
                                 "agitated", "irritated", "switch and focus",
                                 "extremely well", "full of ideas", "racing thoughts")),
         mu = factor(mu, levels = paste0("mu_",1:4), labels = c("neutral","elevated","mixed","lowered"))) %>%
  rename("state" = "mu")

p <- ggplot(data = subj_emiss_nice) +
  geom_line(aes(x = state, y = value, group = patient_id), alpha = 0.2) +
  geom_jitter(aes(x = state, y = value, colour = state), width = 0.15, alpha = 0.5) +
  geom_pointrange(data = group_emiss_nice,
                  aes(x = state, y = group_mean,
                      ymin = group_CCIlwr, ymax = group_CCIupr), size = 0.25) +
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

ggsave(p + theme(text = element_text(size = 14)), filename = "outputs/figures/linked_emiss_a.pdf", dpi = 600)
# ggsave(p, filename = "outputs/figures/linked_emiss_a.jpeg", dpi = 600)
# ggsave(p, filename = "outputs/figures/linked_emiss_a.png", dpi = 600)


# Reordering variables:
# Subject and Group level parameters
subj_emiss_nice <- subj_emiss %>%
  mutate(construct = case_when(dep %in% dep_vars ~ "Depression items",
                               dep %in% man_vars ~ "Mania items"),
         construct = factor(construct, levels = c("Depression items","Mania items")),
         dep = factor(dep,
                      levels = c("bs_diary_5", "bs_diary_13", "bs_diary_22",
                                 "bs_diary_15", "bs_diary_9", "bs_diary_10",
                                 "bs_diary_17",
                                 "bs_diary_8", "bs_diary_14", "bs_diary_16",
                                 "bs_diary_7", "bs_diary_11"),
                      labels = c("down","dread rest of day","worry",
                                 "inadequate", "tired", "content",
                                 "switch and focus",
                                 "extremely well", "full of ideas", "racing thoughts",
                                 "agitated", "irritated")),
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
         construct = factor(construct, levels = c("Depression items","Mania items")),
         dep = factor(dep,
                      levels = c("bs_diary_5", "bs_diary_13", "bs_diary_22",
                                 "bs_diary_15", "bs_diary_9", "bs_diary_10",
                                 "bs_diary_17",
                                 "bs_diary_8", "bs_diary_14", "bs_diary_16",
                                 "bs_diary_7", "bs_diary_11"),
                      labels = c("down","dread rest of day","worry",
                                 "inadequate", "tired", "content",
                                 "switch and focus",
                                 "extremely well", "full of ideas", "racing thoughts",
                                 "agitated", "irritated")),
         mu = factor(mu, levels = paste0("mu_",1:4), labels = c("neutral","elevated","mixed","lowered"))) %>%
  rename("state" = "mu")

p <- ggplot(data = subj_emiss_nice) +
  geom_line(aes(x = state, y = value, group = patient_id), alpha = 0.2) +
  geom_jitter(aes(x = state, y = value, colour = state), width = 0.15, alpha = 0.5) +
  geom_pointrange(data = group_emiss_nice,
                  aes(x = state, y = group_mean,
                      ymin = group_CCIlwr, ymax = group_CCIupr), size = 0.25) +
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

ggsave(p + theme(text = element_text(size = 14)), filename = "outputs/figures/linked_emiss_b.pdf", dpi = 600)
# ggsave(p, filename = "outputs/figures/linked_emiss_b.jpeg", dpi = 600)
# ggsave(p + theme(text = element_text(size = 12)), filename = "outputs/figures/linked_emiss_b.png", dpi = 600)


#------------------------------------------------------------------------------#
# Figure R.4: validated scales (decoding)

# Load weekly data
bs_asrm <- haven::read_sav("bs_asrm_cleaned_UU.sav") %>%
  arrange(patient_id)
bs_qids <- haven::read_sav("bs_qids_cleaned_UU.sav") %>%
  arrange(patient_id)

# Put together
bs_scales <- full_join(bs_asrm, bs_qids) %>%
  dplyr::select(patient_id, bs_asrm_open_from, bs_asrm_tot, bs_qids_tot) %>%
  rename("open_from" = "bs_asrm_open_from") %>%
  group_by(patient_id) %>%
  mutate(from_time_day = format(as.POSIXct(lag(open_from, 1, default = min(open_from)-(3600*24*7))),format = "%Y-%m-%d"),
         to_time_day = format(as.POSIXct(open_from-1),format = "%Y-%m-%d"),
         patient_id = factor(patient_id)) %>%
  ungroup()

# Get decoding
states <- vit_mHMM_cont_mar(object = out, s_data = train_df, burn_in = 1000) %>%
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

data_labelled <- left_join(as.data.frame(train_df) %>% group_by(patient_id) %>%
                             mutate(occasion = row_number()), states) %>%
  group_by(patient_id) %>%
  mutate(occasion = row_number()) %>%
  as.data.frame()

# Add week info
data_labelled <- bind_cols(data_labelled,
                           left_join(extended_data, data) %>%
                             group_by(patient_id) %>%
                             mutate(patient_id = cur_group_id()) %>%
                             ungroup() %>%
                             arrange(patient_id) %>%
                             dplyr::select(open_from, bs_diary_open_from, bs_diary_date)) %>%
  fill(open_from)

for(r in 1:length(data_labelled$bs_diary_open_from)){
  if(is.na(data_labelled$bs_diary_open_from[r])){
    data_labelled$bs_diary_open_from[r] <- data_labelled$bs_diary_open_from[r-1] + 3*60*60
  }
}


# Validated scales
state_data <- data_labelled %>%
  select(patient_id, state, bs_diary_open_from, occasion) %>%
  rename("open_from" = "bs_diary_open_from") %>%
  group_by(patient_id) %>%
  mutate(state = factor(state, levels = 1:4, labels = c("euthymic","manic","mixed","depressive")),
         from_time = lag(open_from, 1, default = min(open_from)-(3600*3)),
         to_time = open_from-1,
         patient_id = factor(patient_id),
         from_occ = occasion-1,
         to_occ = occasion) %>%
  mutate(day = format(as.POSIXct(to_time), format = "%Y-%m-%d")) %>%
  ungroup()

scale_data <- bs_scales %>%
  group_by(patient_id) %>%
  mutate(patient_id = cur_group_id(),
         patient_id = as.factor(patient_id),
         occasion = row_number(),
         day = from_time_day,
         bs_asrm_tot = case_when(is.na(bs_asrm_tot) ~ 99,
                                 !is.na(bs_asrm_tot) ~ bs_asrm_tot),
         bs_qids_tot = case_when(is.na(bs_qids_tot) ~ 99,
                                 !is.na(bs_qids_tot) ~ bs_qids_tot)
  ) %>%
  ungroup() %>%
  select(-occasion, -open_from)

join_data <- full_join(state_data, scale_data) %>%
  tidyr::fill(bs_asrm_tot:bs_qids_tot, .direction = "down") %>%
  gather(variable, value, -patient_id,-state,-open_from,-occasion,-from_time,-to_time,-from_occ,-to_occ,-from_time_day,-to_time_day,-day) %>%
  mutate(value = case_when(value == 99 ~ NA_real_,
                           value != 99 ~ value))

p <- ggplot() +
  geom_line(data = join_data %>%
              group_by(patient_id) %>%
              mutate(variable = factor(variable,
                                       levels = c("bs_asrm_tot","bs_qids_tot"),
                                       labels = c("ASRM","QIDS")),
                     patient_id = factor(patient_id, levels =  c(15,12,19,17,1,3,
                                                                 6,5,7,18,16,13,
                                                                 20,2,14,11,9,8,4,10),
                                         labels = paste0("patient ",c(15,12,19,17,1,3,
                                                                      6,5,7,18,16,13,
                                                                      20,2,14,11,9,8,4,10)))),
            aes(x = from_occ, y = value, colour = variable, linetype=variable)) +
  geom_rect(data = state_data %>%
              mutate(patient_id = factor(patient_id, levels = c(15,12,19,17,1,3,
                                                                6,5,7,18,16,13,
                                                                20,2,14,11,9,8,4,10),
                                         labels = paste0("patient ",c(15,12,19,17,1,3,
                                                                      6,5,7,18,16,13,
                                                                      20,2,14,11,9,8,4,10)))), aes(xmin = from_occ, xmax = to_occ,
                                                                                                   ymin = 25, ymax = 30, fill = state)) +
  scale_fill_viridis(discrete = TRUE, alpha=0.6,direction = -1, option = "plasma") +
  geom_hline(yintercept = 6, linetype = "dashed") +
  facet_wrap(patient_id~., scales = "free", ncol = 3) +
  theme_minimal() +
  labs(colour = "Questionnaire", fill = "Mood state", linetype = "Questionnaire") +
  theme(legend.position = "bottom") +
  # ggtitle(label = "Temporal alignment between mood states and weekly symptom scores") +
  theme(plot.title = element_text(hjust = 0.5)) +
  ylab(label = "Weekly symptom score") +
  xlab(label = "Measurement occasion") 

p

ggsave(plot = p + theme(text = element_text(size = 12)),
       filename = "outputs/figures/decoding_12.pdf", width = 14, height = 11, units = "in", dpi = 300)
ggsave(plot = p + theme(text = element_text(size = 12)),
       filename = "outputs/figures/decoding_12.tiff", width = 14, height = 11, units = "in", dpi = 300)
ggsave(plot = p + theme(text = element_text(size = 14)),
       filename = "outputs/figures/decoding_14.pdf", width = 14, height = 11, units = "in", dpi = 300)
ggsave(plot = p + theme(text = element_text(size = 14)),
       filename = "outputs/figures/decoding_14.tiff", width = 14, height = 11, units = "in", dpi = 300)
ggsave(plot = p + theme(text = element_text(size = 16)),
       filename = "outputs/figures/decoding_16.pdf", width = 14, height = 11, units = "in", dpi = 300)
ggsave(plot = p + theme(text = element_text(size = 16)),
       filename = "outputs/figures/decoding_16.tiff", width = 14, height = 11, units = "in", dpi = 300)
ggsave(plot = p + theme(text = element_text(size = 18)),
       filename = "outputs/figures/decoding_18.pdf", width = 14, height = 11, units = "in", dpi = 300)
ggsave(plot = p + theme(text = element_text(size = 18)),
       filename = "outputs/figures/decoding_18.tiff", width = 14, height = 11, units = "in", dpi = 300)
ggsave(plot = p + theme(text = element_text(size = 20)),
       filename = "outputs/figures/decoding_20.pdf", width = 14, height = 11, units = "in", dpi = 300)
ggsave(plot = p + theme(text = element_text(size = 20)),
       filename = "outputs/figures/decoding_20.tiff", width = 14, height = 11, units = "in", dpi = 300)





#------------------------------------------------------------------------------#
# Table R.4: sojourn times

# Expected duration given transitions
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

# Number of transitions
switches <- full_join(subj_tpm, states %>%
                        drop_na() %>%
                        group_by(patient_id) %>%
                        summarise(switches = length(rle(as.character(state))$values)/length(state)))

switches_table <- states %>%
  drop_na() %>%
  group_by(patient_id) %>%
  summarise(n_occ = length(state),
            n_switch = length(rle(as.character(state))$values),
            rel_switch = length(rle(as.character(state))$values)/length(state))


# Plot empirical duration
p2 <- ggplot(data = patient_state_duration) +
  geom_line(aes(x = decoding, y = median_dur*3, group = patient_id), linetype = "dashed", alpha = 0.2) +
  geom_jitter(aes(x = decoding, y = median_dur*3, colour = decoding), width = 0.05, alpha = 0.75) +
  geom_point(data = group_state_duration,
                  aes(x = decoding, y = mean_dur*3), shape = 2, size = 3) +
  scale_color_viridis(discrete = TRUE,
                      alpha=0.6,
                      direction = -1,
                      option = "plasma") +
  geom_hline(yintercept = 0) +
  theme_minimal() +
  # scale_y_log10() +
  ylab(label = "Duration in hours") +
  xlab(label = "Mood state") +
  ggtitle("Expected state duration") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(legend.position = "none") +
  theme(strip.text.y = element_text(angle = 0)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  coord_cartesian(ylim = c(0,50))

p2

p <- p2

ggsave(plot = p, filename = "outputs/figures/duration_20.pdf", width = 3, height = 6.1, units = "in", dpi = 300)

ggsave(plot = p + theme(text = element_text(size = 12)),
       filename = "outputs/figures/duration_12.pdf", width = 3.5, height = 5.5, units = "in", dpi = 300)
ggsave(plot = p + theme(text = element_text(size = 12)),
       filename = "outputs/figures/duration_12.tiff", width = 3.5, height = 5.5, units = "in", dpi = 300)
ggsave(plot = p + theme(text = element_text(size = 14)),
       filename = "outputs/figures/duration_14.pdf", width = 3.5, height = 5.5, units = "in", dpi = 300)
ggsave(plot = p + theme(text = element_text(size = 14)),
       filename = "outputs/figures/duration_14.tiff", width = 3.5, height = 5.5, units = "in", dpi = 300)
ggsave(plot = p + theme(text = element_text(size = 16)),
       filename = "outputs/figures/duration_16.pdf", width = 3.5, height = 5.5, units = "in", dpi = 300)
ggsave(plot = p + theme(text = element_text(size = 16)),
       filename = "outputs/figures/duration_16.tiff", width = 3.5, height = 5.5, units = "in", dpi = 300)
ggsave(plot = p + theme(text = element_text(size = 18)),
       filename = "outputs/figures/duration_18.pdf", width = 3.5, height = 5.5, units = "in", dpi = 300)
ggsave(plot = p + theme(text = element_text(size = 18)),
       filename = "outputs/figures/duration_18.tiff", width = 3.5, height = 5.5, units = "in", dpi = 300)
ggsave(plot = p + theme(text = element_text(size = 20)),
       filename = "outputs/figures/duration_20.pdf", width = 3.5, height = 5.5, units = "in", dpi = 300)
ggsave(plot = p + theme(text = element_text(size = 20)),
       filename = "outputs/figures/duration_20.tiff", width = 3.5, height = 5.5, units = "in", dpi = 300)




#------------------------------------------------------------------------------#
# Figure R.3: Three switching patterns (backbone for diagrams)

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


#------------------------------------------------------------------------------#
# Linear mixed model for relation between weekly data and decoding

# Load weekly data
bs_asrm <- haven::read_sav("bs_asrm_cleaned_UU.sav") %>%
  arrange(patient_id)
bs_qids <- haven::read_sav("bs_qids_cleaned_UU.sav") %>%
  arrange(patient_id)

# Put together
bs_scales <- full_join(bs_asrm, bs_qids) %>%
  dplyr::select(patient_id, bs_asrm_open_from, bs_asrm_tot, bs_qids_tot) %>%
  rename("open_from" = "bs_asrm_open_from") %>%
  drop_na()

# Get decoding
states <- vit_mHMM_cont_mar(object = out, s_data = train_df, burn_in = 1000) %>%
  as.data.frame() %>%
  gather(subject, state) %>%
  mutate(subject = factor(subject, levels = paste0("Subj_",1:20))) %>%
  group_by(subject) %>%
  mutate(subject = cur_group_id()) %>%
  ungroup() %>%
  rename("patient_id" = "subject") %>%
  group_by(patient_id) %>%
  mutate(occasion = row_number())

data_labelled <- left_join(as.data.frame(train_df) %>% group_by(patient_id) %>%
                             mutate(occasion = row_number()), states) %>%
  group_by(patient_id) %>%
  mutate(occasion = row_number()) %>%
  as.data.frame()

# Add week info
data_labelled <- bind_cols(data_labelled, data %>%
                             # drop_na() %>%
                             group_by(patient_id) %>%
                             mutate(patient_id = cur_group_id()) %>%
                             ungroup() %>%
                             arrange(patient_id) %>%
                             dplyr::select(open_from, bs_diary_open_from, bs_diary_date))


# Validated scales
state_data <- data_labelled %>%
  select(patient_id, state, bs_diary_open_from, occasion) %>%
  rename("open_from" = "bs_diary_open_from") %>%
  group_by(patient_id) %>%
  mutate(state = factor(state, levels = 1:4, labels = c("euthymic","manic","mixed","depressive")),
         from_time = lag(open_from, 1, default = min(open_from)-(3600*3)),
         to_time = open_from-1,
         patient_id = factor(patient_id),
         from_occ = occasion-1,
         to_occ = occasion) %>%
  ungroup()

# Visualize scales
combined_data <- left_join(
  data_labelled %>%
    select(patient_id, state, bs_diary_open_from) %>%
    mutate(week_number = lubridate::isoweek(lubridate::ymd_hms(bs_diary_open_from))),
  bs_scales %>%
    group_by(patient_id) %>%
    mutate(patient_id = cur_group_id(),
           week_number = lubridate::isoweek(lubridate::ymd_hms(open_from))) %>%
    ungroup())

data_lmer <- rbind(combined_data %>%
                     mutate(days_ema = lubridate::wday(lubridate::ymd_hms(bs_diary_open_from), week_start = 1) - 7,
                            days_val = lubridate::wday(lubridate::ymd_hms(open_from))) %>%
                     filter(days_ema == 0) %>%
                     mutate(window = "same day"),
                   combined_data %>%
                     mutate(days_ema = lubridate::wday(lubridate::ymd_hms(bs_diary_open_from), week_start = 1) - 7,
                            days_val = lubridate::wday(lubridate::ymd_hms(open_from))) %>%
                     filter(days_ema >= -3  & days_ema < 0) %>%
                     mutate(window = "3 days before"),
                   combined_data %>%
                     mutate(days_ema = lubridate::wday(lubridate::ymd_hms(bs_diary_open_from), week_start = 1) - 7,
                            days_val = lubridate::wday(lubridate::ymd_hms(open_from))) %>%
                     filter(days_ema >= -6  & days_ema < 0) %>%
                     mutate(window = "6 days before"),
                   combined_data %>%
                     mutate(days_ema = lubridate::wday(lubridate::ymd_hms(bs_diary_open_from), week_start = 1) - 7,
                            days_val = lubridate::wday(lubridate::ymd_hms(open_from))) %>%
                     filter(days_ema >= -7  & days_ema < 0) %>%
                     mutate(window = "7 days before"),
                   combined_data %>%
                     mutate(days_ema = lubridate::wday(lubridate::ymd_hms(bs_diary_open_from), week_start = 1) - 7,
                            days_val = lubridate::wday(lubridate::ymd_hms(open_from))) %>%
                     filter(days_ema >= -6  & days_ema <= 0) %>%
                     mutate(window = "same and 6 days before")) %>%
  mutate(state = factor(state, labels = c("euthymic","manic","mixed","depressive")),
         window = factor(window, levels = c("same day", "3 days before", "6 days before", "7 days before", "same and 6 days before"))) %>%
  ungroup() %>%
  rename("ASRM" = "bs_asrm_tot",
         "QIDS" = "bs_qids_tot") %>%
  dplyr::select(patient_id, state, ASRM, QIDS, window)

# Fit linear mixed models
data_lmer <- data_lmer %>%
  filter(window == "same day") %>%
  filter(window == "3 days before") %>%
  dplyr::select(patient_id, state, ASRM, QIDS)

# Intercept only
model_int_asrm <- lmer(ASRM ~ 1 + (1|patient_id), data = data_lmer, REML = FALSE)
summary(model_int_asrm)
icc(model_int_asrm)
plot(model_int_asrm)

# Fixed effect of state + random intercept
model_state_asrm <- lmer(ASRM ~ 1 + state + (1|patient_id), data = data_lmer, REML = FALSE)
summary(model_state_asrm)
icc(model_state_asrm)
plot(model_state_asrm)

# Fixed effect of state + random slope + random intercept
model_state_slope_asrm <- lmer(ASRM ~ 1 + state + (1 + state | patient_id), data = data_lmer, REML = FALSE)
summary(model_state_slope_asrm)
icc(model_state_slope_asrm)
plot(model_state_slope_asrm)

qqnorm(resid(model_int_asrm))
qqline(resid(model_int_asrm), col = "red")

qqnorm(resid(model_state_slope_asrm))
qqline(resid(model_state_slope_asrm), col = "red")

r2(model_int_asrm)
r2(model_state_asrm)
r2(model_state_slope_asrm)

r.squaredGLMM(model_int_asrm)
r.squaredGLMM(model_state_asrm)
r.squaredGLMM(model_state_slope_asrm)


# Intercept only
model_int_qids <- lmer(QIDS ~ 1 + (1|patient_id), data = data_lmer, REML = FALSE)
summary(model_int_qids)
icc(model_int_qids)
plot(model_int_qids)

# Fixed effect of state
model_state_qids <- lmer(QIDS ~ 1 + state + (1|patient_id), data = data_lmer, REML = FALSE)
summary(model_state_qids)
icc(model_state_qids)
plot(model_state_qids)

# Fixed effect of state + random slope + random intercept
model_state_slope_qids <- lmer(QIDS ~ 1 + state + (1 + state | patient_id), data = data_lmer, REML = FALSE)
summary(model_state_slope_qids)
icc(model_state_slope_qids)
plot(model_state_slope_qids)




