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

# load data
data <- haven::read_sav("data/ESM bipolar cleaned data_UU project.sav")

#------------------------------------------------------------------------------#

# All
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

#------------------------------------------------------------------------------#

# Load output
out <- readRDS(paste0("outputs/results/truncated/out_cont_gamma_prior_emiss_prior_m4_12dv_it4000_c",2,".rds"))


#==============================================================================#
# Main results (figures and tables)
#==============================================================================#

#------------------------------------------------------------------------------#
# Table R.1: patient characteristics


#------------------------------------------------------------------------------#
# Table R.2: AIC states

# Model complexity    AIC
#   2 states        50060.8
#   3 states        48748.7
#   4 states        48231.1


#------------------------------------------------------------------------------#
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
         mu = factor(mu, levels = paste0("S",1:4), labels = c("euthymic","manic","mixed","depressive"))) %>%
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
         mu = factor(mu, levels = paste0("mu_",1:4), labels = c("euthymic","manic","mixed","depressive"))) %>%
  rename("state" = "mu")

p <- ggplot(data = subj_emiss_nice) +
  geom_jitter(aes(x = dep, y = value, colour = dep), width = 0.1, alpha = 0.3) +
  geom_pointrange(data = group_emiss_nice,
                  aes(x = dep, y = group_mean,
                      ymin = group_CCIlwr, ymax = group_CCIupr)) +
  scale_color_manual(values = safe_colorblind_palette) +
  geom_hline(yintercept = 50, linetype = "dashed") +
  geom_hline(yintercept = 0) +
  facet_grid(state~construct, scales = "free_x") +
  theme_minimal() +
  ylab(label = "EMA score") +
  xlab(label = "") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(legend.position = "none") +
  ggtitle(label = "Composition of momentary mood states: group-level and\npatient-specific item emission scores by state") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(strip.text.y = element_text(angle = 0))  

p

ggsave(plot = p + theme(text = element_text(size = 12)),
       filename = "outputs/figures/emissions_12.pdf", width = 9.5, height = 9, units = "in", dpi = 300)
ggsave(plot = p + theme(text = element_text(size = 12)),
       filename = "outputs/figures/emissions_12.tiff", width = 9.5, height = 9, units = "in", dpi = 300)
ggsave(plot = p + theme(text = element_text(size = 14)),
       filename = "outputs/figures/emissions_14.pdf", width = 9.5, height = 9, units = "in", dpi = 300)
ggsave(plot = p + theme(text = element_text(size = 14)),
       filename = "outputs/figures/emissions_14.tiff", width = 9.5, height = 9, units = "in", dpi = 300)
ggsave(plot = p + theme(text = element_text(size = 16)),
       filename = "outputs/figures/emissions_16.pdf", width = 9.5, height = 9, units = "in", dpi = 300)
ggsave(plot = p + theme(text = element_text(size = 16)),
       filename = "outputs/figures/emissions_16.tiff", width = 9.5, height = 9, units = "in", dpi = 300)
ggsave(plot = p + theme(text = element_text(size = 18)),
       filename = "outputs/figures/emissions_18.pdf", width = 9.5, height = 9, units = "in", dpi = 300)
ggsave(plot = p + theme(text = element_text(size = 18)),
       filename = "outputs/figures/emissions_18.tiff", width = 9.5, height = 9, units = "in", dpi = 300)
ggsave(plot = p + theme(text = element_text(size = 20)),
       filename = "outputs/figures/emissions_20.pdf", width = 9.5, height = 9, units = "in", dpi = 300)
ggsave(plot = p + theme(text = element_text(size = 20)),
       filename = "outputs/figures/emissions_20.tiff", width = 9.5, height = 9, units = "in", dpi = 300)


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
                     patient_id = factor(patient_id, levels = 1:20, labels = paste0("patient ",1:20))),
            aes(x = from_occ, y = value, colour = variable)) +
  geom_rect(data = state_data %>%
              mutate(patient_id = factor(patient_id, levels = 1:20, labels = paste0("patient ",1:20))), aes(xmin = from_occ, xmax = to_occ,
                                                                                                            ymin = 25, ymax = 30, fill = state)) +
  scale_fill_viridis(discrete = TRUE, alpha=0.6,direction = -1, option = "plasma") +
  geom_hline(yintercept = 6, linetype = "dashed") +
  facet_wrap(patient_id~., scales = "free", ncol = 4) +
  theme_minimal() +
  labs(colour = "Questionnaire", fill = "Mood state") +
  theme(legend.position = "bottom") +
  ggtitle(label = "Temporal alignment between mood states and weekly symptom scores") +
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



#==============================================================================#
# Tile plots for transitions and emissions
#==============================================================================#

#------------------------------------------------------------------------------#
# Group-level means

# Transition plot
p1 <- out$gamma_prob_bar %>%
  as.data.frame() %>%
  mutate(iter = row_number()) %>%
  gather(state, value, -iter) %>%
  filter(iter > 1000) %>%
  group_by(state) %>%
  summarise(median_value = median(value)) %>%
  separate(state, into = c("from","to"), sep = "to") %>%
  mutate(from = factor(from, levels = paste0("S",4:1) ),
         to = factor(to, levels = paste0("S",1:4) )) %>%
  mutate(from = factor(from, levels = paste0("S",4:1), labels = c("depressive","mixed","manic","euthymic")),
         to = factor(to, levels = paste0("S",1:4), labels = c("euthymic","manic","mixed","depressive") )) %>%
  ggplot(aes(x = to, y = from, fill = median_value)) +
  geom_tile() +
  geom_text(aes(label = format(round(median_value, digits = 2), nsmall = 2)
  )) +
  scale_fill_distiller(palette = "Spectral", limits = c(0,1)) +
  theme_minimal() +
  theme(legend.position = "right") +
  xlab("To mood state") +
  ylab("From mood state") +
  theme(axis.text.x = element_text(angle = 45, hjust=1),
        axis.text.y = element_text(angle = 45)) +
  ggtitle(label = "Group-level switching probabilities") +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(fill = "") +
  theme(legend.key.height= unit(2, 'cm'),
        legend.key.width= unit(0.3, 'cm'))

p1

ggsave(plot = p + theme(text = element_text(size = 12)),
       filename = "outputs/figures/transitions_group_12.pdf", width = 6, height = 6.1, units = "in", dpi = 300)
ggsave(plot = p + theme(text = element_text(size = 12)),
       filename = "outputs/figures/transitions_group_12.tiff", width = 6, height = 6.1, units = "in", dpi = 300)
ggsave(plot = p + theme(text = element_text(size = 14)),
       filename = "outputs/figures/transitions_group_14.pdf", width = 6, height = 6.1, units = "in", dpi = 300)
ggsave(plot = p + theme(text = element_text(size = 14)),
       filename = "outputs/figures/transitions_group_14.tiff", width = 6, height = 6.1, units = "in", dpi = 300)
ggsave(plot = p + theme(text = element_text(size = 16)),
       filename = "outputs/figures/transitions_group_16.pdf", width = 6, height = 6.1, units = "in", dpi = 300)
ggsave(plot = p + theme(text = element_text(size = 16)),
       filename = "outputs/figures/transitions_group_16.tiff", width = 6, height = 6.1, units = "in", dpi = 300)
ggsave(plot = p + theme(text = element_text(size = 18)),
       filename = "outputs/figures/transitions_group_18.pdf", width = 6, height = 6.1, units = "in", dpi = 300)
ggsave(plot = p + theme(text = element_text(size = 18)),
       filename = "outputs/figures/transitions_group_18.tiff", width = 6, height = 6.1, units = "in", dpi = 300)
ggsave(plot = p + theme(text = element_text(size = 20)),
       filename = "outputs/figures/transitions_group_20.pdf", width = 6, height = 6.1, units = "in", dpi = 300)
ggsave(plot = p + theme(text = element_text(size = 20)),
       filename = "outputs/figures/transitions_group_20.tiff", width = 6, height = 6.1, units = "in", dpi = 300)

l <- ggpubr::as_ggplot(ggpubr::get_legend(p +
                                            theme(text = element_text(size = 20)) +
                                            theme(legend.key.height= unit(1, 'in'),
                                                  legend.key.width= unit(0.25, 'in'))))

ggsave(plot = l, filename = "outputs/figures/transitions_legend_12.tiff", width = 1, height = 6, units = "in", dpi = 300)
ggsave(plot = l, filename = "outputs/figures/transitions_legend_14.tiff", width = 1, height = 6, units = "in", dpi = 300)
ggsave(plot = l, filename = "outputs/figures/transitions_legend_16.tiff", width = 1, height = 6, units = "in", dpi = 300)
ggsave(plot = l, filename = "outputs/figures/transitions_legend_18.tiff", width = 1, height = 6, units = "in", dpi = 300)
ggsave(plot = l, filename = "outputs/figures/transitions_legend_20.tiff", width = 1, height = 6, units = "in", dpi = 300)

#------------------------------------------------------------------------------#
# Patient-specific means

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

p3 <- subj_tpm %>%
  gather(state, value, -patient_id) %>%
  mutate(type = case_when(patient_id == 18 ~ "patient 18",
                          patient_id == 8 ~ "patient 8",
                          patient_id == 1 ~ "patient 1",
                          patient_id == 2 ~ "patient 2")) %>%
  drop_na(type) %>%
  separate(state, into = c("from","to"), sep = "to") %>%
  mutate(from = factor(from, levels = paste0("S",4:1), labels = c("depressive","mixed","manic","euthymic")),
         to = factor(to, levels = paste0("S",1:4), labels = c("euthymic","manic","mixed","depressive") ),
         type = factor(type, levels = paste0("patient ",c(18,8,1,2)))) %>%
  ggplot(aes(x = to, y = from, fill = value)) +
  geom_tile() +
  scale_fill_distiller(palette = "Spectral", limits = c(0,1)) +
  geom_text(aes(label=ifelse(format(round(value,2), nsmall = 2) < 0.01, "<0.01", format(round(value,2), nsmall = 2))), size = 3) +
  facet_wrap(type~., nrow = 1) +
  theme_minimal() +
  theme(legend.position = "right") +
  xlab("To mood state") +
  ylab("From mood state") +
  theme(axis.text.x = element_text(angle = 45, hjust=1),
        axis.text.y = element_text(angle = 45)) +
  ggtitle("Patient-specific switching probabilities") +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(fill = "") +
  theme(legend.key.height= unit(1, 'cm'),
        legend.key.width= unit(0.3, 'cm'))

ggsave(plot = p + theme(text = element_text(size = 12)),
       filename = "outputs/figures/transitions_patients_12.pdf", width = 10, height = 4.1, units = "in", dpi = 300)
ggsave(plot = p + theme(text = element_text(size = 14)),
       filename = "outputs/figures/transitions_patients_14.pdf", width = 10, height = 4.1, units = "in", dpi = 300)
ggsave(plot = p + theme(text = element_text(size = 16)),
       filename = "outputs/figures/transitions_patients_16.pdf", width = 10, height = 4.1, units = "in", dpi = 300)
ggsave(plot = p + theme(text = element_text(size = 18)),
       filename = "outputs/figures/transitions_patients_18.pdf", width = 10, height = 4.1, units = "in", dpi = 300)
ggsave(plot = p + theme(text = element_text(size = 20)),
       filename = "outputs/figures/transitions_patients_20.pdf", width = 10, height = 4.1, units = "in", dpi = 300)

ggsave(plot = p + theme(text = element_text(size = 12)),
       filename = "outputs/figures/transitions_patients_12.tiff", width = 10, height = 4.1, units = "in", dpi = 300)
ggsave(plot = p + theme(text = element_text(size = 14)),
       filename = "outputs/figures/transitions_patients_14.tiff", width = 10, height = 4.1, units = "in", dpi = 300)
ggsave(plot = p + theme(text = element_text(size = 16)),
       filename = "outputs/figures/transitions_patients_16.tiff", width = 10, height = 4.1, units = "in", dpi = 300)
ggsave(plot = p + theme(text = element_text(size = 18)),
       filename = "outputs/figures/transitions_patients_18.tiff", width = 10, height = 4.1, units = "in", dpi = 300)
ggsave(plot = p + theme(text = element_text(size = 20)),
       filename = "outputs/figures/transitions_patients_20.tiff", width = 10, height = 4.1, units = "in", dpi = 300)

l <- ggpubr::as_ggplot(ggpubr::get_legend(p + theme(legend.key.height= unit(0.7, 'in'),
                                                    legend.key.width= unit(0.25, 'in'))))
ggsave(plot = l, filename = "outputs/figures/transitions_legend_small_20.pdf", width = 1, height = 4, units = "in", dpi = 300)


# Place plots together with cowplot
library(cowplot)

legend <- cowplot::get_legend(subj_tpm %>%
                                gather(state, value, -patient_id) %>%
                                mutate(type = case_when(patient_id == 18 ~ "patient 18",
                                                        patient_id == 8 ~ "patient 8",
                                                        patient_id == 1 ~ "patient 1",
                                                        patient_id == 2 ~ "patient 2")) %>%
                                drop_na(type) %>%
                                separate(state, into = c("from","to"), sep = "to") %>%
                                mutate(from = factor(from, levels = paste0("S",4:1), labels = c("depressive","mixed","manic","euthymic")),
                                       to = factor(to, levels = paste0("S",1:4), labels = c("euthymic","manic","mixed","depressive") ),
                                       type = factor(type, levels = paste0("patient ",c(18,8,1,2)))) %>%
                                ggplot(aes(x = to, y = from, fill = value)) +
                                geom_tile() +
                                scale_fill_distiller(palette = "Spectral", limits = c(0,1)) +
                                geom_text(aes(label=ifelse(format(round(value,2), nsmall = 2) < 0.01, "<0.01", format(round(value,2), nsmall = 2)) )) +
                                facet_wrap(type~., nrow = 1) +
                                theme_minimal() +
                                theme(legend.position = "bottom") +
                                xlab("To mood state") +
                                ylab("From mood state") +
                                theme(axis.text.x = element_text(angle = 45, hjust=1),
                                      axis.text.y = element_text(angle = 45)) +
                                ggtitle("Patient-specific switching probabilities") +
                                theme(plot.title = element_text(hjust = 0.5)) +
                                labs(fill = "") +
                                theme(legend.key.height= unit(0.1, 'cm'),
                                      legend.key.width= unit(2, 'cm')))

# Create combined figured
p_top <- cowplot::plot_grid(p1, p2, labels=c("a","b"), rel_widths = c(1.4,1), ncol = 2, nrow = 1)

p_bottom <- cowplot::plot_grid(p3, legend, labels=c("c",""),  rel_widths = c(1,0.1), ncol = 2, nrow = 1)

p <- cowplot::plot_grid(p_top, p3, labels=c("","c"),  rel_heights = c(1.5,1), ncol = 1, nrow = 2)

p

# Save figure
ggsave(plot = p, filename = "outputs/figures/transitions_combined.pdf", width = 10, height = 9, units = "in", dpi = 300)


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


#------------------------------------------------------------------------------#
# Figure R.5: validated scales (boxplot)

# Calculate average ASRM & QIDS per state

# Visualize scales
combined_data <- left_join(
  data_labelled %>%
    select(patient_id, state, bs_diary_open_from) %>%
    # rename("open_from" = "bs_diary_date") %>%
    mutate(week_number = lubridate::isoweek(lubridate::ymd_hms(bs_diary_open_from))),
  bs_scales %>%
    group_by(patient_id) %>%
    mutate(patient_id = cur_group_id(),
           week_number = lubridate::isoweek(lubridate::ymd_hms(open_from))) %>%
    ungroup())

# Various days together
retrospective_data <- rbind(combined_data %>%
                              mutate(days_ema = lubridate::wday(lubridate::ymd_hms(bs_diary_open_from), week_start = 1) - 7,
                                     days_val = lubridate::wday(lubridate::ymd_hms(open_from))) %>%
                              filter(days_ema == 0) %>%
                              select(state, days_ema, bs_asrm_tot, bs_qids_tot) %>%
                              gather(variable, value, -state, -days_ema) %>%
                              mutate(window = "same day"),
                            combined_data %>%
                              mutate(days_ema = lubridate::wday(lubridate::ymd_hms(bs_diary_open_from), week_start = 1) - 7,
                                     days_val = lubridate::wday(lubridate::ymd_hms(open_from))) %>%
                              filter(days_ema >= -3  & days_ema < 0) %>%
                              select(state, days_ema, bs_asrm_tot, bs_qids_tot) %>%
                              gather(variable, value, -state, -days_ema) %>%
                              mutate(window = "3 days before"),
                            combined_data %>%
                              mutate(days_ema = lubridate::wday(lubridate::ymd_hms(bs_diary_open_from), week_start = 1) - 7,
                                     days_val = lubridate::wday(lubridate::ymd_hms(open_from))) %>%
                              filter(days_ema >= -6  & days_ema < 0) %>%
                              select(state, days_ema, bs_asrm_tot, bs_qids_tot) %>%
                              gather(variable, value, -state, -days_ema) %>%
                              mutate(window = "6 days before")) %>%
  mutate(state = factor(state, labels = c("euthymic","manic","mixed","depressive")),
         window = factor(window, levels = c("same day", "3 days before", "6 days before"))) %>%
  ungroup() %>%
  rename("questionnaire" = "variable") %>%
  mutate(questionnaire = factor(questionnaire,
                                levels = c("bs_asrm_tot","bs_qids_tot"),
                                labels = c("ASRM","QIDS")))


p <- retrospective_data %>%
  mutate(window = factor(window,
                         levels = c("same day", "3 days before", "6 days before"),
                         labels = c("same day", "three days before", "six days before"))) %>%
  ggplot(aes(x = window, y = value, fill = state)) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_viridis(discrete = TRUE, alpha=0.8, direction = -1, option = "plasma") +
  geom_point(position=position_jitterdodge(jitter.height = 0.5, jitter.width = 0.5), size=0.3, alpha=0.05) +
  # geom_jitter(color="black", size=0.4, alpha=0.4) +
  facet_wrap(questionnaire~., scales = "free_y", nrow = 2) +
  theme_minimal() +
  theme(strip.text.y = element_text(angle = 0)) +
  theme(legend.position = "bottom") +
  ylab(label = "Weekly symptom score") +
  xlab(label = "Reporting window") +
  ggtitle(label = "Effect of mood state on weekly symptom\nscores on three reporting windows") +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(fill = "Mood state")

p

ggsave(plot = p + theme(text = element_text(size = 12)),
       filename = "outputs/figures/boxplot_12.pdf", width = 7, height = 9, units = "in", dpi = 300)
ggsave(plot = p + theme(text = element_text(size = 12)),
       filename = "outputs/figures/boxplot_12.tiff", width = 7, height = 9, units = "in", dpi = 300)
ggsave(plot = p + theme(text = element_text(size = 14)),
       filename = "outputs/figures/boxplot_14.pdf", width = 7, height = 9, units = "in", dpi = 300)
ggsave(plot = p + theme(text = element_text(size = 14)),
       filename = "outputs/figures/boxplot_14.tiff", width = 7, height = 9, units = "in", dpi = 300)
ggsave(plot = p + theme(text = element_text(size = 16)),
       filename = "outputs/figures/boxplot_16.pdf", width = 7, height = 9, units = "in", dpi = 300)
ggsave(plot = p + theme(text = element_text(size = 16)),
       filename = "outputs/figures/boxplot_16.tiff", width = 7, height = 9, units = "in", dpi = 300)
ggsave(plot = p + theme(text = element_text(size = 18)),
       filename = "outputs/figures/boxplot_18.pdf", width = 7, height = 9, units = "in", dpi = 300)
ggsave(plot = p + theme(text = element_text(size = 18)),
       filename = "outputs/figures/boxplot_18.tiff", width = 7, height = 9, units = "in", dpi = 300)
ggsave(plot = p + theme(text = element_text(size = 20)),
       filename = "outputs/figures/boxplot_20.pdf", width = 7, height = 9, units = "in", dpi = 300)
ggsave(plot = p + theme(text = element_text(size = 20)),
       filename = "outputs/figures/boxplot_20.tiff", width = 7, height = 9, units = "in", dpi = 300)



#==============================================================================#
# Supplementary results
#==============================================================================#


#------------------------------------------------------------------------------#
# Table S.1: patient-specific TPMs

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


#------------------------------------------------------------------------------#
# Table S.2: state decoding before and after missing points

# States assignment in sample
data_labelled %>%
  mutate(miss = is.na(bs_diary_date)) %>%
  group_by(state) %>%
  summarise(n_state = n()) %>%
  summarise(prop_state = n_state/sum(n_state)) %>%
  ungroup() %>%
  mutate(state = rep(c("euthymic","manic","mixed","depressive"),),
         state = factor(state, levels = c("euthymic","manic","mixed","depressive"))) %>%
  ungroup()

# Bonus: states by missingness
data_labelled %>%
  mutate(miss = is.na(bs_diary_date)) %>%
  group_by(miss, state) %>%
  summarise(n_state = n()) %>%
  summarise(prop_state = n_state/sum(n_state)) %>%
  ungroup() %>%
  mutate(state = rep(c("euthymic","manic","mixed","depressive"),2),
         state = factor(state, levels = c("euthymic","manic","mixed","depressive"))) %>%
  ungroup() %>%
  spread(miss, prop_state)

# Create ID df with original patient IDs
id_dix <- data.frame("patient_id" = 1:20,
                     "patient_ID" = unique(data$patient_id))
# Get sequence of states:
ldecode <- do.call(rbind, lapply(1:length(out$sample_path), function(s){
  cbind(s,
        1:nrow(out$sample_path[[s]]),
        t(apply(out$sample_path[[s]],1,tabulate, nbins=4))/3999)
}))

ldecode <- cbind(ldecode, apply(ldecode[,3:6],1,which.max))
ldecode <- as.data.frame(ldecode)
names(ldecode) <- c("patient_id", "occasion", paste0("fw_prob_S",1:4), "ldecode")
ldecode <- cbind(ldecode, "miss" = apply(train_df,1,anyNA))

# States before missing observation
ldecode %>%
  group_by(patient_id) %>%
  mutate(one_before = ifelse(miss == FALSE & lead(miss) == TRUE & patient_id == lead(patient_id, n = 1), TRUE, FALSE),
         one_after = ifelse(miss == FALSE & lag(miss, n = 1) == TRUE  & patient_id == lag(patient_id, n = 1), TRUE, FALSE)) %>%
  filter(one_before == TRUE) %>%
  group_by(ldecode) %>%
  summarise(n_state = n()) %>%
  summarise(prop_state = n_state/sum(n_state)) %>%
  ungroup() %>%
  mutate(state = rep(c("euthymic","manic","mixed","depressive"),1),
         state = factor(state, levels = c("euthymic","manic","mixed","depressive"))) %>%
  dplyr::select(state, prop_state)

# States after missing observation
ldecode %>%
  group_by(patient_id) %>%
  mutate(one_before = ifelse(miss == FALSE & lead(miss) == TRUE & patient_id == lead(patient_id, n = 1), TRUE, FALSE),
         one_after = ifelse(miss == FALSE & lag(miss, n = 1) == TRUE  & patient_id == lag(patient_id, n = 1), TRUE, FALSE)) %>%
  filter(one_after == TRUE) %>%
  group_by(ldecode) %>%
  summarise(n_state = n()) %>%
  summarise(prop_state = n_state/sum(n_state)) %>%
  ungroup() %>%
  mutate(state = rep(c("euthymic","manic","mixed","depressive"),1),
         state = factor(state, levels = c("euthymic","manic","mixed","depressive"))) %>%
  dplyr::select(state, prop_state)

#------------------------------------------------------------------------------#
# Figure S.4: stacked dependent variables with overlayed decoding

# Stacking Dbs
data_labelled %>%
  filter(patient_id %in% c(7,8,10,2)) %>%
  gather(variable, value, -patient_id, -occasion, -state) %>%
  mutate(patient_id = factor(patient_id, levels = c(7,8,10,2), labels = paste0("Patient ",c(7,8,10,2))),
         state = factor(state, labels = c("euthymic","manic","mixed","depressive")),
         variable = factor(variable, levels = c("bs_diary_5", "bs_diary_13", "bs_diary_22",
                                                "bs_diary_15", "bs_diary_9", "bs_diary_10",
                                                "bs_diary_7", "bs_diary_11", "bs_diary_17",
                                                "bs_diary_8", "bs_diary_14", "bs_diary_16"))) %>%
  filter(variable %in% c("bs_diary_5", "bs_diary_13", "bs_diary_22",
                         "bs_diary_7", "bs_diary_11", "bs_diary_17")) %>%
  ggplot() +
  geom_line(aes(x = occasion, y = value)) +
  geom_rect(aes(xmin = occasion-1, xmax = occasion,
                ymin = 100, ymax = 140, fill = state)) +
  # scale_fill_viridis(discrete = TRUE, alpha=0.6, direction = -1) +
  scale_fill_viridis(discrete = TRUE, alpha=0.6, direction = -1, option = "plasma") +
  facet_grid(variable~patient_id, scales = "free_x") +
  theme_minimal() +
  ylab(label = "Value") +
  xlab(label = "Occasion")

#------------------------------------------------------------------------------#
# Figure S.3: posterior predictive checks

m <- out$input$m
n_dep <- out$input$n_dep
burn_in <- out$input$burn_in

## Fixed effects:
# Transitions
gamma_ppc <- out$gamma_int_bar %>%
  as.data.frame() %>%
  mutate(iter = row_number()) %>%
  filter(iter > burn_in) %>%
  dplyr::select(-iter) %>%
  summarise(across(.cols = everything(), median)) %>%
  as.numeric() %>%
  matrix(., nrow = m, byrow = TRUE) %>%
  int_to_prob()

gamma_ppc[1,1] <- gamma_ppc[1,1] + 1-apply(gamma_ppc, 1, sum)[1]
gamma_ppc[2,2] <- gamma_ppc[2,2] + 1-apply(gamma_ppc, 1, sum)[2]
gamma_ppc[3,3] <- gamma_ppc[3,3] + 1-apply(gamma_ppc, 1, sum)[3]
gamma_ppc[4,4] <- gamma_ppc[4,4] + 1-apply(gamma_ppc, 1, sum)[4]

# Emissions
emiss_ppc <- lapply(1:length(out$emiss_mu_bar), function(s) {
  
  emiss <- out$emiss_mu_bar[[s]] %>%
    as.data.frame() %>%
    mutate(iter = row_number()) %>%
    filter(iter > burn_in) %>%
    dplyr::select(-iter) %>%
    summarise(across(.cols = everything(), median)) %>%
    gather(mu, value) %>%
    pull(value) %>%
    matrix(., nrow = m)
  
  emiss_var<- out$emiss_var_bar[[s]] %>%
    as.data.frame() %>%
    mutate(iter = row_number()) %>%
    filter(iter > burn_in) %>%
    dplyr::select(-iter) %>%
    summarise(across(.cols = everything(), median)) %>%
    gather(mu, value) %>%
    pull(value) %>%
    matrix(., nrow = m)
  
  cbind(emiss, emiss_var)
  
})


## Random effect: between subject variance
# Transitions
gamma_var_ppc <- apply(out$gamma_V_int_bar[1001:4000,] %>%
                         as.data.frame() %>%
                         dplyr::select(paste0("var_int_S",rep(1:4,each=3),"toS",2:4,"_with_int_S",rep(1:4,each=3),"toS",2:4)),
                       2, median) %>%
  as.numeric()
gamma_var_ppc <- matrix(gamma_var_ppc, nrow = m, ncol = m-1, byrow = TRUE)

# Emissions
emiss_varmu_ppc <- lapply(out$emiss_varmu_bar, function(emiss) {
  emiss %>%
    as.data.frame() %>%
    mutate(iter = row_number()) %>%
    filter(iter > burn_in) %>%
    dplyr::select(-iter) %>%
    summarise(across(.cols = everything(), median)) %>%
    gather(mu, value) %>%
    dplyr::select(value) %>%
    as.matrix()
})


# Simulate data
set.seed(42)
sim_data <- pbapply::pblapply(1:500, function(s) mHMMbayes::sim_mHMM_plnorm(n_t = 642, n = 20,
                                                                            data_distr = "continuous",
                                                                            m = out$input$m, n_dep = out$input$n_dep,
                                                                            gamma = gamma_ppc,
                                                                            emiss_distr = emiss_ppc,
                                                                            var_gamma = gamma_var_ppc,
                                                                            var_emiss = emiss_varmu_ppc,
                                                                            return_ind_par = TRUE))

## PPC1: group-level mean

true_data <- train_df %>%
  as.data.frame() %>%
  group_by(patient_id) %>%
  mutate(occasion = row_number()) %>%
  gather(variable, value, -patient_id, -occasion)

ppc_data <- do.call(rbind, pbapply::pblapply(1:length(sim_data), function(s){
  
  data <- as.data.frame(sim_data[[s]]$obs)
  names(data) <- c("patient_id", out$input$dep_labels)
  data %>%
    group_by(patient_id) %>%
    mutate(occasion = row_number()) %>%
    gather(variable, value, -patient_id, -occasion) %>%
    mutate(rep = s)
  
}))

ppc1_truth <- true_data %>%
  group_by(variable) %>%
  summarise(mean_value = mean(value, na.rm = TRUE), sd_value = sd(value, na.rm = TRUE),
            min_value = min(value, na.rm = TRUE), max_value = max(value, na.rm = TRUE))

dep_vars <- c("bs_diary_5", "bs_diary_9", "bs_diary_10",
              "bs_diary_13", "bs_diary_15", "bs_diary_22")
man_vars <- c("bs_diary_7", "bs_diary_8", "bs_diary_11",
              "bs_diary_14", "bs_diary_16", "bs_diary_17")

# Mean
p <- ppc_data %>%
  group_by(variable, rep) %>%
  summarise(mean_value = mean(value), sd_value = sd(value)) %>%
  mutate(variable = factor(variable, levels = c("bs_diary_5", "bs_diary_13", "bs_diary_22",
                                                "bs_diary_15", "bs_diary_9", "bs_diary_10",
                                                "bs_diary_7", "bs_diary_11", "bs_diary_17",
                                                "bs_diary_8", "bs_diary_14", "bs_diary_16"),
                           labels = c("down","dread nrest of day","worry",
                                      "inadequate", "tired", "content",
                                      "agitated", "irritated", "switch and focus",
                                      "extremely well", "full of ideas", "thoughts are racing"))) %>%
  ggplot(aes(x = mean_value)) +
  geom_histogram() +
  geom_vline(data = ppc1_truth %>%
               mutate(variable = factor(variable, levels = c("bs_diary_5", "bs_diary_13", "bs_diary_22",
                                                             "bs_diary_15", "bs_diary_9", "bs_diary_10",
                                                             "bs_diary_7", "bs_diary_11", "bs_diary_17",
                                                             "bs_diary_8", "bs_diary_14", "bs_diary_16"),
                                        labels = c("down","dread nrest of day","worry",
                                                   "inadequate", "tired", "content",
                                                   "agitated", "irritated", "switch and focus",
                                                   "extremely well", "full of ideas", "thoughts are racing"))), aes(xintercept = mean_value), colour = "dark red") +
  facet_wrap(variable~., nrow = 2) +
  theme_minimal() +
  ylab(label = "Counts") +
  xlab(label = "Value")

p

ggsave(plot = p + theme(text = element_text(size = 12)),
       filename = "outputs/figures/group_mean_ppc_12.pdf", width = 9, height = 6, units = "in", dpi = 300)
ggsave(plot = p + theme(text = element_text(size = 12)),
       filename = "outputs/figures/group_mean_ppc_12.tiff", width = 9, height = 6, units = "in", dpi = 300)
ggsave(plot = p + theme(text = element_text(size = 14)),
       filename = "outputs/figures/group_mean_ppc_14.pdf", width = 9, height = 6, units = "in", dpi = 300)
ggsave(plot = p + theme(text = element_text(size = 14)),
       filename = "outputs/figures/group_mean_ppc_14.tiff", width = 9, height = 6, units = "in", dpi = 300)
ggsave(plot = p + theme(text = element_text(size = 16)),
       filename = "outputs/figures/group_mean_ppc_16.pdf", width = 9, height = 6, units = "in", dpi = 300)
ggsave(plot = p + theme(text = element_text(size = 16)),
       filename = "outputs/figures/group_mean_ppc_16.tiff", width = 9, height = 6, units = "in", dpi = 300)
ggsave(plot = p + theme(text = element_text(size = 18)),
       filename = "outputs/figures/group_mean_ppc_18.pdf", width = 9, height = 6, units = "in", dpi = 300)
ggsave(plot = p + theme(text = element_text(size = 18)),
       filename = "outputs/figures/group_mean_ppc_18.tiff", width = 9, height = 6, units = "in", dpi = 300)
ggsave(plot = p + theme(text = element_text(size = 20)),
       filename = "outputs/figures/group_mean_ppc_20.pdf", width = 9, height = 6, units = "in", dpi = 300)
ggsave(plot = p + theme(text = element_text(size = 20)),
       filename = "outputs/figures/group_mean_ppc_20.tiff", width = 9, height = 6, units = "in", dpi = 300)




## PPC2: between-subject variation

# SD
ppc_data %>%
  group_by(variable, rep) %>%
  summarise(mean_value = mean(value), sd_value = sd(value)) %>%
  mutate(variable = factor(variable, levels = c("bs_diary_5", "bs_diary_13", "bs_diary_22",
                                                "bs_diary_15", "bs_diary_9", "bs_diary_10",
                                                "bs_diary_7", "bs_diary_11", "bs_diary_17",
                                                "bs_diary_8", "bs_diary_14", "bs_diary_16"),
                           labels = c("down","dread nrest of day","worry",
                                      "inadequate", "tired", "content",
                                      "agitated", "irritated", "switch and focus",
                                      "extremely well", "full of ideas", "thoughts are racing"))) %>%
  ggplot(aes(x = sd_value)) +
  geom_histogram() +
  geom_vline(data = ppc1_truth %>% mutate(variable = factor(variable, levels = c("bs_diary_5", "bs_diary_13", "bs_diary_22",
                                                                                 "bs_diary_15", "bs_diary_9", "bs_diary_10",
                                                                                 "bs_diary_7", "bs_diary_11", "bs_diary_17",
                                                                                 "bs_diary_8", "bs_diary_14", "bs_diary_16"),
                                                            labels = c("down","dread nrest of day","worry",
                                                                       "inadequate", "tired", "content",
                                                                       "agitated", "irritated", "switch and focus",
                                                                       "extremely well", "full of ideas", "thoughts are racing"))), aes(xintercept = sd_value), colour = "dark red") +
  facet_wrap(variable~., nrow = 2) +
  theme_minimal() +
  ylab(label = "Counts") +
  xlab(label = "Value")


## PPC3: ranked ids plots

# Patient means
p <- ppc_data %>%
  mutate(variable = factor(variable, levels = c("bs_diary_5", "bs_diary_13", "bs_diary_22",
                                                "bs_diary_15", "bs_diary_9", "bs_diary_10",
                                                "bs_diary_7", "bs_diary_11", "bs_diary_17",
                                                "bs_diary_8", "bs_diary_14", "bs_diary_16"),
                           labels = c("down","dread nrest of day","worry",
                                      "inadequate", "tired", "content",
                                      "agitated", "irritated", "switch and focus",
                                      "extremely well", "full of ideas", "thoughts are racing"))) %>%
  group_by(rep, patient_id, variable) %>%
  summarise(mean_value = mean(value)) %>%
  group_by(rep, variable) %>%
  arrange(desc(mean_value)) %>%
  mutate(id = row_number()) %>%
  group_by(variable, rep) %>%
  arrange(id) %>%
  group_by(variable, id) %>%
  summarise(map_value = median(mean_value),
            cci_lwr = quantile(mean_value, 0.025),
            cci_upr = quantile(mean_value, 0.975)) %>%
  ggplot(aes(x = factor(id), y = map_value )) +
  geom_pointrange(aes(ymin = cci_lwr, ymax = cci_upr), alpha = 1, size = 1, fatten = 1, shape = 3) +
  geom_point(data = true_data %>%
               mutate(variable = factor(variable,
                                        levels = c("bs_diary_5", "bs_diary_13", "bs_diary_22",
                                                   "bs_diary_15", "bs_diary_9", "bs_diary_10",
                                                   "bs_diary_7", "bs_diary_11", "bs_diary_17",
                                                   "bs_diary_8", "bs_diary_14", "bs_diary_16"),
                                        labels = c("down","dread nrest of day","worry",
                                                   "inadequate", "tired", "content",
                                                   "agitated", "irritated", "switch and focus",
                                                   "extremely well", "full of ideas", "thoughts are racing"))
               ) %>% 
               group_by(patient_id, variable) %>%
               summarise(mean_value = mean(value, na.rm = TRUE)) %>%
               group_by(variable) %>%
               arrange(desc(mean_value)) %>%
               mutate(id = row_number()), aes(x = factor(id), y = mean_value), colour = "dark red", shape = 4) +
  facet_wrap(variable~., nrow = 2) +
  coord_flip() +
  theme_minimal() +
  ylab(label = "Value") +
  xlab(label = "Patient no.")

p

ggsave(plot = p + theme(text = element_text(size = 12)),
       filename = "outputs/figures/patient_mean_ppc_12.pdf", width = 9, height = 7, units = "in", dpi = 300)
ggsave(plot = p + theme(text = element_text(size = 12)),
       filename = "outputs/figures/patient_mean_ppc_12.tiff", width = 9, height = 7, units = "in", dpi = 300)
ggsave(plot = p + theme(text = element_text(size = 14)),
       filename = "outputs/figures/patient_mean_ppc_14.pdf", width = 9, height = 7, units = "in", dpi = 300)
ggsave(plot = p + theme(text = element_text(size = 14)),
       filename = "outputs/figures/patient_mean_ppc_14.tiff", width = 9, height = 7, units = "in", dpi = 300)
ggsave(plot = p + theme(text = element_text(size = 16)),
       filename = "outputs/figures/patient_mean_ppc_16.pdf", width = 9, height = 7, units = "in", dpi = 300)
ggsave(plot = p + theme(text = element_text(size = 16)),
       filename = "outputs/figures/patient_mean_ppc_16.tiff", width = 9, height = 7, units = "in", dpi = 300)
ggsave(plot = p + theme(text = element_text(size = 18)),
       filename = "outputs/figures/patient_mean_ppc_18.pdf", width = 9, height = 7, units = "in", dpi = 300)
ggsave(plot = p + theme(text = element_text(size = 18)),
       filename = "outputs/figures/patient_mean_ppc_18.tiff", width = 9, height = 7, units = "in", dpi = 300)
ggsave(plot = p + theme(text = element_text(size = 20)),
       filename = "outputs/figures/patient_mean_ppc_20.pdf", width = 9, height = 7, units = "in", dpi = 300)
ggsave(plot = p + theme(text = element_text(size = 20)),
       filename = "outputs/figures/patient_mean_ppc_20.tiff", width = 9, height = 7, units = "in", dpi = 300)



# Patient sd
ppc_data %>%
  mutate(variable = factor(variable, levels = c("bs_diary_5", "bs_diary_13", "bs_diary_22",
                                                "bs_diary_15", "bs_diary_9", "bs_diary_10",
                                                "bs_diary_7", "bs_diary_11", "bs_diary_17",
                                                "bs_diary_8", "bs_diary_14", "bs_diary_16"),
                           labels = c("down","dread nrest of day","worry",
                                      "inadequate", "tired", "content",
                                      "agitated", "irritated", "switch and focus",
                                      "extremely well", "full of ideas", "thoughts are racing"))) %>%
  group_by(rep, patient_id, variable) %>%
  summarise(sd_value = sd(value)) %>%
  group_by(rep, variable) %>%
  arrange(desc(sd_value)) %>%
  mutate(id = row_number()) %>%
  group_by(variable, rep) %>%
  arrange(id) %>%
  group_by(variable, id) %>%
  summarise(map_value = median(sd_value),
            cci_lwr = quantile(sd_value, 0.025),
            cci_upr = quantile(sd_value, 0.975)) %>%
  ggplot(aes(x = factor(id), y = map_value )) +
  geom_pointrange(aes(ymin = cci_lwr, ymax = cci_upr), alpha = 1, size = 1, fatten = 1, shape = 3) +
  geom_point(data = true_data %>%
               mutate(variable = factor(variable,
                                        levels = c("bs_diary_5", "bs_diary_13", "bs_diary_22",
                                                   "bs_diary_15", "bs_diary_9", "bs_diary_10",
                                                   "bs_diary_7", "bs_diary_11", "bs_diary_17",
                                                   "bs_diary_8", "bs_diary_14", "bs_diary_16"),
                                        labels = c("down","dread nrest of day","worry",
                                                   "inadequate", "tired", "content",
                                                   "agitated", "irritated", "switch and focus",
                                                   "extremely well", "full of ideas", "thoughts are racing"))
               ) %>% 
               group_by(patient_id, variable) %>%
               summarise(sd_value = sd(value, na.rm = TRUE)) %>%
               group_by(variable) %>%
               arrange(desc(sd_value)) %>%
               mutate(id = row_number()), aes(x = factor(id), y = sd_value), colour = "dark red", shape = 4) +
  facet_wrap(variable~., nrow = 2) +
  coord_flip() +
  theme_minimal() +
  ylab(label = "Value") +
  xlab(label = "Patient no.")




