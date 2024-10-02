# This script scores the ChatGPT annotations against human annotated threats
# from the same U.S. Endangered Species Act documentation.

#### SETUP ####
# Load libraries.
library(data.table)
library(tidyverse)

`%!in%` <- negate(`%in%`)

##### Custom Functions #####
# calc.acc - calculates the accuracy, precision, recall, and 
# F1-score of the ChatGPT annotations.
calc.metrics <- function(paired_ann){
  
  trueTrue <- sum(paired_ann$gpt_ann == TRUE & paired_ann$human_ann == TRUE, na.rm = TRUE)
  trueFalse <- sum(paired_ann$gpt_ann == FALSE & paired_ann$human_ann == TRUE, na.rm = TRUE)
  falseTrue <- sum(paired_ann$gpt_ann == TRUE & paired_ann$human_ann == FALSE, na.rm = TRUE)
  falseFalse <- sum(paired_ann$gpt_ann == FALSE & paired_ann$human_ann == FALSE, na.rm = TRUE)
  
  acc <- (trueTrue + falseFalse)/(sum(c(trueTrue, trueFalse, falseTrue, falseFalse)))
  prec <- (trueTrue)/(trueTrue + falseTrue)
  recl <- (trueTrue)/(trueTrue + trueFalse)
  f1 <- (2 * prec * recl)/(prec + recl)
  
  return(c(accuracy = acc, precision = prec, recall = recl, F1 = f1,
           trueTrueRate = trueTrue/(sum(c(trueTrue, trueFalse))),
           trueFalseRate = trueFalse/(sum(c(trueTrue, trueFalse))),
           falseTrueRate = falseTrue/(sum(c(falseTrue, falseFalse))),
           falseFalseRate = falseFalse/(sum(c(falseTrue, falseFalse)))))
}

# Load the human annotations
human_ann <- data.table::fread("../data/FINAL_HUMAN_ANN.csv") %>%
  dplyr::mutate(documentKey = paste0(TaxonID, documentType, urls)) %>%
  dplyr::mutate(across(8:20,
                       ~as.logical(ifelse(. %!in% c(TRUE, FALSE), NA, .))))

# Load the ChatGPT annotations
gpt_ann <- data.table::fread("../data/test_atomized_threats_test.csv") %>%
  dplyr::mutate(documentKey = paste0(TaxonID, documentType, urls)) %>%
  dplyr::mutate(`1 Residential & commercial development` <- ifelse(`1.1 Housing & urban areas`,
                                                                   TRUE, ifelse(
                                                                     `1.2 Commercial & industrial areas`, TRUE,
                                                                     ifelse(`1.3 Tourism & recreation areas`, TRUE, FALSE))))

#### COMPARISONS ####
# Filter the ChatGPT data to those that were also annotated by humans.
gpt_ann <- gpt_ann %>%
  dplyr::filter(documentKey %in% human_ann$documentKey) %>%
  dplyr::filter(duplicated(documentKey) == FALSE) %>%
  unique()
rownames(gpt_ann) <- gpt_ann$documentKey

human_ann <- human_ann %>%
  dplyr::filter(documentKey %in% gpt_ann$documentKey) %>%
  dplyr::filter(duplicated(documentKey) == FALSE) %>%
  unique()
rownames(human_ann) <- human_ann$documentKey

human_ann <- human_ann[match(gpt_ann$documentKey, human_ann$documentKey),]

my_resp <- data.table::fread("../data/subsample_gpt_responses_parseReady.csv") %>%
  dplyr::mutate(documentKey = paste0(TaxonID, documentType, urls)) %>%
  dplyr::right_join(human_ann, by = "documentKey")
write.csv(my_resp, "../data/modelData.csv")

true_freq <- my_resp %>%
  dplyr::summarise(across(where(is.logical), ~sum(.x, na.rm = TRUE))) %>%
  c() %>% unlist()
true_freq <- true_freq[3:length(true_freq)]

##### Residential & Commercial Development #####
t1_df <- data.frame(human_ann = human_ann$`Residential & Commerial Development`,
                    gpt_ann = gpt_ann$`1 Residential & commercial development`)

##### Agriculture and Aquaculture #####
t2_df <- data.frame(human_ann = human_ann$`Agriculture & Aquaculture`,
                    gpt_ann = gpt_ann$`2 Agriculture & aquaculture`)

##### Energy Production and Mining #####
t3_df <- data.frame(human_ann = human_ann$`Energy Producting & Mining`,
                    gpt_ann = gpt_ann$`3 Energy production & mining`)

##### Transportation and Service Corridors #####
t4_df <- data.frame(human_ann = human_ann$`Transportation & Service Corridors`,
                    gpt_ann = gpt_ann$`4 Transportation & service corridors`)

##### Biological Resource Use #####
t5_df <- data.frame(human_ann = human_ann$`Biological Resource Use`,
                    gpt_ann = gpt_ann$`5 Biological resource use`)

##### Human Intrusions and Disturbance #####
t6_df <- data.frame(human_ann = human_ann$`Human Intrusions & Disturbance`,
                    gpt_ann = gpt_ann$`6 Human intrusions & disturbance`)

##### Natural Systems Modifications #####
t7_df <- data.frame(human_ann = human_ann$`Natural System Modifcations`,
                    gpt_ann = gpt_ann$`7 Natural system modifications`)

##### Invasive and Other Problematic Species/Genes #####
t8_df <- data.frame(human_ann = human_ann$`Invasive & Other Problematic Species & Genes`,
                    gpt_ann = gpt_ann$`8 Invasive & other problematic species, genes & diseases`)

##### Pollution #####
t9_df <- data.frame(human_ann = human_ann$Pollution,
                    gpt_ann = gpt_ann$`9. Pollution`)

##### Geological Events #####
t10_df <- data.frame(human_ann = human_ann$`Geological Events`,
                     gpt_ann = gpt_ann$`10 Geological events`)

##### Climate Change and Severe Weather #####
t11_df <- data.frame(human_ann = human_ann$`Climate Change & Severe Weather`,
                     gpt_ann = gpt_ann$`11 Climate change & severe weather`)

#### DATASET AVERAGES ####
scores_df <- sapply(list(t1_df, t2_df, t3_df, t4_df, t5_df,
                         t6_df, t7_df, t8_df, t9_df, t10_df, t11_df),
                    FUN = calc.metrics)
colnames(scores_df) <- c("Residential Development", "Agricultural Development",
                         "Energy Production/Mining", "Transportation", "Biological Resource Use",
                         "Human Intrusions", "Natural Systems Modifications", "Problematic Species",
                         "Pollution", "Geological Events", "Climate Change")

confusion_df <- scores_df[5:8,] %>% t() %>%
  as.data.frame()
confusion_df$threatCategory <- c("Residential Development", "Agriculture", 
                                 "Energy and Mining", "Transportation",
                                 "Biological Resource Use", "Human Intrusions",
                                 "Natural System Mod.", "Problematic Species",
                                 "Pollution", "Geological Events", "Climate Change")

confusion_df <- confusion_df %>%
  tidyr::pivot_longer(-threatCategory, names_to = "outcome", values_to = "metric") %>%
  dplyr::mutate(outcomeCoordX = ifelse(outcome == "trueTrueRate", 1,
                                       ifelse(outcome == "trueFalseRate", 2,
                                              ifelse(outcome == "falseTrueRate", 1, 2))),
                outcomeCoordY = ifelse(outcome == "trueTrueRate", 2,
                                       ifelse(outcome == "trueFalseRate", 2,
                                              ifelse(outcome == "falseTrueRate", 1, 1))))

ggplot(confusion_df)+
  geom_tile(mapping = aes(x = outcomeCoordX, y = outcomeCoordY,
                          fill = metric))+
  geom_label(mapping = aes(x = outcomeCoordX, y = outcomeCoordY,
                           label = round(metric, 2)))+
  colorspace::scale_fill_continuous_divergingx(mid = 0.5, rev = TRUE,
                                               palette = "Zissou 1")+
  scale_x_continuous(name = "CHATGPT PREDICTION",
                     limits = c(0,3), breaks = c(1,2),
                     labels = c("Positive", "Negative"))+
  scale_y_continuous(name = "HUMAN PREDICTION",
                     limits = c(0,3), breaks = c(1,2),
                     labels = c("Negative", "Positive"))+
  facet_wrap(~confusion_df$threatCategory)+
  cowplot::theme_cowplot()+
  theme(legend.position = "none")
ggsave2("../figures/SupplementalFigure_S1.png", dpi = 400,
        height = 10, width = 10)

##### NAIVE BAYES F1 SCORES #####
nb_scores <- data.table::fread("../data/f1_scores.csv") %>%
  pivot_longer(-V1, values_to = "f1", names_to = "cat") %>%
  dplyr::mutate(cat = case_when(
    cat == "t1_f1" ~ "1. Residential & Commercial Dev.",
    cat == "t2_f1" ~ "2. Agriculture & Aquaculture",
    cat == "t3_f1" ~ "3. Energy Prod. & Mining",
    cat == "t4_f1" ~ "4. Transportation & Service Corr.",
    cat == "t5_f1" ~ "5. Biological Resource Use",
    cat == "t6_f1" ~ "6. Human Intrusions & Disturbance",
    cat == "t7_f1" ~ "7. Natural System Modifications",
    cat == "t8_f1" ~ "8. Problematic Biotic Factors",
    cat == "t9_f1" ~ "9. Pollution",
    cat == "t10_f1" ~ "10. Geological Events",
    cat == "t11_f1" ~ "11. Climate Change & Severe Weather",
    cat == "t12_f1" ~ "Misc. Instrinsic Pop. Factors",
    cat == "t13_f1" ~ "Misc. Disease"
  ))
nb_scores$cat <- factor(nb_scores$cat,
                        levels = c("1. Residential & Commercial Dev.",
                                   "2. Agriculture & Aquaculture",
                                   "3. Energy Prod. & Mining",
                                   "4. Transportation & Service Corr.",
                                   "5. Biological Resource Use",
                                   "6. Human Intrusions & Disturbance",
                                   "7. Natural System Modifications",
                                   "8. Problematic Biotic Factors",
                                   "9. Pollution",
                                   "10. Geological Events",
                                   "11. Climate Change & Severe Weather",
                                   "Misc. Instrinsic Pop. Factors",
                                   "Misc. Disease"))

nb_mean <- nb_scores %>%
  dplyr::group_by(cat) %>%
  dplyr::summarise(meanScore = mean(f1))
nb_mean$trueRate <- true_freq

nb_thresh_freq <- nb_scores %>%
  dplyr::group_by(cat) %>%
  dplyr::summarise(numThresh = sum(f1 >= 0.8))

ggplot()+
  geom_point(nb_mean,
             mapping = aes(x = trueRate, y = meanScore),
             size = 3)+
  geom_smooth(nb_mean,
              mapping = aes(x = trueRate, y = meanScore),
              method = "lm")+
  labs(x = "Number of Positive Hits", y = "Mean F1-Score")+
  cowplot::theme_cowplot()+
  theme(axis.text = element_text(size = 16),
        axis.title = element_text(size = 20),
        plot.background = element_rect(fill = "white", color = "white"),
        legend.position = "none")

ggplot(nb_scores)+
  geom_density(mapping = aes(x = f1, fill = cat))+
  geom_boxplot(mapping = aes(x = f1, y = 1))+
  geom_vline(xintercept = 0.8, linetype = 2)+
  facet_wrap(~cat)+
  scale_x_continuous(limits = c(0,1),
                     labels = c(0, 0.25, 0.5, 0.75, 1),
                     name = "F1-Score")+
  scale_y_continuous(name = "Density")+
  cowplot::theme_cowplot()+
  theme(plot.background = element_rect(fill = "white", color = "white"),
        legend.position = "none")
ggsave2("../figures/SupplementalFigure_S2.png", dpi = 400,
        height = 10, width = 10)




