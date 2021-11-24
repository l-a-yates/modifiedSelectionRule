#------------------------------------------------------------
#
# Functions definitions for "Parsimonious Model Selection using 
#    Information Theory: a Modified Selection Rule"
# Authors: Yates, L.A, Richards, S.A, Brook, B.W.
# 
# This script contains function definitions to:
#   - Apply the original and modified one-standard-error rules;
#   - Produce (modified) standard-error plots.
# 
# Date: Original version Dec 2019; this version Sept 2020
#
# Please see manuscript for further details and references
#------------------------------------------------------------

library(tidyverse)
library(ggpubr)

select <- dplyr::select # override MASS::select 


# applies the modified one-standard-error rule
# input: a summary table as a tibble (see main code for format); 
# if multiple = T, then all models satisfying the modified selection condition are returned
# if multiple = F, only the lowest-scoring of the selected models is returned
# output: a tibble with one row per model and columns for model dimension (dim) and model name (model)
modOSErule <- function(summaryTable, multiple = F){
  dimLowest = filter(summaryTable, score == min(score)) %>% pull(dim)
  summaryTable %>%  
    filter(delScore <= summaryTable[["se_mod"]], dim <= dimLowest) %>% 
    filter(dim == min(dim)) %>% 
    {if(!multiple) filter(.,delScore == min(delScore)) else .} %>% 
    select(dim,model)
}

# applies the original one-standard-error rule
# input: a summary table as a tibble (see main code for format)
# returns a one-row tibble with columns for model dimension (dim) and model name (model)
OSErule <- function(summaryTable){
  dimLowest = filter(summaryTable, score == min(score)) %>% pull(dim)
  oseLowest = filter(summaryTable, score == min(score)) %>% pull(se_ose)
  summaryTable %>%  
    filter(delScore <= oseLowest, dim <= dimLowest) %>% 
    filter(dim == min(dim)) %>% 
    filter(delScore == min(delScore)) %>% 
    select(dim,model)
}


# plots model scores, standard errors and displays the ordinary one-standard-error rule
# input: a summary table as a tibble (see main code for format)
# output: ggplot object
plotOSE <- function(summaryTable){
  summaryTable %>% mutate(SD = se_ose, Score = score) %>% 
    ggplot(aes(x = reorder(model,index), y = Score)) + 
    geom_errorbar(aes(ymin = Score - SD, ymax = Score + SD), col = "gray80")  +
    geom_errorbar(aes(ymin = Score - SD, ymax = Score + SD), data = ~ filter(.x, Score == min(Score))) +
    geom_point() +
    geom_point(aes(colour = "Lowest score"), shape =1, size = 4, data = ~ filter(.x, Score == min(Score))) +
    geom_point(aes(colour = "One-standard-error rule"), shape =1, size = 4, data = ~filter(.x,model == OSErule(summaryTable)$model)) +
    geom_blank(aes(colour = "Modified selection rule")) +
    scale_color_manual(name = NULL, values=c("Lowest score" = "blue","One-standard-error rule" = "grey40","Modified selection rule" ="red")) +
    labs(x = NULL, subtitle = "Standard errors", y = "Score") +
    theme(panel.border = element_rect(fill = NA), 
          panel.background = element_blank(), 
          legend.key=element_blank(),
          axis.ticks = element_blank(),
          axis.text.x = element_text(angle=0, hjust= 0.5, vjust = 0.5, size = 8))
}

# plots model scores and adjusted standard errors, and displays the modified one-standard-error rule
# input: a summary table as a tibble (see main code for format); see modOSErule() for 'multiple'
# output: ggplot object
plotModifiedOSE <- function(summaryTable, multiple = F){
  summaryTable %>% mutate(SD = se_mod, Score = delScore) %>% 
    ggplot(aes(x = reorder(model,index), y = Score)) + 
    geom_hline(aes(yintercept = 0), col = "gray30", linetype = "dashed", size = 0.5) +
    geom_errorbar(aes(ymin = Score - SD, ymax = Score + SD), col = "gray80")  +
    geom_errorbar(aes(ymin = Score - SD, ymax = Score + SD), data = ~filter(.x, model %in% modOSErule(summaryTable, multiple)$model)) +
    geom_point() +
    geom_point(aes(colour = "Lowest score"), shape =1, size = 4, data = ~filter(.x, Score == min(Score))) +
    geom_point(aes(colour = "One-standard-error rule"), shape =1, size = 4, data = ~filter(.x,model == OSErule(summaryTable)$model)) +
    geom_point(aes(colour = "Modified selection rule"), shape =1, size = 4, data = ~filter(.x,model %in% modOSErule(summaryTable, F)$model)) +
    scale_color_manual(name = NULL, values=c("Lowest score" = "blue","One-standard-error rule" = "grey40","Modified selection rule" ="red")) +
    labs(x = NULL, subtitle = "Correlation-adjusted errors", y = "Relative score") +
    theme(panel.border = element_rect(fill = NA), 
          panel.background = element_blank(), 
          legend.key=element_blank(),
          axis.ticks = element_blank(),
          axis.text.x = element_text(angle=0, hjust= 0.5, vjust = 0.5, size = 8))
} 

