### Purpose: plot COPRO seq result
### Created: 2024-09-10

# load required library
library(dplyr)
library(tidyr)
library(ggsci)
library(ggplot2)
library(gridExtra)
library(data.table)

setwd("/Volumes/ReyLab_QZ10T/COPROseq/test_BPB/")


# load coproseq result
mapped_reads <- fread("coproseq_out/mapped_read_count.csv", data.table = F)
mapped_stats <- fread("coproseq_out/mapped_read_stats.csv", data.table = F)
rel_abundance <- fread("coproseq_out/relative_abundance.csv", data.table = F)

# barplot
p1 <- mapped_reads %>%
    left_join(mapped_stats %>% select(sample_id, total_read, mapped_read), by = "sample_id") %>%
    mutate(unmapped = total_read - mapped_read) %>%
    select(-total_read, -mapped_read) %>%
    pivot_longer(cols = Anaerobutyricum_soehngenii:unmapped, names_to = "mapped_reads", values_to = "value") %>%
    mutate(sample_id = factor(sample_id)) %>%
    mutate(mapped_reads = factor(mapped_reads, levels = rev(unique(mapped_reads)))) %>%
    ggplot(aes(x = sample_id, y = value, fill = mapped_reads)) + 
    geom_bar(position = "fill", stat = "identity") + 
    theme_classic() +
    theme(axis.title = element_text(size = 15),
          axis.text = element_text(size = 12),
          legend.title = element_text(size = 15),
          legend.text = element_text(size = 15),
          plot.title = element_text(size = 15),
          strip.text.x = element_text(size = 15),
          axis.text.x = element_text(angle = 270, hjust = 0)) +
    scale_y_continuous(expand=c(0,0)) +
    scale_fill_manual(values = c("grey75", pal_npg("nrc", alpha = 0.7)(9))) +
    xlab("sample ID") +
    ylab("proportion of reads mapped") +
    ggtitle("Mapping stats")


p2 <- rel_abundance %>%
    pivot_longer(cols = Anaerobutyricum_soehngenii:Roseburia_intestinalis, names_to = "mapped_reads", values_to = "value") %>%
    mutate(sample_id = factor(sample_id)) %>%
    mutate(mapped_reads = factor(mapped_reads, levels = rev(unique(mapped_reads)))) %>%
    ggplot(aes(x = sample_id, y = value, fill = mapped_reads)) + 
    geom_bar(position = "fill", stat = "identity") + 
    theme_classic() +
    theme(axis.title = element_text(size = 15),
          axis.text = element_text(size = 12),
          legend.title = element_text(size = 15),
          legend.text = element_text(size = 15),
          plot.title = element_text(size = 15),
          strip.text.x = element_text(size = 15),
          axis.text.x = element_text(angle = 270, hjust = 0)) +
    scale_y_continuous(expand=c(0,0)) +
    scale_fill_manual(values = c(pal_npg("nrc", alpha = 0.7)(9))) +
    xlab("sample ID") +
    ylab("proportion of reads mapped") +
    ggtitle("Relative abundance")


grid.arrange(p1, p2, nrow=2)
