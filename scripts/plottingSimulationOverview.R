library(tidyverse)
library(cowplot)
library(ACTdiscover)

    # load("data/20230123.allSimulations.rda")
    ### the data frame is very large, so a subset has been provided with the simulations up to 50000 SNPs
    load("data/20230123.allSimulations.subset.rda")

    phasingTable <- allSimulations %>% 
        group_by(nrSNPs, purity, depth, cpnA, cpnB) %>% 
        tally() %>% 
        ungroup() %>% 
        mutate(region = paste0("Sim", gsub("\\.", "", purity))) %>% 
        mutate(chrom = paste0("chr", as.numeric(factor(nrSNPs)))) %>% 
        select(-n) %>% 
        left_join(allSimulations, by = c("nrSNPs", "purity", "depth", "cpnA", "cpnB")) %>% 
        group_by(region, chrom) %>% 
        mutate(pos = row_number()) %>%
        ungroup() %>%
        mutate(segment = as.numeric(factor(nrSNPs))) %>%
        select(chrom, pos, class, region, BAF, segment)

    newSegTable <- phasingTable %>%
        group_by(region, segment) %>%
        mutate(start = min(pos), end = max(pos), width = end - start) %>%
        filter(class == "bottom") %>%
        mutate(meanBAF = mean(BAF)) %>%
        ungroup() %>%
        group_by(region, chrom, start, end, width, segment, meanBAF) %>%
        tally() %>%
        ungroup() %>%
        select(region, seqnames = chrom, start, end, width, segment, meanBAF) %>%
        mutate(BAF_segmented = ifelse(region == "Sim1", meanBAF, 0.5)) %>%
        select(-meanBAF)

    newPhasingTable <- testingSCNA(phasingTable, newSegTable)

    tmp <- newPhasingTable %>% 
        group_by(chrom, region, segment, type) %>% 
        tally() %>% 
        ungroup() %>% 
        select(-n) %>%
        mutate(type = factor(type, levels = c("not detected", "called", "SCNA rescued"))) %>%
        mutate(region = factor(region, levels = c("Sim0", "Sim1e-04", "Sim5e-04", "Sim0001", "Sim0005", "Sim001", "Sim005", "Sim01", "Sim1"))) %>%
        mutate(chrom = factor(chrom, levels = paste0("chr", 1:length(unique(chrom)))))

    q1 <- ggplot(tmp, aes(region, chrom, fill = type)) + 
        geom_tile(color = "gray") +
        scale_x_discrete(name = "purity", labels = c("0%", "0.01%", "0.05%", "0.1%", "0.5%", "1%", "5%", "10%", "100%")) +
        scale_y_discrete(name = "number of SNPs", labels = c(10, 50, 100, 500, 1000, 5000, 10000, 50000, 100000, 500000, 1000000, 5000000)) +
        scale_fill_manual(name = "Classification", 
                          values = setNames(c("#D95980", "#335120", "#9DCC5F", "#63AAC0", "#bdbdbd"), c("MSAI", "SCNA rescued", "called", "not detected", "not called")),
                          labels = setNames(c("AI not detected", "AI in single region", "AI rescued"), c("not detected", "called", "SCNA rescued"))) +
        theme_cowplot(font_size = 12)

    pdf("plots/20230123.simulationsOverview.pdf")
    q1
    dev.off()

