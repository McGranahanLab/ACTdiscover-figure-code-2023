suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(cowplot))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(lubridate))
suppressPackageStartupMessages(library(grid))
suppressPackageStartupMessages(library(gtable))
suppressPackageStartupMessages(library(gridExtra))
suppressPackageStartupMessages(library(scales))
suppressPackageStartupMessages(library(reshape2))

saveDir      <- getwd()

setwd(saveDir)
plotDir      <- paste0(saveDir, "/plots/")
if (!file.exists(plotDir)) dir.create(plotDir, recursive = TRUE)


cfDNA.overview <- read.table("data/ctDNA.sampleOverview.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
load("data/mutTablePASS.pancreas.20200604properMutFilter.newITHState.newClonality.rda")
load("data/all.snv.23samples.rpdx_full.20200604.RData")
load("data/clonalMutsOverTime.samplePositivity.rda")

mutTablePASS <- subset(mutTablePASS, !mutTablePASS$SampleID %in% c("PAN106", "PAN119", "PAN121"))
mutTablePASS <- mutTablePASS %>% 
    mutate(keep = ifelse(!SampleID %in% c("PAN116", "PAN120"), 
                         "keep",
                         ifelse(grepl("Met:TRUE", is.present),
                                "keep",
                                "remove"))) %>%
    filter(keep == "keep") %>%
    select(-keep)

allPatients  <- unique(mutTablePASS$SampleID)

foo <- lapply(allPatients, function(patient) {
    mutTable <- subset(mutTablePASS, mutTablePASS$SampleID == patient)
    if (nrow(mutTable) == 0) return(NULL)

    regionsOfInterest <- paste0(patient, "_", grep("cfDNA|Met|T1", as.character(sapply(strsplit(mutTable$RegionSum[1], split = ";")[[1]], function(x) as.character(sapply(strsplit(x, split = ":"), function(y) y[1])))), value = TRUE))

    ################################## Getting muttable and snvtables ##################################

    mutTable$key <- gsub("PAN[0-9VH]{3}:", "", mutTable$mutation_id)
    rownames(mutTable) <- mutTable$key

    ### getting number of mutations called using other regions as leverage and the updated filters
    snvTable <- all.snv[[paste0("S_", patient)]]
    snvTable <- subset(snvTable, snvTable$Use.For.Plots == TRUE | snvTable$Use.For.Plots.Indel == TRUE)
    snvTable <- subset(snvTable, snvTable$max.var_count >= 5 | (snvTable$max.var_count < 5 & snvTable$MuTect.any & snvTable$Varscan.any))
    cfDNAmuts <- lapply(regionsOfInterest, function(reg) {
        snv.cfDNA <- snvTable[, c("key", grep(reg, colnames(snvTable), value = TRUE))]
        snv.cfDNA <- subset(snv.cfDNA, snv.cfDNA[, paste0(reg, ".VAF")] > 1.5 &
                                       snv.cfDNA[, paste0(reg, ".var_count")] >= 3)
        snv.cfDNA <- snv.cfDNA %>% 
                        mutate(key = gsub("Y", 24, gsub("X", 23, key))) %>% 
                        mutate(key = gsub(":[-ACGT]+$", "", gsub("chr", "", key))) %>%
                        rename(VAF = grep("VAF", colnames(.), value = TRUE))
        return(snv.cfDNA)
    })
    names(cfDNAmuts) <- regionsOfInterest

    cfDNAmuts.strict <- lapply(regionsOfInterest, function(reg) {
        snv.cfDNA <- snvTable[, c("key", grep(reg, colnames(snvTable), value = TRUE))]
        snv.cfDNA <- subset(snv.cfDNA, snv.cfDNA[, paste0(reg, ".var_count")] >= 5 | (snv.cfDNA[, paste0(reg, ".var_count")] < 5 & snv.cfDNA[, paste0("MuTect.", reg)] & snv.cfDNA[, paste0("Varscan.", reg)]))
        snv.cfDNA <- snv.cfDNA %>% 
                        mutate(key = gsub("Y", 24, gsub("X", 23, key))) %>% 
                        mutate(key = gsub(":[-ACGT]+$", "", gsub("chr", "", key))) %>%
                        rename(VAF = grep("VAF", colnames(.), value = TRUE))
        return(snv.cfDNA)
    })
    names(cfDNAmuts.strict) <- regionsOfInterest

    cfDNA.overview <- subset(cfDNA.overview, cfDNA.overview$Patient == patient)
    cfDNA.overview$nClonalMuts <- sapply(1:nrow(cfDNA.overview), function(x) {
        tmp.key <- cfDNAmuts[[cfDNA.overview$Sample[x]]]$key
        length(which(subset(mutTable, mutTable$key %in% tmp.key)$PyCloneClonal == "C"))
    })
    cfDNA.overview$origCalled <- sapply(1:nrow(cfDNA.overview), function(x) {
        tmp.key <- cfDNAmuts.strict[[cfDNA.overview$Sample[x]]]$key
        length(tmp.key)
    })
    full.cfDNA.overview <- cfDNA.overview %>% 
        mutate(ctDNA2 = ifelse(nClonalMuts > 0, TRUE, FALSE)) %>% 
        mutate(ctDNApos = ifelse(ctDNA & ctDNA2, TRUE, FALSE)) %>%
        select(-ctDNA, -ctDNA2)
    return(full.cfDNA.overview)
})
bar <- Reduce(rbind, foo)


load(file = "data/allPhasingTables.df.rda")

df.cnas <- allPhasingTables %>%
    group_by(patient, region, type) %>%
    tally() %>%
    ungroup() %>%
    filter(grepl("cfDNA", region)) %>%
    complete(nesting(patient, region), type, fill = list(n = 0)) %>%
    mutate(Patient = gsub("^[A-Z]_", "", patient)) %>%
    select(Patient, Sample = region, type, nCNAs = n) %>%
    group_by(Sample) %>%
    mutate(totalCNAs = sum(nCNAs)) %>%
    ungroup()

df.muts <- bar %>% 
    filter(ctDNApos) %>%
    mutate(rescuedMuts = nMuts - origCalled) %>%
    select(Patient, Sample, rescuedMuts, origCalledMuts = origCalled) %>%
    pivot_longer(cols = ends_with("Muts"), values_to = "nMuts", names_to = "type") %>%
    mutate(nMuts = ifelse(nMuts >= 0, nMuts, 0)) %>%
    mutate(type = factor(type, levels = c("rescuedMuts", "origCalledMuts"))) %>%
    group_by(Sample) %>% 
    mutate(totalMuts = sum(nMuts)) %>% 
    ungroup()

df.tmp <- full_join(df.muts %>% select(Patient, Sample, totalMuts) %>% unique(), 
                df.cnas %>% select(Patient, Sample, totalCNAs) %>% unique(),
                by = c("Patient", "Sample")) %>%
    mutate(totalMuts = ifelse(is.na(totalMuts), 0, totalMuts),
           totalCNAs = ifelse(is.na(totalCNAs), 0, totalCNAs)) %>%
    group_by(Patient) %>%
    mutate(medianTotalMuts = median(totalMuts),
           medianTotalCNAs = median(totalCNAs)) %>%
    ungroup() %>%
    arrange(desc(medianTotalMuts), desc(medianTotalCNAs), desc(totalMuts), desc(totalCNAs))

df <- bind_rows(df.muts %>% mutate(class = "muts") %>% rename(nEvents = nMuts, totalEvents = totalMuts),
                df.cnas %>% mutate(class = "cnas") %>% rename(nEvents = nCNAs, totalEvents = totalCNAs)) %>%
    complete(nesting(Patient, Sample), nesting(class, type), fill = list(nEvents = 0, totalEvents = 0)) %>%
    mutate(Patient = factor(Patient, levels = unique(df.tmp$Patient)),
           Sample = factor(Sample, levels = unique(df.tmp$Sample))) %>%
    mutate(type = factor(type, levels = c("rescuedMuts", "origCalledMuts", "SCNA rescued", "called")))


p1 <- ggplot(df %>% filter(class == "muts"), aes(x = Sample, y = nEvents, fill = type)) +
    geom_bar(stat = "identity", position = position_stack()) +
    scale_fill_manual(values = c("darkred", "steelblue"), labels = c("Rescued SNVs", "Called SNVs"), name = "") +
    ylab("# somatic mutations") +
    theme_cowplot(font_size = 12) +
    facet_grid(~ Patient, scales = "free_x", space = "free_x") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    theme(strip.background = element_blank())


p2 <- ggplot(df %>% filter(class == "cnas"), aes(x = Sample, y = nEvents, fill = type)) + 
    geom_bar(stat = "identity", position = position_stack()) +
    scale_fill_manual(values = c("darkred", "steelblue"), labels = c("Rescued SCNAs", "Called SCNAs"), name = "") +
    ylab("# somatic copy number aberrations") +
    theme_cowplot(font_size = 12) +
    facet_grid(~ Patient, scales = "free_x", space = "free_x") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    theme(strip.text = element_blank(), strip.background = element_blank())

### set margins
g.p1 <- ggplotGrob(p1 + theme(plot.margin = unit(c(0.25, 0.25, 0, 0.25), "cm")) + theme(axis.title.x = element_blank(), axis.text.x = element_blank()))
g.p2 <- ggplotGrob(p2 + theme(plot.margin = unit(c(0, 0.25, 0.25, 0.25), "cm")) + theme(axis.title.x = element_blank()))

### align widths
plot_list  <- list(g.p1, 
                   g.p2)
all_widths <- lapply(plot_list, function(x) x$widths)
plot_list_alignedWidths <- lapply(plot_list, function(x) {
    x$widths <- do.call(grid::unit.pmax, all_widths)
    return(x)
})

full_plot <- ggdraw() +
    draw_plot(plot_list_alignedWidths[[2]], x = 0, y = 0, width = 1, height = 0.6) +
    draw_plot(plot_list_alignedWidths[[1]], x = 0, y = 0.6, width = 1, height = 0.4)

pdf("plots/combinedRescuingPlot.pdf", width = 14, height = 7, useDingbats = FALSE)
full_plot
dev.off()



### checking fraction only called with method abstract
df %>% filter(class == "cnas", totalEvents > 0, type == "SCNA rescued") %>% summarize(totalNew = sum(nEvents), total = sum(totalEvents))

### summary of mutations
summary(df %>% filter(class == "muts", totalEvents > 0) %>% select(Sample, totalEvents) %>% unique() %>% pull(totalEvents))
summary(df %>% filter(class == "muts", totalEvents > 0, type == "rescuedMuts") %>% mutate(prop = nEvents / totalEvents) %>% pull(prop))

### summary of SCNAs
summary(df %>% filter(class == "cnas", totalEvents > 0) %>% select(Sample, totalEvents) %>% unique() %>% pull(totalEvents))
summary(df %>% filter(class == "cnas", totalEvents > 0, type == "SCNA rescued") %>% mutate(prop = nEvents / totalEvents) %>% pull(prop))
