####################################################################################################
#########################                    Libraries                    #########################
####################################################################################################

    library(tidyverse)
    library(ggnewscale)


####################################################################################################
#########################                      Paths                      #########################
####################################################################################################

    saveDir      <- getwd()

    setwd(saveDir)
    plotDir      <- paste0(saveDir, "/plots/")
    if (!file.exists(plotDir)) dir.create(plotDir, recursive = TRUE)


####################################################################################################
#########################                    Load data                    #########################
####################################################################################################

    clinicalData <- read.table("data/updatedClinData_fixingTreatment.csv", header = TRUE, sep = ",", stringsAsFactors = FALSE)
    allSamples   <- unique(read.table("data/designFile_full_new_merged103.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)$Label)
    load("data/clonalMutsOverTime.samplePositivity.rda")


####################################################################################################
#########################           Formatting samples and dates           #########################
####################################################################################################

    ### remove germline samples
    allSamples <- gsub("S_", "", grep("GL", allSamples, value = TRUE, invert = TRUE))
    sampleOverview  <- data.frame(patient = sapply(strsplit(allSamples, split = "_"), function(x) x[1]), 
                                  region = allSamples,
                                  stringsAsFactors = FALSE)
    patientOverview <- clinicalData %>%
        select(patient, os_time, cens_os, biopsy_time)
    patientOverview$pdx <- ifelse(patientOverview$patient %in% sapply(strsplit(grep("PDX", allSamples, value = TRUE), split = "_"), function(x) x[1]), patientOverview$biopsy_time, NA)

    cfDNAOverview <- df.clonalMutsOverTime %>% 
        filter(Label == "cfDNA") %>% 
        select(patient, region, TumPresence, meanVAF) %>%
        full_join(sampleOverview %>% filter(grepl("cfDNA", region)), by = c("patient", "region")) %>%
        replace_na(list(TumPresence = FALSE, meanVAF = 0.01)) %>%
        mutate(sample_time = as.numeric(gsub("cfDNA", "", sapply(strsplit(as.character(region), split = "_"), function(x) x[2])))) %>%
        select(patient, TumPresence, sample_time, meanVAF)

    clinicalData2 <- clinicalData %>% 
        select(patient, starts_with("start_time"), starts_with("finish_time")) %>% 
        pivot_longer(cols = starts_with("start_time") | starts_with("finish_time"), 
                     names_to = c(".value", "treatment"), 
                     names_pattern = "([a-z]+_time)\\.([0-9])") %>% 
        select(patient, treatment, start_time, finish_time)


####################################################################################################
#########################              Creating swimmes plot              #########################
####################################################################################################

    patientOverview <- patientOverview %>%
        arrange(desc(os_time)) %>% 
        mutate(patient = forcats::fct_inorder(patient))

    cfDNAOverview <- cfDNAOverview %>%
        mutate(patient = factor(patient, levels = levels(patientOverview$patient))) %>%
        arrange(patient)

    clinicalData2 <- clinicalData2 %>%
        mutate(patient = factor(patient, levels = levels(patientOverview$patient))) %>%
        arrange(patient)

    p.swimmers <- ggplot(patientOverview, aes(os_time, patient)) +
        geom_vline(aes(xintercept = 0), linetype = 2, alpha = 0.25, size = 0.25) +
        geom_segment(aes(x = 0, xend = os_time, y = patient, yend = patient), size = 0.1) +
        geom_segment(data = clinicalData2, aes(x = start_time, xend = finish_time, y = patient, yend = patient), size = 0.5, color = "darkgreen") +
        geom_tile(aes(fill = factor(cens_os)), width = 5) + 
        scale_fill_manual(name = "Event", 
                          values = c("0" = "#000000", "1" = "#969696"),
                          labels = c("0" = "Death", "1" = "Censored")) +
        geom_point(aes(biopsy_time, patient), shape = 21, size = 3, color = "orange", fill = "orange") + 
        # geom_point(aes(pdx, patient), shape = 21, size = 3, color = "red", fill = "red") +
        new_scale_fill() +
        geom_point(data = cfDNAOverview %>% filter(!TumPresence), aes(sample_time, patient, fill = TumPresence), shape = 21, size = 2, color = "#d9d9d9", fill = "#d9d9d9") +
        geom_point(data = cfDNAOverview %>% filter(TumPresence), aes(sample_time, patient, fill = meanVAF), shape = 21, size = 2, color = "black") +
        scale_fill_gradient(name = "ctDNA mean VAF", 
                            low = "#9ecae1", high = "#08306b") +
                            # low = "#A6CEE3", high = "#1F78B4") + 
        scale_y_discrete(drop = FALSE) +
        # scale_x_continuous(expand = c(0,0), limits = c(0, max(patientOverview$os_time) + 5)) +
        xlab("Days post diagnosis") + ylab("") +
        theme_classic() +
        theme(axis.line.y = element_blank())

    pdf("plots/swimmersPlot.pdf", width = 7, height = 5, useDingbats = FALSE)
        p.swimmers
    dev.off()
