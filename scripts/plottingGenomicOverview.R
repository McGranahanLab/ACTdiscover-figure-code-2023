####################################################################################################
######################### Summarizing patient data / mutations / clonality #########################
####################################################################################################

library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Hs.eg.db)
library(reshape2)
library(ggplot2)
library(cowplot)
library(ggpubr)
library(deconstructSigs)
library(RColorBrewer)
library(gtable)
library(grid)
library(gridExtra)
library(plyr)
library(dplyr)
library(scales)
library(lubridate)


g_legend <- function(a.gplot) {
    tmp    <- ggplot_gtable(ggplot_build(a.gplot))
    leg    <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
    if (length(leg) == 0) {
        return(NULL)
    }
    legend <- tmp$grobs[[leg]]
    return(legend)
}


saveDir      <- getwd()

setwd(saveDir)
plotDir      <- paste0(saveDir, "/plots/")
if (!file.exists(plotDir)) dir.create(plotDir, recursive = TRUE)

load("data/mutTablePASS.pancreas.20200604properMutFilter.newITHState.newClonality.rda")
load("data/all.snv.23samples.rpdx_full.20200604.RData")
load("data/all.ASCAT.segs.rda")

mutTablePASS <- subset(mutTablePASS, !mutTablePASS$SampleID %in% c("PAN106", "PAN119", "PAN121"))
mutTablePASS <- mutTablePASS %>% 
    mutate(keep = ifelse(!SampleID %in% c("PAN116", "PAN120"), 
                         "keep",
                         ifelse(grepl("Met:TRUE", is.present),
                                "keep",
                                "remove"))) %>%
    filter(keep == "keep") %>%
    select(-keep)

############################ Collating data from mutTable and reshaping ############################

mutationLoad <- table(mutTablePASS$SampleID, mutTablePASS$PyCloneClonal, useNA = "always")
mutationLoad <- mutationLoad[-dim(mutationLoad)[1],]
mutationLoad <- data.frame(patient = rownames(mutationLoad), 
                           clonal = mutationLoad[, 1], 
                           subclonal = mutationLoad[, 2])#, Unknown = mutationLoad[, 3])
mutationLoadOrig <- mutationLoad


############################# Running deconstructSigs for all patients #############################

fullClonality <- paste(rep(unique(mutTablePASS$SampleID), each = 2), rep(c("C", "S"), length(unique(mutTablePASS$SampleID))), sep = ":")

mutTableCS          <- subset(mutTablePASS, mutTablePASS$PyCloneClonal %in% c("C", "S"))
mutations           <- mutTableCS[, c(1,4,5,7,8)]
mutations$Clonality <- paste(mutTableCS$SampleID, mutTableCS$PyCloneClonal, sep = ":")
mutations$start     <- as.numeric(mutations$start)
use.sigs <- c('SBS1',
              'SBS2', 
              'SBS3', 
              'SBS5', 
              'SBS6', 
              'SBS13')

mutationsoutPerPatient   <- mut.to.sigs.input(mutations, sample.id = "SampleID", pos = "start", alt = "var")
mutationSignaturesPat    <- lapply(unique(mutations$SampleID), function(x) whichSignatures(mutationsoutPerPatient, sample.id = x, contexts.needed = TRUE, tri.counts.method = 'exome2genome', associated = use.sigs, signatures.ref = signatures.genome.cosmic.v3.may2019))

mutationSignaturesAll  <- lapply(mutationSignaturesPat, function(x) x[["weights"]])
mutationSignaturesAll  <- Reduce(rbind, mutationSignaturesAll)
mutationSignaturesAll  <- mutationSignaturesAll[, which(colnames(mutationSignaturesAll) %in% use.sigs)]
mutSigIDs              <- Reduce(rbind, strsplit(rownames(mutationSignaturesAll), split = ":"))
rownames(mutSigIDs)    <- NULL
colnames(mutSigIDs)    <- c("patient")
mutationSignaturesAll  <- cbind(mutSigIDs, mutationSignaturesAll)
mutationSignaturesAll  <- cbind(mutationSignaturesAll, other = 1 - rowSums(mutationSignaturesAll[, -1]))


############################# Getting high confidence driver mutations #############################

highConfDriverTable <- subset(mutTablePASS, mutTablePASS$driverCategory %in% c("1", "1A", "2A") | 
                                            mutTablePASS$indelLungDriver == TRUE | 
                                            mutTablePASS$spliceDriver == TRUE)
driverTable <- data.frame(patient = highConfDriverTable$SampleID, 
                          gene = highConfDriverTable$Hugo_Symbol, 
                          stringsAsFactors = FALSE)
driverTable$patient <- factor(driverTable$patient, levels = unique(mutTablePASS$SampleID))
### TODO update gene list
driverTable$gene <- factor(driverTable$gene, levels = c("TP53", "KRAS", "SMAD4", "BRCA2", "PTEN", "HNF1A", "BAP1", "CREBBP", "KMT2D", "NOTCH2", "SMARCA4", "ASXL1", "POT1", "RNF43", "TSC2", "UBR5"))
driverTable <- as.data.frame.matrix(table(driverTable))

### getting the clonality of drivers
clonalityDrivers <- lapply(rownames(driverTable), function(pat) {
    tmpTable <- subset(highConfDriverTable, highConfDriverTable$SampleID == pat)
    sapply(colnames(driverTable), function(gene) {
        tmpClonal <- tmpTable$PyCloneClonal[which(tmpTable$Hugo_Symbol == gene)]
        if (length(tmpClonal) == 0) return("none")
        if (length(tmpClonal) > 1) return(unique(tmpClonal))
        if (is.na(tmpClonal)) return("Unknown")
        return(tmpClonal)
    })
})
clonalityDrivers <- Reduce(rbind, clonalityDrivers)
rownames(clonalityDrivers) <- rownames(driverTable)
clonalityDrivers <- data.frame(clonalityDrivers)

driverTable <- cbind(patient = rownames(driverTable), clonalityDrivers)

### salvage KRAS mutations from snvtables
krasSalvaged <- lapply(as.character(driverTable$patient[which(driverTable$KRAS == "none")]), function(pat) {
    snvtable <- all.snv[[grep(pat, names(all.snv))]]
    snvtable <- subset(snvtable, snvtable$driverCategory %in% c("1", "1A", "2A") | snvtable$indelLungDriver == TRUE)
    snvtable <- subset(snvtable, snvtable$Gene.refGene == "KRAS")
    if (dim(snvtable)[1] > 0) {
        return("filtered")
    } else {
        return("none")
    }
})
names(krasSalvaged) <- as.character(driverTable$patient[which(driverTable$KRAS == "none")])
driverTable$KRAS <- as.character(driverTable$KRAS)
driverTable$KRAS[which(driverTable$KRAS == "none")] <- unlist(krasSalvaged)

driverTable <- melt(driverTable, id = "patient")
colnames(driverTable)[2] <- "gene"
colnames(driverTable)[3] <- "clonality"
driverTable$clonality <- ifelse(driverTable$clonality == "C", 
                                "clonal", 
                                ifelse(driverTable$clonality == "S", 
                                       "subclonal",
                                       as.character(driverTable$clonality)))


######################### Checking copy number state of genes of interest #########################

bailey <- read.table("data/bailey.driverGenes", fill = TRUE, header = TRUE, stringsAsFactors = FALSE)
bailey <- bailey$Gene[bailey$Cancer%in%c('PANCAN', 'PAAD')]

genesOfInterest <- c("KRAS", "TP53", "CDKN2A", "SMAD4", "TGFBR1", "TGFBR2", "PALB2", "ATM", "BRCA1", "BRCA2")
genesOfInterest <- union(genesOfInterest, unique(driverTable$gene))
genesOfInterest <- genesOfInterest[genesOfInterest %in% bailey]

txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
exons <- exonsBy(txdb, by = "gene")

rangesOfInterest <- lapply(mapIds(org.Hs.eg.db, genesOfInterest, "ENTREZID", "SYMBOL"), function(x) {
    tmp.geneRange <- exons[[x]]
    GRanges(seqnames = unique(seqnames(tmp.geneRange)), ranges = IRanges(start = start(ranges(tmp.geneRange))[1], end = end(ranges(tmp.geneRange))[length(tmp.geneRange)]), strand = unique(strand(tmp.geneRange)))
})

geneCopyNumberStates <- lapply(all.ASCAT.segs, function(pat.segs) {
    lapply(pat.segs, function(region.seg) {
        region.seg$chr <- paste0("chr", region.seg$chr)
        region.seg.GR  <- makeGRangesFromDataFrame(region.seg, ignore.strand = TRUE, seqnames.field = "chr", start.field = "startpos", end.field = "endpos")
        gene.seg <- lapply(genesOfInterest, function(z) {
            region.seg[queryHits(findOverlaps(region.seg.GR, rangesOfInterest[[z]])), ]
        })
        names(gene.seg) <- genesOfInterest
        return(gene.seg)
    })
})

calculate.cpn.classification <- function(sample.segs, use.raw = TRUE, gain.threshold = log2(2.5/2), loss.threshold = log2(1.5/2)) {
    if (use.raw) {
        sample.segs[, "cnTotalNew"] <- sample.segs[, "nAraw"] + sample.segs[, "nBraw"]
    } else {
        sample.segs[, "cnTotalNew"] <- sample.segs[, "cnTotal"]
    }

    sample.segs[, "logValue"] <- log2(sample.segs[, "cnTotalNew"] / sample.segs[, "Ploidy"])
    sample.segs[, "classification"] <- ifelse(sample.segs[, "logValue"] >= gain.threshold, 
                                              "gain",
                                              ifelse(sample.segs[, "logValue"] <= loss.threshold,
                                                     "loss",
                                                     "none"))
    return(sample.segs)
}

tmp <- lapply(geneCopyNumberStates, function(pat.segs) {
    tmp <- lapply(pat.segs, function(region.segs) {
        tmp  <- lapply(region.segs, calculate.cpn.classification)
        tmp2 <- lapply(names(tmp), function(x) cbind(gene = x, tmp[[x]]))
        tmp3 <- Reduce(rbind, tmp2)
        tmp3[, "LOH"] <- ifelse(tmp3[, "nMinor"] == 0, TRUE, FALSE)
        return(tmp3)
    })
    tmp <- Reduce(rbind, tmp)
    return(tmp)
})
geneTableFull <- Reduce(rbind, tmp)

geneTableFull$gene <- factor(as.character(geneTableFull$gene), levels = genesOfInterest)
geneTableFull$classification <- factor(geneTableFull$classification, levels = c("gain", "loss", "none"))
geneTableFull$LOHlabel <- ifelse(geneTableFull$LOH, "*", "")
geneTableFull$combClass <- geneTableFull$classification
geneTableFull$combClass <- ifelse(geneTableFull$classification == "loss" | geneTableFull$LOH, "loss or LOH", as.character(geneTableFull$combClass))
geneTableFull$patient  <- sapply(strsplit(gsub("S_", "", geneTableFull$sample), split = "_"), function(x) x[1])

cosmicGenes <- read.table("data/cancerGeneCensus_ALL_20190206.tsv", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
cosmicGenes <- subset(cosmicGenes, cosmicGenes$Gene.Symbol %in% genesOfInterest)
cosmicGenesLabel <- setNames(cosmicGenes$Molecular.Genetics, cosmicGenes$Gene.Symbol)

geneTable <- geneTableFull[, c("patient", "sample", "gene", "cnTotalNew", "logValue", "classification", "LOH", "LOHlabel", "combClass")]

geneTable <- geneTable %>% group_by(patient, gene) %>% 
                mutate(timing = ifelse(all(combClass == "loss or LOH"),
                                       "clonal loss", 
                                       ifelse(all(combClass == "gain"),
                                              "clonal gain",
                                              ifelse(any(combClass == "loss or LOH") & any(combClass == "gain"),
                                                     ifelse(cosmicGenesLabel[as.character(gene)] == "Rec",
                                                            "subclonal loss",
                                                            ifelse(cosmicGenesLabel[as.character(gene)] == "Dom",
                                                                   "subclonal gain",
                                                                   "mixed")),
                                                     ifelse(any(combClass == "loss or LOH"),
                                                            "subclonal loss",
                                                            ifelse(any(combClass == "gain"),
                                                                   "subclonal gain",
                                                                   "")))))) %>%
                mutate(clonality = ifelse(timing %in% c("clonal loss", "clonal gain"), 
                                          "clonal", 
                                          ifelse(timing %in% c("subclonal loss", "subclonal gain", "mixed"),
                                                 "subclonal",
                                                 "none"))) %>% 
                mutate(eventType = ifelse(timing %in% c("clonal loss", "subclonal loss"),
                                          "loss",
                                          ifelse(timing %in% c("clonal gain", "subclonal gain"),
                                                 "gain",
                                                 ifelse(timing %in% "mixed",
                                                        "mixed",
                                                        "none"))))
geneTable$clonality[which(geneTable$patient == "PAN105" & geneTable$gene == "HNF1A")] <- "clonal"

allPatients    <- unique(driverTable$patient)
remainingGenes <- setdiff(genesOfInterest, unique(driverTable$gene))

combinedGenomeTable <- geneTable %>% 
                            dplyr::select(patient, gene, clonality, eventType) %>% 
                            unique() %>%
                            bind_rows(cbind(rbind(driverTable, data.frame(patient = rep(allPatients, each = length(remainingGenes)), gene = rep(remainingGenes, length(allPatients)), clonality = rep("none", length(allPatients) * length(remainingGenes)), stringsAsFactors = FALSE)), eventType = "mut"))
combinedGenomeTable$clonality[combinedGenomeTable$clonality == "filtered"] <- "clonal"


######################### Combining data sets for plotting summary figure #########################

mutationLoad                  <- melt(mutationLoadOrig)
names(mutationLoad)[2]        <- "clonality"

mutationSignaturesP           <- melt(mutationSignaturesAll)
names(mutationSignaturesP)[2] <- "signatures"
mutationSignaturesP$signatures <- factor(mutationSignaturesP$signatures, levels = c(use.sigs, "other"))

### ordering patients based on mutation load
orderedPat <- mutationLoad %>% group_by(patient) %>% summarize(n = sum(value)) %>% arrange(desc(n)) %>% pull(patient)
mutationLoad$patient          <- factor(as.character(mutationLoad$patient), levels = as.character(orderedPat))
combinedGenomeTable$patient   <- factor(combinedGenomeTable$patient, levels = as.character(orderedPat))
combinedGenomeTable$clonality <- factor(combinedGenomeTable$clonality, levels = c("clonal", "subclonal", "none"))
mutationSignaturesP$patient   <- factor(as.character(mutationSignaturesP$patient), levels = as.character(orderedPat))


######################### Combining data sets for plotting summary figure #########################

fullPurityInfo <- read.table("data/fullPurityInfo.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
fullPurityInfo <- fullPurityInfo %>%
    filter(!is.na(purity)) %>%
    mutate(Patient = factor(Patient, levels = orderedPat)) %>%
    mutate(type = ifelse(grepl("cfDNA", Sample), 
                         "cfDNA", 
                         ifelse(grepl("Met", Sample), 
                                "Met", 
                                ifelse(grepl("PDX", Sample), "PDX", "Tumor"))))


############################################# Plotting #############################################

no.cpn.patients  <- c("PAN106", "PAN119", "PAN121")
no.sign.patients <- c("PAN106", "PAN112", "PAN116", "PAN120", "PAN121", "PAN129", "PANVH3")

col.palette      <- setNames(c("#2b8cbe", "#e34a33", "transparent"), c("clonal", "subclonal", "none"))
col.palette.sigs <- setNames(c("#ed968c", "#88a0dc", "#f9d14a", "#e78429", "#7c4b73", "#ab3329", "#636363"), c(use.sigs, "other"))
shape.palette    <- setNames(c(24, 25, NA, NA), c("gain", "loss", "mixed", "none"))

col.labels       <- setNames(c("clonal", "subclonal", ""), c("clonal", "subclonal", "none"))
shape.labels     <- setNames(c("gain", "loss", "", ""), c("gain", "loss", "mixed", "none"))

plot.mutationLoad <- ggplot(mutationLoad, aes(patient, value, group = clonality, fill = clonality)) + 
                        geom_bar(stat = "identity", position = position_stack(reverse = TRUE)) + 
                        scale_fill_manual(values = col.palette) + 
                        ylab("Number of mutations") +
                        theme_cowplot(font_size = 12) +
                        theme(legend.position = "none") +
                        theme(axis.text.x = element_blank(), axis.title.x = element_blank(), axis.ticks.x = element_blank()) 

plot.genomeOverview <- ggplot(combinedGenomeTable, aes(patient, gene)) +
                        geom_tile(data = combinedGenomeTable %>% filter(eventType == "mut") %>% unique(), aes(fill = clonality), color = "white", linewidth = 0.25) +
                        geom_point(data = combinedGenomeTable %>% filter(eventType != "mut") %>% unique(), aes(fill = clonality, shape = eventType), size = 3) +
                        geom_hline(data = combinedGenomeTable %>% filter(eventType == "mut") %>% unique(), aes(yintercept = as.numeric(gene) + 0.5), color = "gray") +
                        geom_vline(data = combinedGenomeTable %>% filter(eventType == "mut") %>% unique(), aes(xintercept = as.numeric(patient) + 0.5), color = "gray") +
                        scale_fill_manual(name = "Clonality", values = col.palette, labels = col.labels, guide = guide_legend(override.aes = list(shape = c(NA, NA, NA)))) +
                        scale_shape_manual(name = "CN aberration", values = shape.palette, labels = shape.labels) +
                        ylim(rev(levels(geneTableFull$gene))) +
                        xlab("Patients") +
                        ylab("Cancer genes") +
                        theme_cowplot(font_size = 12) +
                        theme(axis.text.x = element_blank(), axis.title.x = element_blank(), axis.ticks.x = element_blank())
                        
plot.mutationSignaturesP <- ggplot(mutationSignaturesP, aes(patient, value, group = signatures, fill = signatures, alpha = patient)) +
                                geom_bar(stat = "identity", position = position_stack()) +
                                scale_fill_manual(name = "Signatures", values = col.palette.sigs) + 
                                scale_alpha_manual(values = ifelse(orderedPat %in% no.sign.patients, 0.5, 1), guide = "none") +
                                ylab("% signatures") +
                                theme_cowplot(font_size = 12) +
                                theme(axis.text.x = element_text(angle = 45, hjust = 1, color = ifelse(orderedPat %in% no.cpn.patients, "gray", "black")))
                                
plot.purityOverview <- ggplot(fullPurityInfo, aes(Patient, purity, shape = type)) + 
    geom_point() +
    scale_y_continuous(name = "Sample purity", limits = c(0,1)) +
    scale_shape_discrete(name = "Sample type") +
    theme_cowplot(font_size = 12) +
    theme(axis.text.x = element_blank(), axis.title.x = element_blank(), axis.ticks.x = element_blank())


g.plot.mutationLoad        <- ggplotGrob(plot.mutationLoad + theme(plot.margin = unit(c(0.5,0.1,0.0,0.1),"cm")))
g.plot.purityOverview      <- ggplotGrob(plot.purityOverview + theme(plot.margin = unit(c(0.0,0.1,0.0,0.1), "cm")))
g.plot.genomeOverview      <- ggplotGrob(plot.genomeOverview + theme(plot.margin = unit(c(0,0.1,0,0.1),"cm")))
g.plot.mutationSignaturesP <- ggplotGrob(plot.mutationSignaturesP + theme(plot.margin = unit(c(0,0.1,0.5,0.1),"cm")))

max.width <- unit.pmax(g.plot.mutationLoad$widths, g.plot.purityOverview$widths, g.plot.genomeOverview$widths, g.plot.mutationSignaturesP$widths)

g.plot.mutationLoad$widths        <- max.width
g.plot.purityOverview$widths      <- max.width
g.plot.genomeOverview$widths      <- max.width
g.plot.mutationSignaturesP$widths <- max.width

plot.all  <- list(g.plot.mutationLoad, g.plot.purityOverview, g.plot.genomeOverview, g.plot.mutationSignaturesP)
plot.all.layout <- matrix(c(rep(1,2), rep(2, 1.5), rep(3,4.5), rep(4,2.5)), ncol = 1)
fullplot.genomicOverview <- arrangeGrob(grobs = plot.all, ncol = 1, nrow = dim(plot.all.layout[1]), layout_matrix = plot.all.layout)

pdf(paste0(plotDir, "genomicOverview.new.pdf"), width = 7, height = 10, useDingbats = FALSE)
grid.draw(fullplot.genomicOverview)
dev.off()


############################# Plotting clonal and subclonal signatures #############################

mutationsoutPerClonality <- mut.to.sigs.input(mutations, sample.id = "Clonality", pos = "start", alt = "var")
mutationSignaturesFull   <- lapply(unique(mutations$Clonality), function(x) whichSignatures(mutationsoutPerClonality, sample.id = x, contexts.needed = TRUE, tri.counts.method = 'exome2genome', associated = use.sigs, signatures.ref = signatures.genome.cosmic.v3.may2019))

mutationSignatures     <- lapply(mutationSignaturesFull, function(x) x[["weights"]])
mutationSignatures     <- Reduce(rbind, mutationSignatures)
mutationSignatures     <- mutationSignatures[, which(colnames(mutationSignatures) %in% use.sigs)]
extraSamples           <- matrix(0, nrow = sum(!fullClonality %in% rownames(mutationSignatures)), ncol = length(use.sigs))
rownames(extraSamples) <- fullClonality[which(!fullClonality %in% rownames(mutationSignatures))]
colnames(extraSamples) <- colnames(mutationSignatures)
mutationSignatures     <- rbind(mutationSignatures, extraSamples)
mutSigIDs              <- Reduce(rbind, strsplit(rownames(mutationSignatures), split = ":"))
rownames(mutSigIDs)    <- NULL
colnames(mutSigIDs)    <- c("patient", "clonality")
mutationSignatures     <- cbind(mutSigIDs, mutationSignatures)
mutationSignatures     <- cbind(mutationSignatures, other = 1 - rowSums(mutationSignatures[, c(-1, -2)]))

mutationSignaturesCS   <- melt(mutationSignatures)
names(mutationSignaturesCS)[3] <- "signatures"
mutationSignaturesCS$signatures <- factor(mutationSignaturesCS$signatures, levels = c(use.sigs, "other"))

orderedPat.mut <- mutationSignaturesCS %>% filter(signatures == "SBS3", clonality == "C") %>% arrange(desc(value)) %>% pull(patient)
orderedPat.mut <- as.character(orderedPat.mut)
mutationSignaturesCS$patient <- factor(mutationSignaturesCS$patient, levels = orderedPat.mut)

no.sign.combination <- c("PAN104:S", "PAN105:S", "PAN106:C", "PAN106:S", "PAN107:S", "PAN108:S", "PAN109:S", "PAN110:S", "PAN112:C", "PAN112:S", "PAN115:C", "PAN115:S", "PAN116:C", "PAN116:S", "PAN117:S", "PAN119:S", "PAN120:C", "PAN120:S", "PAN121:C", "PAN121:S", "PAN125:S", "PAN129:C", "PAN129:S", "PANVH2:S", "PANVH3:C", "PANVH3:S")
combinedPatSign     <- paste(rep(levels(mutationSignaturesCS$patient), each = 2), rep(c("C", "S"), length(unique(mutationSignaturesCS$patient))), sep = ":")

mutationSignaturesCS$combinedPatSign <- paste(mutationSignaturesCS$patient, mutationSignaturesCS$clonality, sep = ":")

plot.mutationSignaturesCS <- ggplot(mutationSignaturesCS, aes(patient, value, group = signatures, fill = signatures, alpha = combinedPatSign)) +
                                geom_bar(stat = "identity", position = position_stack()) +
                                scale_fill_manual(name = "Signatures", values = col.palette.sigs) + 
                                scale_alpha_manual(values = ifelse(combinedPatSign %in% no.sign.combination, 0.5, 1), guide = FALSE) +
                                ylab("% signatures") +
                                facet_wrap(~ clonality, ncol = 1, labeller = labeller(clonality = setNames(c("Clonal", "Subclonal"), c("C", "S")))) +
                                theme_cowplot(font_size = 12) +
                                theme(axis.text.x = element_text(angle = 45, hjust = 1), strip.background = element_rect(fill = "white"))

pdf(paste0(plotDir, "signaturesSplitClonalSubclonal.orderedBySBS3.pdf"), width = 7, height = 4, useDingbats = FALSE)
plot.mutationSignaturesCS
dev.off()


############################## Plotting copy number events per sample ##############################

col.palette.2      <- setNames(c("#2b8cbe", "#e34a33", NA), c("clonal", "subclonal", "none"))
shape.palette.2    <- setNames(c(24, 25, NA), c("gain", "loss", "none"))

col.labels.2       <- setNames(c("clonal", "subclonal", ""), c("clonal", "subclonal", "none"))
shape.labels.2     <- setNames(c("gain", "loss", ""), c("gain", "loss", "none"))

plot.geneTable <- data.frame(geneTable, stringsAsFactors = FALSE)
plot.geneTable$patient   <- gsub("S_", "", plot.geneTable$patient)
plot.geneTable$sample    <- gsub("S_", "", plot.geneTable$sample)
plot.geneTable$eventType <- ifelse(plot.geneTable$combClass == "loss or LOH", "loss", as.character(plot.geneTable$combClass))
plot.geneTable$eventType <- factor(plot.geneTable$eventType, levels = c("gain", "loss", "none"))
plot.geneTable$clonality <- factor(plot.geneTable$clonality, levels = c("clonal", "subclonal", "none"))

cumSumNRegions <- cumsum(plot.geneTable %>% ungroup() %>% select(patient, sample) %>% unique() %>% group_by(patient) %>% tally() %>% unique() %>% pull(n))

plot.copyNumberOverview <- ggplot(plot.geneTable, aes(sample, gene)) +
                                geom_point(aes(fill = clonality, shape = eventType), size = 3) +
                                geom_hline(data = plot.geneTable %>% select(patient, gene) %>% unique(), aes(yintercept = as.numeric(gene) + 0.5), color = "gray") +
                                geom_vline(data = plot.geneTable %>% select(patient, sample, gene) %>% unique(), xintercept = cumSumNRegions + 0.5, color = "black") +
                                scale_fill_manual(name = "Clonality", values = col.palette.2[levels(plot.geneTable$clonality)], labels = col.labels.2[levels(plot.geneTable$clonality)], guide = guide_legend(override.aes = list(shape = c(21, 21, NA)))) +
                                scale_shape_manual(name = "CN aberration", values = shape.palette.2[levels(plot.geneTable$eventType)], labels = shape.labels.2[levels(plot.geneTable$eventType)]) +
                                ylim(rev(levels(geneTableFull$gene))) +
                                xlab("Samples") +
                                ylab("Cancer genes") +
                                theme_cowplot(font_size = 12) +
                                theme(axis.text.x = element_text(angle = 45, hjust = 1))#, color = ifelse(orderedPat %in% no.cpn.patients, "gray", "black")))

pdf(paste0(plotDir, "copynumberClassification.sampleLevel.pdf"), width = 10, height = 7, useDingbats = FALSE)
plot.copyNumberOverview
dev.off()


############################## Creating and plotting BRCA comparisons ##############################

BRCAgl.patients  <- c("PAN103", "PAN119", "PANVH1", "PANVH2", "PANVH3")

### plotting mutationload as boxplot
mutationLoadBRCA <- mutationLoad
mutationLoadBRCA$BRCA <- ifelse(mutationLoadBRCA$patient %in% BRCAgl.patients, "BRCA mutation", "no BRCA mutation")
mutationLoadBRCA$BRCAtype <- ifelse(mutationLoadBRCA$patient %in% BRCAgl.patients, ifelse(mutationLoadBRCA$patient == "PAN103", "GL+somatic BRCA mutation", "GL BRCA mutation"), "no BRCA mutation")
mutationLoadBRCA$value <- sapply(1:dim(mutationLoadBRCA)[1], function(i) if(mutationLoadBRCA$clonality[i] %in% c("clonal", "subclonal") & mutationLoadBRCA$value[i] == 0) {return(NA)}else{return(mutationLoadBRCA$value[i])})

p1 <- ggplot(mutationLoadBRCA %>% group_by(patient) %>% mutate(n = sum(value, na.rm = T), clonality = "total") %>% dplyr::select(patient, BRCA, BRCAtype, n, clonality) %>% unique(), aes(BRCA, n)) + geom_boxplot() + geom_point(aes(color = BRCAtype)) + 
            stat_compare_means(method = "wilcox.test") +
            scale_color_manual(values = c("#e31a1c", "#fdbf6f", "#1f78b4")) +
            xlab("") +
            ylab("Number of mutations") +
            theme_cowplot(font_size = 12) +
            theme(legend.position = "none") +
            theme(strip.background = element_rect(fill = "white"), strip.text.x = element_text(margin = margin(1, 0, 1, 0, "lines"), size = 14)) +
            facet_wrap(~ clonality, labeller = labeller(clonality = setNames("All mutations", "total")))


totalMuts  <- table(mutTablePASS$SampleID)
signature3 <- mutationLoadBRCA[, c("patient", "BRCA", "BRCAtype")] %>% unique()
signature3$weight <- as.numeric(mutationSignaturesAll[as.character(signature3$patient), "SBS3"])
signature3$count  <- as.numeric(signature3$weight * totalMuts[as.character(signature3$patient)])

q1 <- ggplot(signature3, aes(BRCA, count)) + 
            geom_boxplot() + 
            geom_jitter(aes(color = BRCAtype), width = 0.1) + 
            stat_compare_means(method = "wilcox.test") +
            scale_color_manual(values = c("#e31a1c", "#fdbf6f", "#1f78b4")) +
            xlab("") +
            ylab("Number of mutations SBS3") +
            theme_cowplot(font_size = 12) +
            theme(axis.title.x = element_blank())

pdf(paste0(plotDir, "supplementary.sbs3.pdf"), useDingbats = FALSE)
q1
dev.off()



####################################################################################################
#########################                   Random Code                   #########################
####################################################################################################
summary(mutationLoad %>% group_by(patient) %>% summarize(n = sum(value)) %>% pull(n))

mean(mutationLoad %>% filter(clonality == "clonal") %>% pull(value))
sd(mutationLoad %>% filter(clonality == "clonal") %>% pull(value))
mean(mutationLoad %>% filter(clonality == "subclonal") %>% pull(value))
sd(mutationLoad %>% filter(clonality == "subclonal") %>% pull(value))


################################# Calculating wGII of each sample #################################
    source("scripts/calculatewGII.R")
    allPatient.segs <- lapply(all.ASCAT.segs, function(x) {
        if (length(x) == 0) return(NULL)
        tmp <- Reduce(rbind, x)
        names(tmp)[1] <- "region"
        tmp <- cbind(Patient = sapply(strsplit(tmp$region, split = "_"), function(x) x[2]), tmp, stringsAsFactors = FALSE)
    })
    allPatient.segs <- allPatient.segs[!sapply(allPatient.segs, is.null)]
    allPatient.WGII.perSample <- lapply(allPatient.segs, calculate.wGIIscore, pro.Region = TRUE)
    tmp <- data.frame(Reduce(rbind, allPatient.WGII.perSample), stringsAsFactors = FALSE, row.names = NULL)
    summary(as.numeric(as.character(tmp$wGII)))

pdf("plots/correlation.meanVAFpurity.pdf", useDingbats = FALSE)
fullPurityInfo %>% filter(!is.na(purity), !is.na(meanVAF)) %>% ggplot(aes(meanVAF, purity)) + geom_point() + geom_smooth(method = "lm", se = FALSE) + stat_cor(size = 7) + theme_cowplot(font_size = 21) + xlab("mean clonal VAF") + ylab("ASCAT estimated purity")
dev.off()
