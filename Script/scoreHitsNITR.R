# An algorithm predicting nitrofurantoin susceptibility based on genetic presence-absence and alterations.
# This function should be run in an R environment.
# Copyright (C) 2020 Yu Wan <wanyuac@126.com>
# Licensed under the GNU General Public Licence version 3 (GPLv3) <https://www.gnu.org/licenses/>.
# Creation: 25 Sep 2020; the latest modification: 25 Sep 2020

scoreHitsNITR <- function(hit_scores, down_sort = TRUE) {
    scores <- data.frame(Isolate = character(0), nfsA = numeric(0), nfsB = numeric(0),
                         ribE = numeric(0), oqxAB = numeric(0),
                         nfsA_class = character(0), nfsB_class = character(0),
                         ribE_class = character(0), oqxA_class = character(0),
                         oqxB_class = character(0), stringsAsFactors = FALSE)
    isolates <- unique(hit_scores$Isolate)  # This step guarantees that no NULL result from selection subset(hit_scores, Isolate == i).
    
    # Go through each isolate
    for (i in isolates) {
        hits_isolate <- subset(hit_scores, Isolate == i)  # Hits from isolate i
        genes_isolate <- unique(hits_isolate$Gene)  # Genes detected in this isolate
        
        # Scores for intrinsic genes
        scores_intrinsic <- c("nfsA" = NA, "nfsB" = NA, "ribE" = NA)
        classes_intrinsic <- c("nfsA" = NA, "nfsB" = NA, "ribE" = NA)
        for (g in c("nfsA", "nfsB", "ribE")) {
            if (g %in% genes_isolate) {
                hits_gene_isolate <- subset(hits_isolate, Gene == g)  # Hits of the current gene in the current isolate
                if (nrow(hits_gene_isolate) > 1) {  # More than one hit per gene
                    j <- which.min(hits_gene_isolate$Score)    # Full function of each gene renders NIT-S.
                    scores_intrinsic[[g]] <- hits_gene_isolate$Score[j]
                    classes_intrinsic[[g]] <- hits_gene_isolate$Class[j]  # Only record the protein class that gives the minimum risk score
                } else {
                    scores_intrinsic[[g]] <- hits_gene_isolate$Score
                    classes_intrinsic[[g]] <- hits_gene_isolate$Class
                }
            } else {  # Absence of gene g in the current isolate
                scores_intrinsic[[g]] <- 2  # Equivalent to gene knock-out, which is known to confer NIT-R.
                classes_intrinsic[[g]] <- "Absent"
            }
        }
        
        # Scores for acquired resistance genes
        scores_acquired <- c("oqxA" = NA, "oqxB" = NA)
        classes_acquired <- c("oqxA" = NA, "oqxB" = NA)
        for (g in c("oqxA", "oqxB")) {
            if (g %in% genes_isolate) {
                hits_gene_isolate <- subset(hits_isolate, Gene == g)
                if (nrow(hits_gene_isolate) > 1) {
                    j <- which.max(hits_gene_isolate$Score)  # Presence of the gene confers NIT-R.
                    scores_acquired[[g]] <- hits_gene_isolate$Score[j]
                    classes_acquired[[g]] <- hits_gene_isolate$Class[j]  # Only record the protein class giving the largest score
                } else {
                    scores_acquired[[g]] <- hits_gene_isolate$Score
                    classes_acquired[[g]] <- hits_gene_isolate$Class
                }
            } else {
                scores_acquired[[g]] <- 0  # No NIT-R when the gene is absent.
                classes_acquired[[g]] <- "Absent"
            }
        }
        
        # Compile scores of the current isolate
        scores_isolate <- data.frame(Isolate = i, nfsA = scores_intrinsic[["nfsA"]],
                                     nfsB = scores_intrinsic[["nfsB"]],
                                     ribE = scores_intrinsic[["ribE"]],
                                     oqxAB = min(scores_acquired[["oqxA"]], scores_acquired[["oqxB"]]),
                                     nfsA_class = classes_intrinsic[["nfsA"]],
                                     nfsB_class = classes_intrinsic[["nfsB"]],
                                     ribE_class = classes_intrinsic[["ribE"]],
                                     oqxA_class = classes_acquired[["oqxA"]],
                                     oqxB_class = classes_acquired[["oqxB"]],
                                     stringsAsFactors = FALSE)
        scores <- rbind.data.frame(scores, scores_isolate, stringsAsFactors = FALSE)
    }
    
    # Take the sum of per-gene scores
    scores$Sum <- mapply(sum, scores$nfsA, scores$nfsB, scores$ribE, scores$oqxAB)
    if (down_sort) {
        scores <- scores[order(scores$Sum, scores$Isolate, decreasing = TRUE), ]
    }
    scores <- scores[, c("Isolate", "nfsA", "nfsB", "ribE", "oqxAB", "Sum", "nfsA_class", "nfsB_class", "ribE_class", "oqxA_class", "oqxB_class")]
    
    return(scores)
}
