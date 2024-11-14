setwd(
  "~/Dropbox/cary_projects/DODA/reconciliation/strict_reconciliation_extra_dups_removed/ASR_on_no_og_pruned/final_run/fastml_synth_noF_raxml_brlen_opt_brlen_JTT+G/"
  )
df <- read.table("Ancestral_MaxMarginalProb_Char_Indel.txt", header=T)

# nodes <- c("N7", "N8", "N41", "N42", "N9", "N171", "N224", "N232", "N4")
# expanded set for full MS
nodes <- c("N2", "N168", "N169", "N6", "N7", "N3", "N4", "N171", "N224", "N232",
           "N225", "N172", "N8", "N41", "N42", "N9", "N72")

pdf("~/Dropbox/cary_projects/DODA/manuscript/figures/SuppFigX_node_PP.pdf",
    paper="a4")
par(mfrow = c(6,3))
par(mar=rep(2,4))
for (i in nodes) {
  siteProbs <- df[which(df$Node == i),]
  nonGap <- siteProbs[which(siteProbs$Char != "-"),]
  avg <- mean(nonGap$CharProb)
  hist(nonGap$CharProb,
       main = i,
       xlab = NULL,
       ylab = NULL,
       breaks = 10,
       xlim = c(0,1.0),
       ylim = c(0, 250),
       border = FALSE
       )
  text(0.2, 240, paste0("Av PP = ", signif(avg, 3)))
}

dev.off()
