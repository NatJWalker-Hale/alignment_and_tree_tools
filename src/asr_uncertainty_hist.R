setwd(
  "~/Dropbox/cary_projects/DODA/reconciliation/strict_reconciliation_extra_dups_removed/ASR_on_no_og_pruned/final_run/fastml_synth_noF_raxml_brlen_opt_brlen_JTT+G/"
  )
df <- read.table("Ancestral_MaxMarginalProb_Char_Indel.txt", header=T)

nodes <- c("N7", "N8", "N41", "N42", "N9", "N171", "N224", "N232", "N4")

par(mfrow = c(3,3))
for (i in nodes) {
  siteProbs <- df[which(df$Node == i),]
  nonGap <- siteProbs[which(siteProbs$Char != "-"),]
  hist(nonGap$CharProb,
       main = i,
       xlab = NULL,
       ylab = NULL,
       breaks = 10,
       xlim = c(0,1.0),
       ylim = c(0, 250),
       )
}
