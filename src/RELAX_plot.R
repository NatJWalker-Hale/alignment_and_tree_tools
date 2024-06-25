setwd("~/Downloads/")
library(scales)
library(ggplot2)

df <- read.table("tmp.tsv", header = T,
               stringsAsFactors = T, na.strings = "null")

df$gene <- factor(df$gene,
                  levels = c("CHS", "CHI", "F3-H", "F3H", "F3-5-H", "FNS2",
                             "FLS", "DFR", "LAR_1", "LAR_2", "ANS", "ANR",
                             "75C1_1", "75C1_2", "78D2_1", "78D2_2", "79B1",
                             "MATE", "ABCCs", "AHA10", "PAP1", "PAP2", "EGL1",
                             "MYC1", "TT8", "TTG1"))

df <- df[order(df$gene),]

# genelabels <- c("CHS", "CHI", "F3'H", "F3H", "F3'5'H", "FNS2",
#                 "FLS", "DFR", "LAR1", "LAR2", "ANS", "ANR",
#                 "A5GT-1", "A5GT-2", "A3GT-1", "A3GT-2", "F3GT",
#                 "MATE", "ABCCs", "AHA10", "PAP1", "PAP2", "EGL1",
#                 "MYC1", "TT8", "TTG1")

# mytrans <- function(x, base=exp(1)) {
#   ifelse(x <= 0, x, log(x, base=base))
# }

theme_set(theme_bw())
theme_update(text = element_text(size=12),
             panel.grid.major = element_blank(),
             panel.grid.minor = element_blank(),
             strip.background = element_blank()
)

for (i in levels(df$gene)) {
  subdf <- df[ which(df$gene == i), ]
  subdf.plot <- cbind(stack(subdf[,12:14])$values, stack(subdf[,9:11])$values)
  subdf.plot <- data.frame(p = subdf.plot[,1],
                           w = subdf.plot[,2],
                           type = rep(c("ref", "test"), 3)
                           )
  trans <- pseudo_log_trans(sigma = 0.1, base = 10)
  svg(paste("all_RELAX_results/",i,"RELAX_alt.svg"), width = 6, height = 5)
  plot(p ~ trans$transform(w), data = subdf.plot, type = "h", lwd = 2, ylim = c(0,1),
       col = rep(c("#4f71c1", "#a41383"), 3), lend = 2,
       xlab = expression(paste(omega)),
       ylab = "Proportion of sites", xaxt = "n")
  #title(main = i, adj = 0, sub = paste("K = ", round(subdf$k, digits = 3)
                                     #  , "P = ", round(subdf$pval, digits =3)))
  mtext(i, side = 3, line = 1, adj = 0, cex = 1)
  mtext(paste("K = ", round(subdf$k, digits = 3),
              "P = ", round(subdf$pval, digits = 3),
              "AIC-c = ", round(subdf$aicc, digits = 3)),
        side = 3, line = .2, adj = 0, cex = 0.7)
  if (max(subdf.plot$w[!is.na(subdf.plot$w)]) < 10) {
    labvals <- c(0, 0.1, 0.5, 1, 2, 5, 10)
  } else if (max(subdf.plot$w[!is.na(subdf.plot$w)]) < 100) {
    labvals <- c(0, 0.5, 1, 2, 5, 10, 100)
  } else if (max(subdf.plot$w[!is.na(subdf.plot$w)]) < 1000) {
    labvals <- c(0, 1, 10, 100, 1000)
  } else if (max(subdf.plot$w[!is.na(subdf.plot$w)]) < 10000) {
    labvals <- c(0, 1, 10, 100, 1000, 10000)
  } else if (max(subdf.plot$w[!is.na(subdf.plot$w)]) < 100000) {
    labvals <- c(0, 1, 10, 100, 1000, 10000, 100000)
  }
  axis(1, at = trans$transform(seq(0, 1, 0.1)), labels = F)
  axis(1, at = trans$transform(labvals), labels = labvals)
  abline(v=1.004279, lty=2, col = "gray")
  dev.off()
}

df <- read.csv("20230912_RELAX_partitioned_descriptive_results.csv", header = T,
               stringsAsFactors = T)

df$gene <- factor(df$gene,
                  levels = c("CHS", "CHI", "F3-H", "F3H", "F3-5-H", "FNS2",
                             "FLS", "DFR", "LAR_1", "LAR_2", "ANS", "ANR",
                             "75C1_1", "75C1_2", "78D2_1", "78D2_2", "79B1",
                             "MATE", "ABCCs", "AHA10", "PAP1", "PAP2", "EGL1",
                             "MYC1", "TT8", "TTG1"))

df <- df[order(df$gene),]

genelabels <- c("CHS", "CHI", "F3'H", "F3H", "F3'5'H", "FNS2",
                "FLS", "DFR", "LAR1", "LAR2", "ANS", "ANR",
                "A5GT-1", "A5GT-2", "A3GT-1", "A3GT-2", "F3GT",
                "MATE", "ABCCs", "AHA10", "PAP1", "PAP2", "EGL1",
                "MYC1", "TT8", "TTG1")

for (i in levels(df$gene)) {
  subdf <- df[ which(df$gene == i), ]
  subdf.plot <- cbind(stack(subdf[,8:10])$values, stack(subdf[,5:7])$values)
  subdf.plot <- data.frame(p = subdf.plot[,1],
                           w = subdf.plot[,2],
                           type = rep(c("ref", "test"), 3)
  )
  trans <- pseudo_log_trans(sigma = 0.1, base = 10)
  svg(paste("all_RELAX_results/",i,"RELAX_PD.svg"), width = 6, height = 5)
  plot(p ~ trans$transform(w), data = subdf.plot, type = "h", lwd = 2, ylim = c(0,1),
       col = rep(c("#4f71c1", "#a41383"), 3), lend = 2,
       xlab = expression(paste(omega)),
       ylab = "Proportion of sites", xaxt = "n")
  #title(main = i, adj = 0, sub = paste("K = ", round(subdf$k, digits = 3)
  #  , "P = ", round(subdf$pval, digits =3)))
  mtext(i, side = 3, line = 1, adj = 0, cex = 1)
  mtext(paste("AIC-c = ", round(subdf$aicc, digits = 3)),
        side = 3, line = .2, adj = 0, cex = 0.7)
  if (max(subdf.plot$w) < 10) {
    labvals <- c(0, 0.1, 0.5, 1, 2, 5, 10)
  } else if (max(subdf.plot$w) < 100) {
    labvals <- c(0, 0.5, 1, 2, 5, 10, 100)
  } else if (max(subdf.plot$w) < 1000) {
    labvals <- c(0, 1, 10, 100, 1000)
  } else if (max(subdf.plot$w) < 10000) {
    labvals <- c(0, 1, 10, 100, 1000, 10000)
  } else if (max(subdf.plot$w) < 100000) {
    labvals <- c(0, 1, 10, 100, 1000, 10000, 100000)
  }
  axis(1, at = trans$transform(seq(0, 1, 0.1)), labels = F)
  axis(1, at = trans$transform(labvals), labels = labvals)
  abline(v=1.004279, lty=2, col = "gray")
  dev.off()
}
