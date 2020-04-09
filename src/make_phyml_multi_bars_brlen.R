inf <- commandArgs(trailingOnly = TRUE)[1]
intr <- commandArgs(trailingOnly = TRUE)[2]
name_cex <- as.numeric(commandArgs(trailingOnly = TRUE)[3])
num_height <- as.numeric(commandArgs(trailingOnly = TRUE)[4])
library(ape)
paste(inf)
name <- paste0(strsplit(inf,"\\.")[[1]][1])
phyml_multi_bp <- read.csv(inf,header=F)
#phyml_multi_bp <- read.csv("breaks.csv",header=F)
trees <- read.tree(intr)
#trees <- read.tree("trees.txt")

# for (i in 1:length(phyml_multi_bp$V1)) {
#   if (i == 1) {
#     phyml_multi_bp$V3[1] <- 1
#   } else {
#     phyml_multi_bp$V3[i] <- phyml_multi_bp$V1[i-1]+1
#   }
# }

ntree <- length(trees)

for (i in 1:length(phyml_multi_bp$V1)) {
  if (phyml_multi_bp$V3[i] == 0) {
    phyml_multi_bp$V3[i] <- "lightgrey"
  } else if (phyml_multi_bp$V3[i] == 1) {
    phyml_multi_bp$V3[i] <- "lightblue"
  } else {
    phyml_multi_bp$V3[i] <- "lightgreen"
  }
}

layout_matrix <- matrix(c(rep(1,ntree),2:(1+ntree)),2,ntree,byrow = T)
par(mar=c(1,1,1,1),oma=c(3,1,0,0))
layout(layout_matrix)

plot.new()
plot.window(xlim=c(min(phyml_multi_bp$V1),max(phyml_multi_bp$V2)),ylim=c(0,1))
rect(phyml_multi_bp$V1,xright=phyml_multi_bp$V2,ybottom = 0.1,ytop=.2,col = phyml_multi_bp$V3,border = "white")
title(main = paste0(name),adj=0,line=-2)
axis(1)

strlengths <- vector(mode="numeric",length=length(trees))
strheights <- vector(mode="numeric",length=length(trees))
for (i in 1:length(trees)) {
  tr <- ladderize(trees[[i]],right=F)
  plot.phylo(tr,use.edge.length = T,x.lim = max(node.depth.edgelength(tr)),y.lim=max(node.height(tr)),no.margin = F,cex=name_cex)
  add.scale.bar(x=max(node.depth.edgelength(tr))*.75,y=1,cex=name_cex)
  text(x=0,y=num_height,labels=i)
  strlengths[i] <- strwidth(tr$tip.label[nchar(tr$tip.label)==max(nchar(tr$tip.label))])
  strheights[i] <- strheight(tr$tip.label[nchar(tr$tip.label)==max(nchar(tr$tip.label))])
}

svg(paste0("breaks_brlen.svg"))
par(mar=c(1,1,1,1),oma=c(3,1,0,0))
layout(layout_matrix)

plot.new()
plot.window(xlim=c(min(phyml_multi_bp$V1),max(phyml_multi_bp$V2)),ylim=c(0,1))
rect(xleft=phyml_multi_bp$V1,xright=phyml_multi_bp$V2,ybottom = 0.1,ytop=.2,col = phyml_multi_bp$V3,border = "white")
axis(1)
legend(x="topright",title = "Tree",legend = as.character(c(1:ntree)),fill = c("lightgrey","lightblue","lightgreen"))
title(name,adj=0,line=-2)

for (i in 1:length(trees)) {
  tr <- ladderize(trees[[i]],right=F)
  plot.phylo(tr,use.edge.length = T,x.lim = max(node.depth.edgelength(tr))+1.5*strlengths[i],y.lim=max(node.height(tr))+strheights[i],no.margin = F,cex=name_cex)
  add.scale.bar(x=(max(node.depth.edgelength(tr))+1.5*strlengths[i])*.75,y=1,cex=name_cex)
  text(x=0,y=num_height,labels=i)
}
