library(ape)
inf <- commandArgs(trailingOnly = TRUE)[1]
name_cex <- as.numeric(commandArgs(trailingOnly = TRUE)[2])
num_height <- as.numeric(commandArgs(trailingOnly = TRUE)[3])
paste(inf)
name <- paste0(strsplit(inf,"\\.")[[1]][1:3],collapse=".")
gard_bp <- read.csv(inf,sep=" ",header=FALSE)
gard_bp <- gard_bp[order(gard_bp$V3),]
ntree <- length(gard_bp$V3)
nr<- ntree/3
if (nr%%1!=0) {nr <- ceiling(nr)}

mult_of_3 <- function(x) {
  if (x%%3==0) {
    return(x)
  } else {
    mult_of_3(x+1)
  }
}

if (length(c(1,1,1,seq(2,ntree+1,1)))%%3 == 0) {
  layout.matrix <- matrix(c(1,1,1,seq(2,ntree+1,1)),nr+1,3,byrow =T)
} else {
  l <- length(c(1,1,1,seq(2,ntree+1,1)))
  new_l <- mult_of_3(l)
  pad <- new_l - l
  layout.matrix <- matrix(c(1,1,1,seq(2,ntree+1,1),rep(0,pad)),nr+1,3,byrow =T)
}

par(mar=c(1,1,1,1),oma=c(3,1,0,0))
layout(mat = layout.matrix)

plot.new()
plot.window(xlim=c(min(gard_bp$V1),max(gard_bp$V2)),ylim=c(0,1))
rect(xleft=gard_bp$V1,xright=gard_bp$V2,ybottom = 0.1,ytop=.3,col = c("lightblue","lightgrey"),border = "white")
title(main = paste0(name),adj=0,line=-2)
text(x=gard_bp$V1+(gard_bp$V2-gard_bp$V1)/2,y=0.4,labels=gard_bp$V3+1,cex = 0.9)
axis(1)

strlengths <- vector(mode="numeric",length=length(gard_bp$V4))
strheights <- vector(mode="numeric",length=length(gard_bp$V4))
for (i in 1:length(gard_bp$V4)) {
  tr <- ladderize(read.tree(text = as.character(gard_bp$V4[i])),right=F)
  plot.phylo(tr,use.edge.length = T,x.lim = max(node.depth.edgelength(tr)),y.lim=max(node.height(tr)),no.margin = F,cex=name_cex)
  add.scale.bar(x=max(node.depth.edgelength(tr))*.75,y=1,cex=name_cex)
  text(x=0,y=num_height,labels=i)
  strlengths[i] <- strwidth(tr$tip.label[nchar(tr$tip.label)==max(nchar(tr$tip.label))])
  strheights[i] <- strheight(tr$tip.label[nchar(tr$tip.label)==max(nchar(tr$tip.label))])
}

svg(paste0(inf,"breaks_brlen.svg",collapse="_"))
par(mar=c(1,1,1,1),oma=c(3,1,0,0))
layout(mat = layout.matrix)

plot.new()
plot.window(xlim=c(min(gard_bp$V1),max(gard_bp$V2)),ylim=c(0,1))
rect(xleft=gard_bp$V1,xright=gard_bp$V2,ybottom = 0.1,ytop=.3,col = c("lightblue","lightgrey"),border = "white")
title(main = paste0(name),adj=0,line=-2)
text(x=gard_bp$V1+(gard_bp$V2-gard_bp$V1)/2,y=0.4,labels=gard_bp$V3+1,cex = 0.9)
axis(1)

for (i in 1:length(gard_bp$V4)) {
  tr <- ladderize(read.tree(text = as.character(gard_bp$V4[i])),right=F)
  plot.phylo(tr,use.edge.length = T,x.lim = max(node.depth.edgelength(tr))+1.5*strlengths[i],y.lim=max(node.height(tr))+strheights[i],no.margin = F,cex=name_cex)
  add.scale.bar(x=(max(node.depth.edgelength(tr))+1.5*strlengths[i])*.75,y=1,cex=name_cex)
  text(x=0,y=num_height,labels=i)
}

dev.off()




