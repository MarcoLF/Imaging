# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

library(ggplot2)
library(grid)
library(plyr)
library(spatstat)
library(cowplot)
library(Hmisc)
library(gridExtra)


file = c("Image_5709.csv", "Image_5708.csv");


inpt = "~/Desktop/microscopy/"
out = "~/Desktop/microscopy/"


table = function(file, inpt, out, bio, name, exp) {
  
  datalist=list()
  
  if (bio != 1) {
    bio.rep = rep(rep(1:bio, each = bio), length(name));
    tec = rep(c(1:(bio)), bio*length(name));
  }
  else {
    bio.rep = c(rep(bio, length(file)));
    tec = rep(c(1:length(file)), bio);
  }
  
  
  if (length(name) != 1) {
    geno = rep(name, each = length(file)/length(name))
  }
  else{
    print = "Non va bene";
  }
  
  for (i in (1:length(file))) 
    {
    df = read.table(paste0(inpt, file[i]), sep=",", header=T);
    df$genotype = rep(geno[i], length(df[,1]));
    df.n = df[, c(3:4)];
    df$nearest = nndist(df.n);
    df$tec.rep = rep(tec[i], length(df[,1]));
    df$bio.rep = rep(bio.rep[i], length(df[,1]));
    df$source = rep(file[i], length(df[,1]));
    datalist[[i]] = df
  }
  
  all = do.call(rbind, datalist);
  write.table(all, paste0(out, exp, "_", max(tec), "_", max(bio.rep), ".txt"),
              sep = "\t", quote = F, col.names = T, row.names = F);
}


name = c("KF32", "GJV1") 
file = c("Image_5708.csv", "Image_5709.csv", "Image_5710.csv", "Image_5711.csv",
         "Image_5708.csv", "Image_5709.csv", "Image_5710.csv", "Image_5711.csv");
table(file, inpt, out, bio = 2, name, "exp")




