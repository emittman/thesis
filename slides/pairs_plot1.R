load(file = "heterosis_counts.RData")
load("heterosis_design.RData")
library(limma)
library(edgeR)
library(dplyr)
nonzero_id <- which(!(rowSums(my_dat_wide[,-1]) < 8))
my_dat_wide = my_dat_wide[nonzero_id,]
voom_dat <- voom(counts = my_dat_wide[,-1], design = X, plot=TRUE)
str(voom_dat)
fit <- lmFit(voom_dat)
library(GGally)
my_hex <- function(data, mapping, ...) {
ggplot(data = data, mapping=mapping) +
geom_hex(..., bins=35)+
scale_fill_continuous(trans="log",low="white", high="#C8102E")+
    theme_bw()+theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1))
}

my_dens <- function(data, mapping, ...) {
  ggplot(data=data, mapping=mapping) +
    geom_density(..., fill="#C8102E", alpha=.5)+
    theme_bw()+theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1))

}

library(tidyr)
fit.dat <- data.frame(data.matrix(fit$coefficients))
names(fit.dat) <- c("intercept", "parental HD", "hybrid", "hybrid HD", "flow cell")
jpeg(filename = "pairs1.jpg", width=14.3, height=11, units="cm", res = 200)
ggpairs(fit.dat, lower=list(continuous=my_hex),
        upper=list(continuous=my_hex),
        diag=list(continuous=my_dens))
dev.off()
