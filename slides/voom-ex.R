setwd("chapter1/R/data/")

load(file = "heterosis_counts.RData")
load("heterosis_design.RData")
library(limma)
library(edgeR)
library(dplyr)


nonzero_id <- which(!(rowSums(my_dat_wide[,-1]) < 8))

my_dat_wide = my_dat_wide[nonzero_id,]

voom_dat <- voom(counts = my_dat_wide[,-1], design = X, plot=TRUE)

R <- as.numeric(colSums(my_dat_wide[,-1]))

log_cpm <- t(apply(my_dat_wide[,-1], 1, function(row){
  log((row + 0.5)/(R + 1) * 1e6, base=2)
}))

ols_fits <- apply(log_cpm, 1, function(y){
  fit <- lm(y ~ 0 + X)
  list(beta = coef(fit), sigma = sigma(fit))
})

mus <- sapply(ols_fits, function(e) X %*% e$beta)
sigmas <- sapply(ols_fits, function(e) e$sigma)
avg_log_cpm <- drop(apply(log_cpm, 1, mean))
Rtilde <- exp(mean(log(R+1)))
avg_log_count <- avg_log_cpm + log(Rtilde+1, base=2) - log(1e6, base=2)

lo <- lowess(avg_log_count, sqrt(sigmas), f=.5)
lo_m <- approxfun(lo, rule=2)
lo_l <- approxfun(with(lo, list(x=x,y=.5*y)), rule=2)
lo_h <- approxfun(with(lo, list(x=x,y=1.5*y)), rule=2)

#plot loess fit
library(ggplot2)
ord <- order(avg_log_count)
p <- data.frame(avg_log_count = avg_log_count[ord],
           sqrt_sigma = sqrt(sigmas)[ord],
           lo_pred = lo_m(avg_log_count)[ord],
           lo_pred2 = lo_l(avg_log_count)[ord],
           lo_pred3 = lo_h(avg_log_count)[ord]) %>%
  ggplot(aes(x=avg_log_count, y = sqrt_sigma)) + geom_point(pch=".") +
  theme_bw()+xlab(expression(paste(tilde(r)[g]," (average log count)"))) +
  ylab(expression(sqrt(s)[g]))
p
ggsave("voom1.pdf",width=6,height=4)
p1 <- p +   geom_line(aes(y = lo_pred), color = "red") +
annotate(geom="text", x=13, y=.7, label="lowess trend")+
  geom_segment(x=13,y=.6,xend=13,yend=.32,arrow=arrow(length=unit(.2,"cm")))
p1
ggsave("voom2.pdf",width=6,height=4)

p2 <- p1 + geom_line(aes(y = lo_pred2), color = "blue", lty=2) +
  geom_line(aes(y = lo_pred3), color = "blue", lty=2)

library(cowplot)
plot_grid(p,p1,p2,ncol=1)
data.frame(fitted=as.numeric(t(fitted_counts[sample(1:10000,10),])))

fitted_counts <- t(mus + log(rep(R,times=ncol(mus)) + 1, base=2) - log(1e6, base=2))
est_sqsd <- t(apply(fitted_counts, 1, function(x) lo(x)))
prec_weights <- matrix(est_sqsd,
                       nrow(fitted_counts),
                       ncol(fitted_counts)) ^ -4

plot(voom_dat$weights, prec_weights)
setwd("../../../slides/")
save(ols_fits, file="ols_fits.RData")
