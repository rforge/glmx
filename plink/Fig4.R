d <- read.table("http://www.econ.uiuc.edu/~roger/courses/471/data/weco.dat", header=TRUE)
library(gld)
source("Plinks.R")
source("pglm.R")

formula <- kwit ~ sex + dex + poly(lex,2)
f <- pglm(formula , data = d, link = "Pregibon", theta = c(-2.5,.2))
start <- f$f$coef
devstar <-  f$f$deviance
theta <- f$theta

pdf("Fig4.pdf", width = 6, height = 4)
u <- (-100:100)/10
d <- dPregibon(u,theta[1],theta[2])
d0 <- dPregibon(u)
plot(u,d,type="l",xlab="",ylab="")
lines(u,d0,col="grey")
dev.off()
