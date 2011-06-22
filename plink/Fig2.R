# This produces Figure 2 of the Rnews paper

d <- read.table("http://www.econ.uiuc.edu/~roger/courses/471/data/weco.dat", header=TRUE)

source("Plinks.R")
source("pglm.R")

if(TRUE){
  formula <- kwit ~ sex + dex + poly(lex,2)
  f0<- glm(formula, family=binomial(link="probit"),data=d)
  f1 <- pglm(formula,link="Gosset",data=d)
  } 
fv1 <- f1$f$fitted.value
fv0 <- f0$fitted.value
pdf("Fig2.pdf",height=6,width=6)

plot(0:1,0:1,type="n", xlab = "Estimated Probit Probabilities",
	ylab = "Estimated Gosset Probabilities")
points(fv0,fv1)
abline(c(0,1))
dev.off()
