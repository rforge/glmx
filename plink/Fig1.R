#This is figure 1 of the Rnews paper

d <- read.table("http://www.econ.uiuc.edu/~roger/courses/471/data/weco.dat", header=TRUE)

source("Plinks.R")
source("pglm.R")
if(TRUE){
  formula <- kwit ~ sex + dex + poly(lex,2)
  fit <- pglm(formula,data = d, link="Gosset")
  objhat <- fit$f$deviance
  nus <- 4:40/19
  objs <- nus*0
  for(i in 1:length(nus)){
	   objs[i] <-glm(formula, data = d,
		  family=binomial(link=Gosset(nus[i])))$deviance
	}
  } 
pdf("Fig1.pdf", width = 6, height = 4)
plot(nus, -objs, cex = .5, xlab = expression(nu), ylab = "2 log Likelihood")
lines(nus, -objs, lwd = .5, lty=1, col="grey")
cval <- qchisq(.95,1)
segments(fit$nulo, -objhat - cval, fit$nuhi, -objhat - cval)
segments(fit$nulo, -743, fit$nulo, -objhat - cval,2)
segments(fit$nuhi, -743, fit$nuhi, -objhat - cval,2)
dev.off()
