# Fit Finney data with Gosset link function.
# Liu (2004)  gets nuhat = 0.1 but we get something considerable larger

require(forward)
data(vaso)

source("Plinks.R")
source("pglm.R")

#Plot of the profile likelihood
if(TRUE){
  formula <- y ~ volume + rate
  fit <- pglm(formula,data = vaso, link="Gosset")
  objhat <- fit$f$deviance
  nus <- 2:40/29
  objs <- nus*0
  for(i in 1:length(nus)){
           objs[i] <-glm(formula, data = vaso,
                  family=binomial(link=Gosset(nus[i])))$deviance
        }
  }
#pdf("Fig1.pdf", width = 6, height = 4)
plot(nus, -objs, cex = .5, xlab = expression(nu), ylab = "2 log Likelihood")
lines(nus, -objs, lwd = .5, lty=1, col="grey")
cval <- qchisq(.95,1)
#segments(fit$nulo, -objhat - cval, fit$nuhi, -objhat - cval)
#segments(fit$nulo, -743, fit$nulo, -objhat - cval,2)
#segments(fit$nuhi, -743, fit$nuhi, -objhat - cval,2)
#dev.off()



#MLE estimate of nu

if(TRUE){
formula <- y ~ volume + rate
f0 <- glm(formula, family=binomial(link="probit"),data=vaso)
f1 <- pglm(formula,link="Gosset",data=vaso)
}

