d <- read.table("http://www.econ.uiuc.edu/~roger/courses/471/data/weco.dat", header=TRUE)
library(gld)
source("Plinks.R")
source("pglm.R")
source("spiral.R")

formula <- kwit ~ sex + dex + poly(lex,2)
f <- pglm(formula , data = d, link = "Pregibon", theta = c(-2.5,.2))
start <- f$f$coef
devstar <-  f$f$deviance
theta <- f$theta
m  <- 40
as <- theta[1]  + seq(-2.5,3.5,length=m)
bs <- theta[2] + seq(-1.5,1,length=m)

A <- outer(as,bs)
S <- spiral(m)
for(k in 1:(m^2)){ #Spiral indexing to exploit continuity for glm initial values.
	ij <- which(S == k,arr.ind = TRUE)
	i <- ij[1]
	j <- ij[2]
	f <- glm(formula,family=binomial(link=Pregibon(a=as[i],b=bs[j])),
			start = start, data = d)
        A[i,j] <- f$deviance
	start <- f$coef
        }

pdf("Fig3.pdf", width = 6, height = 6)
#rescale contour values to probabilities
B <- pchisq(A - devstar, 2)
#Clip plotting region to hide ugly NW corner where convergence of glm is problematic
plot(as,bs,xlim = c(-5.2,0),ylim = c(-1,1.25), 
	xlab = expression(alpha), ylab = expression(delta), type="n")
Cf <- contour(as[6:40],bs[1:35], B[6:40,1:35], cex = .5, 
	levels = c(.1,.2,.5,.6,.7,.8,.9),add=TRUE)
points(theta[1],theta[2])
dev.off()
