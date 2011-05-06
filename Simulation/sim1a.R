
load("sim1.Rda")


beta <- 3
gamma <- .25 
R <- 1000
taus <- 1:4/5
ns <- c(500,1000,5000)
A <- array(0,c(7,4,length(ns),R))
#B <- rbind(qlogis(taus),rep(1,4))
for(j in 1:length(ns)){
   n <- ns[j]
   x <- runif(n,1,10)
   for(i in 1:R){
      y <- x * beta + (1 + x * gamma) * rlogis(n)
      yc <- quantile(y, taus)
      b0 <- rq(y ~ x, tau = taus)$coef 
      b1 <- rdf(x,y,yc) 
      b2 <- rdx(x,y,yc)
      A[,,j,i] <- rbind(b0,b1,b2)
	print(c(i,j))
      }
   }
