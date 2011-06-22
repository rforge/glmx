# Prentice (Biometrics, 1976, 761-68) Link function

dprentice <- function(x, a=1, b=1)
	exp(x*a) * (1 + exp(x))^(-(a + b))/beta(a,b)

# Base Graphics Version wastes too much space with margins
par(mfrow=c(4,4))
x <- -100:100/10
for(a in 1:4/2){
  for(b in 1:4/2){
	plot(x,dprentice(x,a,b),type="l",ylab="f(x)",ylim=c(0,.4))
	}
   }

#Lattice Version of the Plot

require(lattice)
ab <- as.matrix(expand.grid(a=1:4/2,b=1:4/2))
n <- length(x)
AB <- ab %x% rep(1,n)
Y  <- matrix(0,n,nrow(ab))
X  <- rep(1,nrow(ab)) %x% x
for(k in 1:nrow(ab)){
	Y[,k] <- dprentice(x,ab[k,1],ab[k,2])
	}
D <- cbind(AB,X,c(Y))
dimnames(D)[[2]] <- c("a","b","x","y")
D <- as.data.frame(D)
xyp <- xyplot(y~x|a+b,data=D,type="l")
	

