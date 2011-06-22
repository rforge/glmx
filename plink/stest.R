m <- 10
as <- seq(-1,1,length=m)
bs <- seq(-1,1,length=m)
plot(as,bs,type="n")

A <- outer(as,bs) 
S <- spiral(m)
for(k in 1:(m^2)){
        ij <- which(S == k,arr.ind = TRUE)
        i <- ij[1]
        j <- ij[2]
	text(as[i],bs[j],paste(k))
        }

