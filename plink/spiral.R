spiral <- function(n){
# Rectangular spiral indexing

# Translated from the Matlab 
s <- matrix(0,n,n)
i <- ceiling(n/2)
j <- ceiling(n/2)
s[i,j] <- 1
if(n == 1) return(s)
k <- 1  # Numbering.
d <- 1  # Increasing or decreasing.
for(p in 1:n){
   q = 1:min(p,n-1)  
   j = j+d*q
   k = k+q
   s[i,j] = k
   if(p == n) return(s)
   j = j[p]
   k = k[p]
   i = i+d*t(q)
   k = k+t(q)
   s[i,j] = k
   i = i[p]
   k = k[p]
   d = -d
  }
}
