load("bad.Rda")
plot(x,y)
s <- (y > yc[4])
points(x[s],y[s],col = "red")
f <- hetglm(s ~ x,trace = 10)
z <- sort(x)
o <- order(x)

# Coefficients from trace 
B <- matrix(c(
-5.17324,0.634122,0.00000,
-88.0924,11.5437,0.393495,
-1550.74,212.312,0.773211,
-7608.71,1057.99,0.976548,
-37027.7,5204.26,1.19014,
-196213.,27809.8,1.41494,
-1.09667e+06,156430.,1.65158,
-1.42425e+07,2.04542e+06,2.01629,
-7.47306e+07,1.07710e+07,2.23894
),3,9)

F <- matrix(0,length(z),9)
for(i in 1:9)
    F[,i] <- plogis((cbind(1,z) %*% B[1:2,i])/exp(z*B[3,i]))
matplot(z,F,type = "l")
points(z,.6*s[o])
