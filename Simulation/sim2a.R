
# sim2 analysis
# y <- x * beta + (1 + x * gamma) * rlogis(n)

load("sim2.Rda")
require(Hmisc)

beta <- 3
gamma <- .25 
taus <- 1:4/5
ns <- c(500,1000,5000)
qtaus <- qlogis(taus)
B <- rbind(qtaus, beta + gamma * qtaus)  %o% outer(rep(1,3),rep(1,R))
D <- A[1:2,,,] - B
Bias <- apply(D,1:3,mean)
RMSE <- sqrt(apply(D^2,1:3,mean))

# Now for the hard part -- what about the Fhat estimates?
# Note that I've recentered the intercepts in rdx so at each of the 4 cutoffs
# we are estimating (beta_0, beta_1, gamma) = (0,3,0.25)
D <- A[5:7,,,] - c(0, 3, 0.25)
MeanBias <- matrix(aperm(apply(D,1:3,mean),c(1,3,2)),9,4)
MedBias <- matrix(aperm(apply(D,1:3,median),c(1,3,2)),9,4)
RMSE <- matrix(sqrt(aperm(apply(D^2,1:3,mean),c(1,3,2))),9,4)
MedAE <- matrix(aperm(apply(abs(D),1:3,median),c(1,3,2)),9,4)
Tab1 <- cbind(MedBias, MedAE)
rnames <- rep(c("$\\beta_0$", "$\\beta_1$", "$\\gamma$"),3)
cnames <- rep(paste("$\\tau = $",round(1:4/5,1)),2)
dimnames(Tab1) <- list(rnames, cnames)
rgrp <- paste("n = ",ns)
cgrp <- c("Median Bias","Median Absolute Error")
cap1 <- "Median Bias and Absolute Error Performance for Linear Scale Logit Model"
latex(Tab1, file = "tab1.tex",rowlabel = "", cgroup = cgrp, rgroup = rgrp, dec = 4,
   label = "tab1", caption.loc = "bottom", caption = cap1)
Tab2 <- cbind(MeanBias,  RMSE)
dimnames(Tab2) <- list(rnames, cnames)
cgrp <- c("Mean Bias","Root Mean Squared Error")
cap2 <- "Mean Bias and Root Mean Squared Error Performance for Linear Scale Logit Model"
latex(Tab2, file = "tab2.tex",rowlabel = "", cgroup = cgrp, rgroup = rgrp, dec = 2,
   label = "tab2", caption.loc = "bottom", caption = cap2)

