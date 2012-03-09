#setwd("C:/Documents and Settings/Philipp/Desktop/SimulationenGLMX/sim")

load(file = "results/ao2_sim.Robj")
ao2_sim <- output$output
load(file = "results/gosset_sim.Robj")
gosset_sim <- output$output
load(file = "results/talpha_sim.Robj")
talpha_sim <- output$output
rm(output)

library(xtable)



#### AO2 ####
ao2_sim_x  <- cbind(ao2_sim[,1], log(ao2_sim[,1]), ao2_sim[, c(2:3, 10:15)])

ao2_sim_x[,4] <- as.character(ao2_sim_x[,4]) 
ao2_sim_x[,4] <- gsub("2", "U(-4,4)", ao2_sim_x[,4]) 
ao2_sim_x[,4] <- gsub("1", "U(-2,2)", ao2_sim_x[,4]) 

colnames(ao2_sim_x) <- c("$\\phi$", "$\\log(\\phi)$", "$n$", "$X$",
                         "$  \\beta_0$", "  $\\beta_1$", "  $\\log(\\phi)$",
                         "$   \\beta_0$", "   $\\beta_1$", "   $\\log(\\phi)$")
                         
print(xtable(ao2_sim_x, digits = c(0,2,2,0,rep(2,7)),
             label = "ao2table", caption  = "Parameter Recovery for the Aranda-Ordaz II link"),
      type = "latex", caption.placement = "top",
      file = "latex/ao2table.tex", sanitize.text.function = function(x)x,
      booktabs = TRUE, include.rownames = FALSE)


#### Gosset ####
gosset_sim_x  <- cbind(gosset_sim[,1], log(gosset_sim[,1]), gosset_sim[, c(2:3, 10:15)])

gosset_sim_x[,4] <- as.character(gosset_sim_x[,4]) 
gosset_sim_x[,4] <- gsub("2", "U(-4,4)", gosset_sim_x[,4]) 
gosset_sim_x[,4] <- gsub("1", "U(-2,2)", gosset_sim_x[,4]) 

colnames(gosset_sim_x) <- c("$\\nu$", "$\\log(\\nu)$", "$n$", "$X$",
                         "$  \\beta_0$", "  $\\beta_1$", "  $\\log(\\nu)$",
                         "$   \\beta_0$", "   $\\beta_1$", "   $\\log(\\nu)$")
                         
print(xtable(gosset_sim_x, digits = c(0,2,2,0,rep(2,7)),
             label = "gossettable", caption  = "Parameter Recovery for the Gosset link"),
      type = "latex", caption.placement = "top",
      file = "latex/gossettable.tex", sanitize.text.function = function(x)x,
      booktabs = TRUE, include.rownames = FALSE)

### t_alpha ###
talpha_sim_x  <- cbind(talpha_sim[,1], log(talpha_sim[,1]), talpha_sim[, c(2:3, 10:15)])

talpha_sim_x[,1] <- 2* talpha_sim_x[,1]
talpha_sim_x[,4] <- as.character(talpha_sim_x[,4]) 
talpha_sim_x[,4] <- gsub("1", "U(-2,2)", talpha_sim_x[,4]) 
talpha_sim_x[,4] <- gsub("0.5", "U(-1,1)", talpha_sim_x[,4]) 

colnames(talpha_sim_x) <- c("$\\alpha$", "$\\log(\\alpha/2)$", "$n$", "$X$",
                         "$  \\beta_0$", "  $\\beta_1$", "  $\\log(\\alpha/2)$",
                         "$   \\beta_0$", "   $\\beta_1$", "   $\\log(\\alpha/2)$")
                         
print(xtable(talpha_sim_x, digits = c(0,2,2,0,rep(2,7)),
             label = "talphatable", caption  = "Parameter Recovery for the $t_\\alpha$ link"),
      type = "latex", caption.placement = "top",
      file = "latex/talphatable.tex", sanitize.text.function = function(x)x,
      booktabs = TRUE, include.rownames = FALSE)
