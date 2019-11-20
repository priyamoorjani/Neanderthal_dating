#$ -cwd
options(echo=TRUE) # if you want see commands in output file
args <- commandArgs(trailingOnly = TRUE)
num <- args[1]

# Importance sampling

#num <- sample(1:110000,1)
#num <- 1
set.seed(num)
outlog <- paste("run_", num,".txt", sep="")
print(num)

# variables
sample <- c("CEU", "Clovis", "MA1","Kostenki14","Ust_Ishim_10cM", "Oase") 
radiocarbon_min <- c(0, 12509, 23891, 36262,  43210, 37000)
radiocarbon_max  <- c(1,12722, 24423, 38684, 46880,42000)
radiocarbon_mean <- c(0, 12661,24157,37450,45045,39500)
radiocarbon_error <- (radiocarbon_max  - radiocarbon_mean )/ 1.96
date_neander_mean <- c(1633, 1009,  771,209,92,8)
date_neander_error <- c(41,35,  27,13,1,1)
date_neander_min <- date_neander_mean - 1.96*date_neander_error
date_neander_max <- date_neander_mean + 1.96*date_neander_error

# variables
y <- radiocarbon_mean
x <- date_neander_mean
x.se <- date_neander_error
y.se <- radiocarbon_error

# functions
lik <- function(y,z,beta,theta) {		
	prod <- 1
	sigma <- sqrt((beta*x.se)^2 + y.se^2 + theta^2)
	prod <- prod(((1/(sigma))*exp(-0.5*(((y-z)/sigma)^2))))
	return(prod)	
}

nrun <- 100000
#nrun <- 3
# compute likelihood
beta <- NA; intercept <- NA; theta <- NA; z <- NA; weight <- NA
for (iter in 1:nrun) {
	# sample beta, intercept, theta
	beta[iter] <- runif(1,25,33) # sample from uniform dist [25,33]
	#print(beta[iter])
	intercept[iter] <- rnorm(1,mean=47236, sd=4339 ) # sample from normal dist within the 95% CI of dates of neandertal admixture in CEU
	theta[iter] <- rnorm(1,0,3)  # sample from normal N(0,1) -- change as needed
	
	# computer z = beta*x + intercept
	z <- intercept[iter] - beta[iter]*x  
	
	# compute weight
	weight[iter] <- lik(y,z,beta[iter],theta[iter])	
}

# update weights
weight <- weight /max(weight)

outdata <- data.frame(weight=weight, beta=beta, intercept=intercept, theta=theta)
# write output
write.table(outdata, outlog, sep="\t", col.names=T, row.names=F, append=F, quote=F)

