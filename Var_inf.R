# Please run library(DirichletReg) before running the code 
# Copy and Paste this code into the R console and run the Master function to obtain the desired plots
#install.packages("DirichletReg")

library(DirichletReg)
library(matrixStats)
i <- 1
while(1){
	set.seed(i)
	par(mfrow = c(2,2))

	x <- faithful$waiting
	K <- 2

	sigma <- 0.05
	p1 <- CAVI(x, K, sigma)
	m <- p1$m
	Plot_comparision(x, p1, sigma)

	sigma <- 2
	p2 <- CAVI(x, K, sigma)
	m <- p2$m
	Plot_comparision(x, p2, sigma)

	sigma <- 150
	p3 <- CAVI(x, K, sigma)
	m <- p3$m
	Plot_comparision(x, p3, sigma)

	sigma <- var(x)
	p4 <- CAVI(x, K, sigma)
	m <- p4$m
	Plot_comparision(x, p4, sigma)
	i <- i+1
}
# function to calculate ELBO
my_elbo <- function(x, m, s2, phi, sigma){
	p <- sum(log(s2) - m/(sigma**2))
	q <- (-0.5) *(outer(x**2, m**2 + s2, "+"))
	q <- q + outer(x, m)
	q <- q	 - log(phi)
	q <- sum(q * phi)
	return(p + q)
}

# x: data; K: the desired no. of clusters; sigma: hyperparameter
CAVI <- function(x, K, sigma){ 
	n <- length(x) 
	
	# initialize m, s2 and phi
	m <- runif(K, min(x), max(x)) 
	s2 <- 1
	# runif(K)
	categorical <- rep(1/K, K)
	phi <- rdirichlet(n, categorical)
	
	# calculate the initial elbo
	elbo <- my_elbo(x, m, s2, phi, sigma)
	
	# run the loop up until when the difference in values remain very small
	prev<-1000
	eps <- 1e-6
	flag <- 1
	
	while(flag){
		
		# cavi "phi" update
		e1 <- outer(x, m)
		e2 <- -0.5 * (m**2 + s2)
		e <- t(t(e1) + e2)
		phi <- e - rowMaxs(e) - log(rowSums(exp(e - rowMaxs(e))))
		phi <- exp(phi)
		# phi <- exp(e)/(rowSums(exp(e)))
		
		# cavi "m" update
		m <- colSums(x * phi)
		m <- m / ((1.0/(sigma**2))+ (colSums(phi)))
		
		#cavi "s2" update
		s2 <- 1.0/((1.0/(sigma**2)) + (colSums(phi)))
			
		# check the elbo for termination condition		
		elbo <- my_elbo(x, m, s2, phi, sigma)
		if(abs(prev - elbo)<=eps) flag <- 0
		prev <- elbo
	}
	
	# return the parameters calculated corresponding to the dataset 
	parameters <- list("m" = m,"s2" = s2,"phi" = phi)
	return(parameters)
}

Plot_comparision <- function(x, p, sigma){
	m1 <- p$m[1]
	m2 <- p$m[2]
	s1 <- p$s2[1]
	s2 <- p$s2[2]
	plot(density(x), main = sigma, ylim = c(0, 0.07), xlim = range(m1,m2, x))
	lines(density(rnorm(10000, m1, sqrt(s1))), col = "red")
	lines(density(rnorm(10000, m2, sqrt(s2))), col = "green")
}









