library(DirichletReg)
library(matrixStats)
# library(matlib)

set.seed(1)
MCMC<- function(x, K, lambda){
	n <- length(x)
	
	#Initialize the parameters
	mu <- rep(mean(x), K)
	sigma <- 1
	z <- rep(0, n)
	pis <- rep(1/K, K)
	mu1_stored <- c()
	mu2_stored <- c()
	mu1_stored <- append(mu1_stored, mu[1])
	mu2_stored <- append(mu2_stored, mu[2])


	iter <- 1e4

	while(iter){
		probs <- matrix(1/K, n, K)

		# probs update
		for(i in 1:K){
			probs[, i] <- pis[i] * dnorm(x, mu[i], sigma)
		}

		# normalization
		probs <- probs/rowSums(probs)


		# z update
		for(a in 1:n){
			z[a] <- sample(K, 1, replace = TRUE, probs[a, ])
		}

		# mu update
		for(b in 1:K){

		    Nk <- 0
		    Xk <- 0

		    for(v in 1:n){
		    	if(z[v]==b){
		    		Nk <- Nk + 1
		    		Xk <- Xk + x[v]
		    	}
		    }

			Xk <- Xk/Nk
			muk <- ((Nk/(sigma**2))/(Nk/(sigma**2) + 1/(lambda[b]**2))) * Xk
			lam <- 1/(Nk/(sigma**2) + 1/(lambda[b]**2))
			mu[b] <- rnorm(1, muk, sqrt(lam))
		}
		
		mu1_stored <- append(mu1_stored, mu[1])
		mu2_stored <- append(mu2_stored, mu[2])
		iter <- iter - 1
	}

	# return the parameters calculated corresponding to the dataset 
	parameters <- list("m" = mu,"mu1_stored" = mu1_stored,"mu2_stored" = mu2_stored)
	return(parameters)
}

make_plot <- function(sigma, x, mu1_stored, mu2_stored, mu){
	plot(density(x), main = sigma, ylim = c(0, 6), xlim = range(mu[1],mu[2], x))
	lines(density(mu1_stored), col = "red")
	lines(density(mu2_stored), col = "green")
}

par(mfrow = c(2,2))

x <- faithful$eruptions

K <- 2

sigma <- 0.78
p <- MCMC(x, K, sigma)
mu <- p$mu
mu1_stored <- p$mu1_stored
mu2_stored <- p$mu2_stored
make_plot(sigma, x, mu1_stored, mu2_stored, mu)

sigma <- 1.5
p <- MCMC(x, K, sigma)
mu <- p$mu
mu1_stored <- p$mu1_stored
mu2_stored <- p$mu2_stored
make_plot(sigma, x, mu1_stored, mu2_stored, mu)

sigma <- 2
p <- MCMC(x, K, sigma)
mu <- p$mu
mu1_stored <- p$mu1_stored
mu2_stored <- p$mu2_stored
make_plot(sigma, x, mu1_stored, mu2_stored, mu)

sigma <- var(x)
p <- MCMC(x, K, sigma)
mu <- p$mu
mu1_stored <- p$mu1_stored
mu2_stored <- p$mu2_stored
make_plot(sigma, x, mu1_stored, mu2_stored, mu)


# plot(x, col=z)

# iterations <- c(1:100)

# plot(iterations, mu1_stored, xlab = "iter", ylab = "mu[1]", type='o')
# plot(iterations, mu2_stored)







