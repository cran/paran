`paran` <-
function(x, iterations=0, centile=0, has=FALSE, hasmean=FALSE, quietly=FALSE, status=TRUE) {

# quick validation of centile as an integer value
	centile <- round(centile)
	if (centile > 99 | centile < 0) {
		stop("\nYou must specify a centile value between 1 and 99.\n(Specifying centile 0 will use the mean.)")
		}

# Perform pca.

	pca <- princomp(x, cor = TRUE)

# Get the eigenvalues .  .  .
	Ev = pca[[1]]^2
	P <- length(Ev)

# clean up iteration and determine value
   if (iterations<1) {
		iterations = 30*P
		}
   if (iterations<0) {
		cat("\nInvalid number of iterations! Using default value of ",iterations,"\n",sep="")
		}

# prepare to save the results of each pca
		N <- length(as.matrix(x[1]))
		SimEvs <- matrix(NA,iterations,P)

# Let the user know the program is working if neccesary
	if (status==TRUE) {
		if (iterations >= 10) {
			cat("\nComputing")
			}
		}

		for (i in 1:iterations) {
   
# Yet _more_ letting the user know the program is working!
			if (status == TRUE) {
				if (i %% 10 == 1 & iterations >= 10) {
					disp <- ". "
					}
				else {
					disp <- "" 
					}
				if (i %% 300 != 1) {
					cat(disp)
					}
				else {
					cat("\n",disp)
					}
				}

# initialize previously created random dataset.
			Sim <- matrix(NA,N,P)
			
# Create the random dataset.
		# for normally distributed simulations
     if (has == FALSE & hasmean == FALSE) {
			Sim <- matrix(rnorm(N*P),N,P)
			}

		# as per Hayton, Allen and Scarpello
     if (has == TRUE & hasmean == FALSE) {
     	for (var in 1:P) {
     		MIN <- min(x[[var]])
     		MAX <- max(x[[var]])
     		MIDPOINT <- (MAX-MIN)/2
     		VAR <- var(x[[var]])
     		Sim[,var] <- rnorm(N,MIDPOINT,VAR)
     		passcount <- length(Sim[,var][Sim[,var] > MAX | Sim[,var] < MIN])
     		}
			while (passcount > 0) {
				Sim[,var][Sim[,var] > MAX | Sim[,var] < MIN] <- rnorm(passcount,MIDPOINT,VAR)
     		passcount <- length(Sim[,var][Sim[,var] > MAX | Sim[,var] < MIN])
     		}
			}

		# as per Hayton, Allen and Scarpello, but using the mean, not midpoint
     if (hasmean == TRUE) {
     	for (var in 1:P) {
     		MIN <- min(x[[var]])
     		MAX <- max(x[[var]])
     		MEAN <- mean(x[[var]])
     		VAR <- var(x[[var]])
     		Sim[,var] <- rnorm(N,MEAN,VAR)
     		passcount <- length(Sim[,var][Sim[,var] > MAX | Sim[,var] < MIN])
     		}
			while (passcount > 0) {
				Sim[,var][Sim[,var] > MAX | Sim[,var] < MIN] <- rnorm(passcount,MEAN,VAR)
     		passcount <- length(Sim[,var][Sim[,var] > MAX | Sim[,var] < MIN])
     		}
			}

# Run a principal components on the random dataset (which is 
# the same size and dimension as the user dataset.)

		pca <- princomp(Sim, cor = TRUE)

# Get the eigenvalues .  .  .
		Evs = pca[[1]]^2
# Save eigenvalues
		SimEvs[i,] <- Evs

# end the for i loop
		}


# display if neccesary
	if (quietly == TRUE) {
		cat("\n")
		}
	if (quietly == FALSE) {
	
		cat("\n\nResults of Horn's Parallel Analysis for principal components\n",sep="")

		if (iterations == 1) {
			if (centile == 0) {
				cat("1 iteration, using the mean estimate","\n",sep="")
				}
			if (centile != 0) {
				cat("1 iteration, using the p",centile," estimate","\n",sep="")
				}
			}

		if (iterations > 1) {
			if (centile == 0) {
				cat(iterations," iterations, using the mean estimate","\n",sep="")
				}
			if (centile != 0 & centile != 50) {
				cat(iterations," iterations, using the p",centile," estimate","\n",sep="")
				}
			if (centile == 50) {
				cat(iterations," iterations, using the p",centile," (median) estimate","\n",sep="")
				}		
			}

		cat("\n--------------------------------------------------","\n")
		cat("Component   Adjusted    Unadjusted    Estimated","\n")
		cat("or Factor   Eigenvalue  Eigenvalue    Bias","\n")
		cat("--------------------------------------------------","\n")
		}

	AdjEv = c(1:P)*NA 

	if (centile > 0) {
		for (p in 1:P) {
			AdjEv[[p]] <- quantile(SimEvs[,p],probs=centile/100)[[1]]
			}
		}
	if (centile==0) {
		for (p in 1:P) {
			AdjEv[[p]] <- mean(SimEvs[,p])			}
		}

	if (Ev[[1]] < 1 | AdjEv[[1]] < 1) { 
		if (quietly == FALSE) {
			cat("No components passed.","\n")
			cat("--------------------------------------------------","\n")
			stop
			}
		}

	y <- NA
	for (x in 1:P) {
		y <- x
		if (Ev[x] < 1 | AdjEv[x] < 1) {
			y <- x - 1
			break
			}
		}

	for (x in 1:y) {
		adjusted <- Ev[x]-AdjEv[x]
		if (Ev[x] >= 1) {
			if ( adjusted >=0 ) {
				AdjSpace = " "
				}
			if ( adjusted <0 ) {
				AdjSpace = ""
				}
			if ( Ev[x] >= 0 ) {
				EvSpace = " "
				}
			if ( Ev[x] < 0 ) {
				EvSpace = ""
				}
			if ( AdjEv[x] >= 0 ) {
				AdjEvSpace = " "
				}

# Pad the rear of x in case of single-digits
			if ( x > 9 ) {
				xPad = ""
				}
			if ( x <= 9 ) {
				xPad = " "
				}

# Pad the rear of adjusted in case of right-truncated zeros
			AdjPad = ""
			if ( (adjusted * 10^5) %% 100000 == 0 ) {
				AdjPad = "     "
				}
			if ( (adjusted * 10^5) %% 10000 == 0 ) {
				AdjPad = "    "
				}
			if ( (adjusted * 10^5) %% 1000 == 0 ) {
				AdjPad = "   "
				}
			if ( (adjusted * 10^5) %% 100 == 0 ) {
				AdjPad = "  "
				}
			if ( (adjusted * 10^5) %% 10 == 0 ) {
				AdjPad = " "
				}

# Pad the rear of Ev in case of right-truncated zeros
			EvPad = ""
			if ( (Ev[x] * 10^5) %% 100000 == 0 ) {
				EvPad = "     "
				}
			if ( (Ev[x] * 10^5) %% 10000 == 0 ) {
				EvPad = "    "
				}
			if ( (Ev[x] * 10^5) %% 1000 == 0 ) {
				EvPad = "   "
				}
			if ( (Ev[x] * 10^5) %% 100 == 0 ) {
				EvPad = "  "
				}
			if ( (Ev[x] * 10^5) %% 10 == 0 ) {
				EvPad = " "
				}

			if (quietly == FALSE) {
				cat(x,xPad,"         ",AdjSpace,round(adjusted,digits=6),AdjPad,"   ",EvSpace,round(Ev[x],digits=6),EvPad,"     ",AdjEvSpace,round(AdjEv[x],digits=6),"\n", sep="")
				}
			}
		}
	if (quietly == FALSE) {
		cat("--------------------------------------------------","\n")
		cat("Adjusted eigenvalues >= 1 indicate dimensions to retain.\n\n")
		}
	
	return(AdjEv)
}

