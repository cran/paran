`paran` <-
function(x, iterations=0, centile=0, quietly=FALSE, status=TRUE, all=FALSE, cfa=FALSE, graph=FALSE, color=TRUE, col=c("black","red","blue"), lty=c(1,2,3), lwd=1, file="", width=640, height=640) {

library(MASS)

# x are the data
	x <- data.frame(x)

# quick validation of centile as an integer value
	centile <- round(centile)
	if (centile > 99 | centile < 0) {
		stop("\nYou must specify a centile value between 1 and 99.\n(Specifying centile 0 will use the mean.)")
		}

	P <- length(x)

# Perform pca or cfa
	if (cfa == FALSE) {
		eigenvalues <- eigen(cor(x), only.values = TRUE, EISPACK = FALSE)[[1]]
		}
	if (cfa == TRUE) {
		C <- cor(x)
		C <- C - ginv(diag(diag(ginv(C))))
		eigenvalues <- eigen(C, only.values = TRUE, EISPACK = FALSE)[[1]]
		}
	
# Get the eigenvalues .  .  .
	Ev <- eigenvalues

# note which model
	model <- "component"
	models <- "components"
	Model <- "Component"
	Models <- "Components"
	if (cfa == TRUE) {
		model <- "factor"
		models <- "factors   "
		Model <- "Factor   "
		Models <- "Factors"
		}

# clean up iteration and determine value
   if (iterations<1) {
		iterations <- 30*P
		}
   if (iterations<0) {
		cat("\nInvalid number of iterations! Using default value of ",iterations,"\n",sep="")
		}

# prepare to save the results of each pca
		N <- length(as.matrix(x[1]))
		if ( cfa == FALSE ) {
			SimEvs <- matrix(NA,iterations,P)
			}
		if ( cfa == TRUE ) {
			SimEvs <- matrix(NA,iterations,P)
			}

# Let the user know the program is working if neccesary
	if (status==TRUE) {
		if (iterations >= 10) {
			cat("\nComputing: ")
			}
		}

	set.seed(1000*proc.time()[[1]])
	for (k in 1:iterations) {
   
# Yet _more_ letting the user know the program is working!
		if (status == TRUE) {
			if (k %% (iterations/10) == 1 & iterations >= 10 & k > 1) {
				pct <- (k%/%(iterations/10))*10
				cat(pct,"%  ",sep="")
				}
			if (k == iterations) {
				cat("100%\n")
				}
			}

# initialize previously created random dataset.
		Sim <- matrix(NA,N,P)
			
# Create the random dataset.
		# for normally distributed simulations
		Sim <- matrix(rnorm(N*P),N,P)

# Run a principal components or factor analysis on the random dataset
# (which is the same size and dimension as the user dataset.)

		if (cfa == FALSE) {
			eigenvalues <- eigen(cor(Sim), only.values = TRUE, EISPACK = FALSE)[[1]]
			}
		if (cfa == TRUE) {
			C <- cor(Sim)
			C <- C - ginv(diag(diag(ginv(C))))
			eigenvalues <- eigen(C, only.values = TRUE, EISPACK = FALSE)[[1]]
			}

# Get the eigenvalues .  .  .
		Evs <- eigenvalues

# Save eigenvalues
		SimEvs[k,] <- Evs

# reseed and end the for k loop
	Sys.sleep(.002)
	set.seed(1000*proc.time()[[1]])
	}

# display if neccesary
	if (quietly == TRUE) {
		cat("\n")
		}
	if (quietly == FALSE) {
	
		cat("\n\nResults of Horn's Parallel Analysis for ",model," retention\n",sep="")

		if (iterations == 1) {
			if (centile == 0) {
				cat("1 iteration, using the mean estimate","\n",sep="")
				}
			if (centile != 0) {
				cat("1 iteration, using the ",centile," centile estimate","\n",sep="")
				}
			}

		if (iterations > 1) {
			if (centile == 0) {
				cat(iterations," iterations, using the mean estimate","\n",sep="")
				}
			if (centile != 0 & centile != 50) {
				cat(iterations," iterations, using the ",centile," centile estimate","\n",sep="")
				}
			if (centile == 50) {
				cat(iterations," iterations, using the ",centile," centile (median) estimate","\n",sep="")
				}		
			}

		cat("\n--------------------------------------------------","\n")
		cat(Model,"  Adjusted    Unadjusted    Estimated","\n")
		cat("            Eigenvalue  Eigenvalue    Bias","\n")
		cat("--------------------------------------------------","\n")
		}

	RndEv = c(1:P)*NA 

	if (centile > 0) {
		for (p in 1:P) {
			RndEv[[p]] <- quantile(SimEvs[,p],probs=centile/100)[[1]]
			}
		}
	if (centile==0) {
		for (p in 1:P) {
			RndEv[[p]] <- mean(SimEvs[,p])			}
		}

	if (Ev[[1]] < 1 | RndEv[[1]] < 1) { 
		if (quietly == FALSE) {
			cat("No components passed.","\n")
			cat("--------------------------------------------------","\n")
			stop
			}
		}

	Bias <- rep(0,P)
	AdjEv <- rep(1,P)
	for (p in 1:P) {
		Bias[p] <- RndEv[p] - 1
		AdjEv[p] <- Ev[[p]] - Bias[p]
		}

	# calculate how many components or factors to return by counting those 
	# components or factors with adjusted eigenvalues greater than one
	y <- NA
	for (x in 1:P) {
		y <- x
		if (AdjEv[x] <= 1) {
			y <- x - 1
			retained <- y
			break
			}
		}

	if ( all == TRUE ) {
		y <- P
		}

	for (x in 1:y) {
		if ( AdjEv[x] >=0 ) {
			AdjSpace = " "
			}
		if ( AdjEv[x] < 0 ) {
			AdjSpace = ""
			}
		if ( Ev[[x]] >= 0 ) {
			EvSpace = " "
			}
		if ( Ev[[x]] < 0 ) {
			EvSpace = ""
			}
		if ( Bias[x] >= 0 ) {
			BiasSpace = " "
			}
		if ( Bias[x] < 0 ) {
			BiasSpace = ""
			}

# Pad the rear of x in case of single-digits
		if ( x > 9 ) {
			xPad = ""
			}
		if ( x <= 9 ) {
			xPad = " "
			}

# Pad the front of AdjEv in case of eigenvalues > 10, 100, etc.
		AdjFPad = "   "
		if ( round(AdjEv[x]) >= 10 ) {
			AdjFPad = "  "
			}
		if ( round(AdjEv[x]) >= 100 ) {
			AdjFPad <- " "
			}

# Set the strtrim number SN
		SN <- 8
		if ( abs(AdjEv[x]) >= 10 ) {
			SN <- 9
			}
		if ( abs(AdjEv[x]) >= 100 ) {
			SN >= 10
			}
		if ( AdjEv[x] < 0 ) {
			SN <- SN + 1
			}

# Pad the front of Ev in case of eigenvalues > 10, 100, etc.
		EvFPad = "   "
		if ( round(Ev[[x]]) >= 10 ) {
			EvFPad = "  "
			}
		if ( round(Ev[[x]]) >= 100 ) {
			EvFPad = " "
			}

# Set the strtrim number SN
		EvSN <- 8
		if ( abs(Ev[[x]]) >= 10 ) {
			EvSN <- 9
			}
		if ( abs(Ev[[x]]) >= 100 ) {
			EvSN <- 10
			}
		if (abs(Ev[[x]]) >= .0000005) {
			EvZPad <- ""
			}
		if (abs(Ev[[x]]) < .0000005) {
			Ev[[x]] <- 0
			EvZPad <- ".000000"
			}

# Set the strtrim number SN
		BiasSN <- 8
		if ( Bias[x] >= 10 ) {
			BiasSN <- 9
			}
		if ( Bias[x] >= 100 ) {
			BiasSN >= 10
			}

		if (quietly == FALSE) {
			cat(x,xPad,"      ",AdjFPad,AdjSpace,strtrim(AdjEv[x],SN),EvFPad,EvSpace,strtrim(Ev[[x]],EvSN),EvZPad,"     ",BiasSpace,strtrim(Bias[x],BiasSN),"\n", sep="")
			}
		}
	if (quietly == FALSE) {
		cat("--------------------------------------------------","\n")
		cat("\nAdjusted eigenvalues > 1 indicate dimensions to retain.\n(",retained," ",models," retained)\n\n",sep="")
		}

# Graph it if needed
	if (graph == TRUE) {
		AdjEvCol = col[1]
		EvCol = col[2]
		RndEvCol = col[3]
		AdjEvLty = 1
		EvLty = 1
		RndEvLty = 1
		if (color == FALSE) {
			EvCol = "black"
			RndEvCol = "black"
			EvLty = lty[2]
			RndEvLty = lty[3]
			}
		if (cfa==FALSE) {
			par(yaxs='i', xaxs='i', lab=c(P,ceiling(Ev[1]),2))
			plot.default(c(1:P), RndEv, type='o', main='Parallel Analysis', xlab='Components', ylab='Eigenvalues', pch=20, col=RndEvCol, lty=RndEvLty, lwd=lwd, xlim=c(.5,P+.5), ylim=c(-.5,ceiling(Ev[[1]])))
			}
		if (cfa==TRUE) {
			par(xaxp=c(1,P,1))
			plot.default(c(1:P), RndEv, type='o', main='Parallel Analysis', xlab='Factors', ylab='Eigenvalues', pch=20, col=RndEvCol, lty=RndEvLty, lwd=lwd, xlim=c(.5,P+.5), ylim=c(-.5,ceiling(Ev[[1]])))
			}
		abline(h=1, col='lightgrey',lwd=.5)
		points(c(1:P),AdjEv, type='o', col=AdjEvCol, lty=AdjEvLty, pch=21, bg='white', lwd=lwd)
		points(c(1:P),Ev, type='o', col=EvCol, lty=EvLty, pch=20, lwd=lwd)
		points(c(1:retained), AdjEv[1:retained], type='p', pch=19, col=AdjEvCol, lty=AdjEvLty, lwd=lwd)
		if (file != "" & typeof(file) == "character") {
			dev.copy(png, height=height, width=width, file=file)
			dev.off()
			}
		}

	names(retained) <- c(paste("Retained ",Models,sep=""))
	names(AdjEv) <- c("Adjusted Eigenvalues")
	names(Ev) <- c("Unadjusted Eigenvalues")
	names(RndEv) <- c("Random Eigenvalues")
	names(Bias) <- c("Estimated Bias")
	return(list(retained,AdjEv,Ev,RndEv,Bias,SimEvs))
}

