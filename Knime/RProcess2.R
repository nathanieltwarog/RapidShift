# Node 12: Table to R
# Read in raw data that has been selected for refitting by user
rawdata <- knime.in

# Node 13: Add Table to R
# Read in upper and lower temperature bounds provided by user
welldata <- knime.in

# Node 23: R to R
# Again, the main processing block.  Fits data to Boltzmann equation
# using non-linear optimization, excluding supramaximum measurements.
# If fit is poor or invalid, refits data, also exlcuding subminimum
# measurements.  Fit parameters, metrics for goodness of fit, and fit
# flags are stored in 'welldata' data frame.
fitBoltzmann <- function(temp,effect,def=NULL) {
	if (length(temp)!=length(effect)) {
		stop("Parameters 'temp' and 'effect' must be of the same length!")
	}
	nls <- nlinBoltzmannFit(temp,effect,def)
	if (class(nls)=="try-error" || nls$convergence>1) {
		class(nls) <- "try-error"
		return(nls)
	}
	pred_effect <- evalBoltzmann(temp,nls$par)
	bfit <- list(temp=temp,effect=effect,fitted.values=pred_effect,residuals=effect-pred_effect,
					coefficients=nls$par,convergence=nls$convergence,message=nls$message,call=match.call())
	names(bfit$coefficients) <- c("E0","Ef","T50","s")
	class(bfit) <- "boltzfit"
	return(bfit)
}

evalBoltzmann <- function(temp,par,calcderivs=FALSE) {
	den <- 1+exp(-(temp-par[3])/par[4])
	eff <- par[1]+(par[2]-par[1])/den
	eff[is.infinite(den)] <- par[1]
	if (!calcderivs) { return(eff) }
	ederivs <- array(0,dim=c(length(temp),5))
	ederivs[,1] <- eff
	ederivs[,3] <- 1/den
	ederivs[,2] <- 1-ederivs[,3]
	ederivs[,4] <- -((den-1)/par[4])*(eff-par[1])/den
	ederivs[,5] <- ((temp-par[3])/par[4])*ederivs[,4]
	ederivs[is.infinite(den),2:5] <- 0
	ederivs[is.nan(ederivs)] <- 0
	return(ederivs)
}

nlinBoltzmannFit <- function(temp,effect,def=NULL) {
	parfunc <- function(par) {
		est <- evalBoltzmann(temp,c(par[1:3],exp(par[4])))
		err <- est-effect
		return(sum(err^2))
	}
	derivfunc <- function(par) {
		ederivs <- evalBoltzmann(temp,c(par[1:3],exp(par[4])),calcderivs=TRUE)
		err <- ederivs[,1]-effect
		pderivs <- apply(ederivs[,2:5],2,function(v) 2*sum(v*err))
		pderivs[4] <- exp(par[4])*pderivs[4]
		return(pderivs)
	}
	ord <- order(temp)
	if (is.null(def)) {
		spar <- c(min(effect),max(effect),mean(temp),log(2))
		llims <- c(rep(min(effect),2),min(temp),-log(10))
		ulims <- c(rep(max(effect),2),max(temp),log(10))
	} else {
		spar <- c(def[1:3],log(def[4]))
		if (def[1]<def[2]) {
			llims <- c(min(effect),(def[1]+def[2])/2,min(temp),-log(10))
			ulims <- c((def[1]+def[2])/2,max(effect),max(temp),log(10))
		} else {
			llims <- c((def[1]+def[2])/2,min(effect),min(temp),-log(10))
			ulims <- c(max(effect),(def[1]+def[2])/2,max(temp),log(10))
		}
	}
	nls <- try(optim(spar,parfunc,derivfunc,method="L-BFGS-B",lower=llims,upper=ulims,hessian=FALSE),silent=TRUE)
	if (class(nls)!="try-error") { nls$par[4] <- exp(nls$par[4]) }
	return(nls)
}

welldata$MinTemp <- welldata$StartTemp
welldata$MaxTemp <- welldata$StopTemp
pvec <- c("E0","Ef","T50","s","R2","refit")
welldata[,pvec] <- 0
welldata$flags <- ""
for (wi in 1:nrow(welldata)) {
	wrel <- rawdata$Well==welldata$Well[wi]&rawdata$Temperature>=welldata$StartTemp[wi]&
			rawdata$Temperature<=welldata$StopTemp[wi]
	reldata <- rawdata[wrel,]
	nrel <- nrow(reldata)
	reldata <- reldata[order(reldata$Temperature),]
	dtemp <- reldata$Value[2:nrel]-reldata$Value[1:(nrel-1)]
	maxd <- which.max(dtemp)
	maxdT <- mean(reldata$Temperature[maxd:(maxd+1)])
	subtemp <- reldata[reldata$Temperature<maxdT,]
	welldata$MinTemp[wi] <- subtemp$Temperature[which.min(subtemp$Value)]
	suptemp <- reldata[reldata$Temperature>maxdT,]
	welldata$MaxTemp[wi] <- suptemp$Temperature[which.max(suptemp$Value)]
	
	
	frel <- reldata$Temperature<=welldata$MaxTemp[wi]
	if (length(which(frel))<=2) {
		nofit <- TRUE
		welldata[wi,pvec] <- c(min(reldata$Value[frel]),max(reldata$Value[frel]),mean(reldata$Temperature[frel]),2,0,0)
		welldata$flags[wi] <- "Insufficient data"
	} else {
		nofit <- FALSE
		flags <- c()
		bfit <- fitBoltzmann(reldata$Temperature[frel],reldata$Value[frel])
		if (class(bfit)=="try-error" || var(fitted(bfit))==0) { bcor <- 0 }
		else { bcor <- cor(bfit$effect,fitted(bfit))^2 }
		refit <- 1
		if (any(reldata$Temperature<welldata$MinTemp[wi]) || bcor==0 || coef(bfit)[1]>coef(bfit)[2]) {
			cthresh <- 0.8
			if (bcor < cthresh || coef(bfit)[1]>coef(bfit)[2]) {
				frel <- reldata$Temperature>=welldata$MinTemp[wi] & reldata$Temperature<=welldata$MaxTemp[wi]
				if (length(which(frel))>2) {
					wdef <- c(min(reldata$Value[frel]),max(reldata$Value[frel]),maxdT,2)
					bfit2 <- fitBoltzmann(reldata$Temperature[frel],reldata$Value[frel],def=wdef)
					if (class(bfit2)=="try-error" || var(fitted(bfit2))==0) { bcor2 <- 0 }
					else { bcor2 <- cor(bfit2$effect,fitted(bfit2))^2 }
					if (bcor2>=0.9 || (bcor2>0 && (bcor==0 ||coef(bfit)[1]>coef(bfit)[2]))) {
						cthresh <- 0.9
						bfit <- bfit2
						bcor <- bcor2
						refit <- 2
					}
				}
			}
		} else { cthresh <- 0.9 }
		if (bcor<cthresh) { flags <- c(flags,"Poor fit") }
		if (bcor>0) { welldata[wi,pvec] <- c(coef(bfit),bcor,refit) }
		else {
			nofit <- TRUE
			wrel <- tab$WELL==w & tab$HighCut==0 & tab$LowCut==0
			welldata[wi,pvec] <- c(min(reldata$Value[frel]),max(reldata$Value[frel]),mean(reldata$Temperature[frel]),2,0,0)
		}
		welldata$flags[wi] <- paste(flags,collapse="; ")
	}
}

# Node 9: R to Table
# Save 'welldata' data frame to outgoing Knime data stream
knime.out <- welldata
