# Node 12: Table to R
# Read in raw temperature-fluorescence data
rawdata <- knime.in

# Node 19: Add Table to R
# Read in fitted Boltzmann parameters for each well
welldata <- knime.in

# Node 20: R to Table
# Evaluate best-fit Boltzmann curve for each well and
# add fitted data to raw data data frame.  Send 'rawdata'
# to outgoing Knime datastream
evalBoltzmann <- function(temp,par) {
	den <- 1+exp(-(temp-par[3])/par[4])
	eff <- par[1]+(par[2]-par[1])/den
	eff[is.infinite(den)] <- par[1]
	return(eff)
}

rawdata$Fit <- 0
for (wi in 1:nrow(welldata)) {
	rel <- rawdata$Well==welldata$Well[wi]
	par <- as.numeric(welldata[wi,c("E0","Ef","T50","s")])
	rawdata$Fit[rel] <- evalBoltzmann(rawdata$Temperature[rel],par)
}

knime.out <- rawdata
