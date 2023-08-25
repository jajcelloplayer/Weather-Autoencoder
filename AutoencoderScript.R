################################################################################
################################################################################
####################### General Melding Autoencoder ############################
################################################################################
################################################################################

# Libraries ---------------------------------------------------------------
library(tidyverse)
library(truncnorm)
library(mvtnorm)
library(viridis)
library(colorspace)

# Configuration Options ---------------------------------------------------
#### Neural Network Shape
n.neurons <- c(4,3,2,1,2,3,4)

#### Network Activation Functions
act.fns <- c(NA, rep("relu", length(n.neurons)-2),"ident") 

#### Data Options
scaledInput <- TRUE
tied_down <- 'TRMM'
month <- '201507'

#### Basis Function Setup
basis_rows <- 6 
basis_cols <- 6

#### MCMC Options
n.draws <- 50000
n.burn <-  100000
n.thin <- 50

#### Output Filename
outfilename <- 'outfile.Rdata'


# Data Preparation --------------------------------------------------------
#### Load in the Key for which Data is available which month
load("month_key.Rdata")
#### Load in the data files
load("aligned_era5.Rdata")
load("aligned_merra2.Rdata")
load("trmm.Rdata")
load('aphro.Rdata')

#### filter the data for the specific month
prods <- c('aphro', 'trmm', 'era5', 'merra2')
to_use <- prods[as.logical(month_key[which(month_key$month == month),2:5])]
dat <- matrix(NA, nrow=length(c(get(to_use[1])[,,month])), ncol=length(to_use))
for(i in 1:length(to_use)){
  dat[,i] <- c(get(to_use[i])[,,month])
}
dat <- dat %>%
  as_tibble(.name_repair="unique")
names(dat) <- colnames(month_key[,2:5])[as.logical(month_key[which(month_key$month == month),2:5])]
models <- names(dat)[-which(names(dat)==tied_down)] ## Tie down a model to be constrained
dat <- select(dat, all_of(tied_down), all_of(models))
N <- prod(dim(dat))

#### Change the NN structure based on available data products
l <- length(to_use)
for(i in 1:length(n.neurons)){
  if(n.neurons[i] > l){
    n.neurons[i] <- l
  }
}

#### Latitude and Longitude values
dat$lon <- rep(as.numeric(rownames(aphro)), ncol(aphro))
dat$lat <- rep(as.numeric(colnames(aphro)), each=nrow(aphro))


# Defining Functions ------------------------------------------------------
#### Activation Functions
act <- function(input, fx="ident"){
  if(fx%in%c("relu", "tanh", "sigmoid", "ident")){
    output <- switch(fx,
                     "relu" = pmax(matrix(0, nrow=nrow(input), ncol=ncol(input)), input),
                     "tanh" = tanh(input),
                     "sigmoid" = plogis(input),
                     "ident" = input)
    output <- 
      return(output)
  } else {
    stop(paste(fx, 'not a valid activation function'))
  }
}
  
#### Negation Function
`%!in%` = Negate(`%in%`)

#### Matrix Multiplication Functions
WMultiply1 <- function(input, weight, which="in"){
  if(is.matrix(input)){
    if(is.matrix(weight)){
      
      w <- ncol(input) # How many nodes do we have on the previous layer
      j <- ncol(weight)/w # How many nodes we will have on the next layer
      
      in_function <- function(input=input, weight=weight){
        R <- matrix(NA, nrow=nrow(input), ncol=j)
        for(u in 1:j){
          A <- as.matrix(rowSums(input * weight[,((u-1)*w+1):(u*w)]))
          R <- cbind(R, A)
          R <- as.matrix(R[,-1])
        }
        return(R)
      }
      
      out_function <- function(input=input, weight=weight){
        R <- matrix(0, nrow=nrow(input), ncol=j)
        for(o in 1:ncol(input)){
          R <- R + (as.numeric(input[,o]) * weight[,((o-1)*j+1):(o*j)])#/ncol(input)  # Do we need the divide?
        }
        R
      }
      
      if(j == round(j)){   ## Or we can use  j %% 1 == 0
        
        output <- switch(which,
                         "in" = in_function(input, weight),
                         "out" = out_function(input, weight))
        output <- 
          return(output)
        
      } else {
        stop("The number of columnns in the weight matrix is not a multiple of the number of columns in the input matrix")
      }
    } else {
      stop("invalid weight matrix")
    }
  } else {
    stop("invalid input matrix")
  }
}

#### AMCMC Function
source("AMCMCUpdate.R")
#### AMCMC Options
amcmc.draws <- 250
lr <- 0.1^2

# Basis Function Setup ----------------------------------------------------
x <- basis_cols+1 ; y <- basis_rows+1 
m <- x * y  
n <- length(dat$lon) 
x_support <- seq(min(dat$lon)-0.25, max(dat$lon)+0.25, length=x)
y_support <- seq(min(dat$lat)-0.25, max(dat$lat)+0.25, length=y)
#### Create a matrix for all the support points
M <- expand.grid(x_support, y_support)
#### Create a Data Location Matrix
S <- matrix(c(dat$lon,dat$lat), ncol=2, nrow=n)

#### Create the Sigma Matrix & Initialize values for x_sd, y_sd, r
x_sd <- (max(M[,1])-min(M[,1]))/x/((max(M[,1])-min(M[,1]))/(max(M[,2])-min(M[,2])))
y_sd <- (max(M[,2])-min(M[,2]))/y/((max(M[,1])-min(M[,1]))/(max(M[,2])-min(M[,2])))
r <-  0
Sigma <- matrix(c(x_sd^2,r*x_sd*y_sd,
                  r*x_sd*y_sd,y_sd^2),
                  nrow=2, ncol=2)

#### Create our Basis (B(s)) Matrix
x_support2 <- numeric(x-1)
for(i in 1:(x-1))  x_support2[i] <- x_support[i] + (x_support[2]-x_support[1])/2 
y_support2 <- numeric(y-1)
for(i in 1:(y-1))  y_support2[i] <- y_support[i] + (y_support[2]-y_support[1])/2 
M2 <- expand.grid(x_support2, y_support2)
B <- matrix(1,nrow=n,ncol=(x-1)*(y-1)+1)
if(m > 0){
  for(i in 1:((x-1)*(y-1))){
    B[,i+1]=dmvnorm(x = S,
                    mean = as.numeric(M2[i,]),
                    sigma = Sigma)
  }
}
k <- ncol(B) 


# NN Initialization -------------------------------------------------------
#### Neural Network Setup
mid.node <- which(n.neurons==1)
L <- length(n.neurons)

#### Initialize the NN weights
dat <- dat %>% select(-lat, -lon) # get rid of the lat and lon columns
nn <- vector("list", length=L-1)
for(l in 1:length(nn)){
  if(l==1){ 
    ## Initialize input variables
    if(scaledInput){
      nn[[l]]$I <- apply(as.matrix(dat), 2, scales::rescale)
    } else {
      nn[[l]]$I <- as.matrix(dat)
    }
    
    ## Initialize the coefficients for our Theta Matrices 
    nn[[l]]$Theta <- matrix(1/n.neurons[l] + rnorm(n.neurons[l]*n.neurons[l+1], 0, .01),
                            nrow=k, ncol=n.neurons[l]*n.neurons[l+1]) 
    
    ## Create the Weight matrix from B %*% Theta 
    nn[[l]]$W <- B %*% nn[[l]]$Theta 
    
    ## Look at which parameters we are estimating
    nn[[l]]$ep <- which(row(nn[[l]]$Theta)>0)   # This is the index of Theta parameters that we are estimating on this layer
    nn[[l]]$n.ep <- length(nn[[l]]$ep)          # This is the number of Theta parameters that we are estimating on this layer
    nn[[l]]$ew <- which(nn[[l]]$W %!in% c(0,1)) # This is the index of which weights are not fixed on this layer
    nn[[l]]$n.ew <- length(nn[[l]]$ew)          # This is the number of weights that are not fixed on this layer
    
  } else {
    if(l < mid.node){
      ## Initialize input variables from the previous layer
      nn[[l]]$I <- act(WMultiply1(nn[[l-1]]$I, nn[[l-1]]$W, "in"), fx=act.fns[l]) 
      
      ## Initialize the coefficients for our Theta Matrices 
      nn[[l]]$Theta <- matrix(1/n.neurons[l] + rnorm(n.neurons[l]*n.neurons[l+1], 0, .01),
                              nrow=k, ncol=n.neurons[l]*n.neurons[l+1]) 
      
      ## Create the Weight matrix from B %*% Theta 
      nn[[l]]$W <- B %*% nn[[l]]$Theta 
      
      ## Look at which parameters we are estimating
      nn[[l]]$ep <- which(row(nn[[l]]$Theta)>0)   # This is the index of Theta parameters that we are estimating on this layer
      nn[[l]]$n.ep <- length(nn[[l]]$ep)          # This is the number of Theta parameters that we are estimating on this layer
      nn[[l]]$ew <- which(nn[[l]]$W %!in% c(0,1)) # This is the index of which weights are not fixed on this layer
      nn[[l]]$n.ew <- length(nn[[l]]$ew)          # This is the number of weights that are not fixed on this layer
      
    } else {
      ## Initialize input variables from the previous layer
      nn[[l]]$I <- act(WMultiply1(nn[[l-1]]$I, nn[[l-1]]$W, "out"), fx=act.fns[l]) 
      
      ## Initialize the coefficients for our Theta Matrices 
      nn[[l]]$Theta <- matrix(1/n.neurons[l] + rnorm(((n.neurons[l]*n.neurons[l+1])-sum(1:n.neurons[l])), 0, .01),
                              nrow=k, ncol=(n.neurons[l]*n.neurons[l+1])-sum(1:n.neurons[l])) 
      
      ## Create the Weight matrix from B %*% Theta 
      W <- vector("list", n.neurons[l+1]) 
      W[[n.neurons[l]+1]] <- matrix(NA, nrow=n)  
      names(W) <- c(rep(NA,n.neurons[l]), "Complete") 
      
      ## Adding in CFA Constraints
      f <- 1
      for(v in 1:n.neurons[l]){  
        W[[v]] <- matrix(0, nrow = n, ncol = n.neurons[l+1])
        W[[v]][,v] <- 1
        z <- n.neurons[l+1]-v
        if(f <= ncol(nn[[l]]$Theta)){
          W[[v]][,((v+1):(n.neurons[l+1]))] <- B %*% nn[[l]]$Theta[,(f:(f+z-1))]
        }
        W$Complete <- cbind(W$Complete, W[[v]])
        f <- f+z
      }
      nn[[l]]$W <- as.matrix(W$Complete[,-1])  
      
      ## Look at which parameters we are estimating
      nn[[l]]$ep <- which(row(nn[[l]]$Theta)>0) # This is the index of Theta parameters that we are estimating on this layer
      nn[[l]]$n.ep <- length(nn[[l]]$ep)        # This is the number of Theta parameters that we are estimating on this layer
      nn[[l]]$ew <- which(nn[[l]]$W %!in% c(0,1)) # This is the index of which weights are not fixed on this layer
      nn[[l]]$n.ew <- length(nn[[l]]$ew)        # This is the number of weights that are not fixed on this layer
    }
  }
}
nn.out <- act(WMultiply1(nn[[L-1]]$I, nn[[L-1]]$W, "out"), fx=act.fns[L])

#### Initialize our Latent Variable (LV)
LV <- as.matrix(dat)
LV[which(dat==0)] <- rtruncnorm(n=length(which(dat==0)),a=-Inf,b=0, mean=0,sd=.01)

#### Initialize the variance parameters   
sig2 <- max(colMeans((LV-nn.out)^2))

#### Calculate the current log-likelihood 
cur.llike <- -0.5*log(2*pi*sig2)-sum((LV-nn.out)^2)/(2*sig2)                   


# Priors ------------------------------------------------------------------
#### Prior on Thetas
theta.mn <- 0
theta.var <- 5
#### Prior on sigma^2
sig2.a <- 2.1
sig2.b <- 1.1


# MCMC Setup --------------------------------------------------------------
#### Set up Adaptive-MCMC
tot.it <- n.burn+n.thin*n.draws
kp.seq <- seq(n.burn+n.thin, tot.it, by=n.thin)
kp <- 0
amcmc <- vector("list", length=length(nn))
for(l in 1:length(amcmc)){
  amcmc[[l]]$mn <- matrix(0, nrow=nn[[l]]$n.ep, ncol=1)
  amcmc[[l]]$var <- matrix(0, nrow=nn[[l]]$n.ep, ncol=nn[[l]]$n.ep)
}
eps <- 0.0001

#### Matrices to hold draws
Theta.draws <- vector("list", length=length(nn))
for(l in 1:length(nn)){
  Theta.draws[[l]] <- matrix(0, nrow=n.draws, ncol=length(nn[[l]]$Theta))
}
sig2.draws <- matrix(0, nrow=n.draws, ncol=length(sig2))
fit.vals <- matrix(0, nrow=nrow(nn.out), ncol=ncol(nn.out))
consensus <- matrix(0, nrow=nrow(nn.out), ncol=1)
nrconsensus <- matrix(0, nrow=nrow(nn.out), ncol=1)
all_mid_nodes <- matrix(NA, nrow=n.draws, ncol=nrow(dat))    


# MCMC Sampling -----------------------------------------------------------
#### Run the Sampler
pb <- txtProgressBar(min = 0, max = tot.it, style = 3)
for(it in 1:tot.it){
  for(l in 1:length(nn)){
    np <- nrow(amcmc[[l]]$mn)
    prop.nn <- nn
    if(it > amcmc.draws){
      prop.var <- lr*(amcmc[[l]]$var+eps*diag(np))
    } else {
      prop.var <- lr*eps*diag(np)
    }
    
    ## Propose new values for all of the estimated parameters
    prop.nn[[l]]$Theta[prop.nn[[l]]$ep] <- (matrix(nn[[l]]$Theta, ncol=1) + t(chol(prop.var))%*%rnorm(np))
    
    ## Create the Weight matrix from B %*% Theta (not needed to be done and saved)
    prop.nn[[l]]$W[nn[[l]]$ew] <- B %*% prop.nn[[l]]$Theta
    
    ## Do the multiplication across the different nodes
    if(l!=length(nn)){
      for(i in (l+1):length(nn)){ 
        if(i < mid.node){
          prop.nn[[i]]$I <- act(WMultiply1(prop.nn[[i-1]]$I, prop.nn[[i-1]]$W, "in"), fx=act.fns[i])
        }else{
          prop.nn[[i]]$I <- act(WMultiply1(prop.nn[[i-1]]$I, prop.nn[[i-1]]$W, "out"), fx=act.fns[i])
        }
      }
    }
    prop.nn.out <- act(WMultiply1(prop.nn[[L-1]]$I, prop.nn[[L-1]]$W, "out"), fx=act.fns[L])
    
    ## Calculate our Metropolis Hastings Value
    prop.llike <- -0.5*log(2*pi*sig2)-sum((LV-prop.nn.out)^2)/(2*sig2)
    MH <- prop.llike - cur.llike + 
      sum(dnorm(prop.nn[[l]]$Theta, theta.mn, sqrt(theta.var), log=TRUE)) - 
      sum(dnorm(nn[[l]]$Theta, theta.mn, sqrt(theta.var), log=TRUE))
    
    ## Decide to accept or reject
    if(log(runif(1)) < MH){
      nn <- prop.nn
      nn.out <- prop.nn.out
      cur.llike <- prop.llike
    }
    
    amcmc[[l]] <- AMCMC.update(matrix(nn[[l]]$Theta, ncol=1), 
                               amcmc[[l]]$mn, amcmc[[l]]$var,
                               it)
  }
  
  ## Sample the variances
  a.star <- sig2.a + N/2 
  b.star <- 0.5*sum((LV-nn.out)^2)+sig2.b
  sig2 <- 1/rgamma(length(sig2), shape=a.star, rate=b.star)
  
  LV[which(dat==0)] <- rtruncnorm(n=length(which(dat==0)),a=-Inf,b=0, mean=nn.out[which(dat==0)],sd=sqrt(sig2))
  cur.llike <- -0.5*log(2*pi*sig2)-sum((LV-nn.out)^2)/(2*sig2)
  
  ## Store the draws if necessary
  if(it %in% kp.seq){
    kp <- kp + 1
    for(l in 1:length(nn)){
      Theta.draws[[l]][kp,] <- c(nn[[l]]$Theta)
    }
    sig2.draws[kp,] <- sig2
    fit.vals <- fit.vals + (1/n.draws)*nn.out
    consensus <- consensus + (1/n.draws)*nn[[mid.node]]$I
    all_mid_nodes[kp,] <-  nn[[mid.node]]$I
  }
  
  ## Increment progress bar
  setTxtProgressBar(pb, it)
}
close(pb)


# Saving Sampled Values ---------------------------------------------------
#### Save the Latitude and Longitude Values
y <- dimnames(aphro)[[1]] %>% as.numeric()
x <- dimnames(aphro)[[2]] %>% as.numeric()

#### Save the quantiles for our consensus
mid_quantiles <- t(apply(all_mid_nodes, 2, quantile, probs=c(.025, 0.5, .975)))
mid_quantiles <- cbind(mid_quantiles, consensus)
colnames(mid_quantiles) <- c("Lower", "Median", "Upper", "Mean")

#### Save our selected values
save(file=outfilename, 
       list=c("Theta.draws", "sig2.draws", "dat", "fit.vals", "mid_quantiles", "B", "LV",
            "x", "y", "n.neurons", "act.fns", "act", "consensus"))


# Plot Values of Interest -------------------------------------------------
#### Plot the ratio plots for every month
Rdat <- matrix(0, 4346,length(to_use))
rfit.vals <- act(fit.vals, "relu")
for(i in 1:length(to_use)){
  Rdat[,i] <- rfit.vals[,i]/(consensus)
}
for(i in 1:4346){
  for(j in 1:length(to_use)){
    if(is.na(Rdat[i,j])) Rdat[i,j] <- 0
  }
}
R1dat <- matrix(0, 4346,length(to_use))
for(i in 1:4346){
  R1dat[i,] <- Rdat[i,]/rowSums(Rdat)[i]
}
for(i in 1:4346){
  for(j in 1:length(to_use)){
    if(is.na(R1dat[i,j])) R1dat[i,j] <- .25
  }
}
R2dat <- R1dat -.25
plot2dat <- cbind(expand.grid(y,x),R2dat)
colnames(plot2dat) <- c("x", "y", names(dat))
plot2dat <- as_tibble(plot2dat)

null.plot <- ggplot(plot2dat, aes(x,y,fill=0)) + geom_raster() + scale_fill_continuous_diverging(palette="Blue-Red 3", limits=c(-.75,.75))
for(ratio in c("TRMM", "APHRO", "MERRA", "ERA")) assign(paste0(ratio,"plot"), null.plot)

for(ratio in names(dat)){
  fit.plot <- ggplot(plot2dat,aes_string(x="x", y="y", fill=ratio)) + geom_raster() + scale_fill_continuous_diverging(palette="Blue-Red 3", limits=c(-.75,.75))
  assign(paste0(ratio,"plot"), fit.plot)
}

jpeg(file=paste0("./Plots/Ratio_Plots/Ratio", paste(month),".jpeg"), width = 1200, height = 800,)
gridExtra::grid.arrange(APHROplot, TRMMplot, ERAplot, MERRAplot)  
dev.off()


#### Fitted vs Actual plots
plot.df <- bind_cols(expand.grid(y,x), as_tibble(act(fit.vals, "relu")))
names(plot.df) <- c("x","y",names(dat))
plotdat <- bind_cols(expand.grid(y,x), dat)
names(plotdat)[1:2] <- c("x", "y")

null.plot <- ggplot(data=plot.df, aes_string(x="x", y="y", fill=0)) + 
  geom_raster() + scale_fill_viridis(option="E") 
for(ratio in c("TRMM", "APHRO", "MERRA", "ERA")){
  assign(paste0(ratio,"fitplot"), null.plot)
  assign(paste0(ratio,"obsplot"), null.plot)
}

for(product in names(dat)){
  
  fitplot <- ggplot(data=plot.df, aes_string(x="x", y="y", fill=product)) + 
    geom_raster() + scale_fill_viridis(option="B", 
                                       limits=c(0,max(c(fit.vals[,which(names(dat)==product)],dat[[product]])))) +
    labs(fill=paste("Fitted",product))
  
  obsplot <- ggplot(data=plotdat, aes_string(x="x", y="y", fill=product)) + 
    geom_raster() + scale_fill_viridis(option="B", 
                                       limits=c(0,max(c(fit.vals[,which(names(dat)==product)],dat[[product]])))) 
  
  assign(paste0(product,"fitplot"), fitplot)
  assign(paste0(product,"obsplot"), obsplot)
  
}

## Plot the middle node
mid.df <- bind_cols(expand.grid(y,x), consensus)
names(mid.df) <- c("x", "y", "Consensus")
mid.GS <- ggplot(data=mid.df, aes(x=x, y=y, fill=Consensus)) + geom_raster() +
  scale_fill_viridis(option="B")+labs(x="Longitude", y="Lattitude", fill="Precipitation")+ggtitle("Consensus") 

## Save the picture
jpeg(file=paste0("./Plots/Fitted_Plots/Fitted", paste(month),".jpeg"), width = 1200, height = 800,)
gridExtra::grid.arrange(APHROfitplot, APHROobsplot,
                        TRMMfitplot,  TRMMobsplot,
                        MERRAfitplot, MERRAobsplot,
                        ERAfitplot,   ERAobsplot, 
                        mid.GS)
dev.off()


#### Save the middle node as its own picture
jpeg(file=paste0("./Plots/Mid/Mid", paste(month),".jpeg"), width = 1200, height = 800,)
gridExtra::grid.arrange(mid.GS)
dev.off()


# Garbage Collection ------------------------------------------------------
#### Remove everything that we don't need 
rm(list=ls()[ls() %!in% c("month", "month_key", "nn.train", "%!in%")])
gc()
