## Importance
#### Code for calculating the importance of each data product to the consensus

## Working Directory
setwd("~/R/Research/TruncNorm/Movie2")

## Libraries
library(tidyverse)
library(mvtnorm)
library(viridis)
library(colorspace)

## Functions ------------------------------------------------------
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


## Data
load('Results/Movie2015074.RData')

## Average Weights
load('aphro.Rdata')
lon <- rep(as.numeric(rownames(aphro)), ncol(aphro))
lat <- rep(as.numeric(colnames(aphro)), each=nrow(aphro))

# Basis Function Setup ----------------------------------------------------
x <- 7; y <- 7 
m <- x * y  
n <- length(lon) 
x_support <- seq(min(lon)-0.25, max(lon)+0.25, length=x)
y_support <- seq(min(lat)-0.25, max(lat)+0.25, length=y)
#### Create a matrix for all the support points
M <- expand.grid(x_support, y_support)
#### Create a Data Location Matrix
S <- matrix(c(lon, lat), ncol=2, nrow=n)

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

#### Extract the average NN weights
mid.node <- which(n.neurons==1)
L <- length(n.neurons)
nn <- vector("list", length=L-1)
for(l in 1:length(nn)){
  if(l==1){ 
    nn[[l]]$I <- apply(as.matrix(dat), 2, scales::rescale)
    
    ## Initialize the coefficients for our Theta Matrices 
    nn[[l]]$Theta <- matrix(colMeans(Theta.draws[[l]]),
                            nrow=k, ncol=n.neurons[l]*n.neurons[l+1]) 
    
    ## Create the Weight matrix from B %*% Theta 
    nn[[l]]$W <- B %*% nn[[l]]$Theta 
    
  } else {
    if(l < mid.node){
      
      ## Initialize the coefficients for our Theta Matrices 
      nn[[l]]$Theta <- matrix(colMeans(Theta.draws[[l]]), 
                              nrow=k, ncol=n.neurons[l]*n.neurons[l+1])
      
      ## Create the Weight matrix from B %*% Theta 
      nn[[l]]$W <- B %*% nn[[l]]$Theta 
      
    } else {
      
      ## Initialize the coefficients for our Theta Matrices 
      nn[[l]]$Theta <- matrix(colMeans(Theta.draws[[l]]), 
                              nrow=k, ncol=n.neurons[l]*n.neurons[l+1])
      
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
    }
  }
}

# Jitter Each Product -----------------------------------------------------
## Jitter the data and see how it changes the consensus

cj <- function(reps){
  new_consensus <- array(NA, c(reps, length(consensus), 4))
  for(m in 1:reps){
    scale <- runif(nrow(dat), 0, 2)
    for(p in 1:4){
      new_dat <- as.matrix(dat)
      new_dat[,p] <- new_dat[sample(1:nrow(dat)),p]
      
      new_scaled <- apply(as.matrix(new_dat), 2, scales::rescale)
      
      nn[[1]]$I <- new_scaled
      
      ## Multiply to get the consensus
      for(l in 1:length(nn)){
        ## Do the multiplication across the different nodes
        if(l!=length(nn)){
          for(i in (l+1):length(nn)){ 
            if(i < mid.node){
              nn[[i]]$I <- act(WMultiply1(nn[[i-1]]$I, nn[[i-1]]$W, "in"), fx=act.fns[i])
            }else{
              nn[[i]]$I <- act(WMultiply1(nn[[i-1]]$I, nn[[i-1]]$W, "out"), fx=act.fns[i])
            }
          }
        }
      }
      new_consensus[m,,p] <- nn[[mid.node]]$I
    }
  }
  return(new_consensus)
}

set.seed(27)
results <- apply(cj(5), c(2,3),mean)
summary(results)


scaled_res <- results
for(i in 1:4){
  scaled_res[,i] <- results[,i] / consensus
}

cms <- colMeans(scaled_res[-which(consensus == 0),] ) 
cms <- cms/sum(cms)
cms

library(xtable)
xtable(t(cms), )

# Maps --------------------------------------------------------------------
dim(results)
scaled_res1 <- scaled_res - 1
ggplot() +
  geom_raster(aes(lon, lat, fill=scaled_res1[,1])) +
  scale_fill_continuous_diverging(palette="Blue-Red", limits=c(-.35,.35))

ggplot() +
  geom_raster(aes(lon, lat, fill=scaled_res1[,2])) +
  scale_fill_continuous_diverging(palette="Blue-Red", limits=c(-.35,.35))
ggplot() +
  geom_raster(aes(lon, lat, fill=scaled_res1[,3])) +
  scale_fill_continuous_diverging(palette="Blue-Red", limits=c(-.35,.35))
ggplot() +
  geom_raster(aes(lon, lat, fill=scaled_res1[,4])) +
  scale_fill_continuous_diverging(palette="Blue-Red", limits=c(-.35,.35))


part_vals <- scaled_res
for(i in 1:4){
  part_vals[,i] <- scaled_res[,i] / rowSums(scaled_res)
}
part_vals <- part_vals - 0.25

summary(part_vals)
colnames(part_vals) <- colnames(dat)

r1 <- ggplot() +
  geom_raster(aes(lon, lat, fill=part_vals[,1]), show.legend=F) +
  scale_fill_continuous_diverging(palette="Blue-Red", limits=c(-.28,.28)) +
  labs(x='Longitude', y='Lattitude') + ggtitle('d) TRMM') +
  theme(text=element_text(size=16),
        legend.title=element_text(size=14),
        legend.text=element_text(size=12))
r2 <- ggplot() +
  geom_raster(aes(lon, lat, fill=part_vals[,2]), show.legend=F)  +
  scale_fill_continuous_diverging(palette="Blue-Red", limits=c(-.28,.28)) +
  labs(x='Longitude', y='Lattitude') + ggtitle('a) APHRODITE') +
  theme(text=element_text(size=16),
        legend.title=element_text(size=14),
        legend.text=element_text(size=12))
r3 <- ggplot() +
  geom_raster(aes(lon, lat, fill=part_vals[,3]), show.legend=F)  +
  scale_fill_continuous_diverging(palette="Blue-Red", limits=c(-.28,.28)) +
  labs(x='Longitude', y='Lattitude') + ggtitle('c) ERA5') +
  theme(text=element_text(size=16),
        legend.title=element_text(size=14),
        legend.text=element_text(size=12))
r4 <- ggplot() +
  geom_raster(aes(lon, lat, fill=part_vals[,4]), show.legend=F)  +
  scale_fill_continuous_diverging(palette="Blue-Red", limits=c(-.28,.28)) +
  labs(x='Longitude', y='Lattitude') + ggtitle('b) MERRA2') +
  theme(text=element_text(size=16),
        legend.title=element_text(size=14),
        legend.text=element_text(size=12))

Contribution <- 0
rlegend <- ggplot() +
  geom_raster(aes(lon, lat, fill=Contribution))  +
  scale_fill_continuous_diverging(palette="Blue-Red", limits=c(-.28,.28)) +
  theme(axis.ticks = element_blank(), 
        axis.text = element_blank()) + 
  xlim(c(0,0)) +
  labs(x=element_blank(),y=element_blank(), fill='Contribution')

library(gridExtra)
library(ggplotify)

pdf(file="./PaperPlots/Contributions_Shuffle.pdf", width = 12, height = 8)
grid.arrange(r2, r4, r3, r1, right=as.grob(rlegend))

dev.off()
