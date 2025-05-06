#-----------Functions for simulating Ornstein-Uhlenbeck process-------------------

#Random effects for Ornstein-Uhlenbeck process
random_effects_OU <- function(n, par) {
  rnorm(n, mean = par[1], sd = par[2])
}


#Here we define the functions $F$ ang $G$ for the OU-process. These are the drift and diffusion functions.
#Drift-function
F_OU <- function(x, U){
  return(U*x)
}

#Diffusion-function, here it is simply constant 1:
G_OU <- function(x, U){
  return(1)
}


#-------------Functions for estimating the parameters of the OU-process-------------------

#For the loglikelihood approach from Delattre et al. (2013)
#U-function
U_in <- function(path){
  sum(path[-length(path)] * (path[-1] - path[-length(path)]))
}

#V-function
V_in <- function(path, sampleinterval){
  sampleinterval * sum(path[-length(path)]^2)
}



#The actual estimatorfunction
OU_estimator <- function(data){
  #applying the U and V functions to the data
  U <- apply(data, 2, U_in)
  V <- apply(data, 2, V_in, sampleinterval = 0.01)
  
  #negative contrast function
  l_Nn <- function(par){
    1/2 * sum(log(1 + V * par[2])) + 1/2 * sum((par[1] - 1/V*U)^2 * (1 + V * par[2])^(-1))
  }
  
  #Optimizing the contrast function using optim
  optim_result <- optim(par = c(0, 1), 
                        fn = l_Nn
  )
  
  return(optim_result$par)
}