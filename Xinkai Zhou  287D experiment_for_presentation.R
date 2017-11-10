# The experiment here are based on the functions below.
# No need to run the functions in "Xinkai Zhou  287D final_project.R"

set.seed(24500418)

delta_alpha = function(x_i, y_i, w_t_1, alpha_i_t_1, lambda, row_number,loss,gamma)
{
  if (loss==1)
  { # Smoothed Hinge Loss Ref. p578 top
    numerator=1-(t(x_i)%*%w_t_1)*y_i - gamma*alpha_i_t_1*y_i
    denominator=(t(x_i)%*%x_i)/(lambda*row_number)+gamma
    inside=y_i * alpha_i_t_1 + numerator/denominator
    delta_alpha_value=y_i * max(0,min(1,inside)) - alpha_i_t_1
  }
  
  if (loss==2)
  { # Non-Smoothed Hinge Loss Ref. p577 top equ
    inside=y_i*alpha_i_t_1+(lambda*row_number)*(1-(t(x_i)%*%w_t_1)*y_i)/(t(x_i)%*%x_i)
    delta_alpha_value=y_i*max(0,min(1,inside)) - alpha_i_t_1
  }
  
  if (loss==3)
  { # Squared Loss Ref. p577
    numerator=y_i-t(x_i)%*%w_t_1-0.5*alpha_i_t_1
    denominator=0.5+(t(x_i)%*%x_i)/(lambda*row_number)
    delta_alpha_value=numerator/denominator
  }
  
  return(delta_alpha_value)
}

primal_solution = function(w_t,X,y,lambda,loss,gamma)
{
  summation=0
  
  if (loss==1) # Smooth Hinge
  {
    for (i in 1:nrow(X))
    {
      argument=t(w_t)%*%X[i,]
      
      if (argument>1)
      {
        summation=summation
      }
      
      if (argument<1-gamma)
      {
        summation=summation+1-argument-0.5*gamma
      }
      
      if (argument>=(1-gamma) && argument<=1)
      {
        summation=summation+ ((1-argument)^2)/(2*gamma)
      }
    }
  }
  
  if (loss==2) # Non Smooth Hinge
  {
    for (i in 1:nrow(X))
    {
      summation=summation+max(0,1-y[i]*(t(w_t)%*%X[i,]))
    }
  }
  
  if (loss==3) # Squared Loss
  {
    # Use matrix manipulation instead of for loop: summation=
    for (i in 1:nrow(X))
    {
      summation=summation+(t(w_t)%*%X[i,]-y[i])^2
    }
  }
  
  primal_solution_value=summation/nrow(X)+0.5*lambda*(t(w_t)%*%w_t)
  return(primal_solution_value[1,1])
}

dual_solution = function(alpha_t,X,y,lambda,loss,gamma)
{  
  if (loss==1) # Smooth Hinge
  {
    sum_1=t(alpha_t)%*%y - 0.5*gamma*(t(alpha_t)%*%alpha_t)
  }
  
  if (loss==2) #Non-Smooth Hinge
  {
    sum_1=t(alpha_t)%*%y
  }
  
  if (loss==3) #Squared Loss
  {
    sum_1=t(alpha_t)%*%y - 0.25*t(alpha_t)%*%alpha_t
    #for (i in 1:nrow(X))
    #{
    #  sum_1=sum_1 + alpha_t[i]*y[i] - 0.5*gamma*(alpha_t[i])^2
    #}
  }
  
  sum_2=rep(0, ncol(X))
  
  for (i in 1:nrow(X))
  {
    sum_2=sum_2 + alpha_t[i]*X[i,]
  }
  
  dual_solution_value=sum_1/nrow(X) - (t(sum_2)%*%sum_2)/(2*lambda*(nrow(X))^2)
  return(dual_solution_value[1,1])
}

sdca = function(X,y,lambda,epsilon,stochastic,loss,gamma)
{
  t=1; duality_gap=100
  e=diag(nrow(X)); record=matrix(nrow=4,ncol=1)
  w=matrix(rep(0,ncol(X)),1,ncol(X)) 
  alpha=matrix(rep(0,nrow(X)),1,nrow(X))    
  
  if (stochastic==1)
  {
    while (duality_gap>epsilon)
    {
      sequence=sample(1:nrow(X))
      
      for (i in sequence) 
        # i: Indicating which dual coordinate to optimize
      {
        t=t+1 
        delta_alpha_value=delta_alpha(X[i,], y[i], w[t-1,], alpha[t-1,i], lambda, nrow(X),loss,gamma) # Compute delta_alpha 
        alpha=rbind(alpha, alpha[t-1,] + delta_alpha_value*e[i,]) # Update alpha
        w=rbind(w, w[t-1,] + (delta_alpha_value/(lambda*nrow(X)))*X[i,]) # Update w
        cat("When epoch is", floor(t/nrow(X)), "\t","delta_alpha_value is", delta_alpha_value, "\n")
      }
      
      primal_solution_value=primal_solution(w[t,],X,y,lambda,loss,gamma)
      dual_solution_value=dual_solution(alpha[t,],X,y,lambda,loss,gamma)
      duality_gap=primal_solution_value - dual_solution_value
      
      record=cbind(record, matrix(c(floor(t/nrow(X)),primal_solution_value,
                                    dual_solution_value, duality_gap),4,1))
      #if (floor(t/nrow(X))>25) 
      #{
       # return(record[,2:ncol(record)])
        #stop("XXX")
      #}
      #cat("Print:",record[,floor(t/nrow(X))],"\n")
    }
  }
  
  if (stochastic==2)
  {
    while (duality_gap>epsilon)
    {
      sequence=sample(seq(1:nrow(X)),replace=TRUE)
      
      for (i in sequence) 
        # i: Indicating which dual coordinate to optimize
      {
        t=t+1 
        delta_alpha_value=delta_alpha(X[i,], y[i], w[t-1,], alpha[t-1,i], lambda, nrow(X),loss) # Compute delta_alpha 
        alpha=rbind(alpha, alpha[t-1,] + delta_alpha_value*e[i,]) # Update alpha
        w=rbind(w, w[t-1,] + (delta_alpha_value/(lambda*nrow(X)))*X[i,]) # Update w
      }
      
      primal_solution_value=primal_solution(w[t,],X,y,lambda,loss)
      dual_solution_value=dual_solution(alpha[t,],X,y,lambda,loss)
      duality_gap=primal_solution_value - dual_solution_value
      
      record=cbind(record, matrix(c(floor(t/nrow(X)),primal_solution_value,
                                    dual_solution_value, duality_gap),4,1))
      
      #cat("Print:",record[,floor(t/nrow(X))],"\n")
    }
  }
  
  if (stochastic==3)
  {
    sequence=sample(1:nrow(X))
    while (duality_gap>epsilon)
    { 
      for (i in sequence) 
        # i: Indicating which dual coordinate to optimize
      {
        t=t+1 
        delta_alpha_value=delta_alpha(X[i,], y[i], w[t-1,], alpha[t-1,i], lambda, nrow(X),loss) # Compute delta_alpha 
        alpha=rbind(alpha, alpha[t-1,] + delta_alpha_value*e[i,]) # Update alpha
        w=rbind(w, w[t-1,] + (delta_alpha_value/(lambda*nrow(X)))*X[i,]) # Update w
      }
      
      primal_solution_value=primal_solution(w[t,],X,y,lambda,loss)
      dual_solution_value=dual_solution(alpha[t,],X,y,lambda,loss)
      duality_gap=primal_solution_value - dual_solution_value
      
      record=cbind(record, matrix(c(floor(t/nrow(X)),primal_solution_value,
                                    dual_solution_value, duality_gap),4,1))
      
      #cat("Print:",record[,floor(t/nrow(X))],"\n")
    }
  }
  
  return(record[,2:ncol(record)])
}


# Generate Data
X <- rbind(matrix(rnorm(120), , 2), matrix(rnorm(120, mean = 3), , 2))
y <- matrix(c(rep(1, 60), rep(-1, 60)))
epsilon=10^-1
stochastic=1
loss=1
gamma=10^-1
lambda=10

# Experiment
expe_presentation=sdca(X,y,lambda,epsilon,stochastic,loss,gamma)
