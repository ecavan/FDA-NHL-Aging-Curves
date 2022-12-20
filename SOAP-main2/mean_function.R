#Packages in use
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(fda))
suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(locfit))
suppressPackageStartupMessages(library(fdapace))

#define b(t): basis function vector
#I simply copied the parameters used in the demo with the corresponding meanings
#spline_basis = fda::create.bspline.basis(rangeval = c(0,6), #total number of basis functions
#                                         nbasis = 10, #L: the length of each basis functions
#                                         norder = 4) #cubic splines: 4 as default

########### MEAN FUNCTION ###############

mean_function = function(beta1, #initial beta
                         observed, #y_ij
                         timepoints, #t_ij
                         threshold = 1e-3, #threshold for iteration stop
                         minit = 3) #minimum iteration round
{
  thresh = 1 #define the computed thresh for combination with threshold set
  it = 1 #round checking for iteration, starting from 1
  
  beta1_before = rep(0, length(beta1)) #define vector for initial beta with fixed length
  value = -1e14 #loss function default value
  
  #Fit the initial mean function
  mean_fit = fda::fd(beta1, spline_basis) #mean fitted as a linear combination of basis functions
  #mean_fit = (1 / sqrt(inprod(mean_fit, mean_fit)) * mean_fit) #normalization
  beta1 = coef(mean_fit) %>% as.numeric #coefs for beta1
  #print(beta1)
  #When: computed thresh greater than the threshold set 
  #Or, the iteration time less than the minimum round 
  # --------->> Run the iteration process as following: 
  while (thresh > threshold || it < minit){
    beta1_before = beta1
    value_before = value
    #print(beta1_before)
    #Calculate the residuals
    residual_fit = function(subj_index){
      timei = timepoints[[subj_index]]
      xmati = eval.basis(timei, spline_basis) #b(t_ij) 
      lm(observed[[subj_index]] ~ xmati + 0) %>% residuals
    }
    rfit0 = lapply(1:length(observed), residual_fit)
    rfit = lapply(rfit0, function(x)
    {x^2 %>% mean}) %>% do.call(c,.)
    
    value = sum(rfit^2) #2.1 there is no 1/n mean, replaced with sum
    
    # print(value)
    #Converge judgment
    if(abs(value_before - value/value_before) < threshold) break
    
    yem = do.call(c,observed) #convert into one vector
    
    #Design matrix in use
    xalpha = lapply(1:length(observed), function(subj_index){
      timei = timepoints[[subj_index]]
      xmati = eval.basis(timei, spline_basis) #b(t_ij)
      xmati
      }) %>% do.call(rbind,.)
    
    #print(xalpha)
    #print(yem)
    #ax = b, x = solve(a,b)
    #Solve the new beta(l+1), LSE method
    #beta1 = 
    yem[is.na(yem)] <- 0
    beta1 = solve(t(xalpha%>%as.matrix)%*%(as.matrix(xalpha)), t(xalpha%>%as.matrix)%*%as.matrix(yem))
    #print(beta1)
    #Update new mean function and beta
    mean_fit = fd(beta1, spline_basis)
    #mean_fit = (1/sqrt(inprod(mean_fit, mean_fit)) * mean_fit) #normalization
    beta1 = coef(mean_fit) %>% as.numeric
    #print(beta1)
    #Calculate the thresh for iteration judgment
    thresh =  max(abs(beta1_before - beta1))
    #print(thresh)
    it = it + 1
    #Record the iteration round every 30 times
    if(it%%30 == 0) {print(it); print(as.numeric(thresh))}
  }
  
  #Satisfied beta and corresponding residuals
  mean_fit = fd(beta1, spline_basis)
  #mean_fit = (1 / sqrt(inprod(pc_fit, pc_fit)) * pc_fit) #normalization
  beta1 = coef(mean_fit) %>% as.numeric
  rfit = lapply(1:length(observed), residual_fit) %>% do.call(c,.)
  
  print("Done!")
  cat(sprintf('The threshold is %s . \n', thresh))
  
  #Outputs
  return(list(beta = beta1,
              mean_fit = mean_fit,
              thresh = thresh,
              it = it,
              value = value,
              residual = rfit0)) 
}