###############################################################################
## This script is for mJAM (multi-ethnic Joint Analysis of Marginal statistics)
## with SuSiE (sufficient statistic version) as the model selection method.
##
## Section 1: Functions to assist intermediate steps in the main mJAM-SuSiE function
## Section 2: The main function of mJAM-SuSiE
## Section 3: Functions to draw inferences from a fitted mJAM-SuSiE model.
##

##########################  Section 1  #################################
## Function 1.1 -- Get transformed statistics: z, or Xty
get_z <- function(maf, betas, N_outcome){

  ## follows JAM supplementary material section 1
  n0 = N_outcome*(1-maf)^2
  n1 = N_outcome*2*maf*(1-maf)
  n2 = N_outcome*maf^2

  y0 = -(n1*betas+2*n2*betas)/(n0+n1+n2)
  y1 = y0+betas
  y2 = y0+2*betas
  z = n1*y1 + 2*n2*y2

  return(z)
}


get_Xty <- function(maf, betas, N_outcome){

  Dj <- 2*maf*(1-maf)*N_outcome
  Xty <- diag(Dj) %*% betas

  return(Xty)
}

