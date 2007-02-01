degreeNorm <- function(unrecipIn, unrecipOut){

  sumUnrecip <- unrecipIn+unrecipOut
  diffUnrecip <- unrecipIn-unrecipOut

  normDeg <- (diffUnrecip)/(sqrt(sumUnrecip))
  return(normDeg)

}

  
