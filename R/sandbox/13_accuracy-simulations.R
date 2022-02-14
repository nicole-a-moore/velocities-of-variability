### trying to recreate mathematica script in R



Levy = function (alpha) {
  phi = runif(1, -pi/2, pi/2)
  phi0 = -(pi/2)*((1-abs(1-alpha))/alpha)
  result = sign(alpha-1)*sin(alpha*(phi-phi0))*(cos(phi)*abs(alpha-1))^-(1/alpha)
  result = result*as.complex((cos(phi-alpha*(phi-phi0))/rexp(1,1)))^((1-alpha)/alpha)
  
  return(result)
}

as.numeric(sapply(rep(2, 10), Levy))

