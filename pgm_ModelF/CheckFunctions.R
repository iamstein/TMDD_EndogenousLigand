
read.param.file = function(filename) {
  d                      = read_excel(filename, 1)
  param.as.double        = as.numeric(d$Value)
  names(param.as.double) = d$Parameter
  param.as.double        = param.as.double[model$pin] #keep only parameters used in ODE
}


compare_eigenvalues = function(param.as.double){
  # This function compares eigenvalues of matrix A
  # and \alpha, \beta, \gamma

  pars = as.data.frame(t(param.as.double))
  
  # Matrix A is defined below
  a11 = with(pars, -(keD1 + k12D + k13D))
  a12 = with(pars, (VD2/VD1)*k21D)
  a13 = with(pars, (VD3/VD1)*k31D)
  a21 = with(pars, (VD1/VD2)*k12D)
  a22 = with(pars, -k21D)
  a23 = 0
  a31 = with(pars, (VD1/VD3)*k13D)
  a32 = 0
  a33 = with(pars, -(keD3 + k31D))
  
  M = matrix(c(a11,a12,a13,a21,a22,a23,a31,a32,a33), 3, 3, byrow=TRUE)
  
  # Eigenvalues of M is computed below
  ev = eigen(M)$values
  
  # Implement the formula 1.21 from http://www.pfim.biostat.fr/PFIM_PKPD_library.pdf
  
  #CL = with(pars, (keD1*VD1))
  #Q2 = with(pars, k12D * VD1)
  #Q3 = with(pars, k13D * VD1)
  #a0 = with(pars, (CL/VD1)*(Q2/VD2)*(Q3/VD3))
  #a1 = with(pars, (CL/VD1)*(Q3/VD3) + (Q2/VD2)*(Q3/VD3) + (Q2/VD2)*(Q3/VD1) + (CL/VD1)*(Q2/VD2) + (Q3/VD3)*(Q2/VD1))
  #a2 = with(pars, (CL/VD1)+ (Q2/VD1) + (Q3/VD1) + (Q2/VD2) + (Q3/VD3))

  a0 = with(pars, keD1*k21D*k31D)
  a1 = with(pars, keD1*k31D + k21D*k31D + k21D*k13D + keD1*k21D + k31D*k12D)
  a2 = with(pars, keD1 + k12D + k13D + k21D + k31D)
  
  p  = a1 - (a2^2)/3
  q  = 2*(a2^3)/27 - a1*a2/3 + a0 
  r1 = (-(p^3)/27)^0.5
  r2 = 2*(r1^(1/3))

  phi = acos(-q/(2*r1))/3


  alpha = -(cos(phi)         *r2 - a2/3)
  beta  = -(cos(phi + 2*pi/3)*r2 - a2/3)
  gamma = -(cos(phi + 4*pi/3)*r2 - a2/3)
  
  print("Eigenvalues of the matrix A are")
  print(ev)
  print("alpha, beta, gamma from the formula")
  print(c(alpha, beta, gamma))
}

