#
# Library of MYfunctions (Chidambaran)
#
# Unifrom Random Number Generator
# Normal Random Number Generator
# Poisson Random Number Generator
# Black-Scholes Call Value
# Black-Scholes Put Value
#
MYUNIFORM <- function(inputvar) {
  #
  # Initialize Constants
  #
  IM1<-2147483563
  IM2<-2147483399
  IA1<-40014
  IA2<-40692
  IQ1<-53668
  IQ2<-52774
  IR1<-12211
  IR2<-3791
  NTAB<-32
  EPS<-1.2e-7
  RNMX<-1.-EPS
  #
  # Transform Variables
  #
  IMM1<-IM1-1
  NDIV<-as.integer(1+IMM1/NTAB)
  AM<-1.0/IM1
  #
  # Initialize variables and arrays
  #
  idum<-inputvar[1]
  idum2<-123456789
  numran<-inputvar[2]
  ran2<-0
  iy<-0
  iv<-rep(0,NTAB)
  rand_uniform_c<-rep(0,numran)
  #
  # Run the random number loop
  #  
  icount<-1
  for (icount in 1:numran) {
    if (idum <= 0) {
      idum<-max(-idum,1)
      idum2<-idum
      j<-NTAB+8
      while (j > 0) {
        k=as.integer(idum/IQ1)
        idum<-IA1*(idum-k*IQ1)-k*IR1
        if (idum < 0) {idum=idum+IM1}
        if (j <= NTAB) {iv[j]<-idum}
        j<-j-1
      }
      iy<-iv[1]
    }
    k<-as.integer(idum/IQ1)
    idum<-IA1*(idum-k*IQ1)-k*IR1
    if(idum < 0) {idum=idum+IM1}
    k=as.integer(idum2/IQ2)
    idum2<-IA2*(idum2-k*IQ2)-k*IR2 
    if (idum2 < 0) {idum2<-idum2+IM2}
    j<-as.integer(iy/NDIV)+1
    iy<-iv[j]-idum2
    iv[j]<-idum
    if(iy < 1) {iy<-iy+IMM1}
    ran2<-min(AM*iy,RNMX)
    rand_uniform_c[icount]<-ran2
  }
  return(rand_uniform_c)
}
#
# Inverse Normal Generator
# Calls Uniform Random Number Generator
#
# Input is SEED and Number of Random Numbers
#
MYNORM <- function(seed,numran) {
  inputvar<-rep(0,2)
  inputvar[1]<-seed
  inputvar[2]<-numran
  #
  # Call Uniform Random Number Generator
  #
  rand_uniform_c<-MYUNIFORM(inputvar)
  # Initialize Constants
  a0<-2.50662823884
  a1<--18.61500062529
  a2<-41.39119773534
  a3<--25.44106049637
  b0<--8.47351093090
  b1<-23.08336743743
  b2<--21.06224101826
  b3<-3.13082909833
  c0<-0.3374754822726147
  c1<-0.9761690190917186
  c2<-0.1607979714918209
  c3<-0.0276438810333863
  c4<-0.0038405729373609
  c5<-0.0003951896511919
  c6<-0.0000321767881768
  c7<-0.0000002888167364
  c8<-0.0000003960315187
  #
  # Call Uniform Random Number Generator
  #
  inputvar<-c(seed,numran)
  rand_uniform_c<-MYUNIFORM(inputvar)
  #
  # Loop over set of uniform random numbers and transform
  #
  jcount<-1
  rand_norm_c<-rep(0,numran)
  while(jcount <= numran) {
    u<-rand_uniform_c[jcount]
    y<-u-0.5
    if(abs(y) < 0.42) {
      r<-y*y
      x<-y*(((a3*r+a2)*r+a1)*r+a0)/((((b3*r+b2)*r+b1)*r+b0)*r+1)
    } else {
      r<-u
      if(y>0){r<-1-u}
      r<-log(-log(r))
      x<-c0+r*(c1+r*(c2+r*(c3+r*(c4+r*(c5+r*(c6+r*(c7+r*c8)))))))
      if(y<0){x<--x}
    }
    #  cat("JCOUNT",jcount,"",u,"",x,"\n")
    rand_norm_c[jcount]<-x
    jcount=jcount+1
  }
  return(rand_norm_c)
}
#
# Inverse Poisson Generator
# Input is theta, seed and numran
#
# Calls uniform random number generator
#
MYPOISSON <- function(theta,seed,numran) {
  Njumps<-rep(0,numran)
  inputvar<-rep(0,2)
  inputvar[1]<-seed
  inputvar[2]<-numran
  #
  # Call Uniform Random Number Generator
  #
  rand_uniform_c<-MYUNIFORM(inputvar)
  #
  # Loop over set of uniform random numbers and transform
  #
  for (jcount in 1:numran) {
    p_val<-exp(-theta)
    Njumps[jcount]<-0
    Fp<-p_val
    random_u=rand_uniform_c[jcount]
    while(random_u > Fp) {
      Njumps[jcount]<-Njumps[jcount]+1
      p_val<-p_val*theta/Njumps[jcount]
      Fp<-Fp+p_val
    }
  }
  return(Njumps)
}
#
# Black-Scholes Call Value
#
MYBSCall<- function(S0,Strike,rf,T,Volatility) {
  d1<-(log(S0/Strike)+(rf+0.5*Volatility^2)*T)/(Volatility*sqrt(T))
  d2<-d1-Volatility*sqrt(T)
  Nd1<-pnorm(d1, mean=0, sd=1)
  Nd2<-pnorm(d2, mean=0, sd=1)
  C<-S0*Nd1-Strike*exp(-rf*T)*Nd2
  return (C)
}
#
# Black-Scholes Put Value
#
MYBSPut<- function(S0,Strike,rf,T,Volatility) {
  d1<-(log(S0/Strike)+(rf+0.5*Volatility^2)*T)/(Volatility*sqrt(T))
  d2<-d1-Volatility*sqrt(T)
  Nmd1<-pnorm(-d1, mean=0, sd=1)
  Nmd2<-pnorm(-d2, mean=0, sd=1)
  P<-Strike*exp(-rf*T)*Nmd2-S0*Nmd1
  return (P)
}
#
# Black-Scholes Call Delta
#
MYBSCall_Delta<- function(S0,Strike,rf,T,Volatility) {
  d1<-(log(S0/Strike)+(rf+0.5*Volatility^2)*T)/(Volatility*sqrt(T))
  Nd1<-pnorm(d1, mean=0, sd=1)
  return (Nd1)
}
#
# Black-Scholes Put Delta
#
MYBSPut_Delta<- function(S0,Strike,rf,T,Volatility) {
  d1<-(log(S0/Strike)+(rf+0.5*Volatility^2)*T)/(Volatility*sqrt(T))
  Nmd1<-pnorm(-d1, mean=0, sd=1)
  return (Nmd1)
}
