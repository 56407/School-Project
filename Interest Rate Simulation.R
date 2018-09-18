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
  #v
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
# Input is vector of uniform random numbers
#
MYNORM <- function(rand_c) {
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
  # Loop over set of uniform random numbers and transform
  #
  jcount<-1
  numran<-length(rand_c)
  rand_norm_c<-rep(0,numran)
  while(jcount <= numran) {
    u<-rand_c[jcount]
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
# Call Uniform Random Number Generator (INPUT: SEED and NUMBER)
#
seed<-0
numran<-1000000
inputvar<-c(seed,numran)
rand_uniform_c<-MYUNIFORM(inputvar)
rand_uniform_c[1:10]
#hist(rand_uniform_c)
#
# Call Normal Random Number Generator (INPUT: UNIFORM RANDOM NUMBER VECTOR)
#
rand_norm_c<-MYNORM(rand_uniform_c)

#generate stock price path using Ornstein-Uhlenbeck process
OUstockprice<-function(S0,mu,alpha,Volatility,T,numpath){
  numstep<-T*12
  dt<-1
  ST_Vals<-matrix(rep(0),nrow=numstep, ncol=numpath)
  jcount<-1
  lcount<-1
  while(jcount <= numpath) {
    St<-S0
    rt<-mu
    for (kcount in 1:numstep) {
      rt<-rt*exp(-alpha*dt)+mu*(1-exp(-alpha*dt))+Volatility*sqrt((1-exp(-2*alpha*dt))/(2*alpha))*rand_norm_c[lcount]
      ST_Vals[kcount,jcount]<-rt
      lcount<-lcount+1
    }
    jcount=jcount+1
  }
  return (ST_Vals)
}



#Moving average strategy simulation
S0<-0.01003
rf<-0.05
Volatility<-0.0025
numpath<-10
T<-122/12
alpha<-1.0082
mu<-0.01502747
ST_Vals = OUstockprice(S0,mu,alpha,Volatility,T,numpath)

write.csv(ST_Vals,'Interest_Rate_Simulation.csv')