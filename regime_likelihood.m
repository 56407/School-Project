function LL= likelihood(p,q,alpha0,beta0,alpha1,beta1,sigma0,sigma1,gamma,R,deltareturn)
P1 = (1-p)/(2-p-q);
P2 = (1-q)/(2-p-q);
L = 0;
for i = (1: num+1)
    PP1 = q*P1 + (1-p)*P2;
    PP2 = (1-q)*P1 + p*P2;
    PP1=(normpdf(deltareturn(t),alpha0+beta0*R(t-1),sigma0*R(t-1)^(gamma))*pp1)/(normpdf(deltareturn(t),alpha0+beta0*R(t-1),sigma0*R(t-1)^(gamma))*pp1+normpdf(deltareturn(t),alpha1+beta1*R(t-1),sigma1*R(t-1)^(gamma))*pp2);
    PP2=(normpdf(deltareturn(t),alpha1+beta1*R(t-1),sigma1*R(t-1)^(gamma))*pp2)/(normpdf(deltareturn(t),alpha0+beta0*R(t-1),sigma0*R(t-1)^(gamma))*pp1+normpdf(deltareturn(t),alpha1+beta1*R(t-1),sigma1*R(t-1)^(gamma))*pp2);
    L=L+log(PP1*normpdf(deltareturn(t),alpha0+beta0*R(t-1),sigma0*R(t-1)^(gamma))+PP2*normpdf(deltareturn(t),alpha1+beta1*R(t-1),sigma1*R(t-1)^(gamma)));
    P1=PP1;
    P2=PP2;
end
LL = -L;
end
