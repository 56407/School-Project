%% Part1
syms a  b  c
syms xs ys zs
syms x  y  z
syms h t s 
K=5;
J=K+1;
%% CIR
muX=b*(a-x);
sigmaX=c*sqrt(x);
%% Xt to Yt
x2y=int(1/sigmaX,x);
y2x=subs((finverse(x2y)),x,y);
muY=muX/sigmaX-sym('1')/sym('2')*diff(sigmaX,x,1);
muY=subs(muY,x,y2x);
muY=simplify(muY);
sigmaY=sym('1');
%% Yt to Zt
y2z=(y-ys)*h^(-1/2);
z2y=subs(finverse(y2z),y,zs);
%% Calculate Beta
syms Htemp Expe Beta
Beta=cell(K,1);
for i=1:K
    Htemp=subs(Hermite(i),z,y2z);
    Expe=Htemp;
    for j=1:J
        Htemp=muY*diff(Htemp,y,1)+1/2*diff(Htemp,y,2);
        Expe=Expe+h^j/factorial(j)*Htemp;
    end
    Beta{i}=1/factorial(i-1)*subs(Expe,y,ys);
end
%% Calculate density
pZ=0;
for m=1:K
    pZ=pZ+Beta{m}*Hermite(m);
end
pZ=exp(-z^2/2)/sqrt(2*pi)*pZ;
pY=h^(-1/2)*subs(pZ,z,y2z);
pX=1/sigmaX*subs(pY,y,x2y);
temp=subs(x2y,x,xs);
pX=subs(pX,ys,temp);
pX=simplify(pX);
%% Plot density function
cir=subs(pX,{a,b,c,h,xs},{1,1,2,1/250,1});
ezplot(cir,[-1,5]);
%parameters for real density function
cc=2*a/((1-exp(-a*h))*c^2);
q=2*a*b/c^2-1;
density1=cc*exp(-cc*(x+exp(-a*h)*xs))*(x/(exp(-a*h)*xs))^(q/2)*besseli(q,2*cc*sqrt(x*exp(-a*h)*xs));
v3=subs(density1,{a,b,c,h,xs},{1,1,2,1/250,1});
ezplot(v3,[-1,5]);
% Error term
ezplot(v3-cir,[-1,5]);
