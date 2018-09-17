%Financial Econometrics HW 4

%Easy way for Hermite
global z;
syms z;
H0 = sym('1');
H1 = simplify(z*H0 - diff(H0,z));
H2 = simplify(z*H1 - diff(H1,z));
H3 = simplify(z*H2 - diff(H2,z));
H4 = simplify(z*H3 - diff(H3,z));
H5 = simplify(z*H4 - diff(H4,z));
H6 = simplify(z*H5 - diff(H5,z));
H7 = simplify(z*H6 - diff(H6,z));

%The main program
%Define variables
syms a b c;
syms xs ys zs;
syms x y z;
syms h t s;

K = 5;
J = K + 1;

%Drift and diffusion for X(t)
%Vasicek
muX = a*x;
sigmaX = c;

%Black-Scholes 
muX = a*x;
sigmaX = b*x;

%CIR 
muX = a*x;
sigmaX = c*sqrt(x);

%Transformation X(t) to Y(t)
fX2Y = int(1/sigmaX,x);
%simplify(finverse(fX2Y));
fY2X = subs((finverse(fX2Y)),x,y);

%Drift and diffusion for Y(t)
muY_temp = muX/sigmaX - sym('1')/sym('2')*diff(sigmaX,x,1);
muY = subs(muY_temp,x,fY2X);
muY = simplify(muY);
sigmaY = sym('1');

%Transform Y(t) to Z(t)
fY2Z = h^ (-1/2)*(y-ys);
% subs(fY2Z,y,ys);
fZ2Y = h^ (1/2)* z + ys;
% subs(fZ2Y,z,zs);

%Testing symbolic loop
syms z;
Temp = sym('1');
for n = 1:10
    Temp = Temp + z^n;
end

%Generating beta
syms Htemp Exceptation Beta;
clear Beta Htemp Expectation;

for n= 1:K
    HTemp = subs(Hermite(n),z,fY2Z);
    Expectation = HTemp;
    
    for k = 1:J
        HTemp = muY*diff(HTemp,y,1) + sym('1')/sym('2')*diff(HTemp,y,2);
        Expectation = Expectation + h ^ k /factorial(k)*HTemp;
    end
    Beta{n} = sym('1')/factorial(n-1) * subs(Expectation,y,ys);
end

%subs(Beta{3},{a,b,c,h,ys},{1,1,2,1/250,1});

%Generating pZ with Loop
pZ = sym('0');

for m= 1:K
    pZ = pZ + Beta{m}*Hermite(m);
end
%symvar(pZ);

%Generating pZ without loop

pZ = Beta{1}*Hermite(1) + Beta{2}*Hermite(2)+Beta{3}*Hermite(3) + Beta{4}*Hermite(4) + Beta{5}*Hermite(5);
subs(Beta{5}, {a,b,c,h,ys},{1,1,2,1/250,1});

%Generating pY pX
pZ = exp(-z^2/2)/sqrt(2*pi)*pZ;
pY = (h^(-1/2))*subs(pZ,z,fY2Z);
pX = (sigmaX ^ (-1))*subs(pY,y,fX2Y);
pX = subs(pX,ys,subs(fX2Y,x,xs));
pX = simplify(pX);

%Plotting pX
%Plotting Black Scholes
g1 = subs(pX, {a,b,h,xs},{1,1,3,1});
ezplot(g1,[0.01,5]);

g1 = subs(pX, {a,b,h,xs},{1,1,0.25,1});
ezplot(g1,[0.01,2]);

g1 = subs(pX,{a,b,h,xs},{1,1,0.25,1});
ezplot(g1,[0.001,5]);
%explot(log(g1),[0.5,2]);

%Plotting Exact Black-Scholes Density
density = (2*pi*b^2*h)^(-1/2)/x * exp( - (log(x) - (log(xs)+ (a - b^2/2)*h))^2/(2*b^2*h));
g2 = subs(density,{a,b,h,xs},{0.06,0.2,7,1});
g2 = simplify(g2);
ezplot(g2,[0.01,2]);

%Ploting exact density for Vasicek interest rate model
gamm = sigmaX * sqrt(1 - exp(-2*a*h));
density = (pi*gamm^2/a)^(-1/2)*exp(-(x-b-(xs-b)*exp(-a*h))^2*a/(gamm^2));
g2 = subs(density, {a,b,c,h,xs}, {1,1,2,1/250,1});
g2 = simplify(g2);
ezplot(g2);

%Ploting density difference
gDiff = g1-g2;
ezplot(gDiff,[0.001,2]);
ezplot(log(gDiff));
int(gDiff, [-1,2]);

%Exact density for Black Scholes
muX = a*x;
sigmaX = b * x;

densityBS = (sym('1')/sqrt(2*pi*b^2*h)/x)*exp(-(log(x)-log(xs)-(a-b^2)/2*h)^2/(2*b^2*h));
findsym(densityBS);
g3 = subs(densityBS, {a,b,h,xs},{1,1,1/250,2});
ezplot(g3,[0.5,3]);
solve(g3);

disp(pX);
pdfX = matlabfunction(pX);
disp(pdfX);

%Hermite function 
function [temp] = Hermite(k)
syms z;
H{1} = sym('1');
for n =2:k
    H{n} = simplify(z*H{n-1} - diff(H{n-1},z));
end
temp = H{k};
end
