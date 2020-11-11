function [result]=mmWavechannel(M,N,L,d)
tempa=0:M-1;
tempb=0:N-1;

%Los=================
%Los path loss 
alpha=72;
belta=2.92;
sigma=10^(0.87);
pl=alpha+belta*10*log10(d)+sigma*randn(1);
pl=-pl;

%channel matrix 
a=1/M*exp(j*pi*tempa*sin(pi/2));
a=a';
b=1/N*exp(j*pi*tempb*sin(pi/2));
result=10^(pl/20)*a*b;


%Nlos==================
%path loss
for k=1:L-1    
    %nLos path loss
    alpha=61.4;
    belta=2;
    sigma=10^(0.58);
    pl=alpha+belta*10*log10(d)+sigma*randn(1);
    pl=-pl;
    
    %channel matrix
    bias=pi*(rand(1)-0.5); %angle voilation
    a=1/M*exp(j*pi*tempa*sin(pi/2+bias));
    a=a';
    b=1/N*exp(j*pi*tempb*sin(pi/2+bias));
   
    result=result+10^(pl/20)*a*b;
end
result=sqrt(M*N/L)*result;
