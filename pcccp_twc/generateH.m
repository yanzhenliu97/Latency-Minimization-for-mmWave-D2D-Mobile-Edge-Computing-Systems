function [result]=generateH(M,N,L)
tempa=0:M-1;
tempb=0:N-1;
result=zeros(M,N);
for k=1:L
    a=1/sqrt(M)*exp(j*pi*tempa*sin(2*pi/L*k));
    a=a';
    b=1/sqrt(N)*exp(j*pi*tempb*sin(2*pi/L*k));
    c=randn;
    d=randn;
    result=result+(c+j*d)*a*b;
end
result=sqrt(M*N/L)*result;
