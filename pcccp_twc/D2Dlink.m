% function [D2D_rate,D2D_rate2] = D2Dlink(Na,Nb,H3,sigma)
%这是调试部分＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝
%矩阵维度参数
sigma=10^(-6); 
L=15;
Na=8;
Nb=8;
H3=mmWavechannel(Nb,Na,L,20);
%调试部分结束====================================
%归一化
H3=H3*1e6;
sigma=sigma*1e6;

Ua=exp(2*pi*j*(rand(Na,1)-0.5));
Ub=exp(2*pi*j*(rand(Nb,1)-0.5));

Nmax=100;
Fx=zeros(1,Nmax);
D2D_rate0=log2(1+norm(Ub'*H3*Ua,'fro')^2/(Nb*sigma^2));
for i=1:Nmax
    Fx(i)=norm(Ub'*H3*Ua,'fro')^2;
    temp=H3*Ua;
    Ub=temp./abs(temp);
    
    temp1=H3'*Ub;
    Ua=temp1./abs(temp1);
   
end
D2D_rate=log2(1+Fx(Nmax)/(Nb*sigma^2));

plot(log2(1+Fx/(Nb*sigma^2)))
x=1:8;
plot(x,log2(1+Fx(1:8)/(Nb*sigma^2)),'r','linewidth',2);
xlabel('Number of iterations')
ylabel('Objective function (dB)')
set(get(gca,'XLabel'),'Fontsize',16)
set(get(gca,'YLabel'),'Fontsize',16)
grid on

 [U,diag,V]=svd(H3);
Ub=U(1,:)';
Ub=Ub.\(abs(Ub));
Ua=V(:,1);
Ua=Ua.\(abs(Ua));
D2D_rate2=log2(1+norm(Ub'*H3*Ua,'fro')^2/(Nb*sigma^2));