%function [uplink_rate,uplink_rate2] = pcccp_uplink(Nu,Na,Nrf1,H1,sigma)
%这是调试部分＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝
%矩阵维度参数
Nu=32;
Na=8;
Nb=8;
Nrf1=4;
Nrf2=4;

%L是路径数
L=15; 
%噪声方差
sigma=10^(-6); 
%功率p
p=1;
Pmax=1000;
%产生信道
H1=mmWavechannel(Nu,Na,L,120);
H2=mmWavechannel(Nb,Nu,L,120);
%调试部分结束＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝

%归一化
H1=H1*1e6;
sigma=sigma*1e6;

%矩阵初始化
V1=rand(Nrf1,1)-0.5;
U1=exp(2*pi*j*(rand(Nrf1,Nu)-0.5));
Ua=exp(2*pi*j*(rand(Na,1)-0.5));
V1power=norm(V1,'fro')^2;
outer_Nmax=100;

inner_Nmax=20;
opX=V1'*U1;
barU=U1;
barUa=Ua;
stop_cre=zeros(inner_Nmax,outer_Nmax);
object=zeros(inner_Nmax,outer_Nmax);

%参数设置
rho=1e-7;c=2;delta=1e-8;


for outer_i=1:(outer_Nmax)
   
    rho=c*rho;
    
    for i=1:inner_Nmax
    object(i,outer_i)=log2(1+abs(V1'*U1*H1*Ua)^2/(norm(V1'*U1,'fro')^2*sigma^2));
    %Ua==================
    temp=H1'*U1'*V1;
    Ua=temp./abs(temp);
    barUa=Ua;
    %z===================
    z=V1'*barU*H1*barUa/(abs(V1'*barU*H1*barUa)^2+norm(V1'*barU,'fro')^2*sigma^2);
    
    %V1===================
    V1=inv(abs(z)^2*barU*H1*barUa*barUa'*H1'*barU'+abs(z)^2*sigma^2*barU*barU')*z'*barU*H1*barUa;
    
    %barU=================
    A=abs(z)^2*V1*V1';
    B=inv(H1*barUa*barUa'*H1'+sigma^2*eye(size(H1*H1')));
    C=(z*V1*barUa'*H1'+rho*U1)*B;
    barU=sylvester(A,rho*B,C);
    
    %barUa================
    %inv(rho*eye(size(H1'*H1))+abs(z)^2*H1'*barU'*V1*V1'*barU*H1)*(rho*Ua+z*H1'*barU'*V1);
    
    %U1===================
    U1=barU./abs(barU);
    
%     %Ua===================
%     Ua=barUa./abs(barUa);
 
    %V1 normalization===================
    %V1=V1*(sqrt(V1power)/norm(V1,'fro'));
    constrain1=max(max(U1-barU));
    constrain2=max(max(barUa-Ua));
     
    stop_cre(i,outer_i)=norm([constrain1,constrain2],inf);
    end
    
     if stop_cre(i,outer_i)<delta && (object(i,outer_i)-object(i-1,outer_i)<delta)
         outer_i;
         break;
    end


    
end
uplink_rate=log2(1+abs(V1'*U1*H1*Ua)^2/(norm(V1'*U1,'fro')^2*sigma^2));
%object=reshape(object,1,[]);
%stop_cre=reshape(stop_cre,1,[]);

figure
%plot(object(1:outer_i*inner_Nmax),'r','linewidth',2)
plot(object(1,(1:outer_i)),'r','linewidth',2)
% bound=ceil(outer_i*inner_Nmax/200)*200;
% set(gca,'XLim',[0 bound]);
% set(gca,'xTick',[0:200:bound])
% set(gca,'YLim',[0 10]);
% set(gca,'yTick',[0:1:10])
xlabel('Number of iterations')
ylabel('Objective function (dB)')
set(get(gca,'XLabel'),'Fontsize',14)
set(get(gca,'YLabel'),'Fontsize',14)
grid on


figure
%semilogy(stop_cre(1:outer_i*inner_Nmax),'b','linewidth',2)
semilogy(stop_cre(1,(1:outer_i)),'b','linewidth',2)
% set(gca,'XLim',[0 bound]);
% set(gca,'xTick',[0:200:bound])
% set(gca,'YLim',[1e-8 1e0 ]);

xlabel('Number of iterations')
ylabel('Penalty terms')
set(get(gca,'XLabel'),'Fontsize',14)
set(get(gca,'YLabel'),'Fontsize',14)
grid on

[U,diag,V]=svd(H1);
U1=U(:,1:Nrf1)';
U1=U1.\(abs(U1));
Ua=V(:,1);
Ua=Ua.\(abs(Ua));
V1=inv(U1*H1*Ua*Ua'*H1'*U1'+sigma^2*U1*U1')*U1*H1*Ua;
uplink_rate2=log2(1+abs(V1'*U1*H1*Ua)^2/(norm(V1'*U1,'fro')^2*sigma^2));