%function [uplink_rate,uplink_rate2] = uplink_new_version(Nu,Na,Nrf1,H1,sigma)
% %这是调试部分＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝
%矩阵维度参数
Nu=32;
Na=8;
Nb=8;
Nrf1=4;
Nrf2=4;

%L是路径数
L=15; 
%噪声方差
sigma=0.1; 
%功率p
p=1;
Pmax=1000;
%产生信道
H1=generateH(Nu,Na,L);
H2=generateH(Nb,Nu,L);
%调试部分结束＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝

%矩阵初始化
V1=rand(Nrf1,1)-0.5
U1=exp(2*pi*j*(rand(Nrf1,Nu)-0.5))
Ua=exp(2*pi*j*(rand(Na,1)-0.5));
V1power=norm(V1,'fro')^2
outer_max=60;
middle_Nmax=100;
inner_Nmax=50;
opX=V1'*U1;

b=V1'*U1
barU=U1;
barUa=Ua;


%参数设置
rho=1e-10;c=2;delta=1e-8;
stop_cre=zeros(1,outer_max);
object=zeros(1,outer_max);

for outer_i=1:(outer_max)
    outer_i
    rho=c*rho;
    object(outer_i)=log2(1+abs(V1'*U1*H1*Ua)^2/(norm(V1'*U1,'fro')^2*sigma^2));
    inner_Nmax=ceil(1000/outer_i)+20;
    for i=1:inner_Nmax    
    %z===================
    z=b*H1*barUa/(abs(b*H1*barUa)^2+norm(b,'fro')^2*sigma^2);
    
    %b===================
    b=(rho*V1'*barU+z*barUa'*H1')*pinv(abs(z)^2*H1*Ua*Ua'*H1'+(abs(z)^2*sigma^2+rho)*eye(size(H1*H1')));

    
    %V1===================
    V1=pinv(barU*barU')*barU*b';
    
    %barU=================
    barU=pinv(V1*V1'+eye(size(V1*V1')))*(U1+V1*b);
    
    %barUa================
    barUa=pinv(rho*eye(size(H1'*H1))+abs(z)^2*H1'*b'*b*H1)*(rho*Ua+z*H1'*b');
    
   
    
    %Ua===================
    Ua=barUa./abs(barUa);
    
        
    %U1===================
    U1=barU./abs(barU);
    %V1 normalization===================
    %V1=V1*(sqrt(V1power)/norm(V1,'fro'));

    end
    
    constrain1=max(max(U1-barU));
    constrain2=max(max(barUa-Ua));
    constrain3=max(max(b-V1'*barU));
    stop_cre(outer_i)=norm([constrain1,constrain2,constrain3],inf);
    
    if outer_i<=1
        difference= object(1);
    else
        difference= object(outer_i)-object(outer_i-1);
    end

    
end
subplot(2,1,1)
plot(stop_cre,'r');
subplot(2,1,2)
plot(object,'g');