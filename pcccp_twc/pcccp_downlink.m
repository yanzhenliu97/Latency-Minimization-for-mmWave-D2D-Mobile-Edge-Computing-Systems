%function [downlink_rate,downlink_rate2] = pcccp_downlink(Nu,Nb,Nrf2,H2,Pmax,sigma)
% % 这是调试部分＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝
%矩阵维度参数
Nu=32;
Na=8
Nb=8;
Nrf1=4;
Nrf2=4;
level=4
% L是路径数
L=15; 
% 噪声方差
sigma=10^(-6); 
% 功率p
p=1;
Pmax=1000;
% 产生信道
H1=mmWavechannel(Nu,Na,L,120);
H2=mmWavechannel(Nb,Nu,L,120);


% % 调试部分结束＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝
%归一化
H2=H2*1e7;
sigma=sigma*1e7;


%矩阵初始化
V2=rand(Nrf2,1)-0.5;
U2=exp(2*pi*j*(rand(Nu,Nrf2)-0.5));
V2=0.6*V2*(sqrt(Pmax/norm(U2*V2,'fro')^2));
Ub=exp(2*pi*j*(rand(Nb,1)-0.5));
start_power=norm(U2*V2,'fro')^2;
outer_Nmax=100;
middle_Nmax=20;

opX=U2*V2;
barU=U2;    
%初始化B
B=U2*V2;
%参数设置
rho=0.01;c=2;delta=1e-10;
stop_cre=zeros(middle_Nmax,outer_Nmax);
object=zeros(middle_Nmax,outer_Nmax);
inner_object=zeros(1,middle_Nmax);
inner_penalty=zeros(1,middle_Nmax);

for outer_i=1:outer_Nmax
    Bi=B;
    %pcccp==================
    rho=c*rho;
    for i=1:middle_Nmax 
        object(i,outer_i)=log2(1+norm(Ub'*H2*U2*V2,'fro')^2/(Nb*sigma^2));
        
        inner_object(i)=log2(1+norm(Ub'*H2*U2*V2,'fro')^2/(Nb*sigma^2));
        %Ub================
        temp=H2*U2*V2;
        Ub=temp./abs(temp);
        %B========================
        B0=barU*V2+1/rho*(H2'*Ub*Ub'*H2*Bi);
        B=sqrt(Pmax)*B0/(max(0,sqrt(Pmax)-norm(B0,'fro'))+norm(B0,'fro'));
        
        %barU=====================
        barU=(B*V2'+U2)*inv(V2*V2'+eye(size(V2*V2')));
        
        %U2=======================
        [UM,UN]=size(U2);
        for Ui=1:UM
            for Uj=1:UN
                   U2(Ui,Uj)=barU(Ui,Uj)/abs(barU(Ui,Uj));
            end
        end
        %V2=======================
        V2=pinv(barU)*B;
            
        
        
        constrain1=max(max(B-barU*V2));
        constrain2=max(max(barU-U2));
        inner_penalty(i)=norm([constrain1,constrain2],inf);   
        stop_cre(i,outer_i)=norm([constrain1,constrain2],inf);
    end
 
     if stop_cre(i,outer_i)<delta && (object(i,outer_i)-object(i-1,outer_i)<delta)
         outer_i;
         break;
     end
    
    
    
    end
    
    

for i=1:20
   temp=H2*opX;
    Ub=temp./abs(temp);

       tempX=Ub'*H2;
       for xi=1:length(tempX)
           opX(xi)=sqrt(Pmax)*tempX(xi)'/norm(tempX,'fro');
       end
      
end
log2(1+norm(Ub'*H2*opX,'fro')^2/(Nb*sigma^2));
%object=reshape(object,1,[]);
%stop_cre=reshape(stop_cre,1,[]);
downlink_rate=log2(1+norm(Ub'*H2*U2*V2,'fro')^2/(Nb*sigma^2));

figure
%plot(object(1:(outer_i-1)*middle_Nmax),'r','linewidth',2)
plot(object(1,(1:outer_i)),'r','linewidth',2)
% set(gca,'XLim',[0 1000]);
% set(gca,'xTick',[0:200:1000])
xlabel('Number of iterations')
ylabel('Objective function (dB)')
set(get(gca,'XLabel'),'Fontsize',14)
set(get(gca,'YLabel'),'Fontsize',14)
grid on

figure
%semilogy(stop_cre(1:(outer_i-1)*middle_Nmax),'b','linewidth',2)
semilogy(stop_cre(1,(1:outer_i)),'b','linewidth',2)
% set(gca,'XLim',[0 1000]);
% set(gca,'xTick',[0:200:1200])
set(gca,'YLim',[1e-9 1e0 ]);

xlabel('Number of iterations')
ylabel('Penalty terms')
set(get(gca,'XLabel'),'Fontsize',14)
set(get(gca,'YLabel'),'Fontsize',14)
grid on

 [U,diag,V]=svd(H2);
Ub=U(1,:)';
Ub=Ub.\(abs(Ub));
U2=V(:,1:Nrf2);
U2=U2.\(abs(U2));

V2=sqrt(Pmax)*(U2'*H2'*Ub)/norm(U2*U2'*H2'*Ub,'fro');
downlink_rate2=log2(1+norm(Ub'*H2*U2*V2,'fro')^2/(Nb*sigma^2));