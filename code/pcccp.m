function [proposed_rate,heuristic_rate,Wa] = pcccp(Nrfa,Na,N,Nrf,P,Pmax,sigma,Fa,U,H)

%%tunning part ================

% N = 64  %BS attenna number
% Nrf = 6  %BS RF chain number
% Na = 8  %user A attenna number
% Nrfa = 2  %user A RF chain number
% Pmax = 20  %PA max power output (dBm)
% Pmax = 10^((Pmax-30)/10);
% P = 30;
% P = 10^((P-30)/10);
% sigma = 1;
% L=15;
% Theta_Fa = 2*pi*rand(Na,Nrfa);
% Theta_U1 = 2*pi*rand(N,Nrf);
% Fa = exp(1j*Theta_Fa);
% U = exp(1j*Theta_U1);
% dist1 = 50;
% AoA1 = 0.1;
% scale_factor = sqrt(10^11);
% tau = 1e-3;
% tts_tau = tau*(Nrf*Nrfa)/(N*Na);
% [H,H_tau,H_tau_tts,Ar,At] = mmWavechannel(N,Na,L,dist1,AoA1,tau,tts_tau,scale_factor);

%%tuning part over ==================================



%% short term variables initialization
D = min(Nrfa,Nrf); %send data number
if 4*Pmax/pi*Na<=P
    P = 4*Pmax/pi*Na;
end
V = randn(Nrf,D);
%proposed heuristic method for digital beamformer design
[Wa,heuristic_rate] = heuristic(H,Fa,U,sigma,P,D,Pmax); 
Wa = Wa*0.8;

Z = randn(D,D); %weight matrix

Ppa = zeros(Na,1);
Vout = zeros(Na,1);
for i = 1:Na
    Vout(i) = norm(Fa(i,:)*Wa,'fro');
    Ppa(i) = f(Vout(i),Pmax);
end
sum(Ppa);
%% coefficient
outer_max = 100; %outer layer max iteration
inner_max = 100;
rho=0.1;c=0.8;
delta1 = 1e-6;
delta2 = 1e-10;
delta3 = 1e-6;
rate=zeros(1,outer_max);
constraint = zeros(1,outer_max);
object = zeros(inner_max,5);
penalty =zeros(outer_max,3);

%% start iteration
for outer_i =1:outer_max
    outer_i;
    rate(outer_i) = compute_rate(H,U,Fa,Wa,sigma);
    for inner_i =1:inner_max
        %V ==================================
        V = pinv(U'*(sigma^2*eye(N)+H*Fa*Wa*Wa'*Fa'*H')*U)*U'*H*Fa*Wa;
        %[object(inner_i,1),penalty(outer_i,1),penalty(outer_i,2),penalty(outer_i,3)] = compute_object(H,U,V,Fa,Wa,Z,sigma,rho,Ppa,Vout,P,Pmax);    
        
        %Z ==================================
        E = eye(D) - V'*U'*H*Fa*Wa;
        Z = inv(E);
        %Z = eye(size(Z));
        %[object(inner_i,2),penalty(outer_i,1),penalty(outer_i,2),penalty(outer_i,3)] = compute_object(H,U,V,Fa,Wa,Z,sigma,rho,Ppa,Vout,P,Pmax);
        
        %Wa =================================
        Wabar = Wa;
        temp =zeros(Nrfa,D);
        for i = 1:Na
           Voutbar = norm(Fa(i,:)*Wabar,'fro');
           temp = temp+Vout(i)/(2*rho*Voutbar)*Fa(i,:)'*Fa(i,:)*Wabar;
        end
        Wa = pinv(Fa'*H'*U*V*Z*V'*U'*H*Fa+1/(2*rho)*Fa'*Fa)*(Fa'*H'*U*V*Z+temp);
        %[object(inner_i,3),penalty(outer_i,1),penalty(outer_i,2),penalty(outer_i,3)] = compute_object(H,U,V,Fa,Wa,Z,sigma,rho,Ppa,Vout,P,Pmax);        
        
        %Vout ===============================
        for i =1:Na
            a = Ppa(i);
            b = norm(Fa(i,:)*Wa,'fro');
            Vout(i) = compute_Vout(a,b,Pmax);
        end
        %[object(inner_i,4),penalty(outer_i,1),penalty(outer_i,2),penalty(outer_i,3)] = compute_object(H,U,V,Fa,Wa,Z,sigma,rho,Ppa,Vout,P,Pmax);
        
        %Ppa ================================
        sumPpa = sum(Ppa);
        for i =1:Na
            sumPpa_i = sumPpa - Ppa(i);
            P1 = P - sumPpa_i;
            P2 = f(Vout(i),Pmax);
            Ppa(i) = max(0,min(4*Pmax/pi,(P1+P2)/2));
            sumPpa = sumPpa_i+Ppa(i);
        end
        [object(inner_i,5),penalty(outer_i,1),penalty(outer_i,2),penalty(outer_i,3)] = compute_object(H,U,V,Fa,Wa,Z,sigma,rho,Ppa,Vout,P,Pmax);
        
        if inner_i >5 & abs(object(inner_i,5)-object(inner_i,5)) < delta1
            break;
        end
    end
    constraint(outer_i) = sum(penalty(outer_i,:),2); 
    rho=c*rho;
    if constraint(outer_i)<delta2 & outer_i>5 & abs(rate(outer_i)-rate(outer_i-1))<delta3
        break;
    end
    
end
proposed_rate = rate(outer_i);

%tuning part ========================================

%% plot figure
% figure()
% plot(real(rate(1:outer_i)),'r','linewidth',2)
% xlabel('Number of iterations')
% ylabel('Objective function (bps/Hz)')
% set(get(gca,'XLabel'),'Fontsize',12)
% set(get(gca,'YLabel'),'Fontsize',12)
% grid on
% 
% figure()
% semilogy(constraint(1:outer_i),'b','linewidth',2)
% xlabel('Number of iterations')
% ylabel('Penalty terms')
% set(get(gca,'XLabel'),'Fontsize',12)
% set(get(gca,'YLabel'),'Fontsize',12)
% grid on

%%tuning part over===========================================

        
        
