function [H,H_tau,H_tau_tts,Ar,At]= mmWavechannel(M,N,nLOS_AoA_angle,nLOS_AoD_angle,d,belta,AoA,tau,tau_tts,scale_factor)
%M receiver antenna dim
%N transmit antenna dim
n_clu = size(nLOS_AoA_angle,1);
n_ray = size(nLOS_AoA_angle,2);
L = n_clu*n_ray+1;

fd = 70; %maximum doppler shift
C0 = 1e-3; %pathloss at d=1m
pl = C0*d^(-belta);

tempa=0:M-1;
tempb=0:N-1;
Ar = zeros(M,L);
At = zeros(N,L);

%Los=================
%Los complex gain 
alpha = randn(1)+1j*randn(1);

%channel matrix 
ar=1/M*exp(1j*pi*tempa*sin(AoA));
ar=ar.';
ar_tau = ar*exp(1j*2*pi*fd*tau*cos(AoA));
ar_tau_tts = ar*exp(1j*2*pi*fd*tau_tts*cos(AoA));
Ar(:,1) = ar_tau;

at=1/N*exp(1j*pi*tempb*sin(pi-AoA));
at=at.';
At(:,1) = at;

H = alpha*ar*at';
H_tau = alpha*ar_tau*at';
H_tau_tts = alpha*ar_tau_tts*at'; 
%Nlos==================
%path loss
for i = 1:n_clu
    for j=1:n_ray
    %complex_gain
    alpha = sqrt(0.1)*(randn(1)+1j*randn(1));
    
    %channel matrix
    bias=nLOS_AoA_angle(i,j); %NLOS AoA
    ar=1/M*exp(1j*pi*tempa*sin(bias));
    ar=ar';
    ar_tau = ar*exp(1j*2*pi*fd*tau*cos(bias));
    ar_tau_tts = ar*exp(1j*2*pi*fd*tau_tts*cos(bias));
    Ar(:,j+(i-1)*n_ray) = ar_tau;
    
    bias=nLOS_AoD_angle(i,j); %NLOS AoD
    at=1/N*exp(1j*pi*tempb*sin(bias));
    at=at.';
    At(:,j+(i-1)*n_ray) = at;
        
    H=H+alpha*ar*at';
    H_tau = H_tau + alpha*ar_tau*at';
    H_tau_tts = H_tau_tts + alpha*ar_tau_tts*at';
    end
end

H=sqrt(M*N/L)*H*sqrt(pl)*scale_factor;
H_tau = sqrt(M*N/L)*H_tau*sqrt(pl)*scale_factor;
H_tau_tts = sqrt(M*N/L)*H_tau_tts*sqrt(pl)*scale_factor;
