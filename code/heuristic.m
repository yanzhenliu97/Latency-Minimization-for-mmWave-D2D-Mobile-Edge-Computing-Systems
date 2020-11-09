function [Wabar,heuristic_non_linear_power] = heuristic(H,Fa,U,sigma,P,D,Pmax)
%% a heuristic method for transmitter design
Na = size(H,2);
Hbar = H*Fa;
QU = U*inv(U'*U)*U';
[U_QU,Sigma_QU,V_QU] =svd(QU);
QU_eff = U_QU*sqrt(Sigma_QU);

QF = Fa'*Fa;

Heff = QU_eff'*Hbar*inv(sqrtm(QF));
%Heff = Hbar*inv(sqrtm(QF));
[Wabar,linear_power_rate] = waterfilling(Heff,P,sigma,D);
Wabar = inv(sqrtm(QF))*Wabar;
V = pinv(U'*H*Fa*Wabar*Wabar'*Fa'*H'*U+sigma^2*U'*U)*U'*H*Fa*Wabar;
rate_not_scaled = compute_rate(H,U,Fa,Wabar,sigma);

%scale the power to satisfy the power constraint
Vout_list = sqrt(sum(abs(Fa*Wabar).^2,2));
scale_factor1 = 1;
if max(Vout_list)>sqrt(Pmax)
scale_factor1 = sqrt(Pmax)/max(Vout_list);
end
Vout_list = scale_factor1*Vout_list;

%bisection search the scale factor
totalPApower = 0;
for i =1:Na
    totalPApower = totalPApower + f(Vout_list(i),Pmax);
end

if totalPApower <P
    %this satisfy the constraint and can stop
    Wabar = Wabar*scale_factor1;    
    V = pinv(U'*H*Fa*Wabar*Wabar'*Fa'*H'*U+sigma^2*U'*U)*U'*H*Fa*Wabar;
    heuristic_non_linear_power = compute_rate(H,U,Fa,Wabar,sigma);


else
    %search the scale factor
scale_factor = 1;
scale_delta = 10^-5; %accuracy
while totalPApower > P
    scale_factor = 0.5*scale_factor;
    totalPApower = 0;
    for i =1:Na
    totalPApower = totalPApower + f(Vout_list(i)*scale_factor,Pmax);
    end 
end

min_scale_factor = scale_factor;
max_scale_factor = 2*scale_factor;
while max_scale_factor-min_scale_factor > scale_delta
    scale_factor = (min_scale_factor+max_scale_factor)/2;
    totalPApower = 0;
    for i =1:Na
    totalPApower = totalPApower + f(Vout_list(i)*scale_factor,Pmax);
    end 
    
    if totalPApower>P
        max_scale_factor = scale_factor;
    else
        min_scale_factor = scale_factor;
    end
end
Wabar = Wabar*scale_factor1*scale_factor;
V = pinv(U'*H*Fa*Wabar*Wabar'*Fa'*H'*U+sigma^2*U'*U)*U'*H*Fa*Wabar;
heuristic_non_linear_power = compute_rate(H,U,Fa,Wabar,sigma);
% disp('heuristic+receiver rate is')
% compute_rate1(H,U,Fa,Wabar,V,sigma)
end