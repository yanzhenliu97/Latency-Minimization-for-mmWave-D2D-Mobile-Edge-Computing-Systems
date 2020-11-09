function result = scale_digital_transmitter(Fa,Wabar,Pmax,P)
Na = size(Fa,1);

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
end
result = Wabar;