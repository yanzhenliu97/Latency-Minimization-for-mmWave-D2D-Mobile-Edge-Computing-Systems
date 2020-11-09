function result = f(Vout,Pmax)

if Vout <= sqrt(Pmax)/2 & Vout>0
    result = 2*Vout*sqrt(Pmax)/pi;
elseif Vout >sqrt(Pmax)/2 & Vout<=sqrt(Pmax)
        result = 6*Vout*sqrt(Pmax)/pi-2*Pmax/pi;
elseif Vout <0
    disp("invalid input of Vout");
else  
    disp("invalid input of Vout");
    result = sqrt(Pmax);
end