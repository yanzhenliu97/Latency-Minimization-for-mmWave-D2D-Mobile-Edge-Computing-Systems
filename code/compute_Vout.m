function result = compute_Vout(a,b,Pmax)

%a heuristic solution
% if a>0 & a<=Pmax/pi
%     result = pi^2*(2*sqrt(Pmax)*a + b)/(pi^2+4*Pmax);
%     result = max(0,result);
% elseif a>Pmax/pi & a<=4*Pmax/pi
%     result = pi^2*(b + 6*sqrt(Pmax)*a/pi +12*Pmax*sqrt(Pmax)/pi^2)/(pi^2+36*Pmax);
%     result = min(result,Pmax);
% else
%     disp("invaild input when compute Vout")
%     result = sqrt(Pmax);
% end

%optimal solution
point1 = pi^2*(2*sqrt(Pmax)*a/pi + b)/(pi^2+4*Pmax);
point1 = tail(point1,0,sqrt(Pmax)/2);
min1 = (a-f(point1,Pmax))^2+(b-point1)^2;

point2 = pi^2*(b + 6*sqrt(Pmax)*a/pi +12*Pmax*sqrt(Pmax)/pi^2)/(pi^2+36*Pmax);
point2 = tail(point2,sqrt(Pmax)/2,sqrt(Pmax));
min2 = (a-f(point2,Pmax))^2+(b-point2)^2;

if min1<min2
    result = point1;
else
    result = point2;
end





end


function result = tail(x,a,b)
if a>b
    disp("invalid input when using tail")
end

if x<a
    result = a;
elseif x<=b
    result = x;
else
    result = b;
end

end


