function [result]= compute_total_delay(k1,ke,k2,kL,k3,rho)
tlc=(1-rho)*kL;
tec=rho*ke;
tupt=rho*k1;
tdownt=rho*k2;
td2dt=(1-rho)*k3;

if tupt>=tlc
    result=tupt+max(tec,td2dt)+tdownt;
else
    result=max(tupt+tec,td2dt+tlc)+tdownt;
end