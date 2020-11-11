function [optimal_rho,optimal_time,condition]=best_partition(k1,ke,k2,kL,k3)
if k1*k3>=kL*ke
    if (kL+k3)>=k2
        if k1+k2<=k3
            optimal_rho=k3/(ke+k3);
            optimal_time=k3*(k1+ke+k2)/(ke+k3);
            condition=1;
        end
        if k1+k2>k3
             optimal_rho=kL/(k1+kL);
             optimal_time=(k1*kL+k1*k3+kL*k2)/(k1+kL);
             condition=2;
        end
    end
    
    if (kL+k3)<k2
        if k1+k2>k3
            optimal_rho=0;
            optimal_time=kL+k3;
            condition=3;
        end
    end
end
if  k1*k3<kL*ke
    if kL+k3>=k2
        optimal_rho=(kL+k3)/(kL+k3+k1+ke);
        optimal_time=(kL+k3)*(k1+ke+k2)/(kL+k3+k1+ke);
        condition=4;
    end
    if kL+k3<k2
        optimal_rho=0;
        optimal_time=kL+k3;
        condition=5;
    end
end
        
    
end