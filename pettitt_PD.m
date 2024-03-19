
% https://pure.tudelft.nl/ws/portalfiles/portal/62790278/Comparative_analysis_of_nonparametric_change_point_detectors_commonly_used_in_hydrology.pdf


clear; clc

% x = [16 18 23 18 14 ]; n = length(x);
x = [12 14 16 18 23 29 25 21 18 11 6 5 4 2]; n = length(x);


Ut = [];


for t = 1:n
    
    
    sut = 0;
    for i = 1:t
        
        
        sui = 0;
        for j = t+1:n
            uj = sign(x(i) - x(j));
            sui = sui+uj;
        end
        
        sut = sut+sui;
        
        
        
    end
    
    
    
    Ut = [Ut;sut];
end

loc=find(abs(Ut)==max(abs(Ut)));
K=max(abs(Ut));
pvalue=2*exp((-6*K^2)/(n^3+n^2));
a=[loc; K ;pvalue];

close all; figure
plot(x)
hold on
scatter(loc-1,0, 'o', 'filled')
