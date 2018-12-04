function [peind] = pectoral(I0,roih)
%  function to find the pectoral muscle
%  peind returns the index of point along the pectoral muscle
%
%  Version info:
%  $Rev: 450 $:     Revision of last commit
%  $Author: hsiehm@UPHS.PENNHEALTH.PRV $:  Author of last commit
%  $Date: 2015-11-20 11:34:20 -0500 (Fri, 20 Nov 2015) $:    Date of last commit
%
%  Contact:
%     CBIG Group <software at cbica.upenn.edu>

[dimx, ~] =  size(I0);
I = roih;
emap = edge(I, 'canny', [], 2.0);

theta=linspace(40*pi/180, 80*pi/180, 128);
rho_max= sqrt(sum(size(I).^2));
rho=linspace(-rho_max,-0.25*rho_max,128);
[x, y] = ind2sub(size(emap), find(emap > 0));
x0 = x - size(I,1);
y0 = y - size(I,2);

rho_exact = x0 * cos(theta) + y0 * sin(theta);
A = histc(rho_exact, rho);

[j, i] = ind2sub(size(A), find(A >= .98* max(A(:))));  %0.98

for k = 1:length(i)
    tk = theta(i(k));
    rk = rho(j(k));
    for xa = 1:dimx;
        x0a = xa - size(I,1);
    end 
end
for l = 1:dimx
    l0 = l - size(I,1);
    peind(l)= 0;
    for k = 1:length(i)
        tk = theta(i(k));
        rk = rho(j(k));
        if (rk - l0 * cos(tk))/sin(tk) + size(I,2)> peind(l)
            peind(l) = (rk - l0 * cos(tk))/sin(tk) + size(I,2);
        end
    end
end

%clear temporary variables
clear emap y y0 x theta rho_exact rho x0 xa x0a tk rk A l0
