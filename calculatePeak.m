function [c_ratio, t_peak_calc, num_iter, all_var] = calculatePeak(volume, clearance, ke0, t_peak)
% Calculates R and time to peak effrct from model parameters


k1 = clearance./volume(1);
k2 = clearance(2)/volume(2);

if length(k1)<3
    k1 = [k1 0];
end

r1 = ((sum(k1) + k2) + sqrt((sum(k1) + k2)^2-4*((sum(k1)-k1(2))*k2)))/2;
r2 = ((sum(k1) + k2) - sqrt((sum(k1) + k2)^2-4*((sum(k1)-k1(2))*k2)))/2;

%Hack #2:
if abs((r1-ke0)/ke0)< 0.025
    if r1 > ke0
        ke0 = ke0 *.975;
    else
        ke0 = ke0 *1.025;
    end
    
end
if abs((r2-ke0)/ke0)< 0.025
    if r2 > ke0
        ke0 = ke0 *.975;
    else
        ke0 = ke0 *1.025;
    end
end

try
    [t_peak_calc, ~, ~, output] = fzero(@(t) -r2 *(ke0 - r1)*(k2 - r2)*exp(-r2 *t) + ...
        r1*(ke0 - r2)*(k2 - r1)*exp(-r1 *t)+ ...
        -ke0 *(k2 - ke0)*(r1 - r2)*exp(-ke0 *t), [(log(ke0)-log(r1))/(ke0-r1) (log(ke0)-log(r2))/(ke0-r2)]);
catch
    t_peak_calc = -1;
    output.iterations =0;
end
if nargin <4
    t_peak = t_peak_calc;
end

if t_peak_calc>0
    sol = (ke0)*((k2 - r2)/((ke0 - r2)*(r1 - r2))*exp(-r2 *t_peak) + ...
        (k2 - r1)/((ke0 - r1)*(r2 - r1))*exp(-r1 *t_peak)+ ...
        (k2 - ke0)/((r1 - ke0)*(r2 - ke0))*exp(-ke0 *t_peak));
    
    c_ratio = 1/sol;
    num_iter = output.iterations;
    
    all_var.k1 =k1;
    all_var.k2 =k2;
    all_var.r1 =r1;
    all_var.r2 =r2;
    
else
    c_ratio = -1;
    num_iter =-1;
    all_var= [];
end




