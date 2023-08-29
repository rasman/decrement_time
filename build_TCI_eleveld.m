function input =  build_TCI_eleveld (mPatient, plan, end_time, display_out)

cycle_duration = 2.5; % each cycle is 10 minutes
timeMax = end_time/cycle_duration; % Number of cycles 


[sys, V, Cl] = eleveld18(mPatient); % builds model for demostration
sys.C = [sys.C; 1/V(1)/1000 0 0 0]; % to see the circulating plasma concentration

[c_ratio, t_peak_calc] = calculatePeak(V, Cl, sys.A(1,4)); % this is from my manuscript, returns gain ratio and time to peak effect

tciformula= @(t) (Cl(1) + Cl(2)*exp(-(Cl(2)*t)/V(2)) + Cl(3)*exp(-(Cl(3)*t)/V(3))); % Replace what gets lost
%tciformulaAv= @(t1,t2) Cl(1) - (V(2)*(exp(-(Cl(2)*t2)/V(2))-exp(-(Cl(2)*t1)/V(2))) + V(3)*(exp(-(Cl(3)*t2)/V(3))-exp(-(Cl(3)*t1)/V(3))))/(t2-t1); % Average of TCI formula

time = 0:1/60:end_time;
input = zeros(size(time));

plan(plan<0) = 0;

goal = 0;
for stepper  = 1:size(plan,1)
    f_peak = 0;
    if(plan(stepper,2)>cycle_duration*timeMax)
        continue
    end
    f_p =floor(plan(stepper,2)*60); % index of next transition   
    ceil_val = floor((f_peak +f_p + 1)/(60*cycle_duration))+1; % index of next cycle transition
    ceil_val_index = ceil_val*(60*cycle_duration);

    goaldiff = plan(stepper,1)- goal;
    goal = plan(stepper,1);
    
    if goaldiff >= 0
        f_peak = floor(t_peak_calc*60); % amount of time to set aside following bolus
        
        input(f_p+1) = c_ratio*V(1)*goaldiff*60 + input(f_p+1); %60 is to bolus amount to instant rate in mg/min.
        
        %calculating maintenance for remainder of cycle_duration
        input(f_peak +f_p + 1: ceil_val_index) = goaldiff* mean(tciformula(time(1:ceil_val_index-(f_peak +f_p + 1)))) +...
            input(f_peak +f_p + 1:ceil_val_index);

        %calculating maintenance for all cycle_durations in iteration
        for i = ceil_val:(timeMax-1)
            input((1:60*cycle_duration)+i*60*cycle_duration) = goaldiff* mean(tciformula(time((1:60*cycle_duration)+ i*60*cycle_duration-f_p)))+...
                input((1:60*cycle_duration)+i*60*cycle_duration);
        end   
    else
        remove_drug = -c_ratio*V(1)*goaldiff; %calculate drug to remove
        T_m = find(cumsum(input(f_p+1:end)/60) - remove_drug>0,1,'first'); % replace negative drug amount by infusion pause of equivalent duration

        input(f_p + 1: ceil_val_index) = goaldiff* mean(tciformula(time(1:ceil_val_index-(f_p + 1)))) + input(f_p + 1:ceil_val_index);
        for i = ceil_val:(timeMax-1)
            input((1:60*cycle_duration)+i*60*cycle_duration) = goaldiff* mean(tciformula(time((1:60*cycle_duration)+i*60*cycle_duration-f_p)))+...
                input((1:60*cycle_duration)+i*60*cycle_duration);
        end       
        if isempty(T_m)
            input(f_p+1:end) = 0;
        else
            input(f_p+1:f_p+T_m) = 0;
        end
        
    end
    
end
output = lsim(sys,input,time);
if nargin == 4
if (display_out)
plot(time, output*1e3) % to show plot in mcg/mL
xlabel ('time (min)')
ylabel ('drug concentration (mcg/mL)')
ylim([0 max(plan(:,1))*1.2])

legend({'Effect Site Concentration','Plasma Concentration'})

planVal = zeros(size(time));
for i = 1:size(plan,1)
planVal(find(time>plan(i,2),1,'first'):end) = plan(i,1);
end
figure
plot(time, output(:,1)*1e3- planVal')
ylim([-1 1])
title('Error graph')
xlabel ('time (min)')
ylabel ('absolute error (mcg/mL)')
end
end
end