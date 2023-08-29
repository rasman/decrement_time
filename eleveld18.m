function [sys, volume, clearance] = eleveld18(mPatient, stepSize)
% Implements the model of Eleveld as described in page 942 of :
% Eleveld DJ, Colin P, Absalom AR, Struys MMRF. 
% Pharmacokinetic-pharmacodynamic model for propofol for broad application
% in anaesthesia and sedation. Br J Anaesth. 2018 May;120(5):942-959

% 70 kg male, 35 yr of age, and 170 cm tall 


mPatientRef = patBuilder(35,70,170,1); %70-kg, 35-year-old, 170-cm male patient.

wgt = mPatient.Weight;
age = mPatient.Age;
isOpiates = mPatient.Opiates;

th= [6.28; 25.5; 273; 1.79; 1.83; 
    1.11; 0.191; 42.3; 9.06; -0.0156; 
    -0.00286; 33.6; -0.0138; 68.3 ; 2.10; 
    1.30; 1.42; 0.68];
thke0 = exp(-1.921620);
%omega=[0.610; 0.565; 0.597; 0.265; 0.346; 0.209; 0.463];
eta = zeros(1,6);

    function f = alsallami(mPatient)
        BMI = mPatient.Weight/((mPatient.Height/100)^2);
        if mPatient.gender
            f = (0.88 + (1-0.88)/(1+((mPatient.Age/13.4)^-12.7)))*((9270*mPatient.Weight)/(6680 + 216*BMI));
        else
            f = (1.11 + (1-1.11)/(1+((mPatient.Age/7.1)^-1.1)))*((9270*mPatient.Weight)/(8780 + 244*BMI));
        end
    end
    
    function f = sigmoid(x, E50, lambda)
        f = (x.^lambda)./(x.^lambda + E50.^lambda);
    end

    function f = aging(x, age)
        f = exp(x*(age - 35));
    end

    function f = central(x)
        f = sigmoid(x, th(12), 1);
    end

    function f = CLmaturation (PMA)
        f =sigmoid(PMA*52,th(8),th(9));
    end

    function f = Q3maturation (Age)
        f =sigmoid(Age*52+40,th(14),1);
    end

    function f = opiates (x, age, isOpiate)
        if isOpiate
           f =exp(x*age);
        else
            f = 1;
        end
    end

V1art= th(1)*(central(wgt)/central(mPatientRef.Weight))*exp(eta(1));
%V1ven = V1art * (1 + th(17)*(1-central(wgt)));
V1 = V1art;
V2 = th(2)*(wgt/mPatientRef.Weight)*aging(th(10), age)*exp(eta(2));
V3 = th(3)*(alsallami(mPatient)/alsallami(mPatientRef))*opiates(th(13), age, isOpiates)*exp(eta(3));

volume = [V1 V2 V3];

if mPatient.gender
    var1 = th(4);
else
   var1 = th(15);
end
CL1 = var1*((wgt/mPatientRef.Weight)^0.75)*...
    CLmaturation(age+40/52)/CLmaturation(35+40/52) *...
    opiates(th(11), age, isOpiates)*exp(eta(4));

CL2art = th(5)*((V2/th(2))^0.75)*(1+th(16)*(1-Q3maturation(age)))*exp(eta(5));
%CL2ven = CL2art*th(18);
CL2 = CL2art;
CL3 = th(6)*((V3/th(3))^0.75)*(Q3maturation(age)/Q3maturation(35))*exp(eta(6));
clearance = [CL1 CL2 CL3];

ke0 = thke0*((wgt/70)^-0.25)*exp(eta(2));

sys = mam2ss(clearance,volume, ke0);
if (nargin==3)
    sys=c2d(sys, stepSize,'zoh');
end
end