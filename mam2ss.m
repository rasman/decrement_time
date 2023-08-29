function [sys] = mam2ss(clearance,volume, ke0)
% returns a state space model from vectors of clearance
% and compartment volumes. All states are drug amount
% in the compartment; to get concentration, divide by compartment
% volume in the C matrix. Note that this code can handle an
% arbitrary number of compartments, but assumes the first compartment is
% the central compartment (mam is short for mamillary).
% The effect site compartment is represented as a compartment of negligible
% volume adjoined to the central compartment.

n = length(clearance);
if (nargin == 2)
    k1 = clearance./volume(1);
    k2 = clearance(2:n)./volume(2:n);
    A = [-sum(k1) k2;k1(2:n)' -diag(k2)];
    B = [1;zeros(n-1,1)];                   % Direct addition of drug into the central compartment
    C = [1/volume(1) zeros(1,n-1)];     % Observation of the central compartment concentration.
else
    EFFECT_VOL_FACTOR=10000; % ratio of central compartment vol to effect
    k1 = [clearance./volume(1) ke0/EFFECT_VOL_FACTOR];
    k2 = [clearance(2:n)./volume(2:n) ke0];
    A = [-sum(k1) k2;k1(2:n+1)' -diag(k2)];
    B = [1;zeros(n,1)];                   % Direct addition of drug into the central compartment
    C = [zeros(1,n) EFFECT_VOL_FACTOR/volume(1)];     % Observation of the effect site compartment concentration.
end


C = C./1000; % All compartment volumes are in L, convert to ml
D=0;

sys = ss(A,B,C,D);
end