function [ mPatient ] = patBuilder( Age,Weight,Height,gender, isPatient, opiates )
%PATBUILDER Summary of this function goes here
%   Detailed explanation goes here

mPatient.Age = Age;
mPatient.Weight = Weight;
if nargin>2
    mPatient.Height = Height;
end
if nargin>3
    mPatient.gender = gender;
end
if nargin>4
    mPatient.isPatient = isPatient;
else
    mPatient.isPatient = false;
end
if nargin>5
    mPatient.Opiates = opiates;
else
    mPatient.Opiates = false;
end
end

