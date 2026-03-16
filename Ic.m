function I = Ic(T,Ic0,Tc)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
I=Ic0*(1-(T/Tc)^2)^2;
if T>Tc
    I=0;
end
end