function y = ct(T,Tc,I_0,Ic0)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
y=zeros(size(T));
kB=1.380649*10^(-23);%m2kg/s2K
% Tc=10.5;%K
%ces=2400J/m3K @T=10K
Delta_E=2.25*10^(-3)*1.60218*10^(-19);%J
A=2400/exp(-Delta_E/kB/10);

Delta_T=@(T1)1.76*kB*Tc*tanh(pi/1.76*sqrt(2/3*1.43*(Tc/T1-1)));
% Delta_T=@(T1)2.15*kB*Tc*(1-(T1/Tc)^2);
% A2=2400/exp(-Delta_T(10)/kB/10);
% %ces(Tc)=2.43cen(Tc)
% gamma=A2*exp(-Delta_T(Tc)/kB/Tc)/Tc/2.43;

gamma=1303/8.6;%2400/10;
A2=2.43*gamma*Tc;

for j=1:size(T,1)
    Ic_0=Ic(T(j),Ic0,Tc);%uA
    if I_0>=Ic_0
        Ce=gamma*T(j) ;
    else
        Ce=A2*exp(-Delta_T(T(j))/kB/T(j));
    end
    Cp=12.849*T(j)^3;

    y(j)=Ce+Cp;
end
%c phonon

end







































































































































