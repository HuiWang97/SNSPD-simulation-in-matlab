function y = kappa(T,Tc,I_0,Ic0,Rd,d)
%KAPPA Summary of this function goes here
%   Detailed explanation goes here
y=zeros(size(T));
L=2.45*10^(-8);%W omega/K2
for j=1:size(T,1)
    Ic_0=Ic(T(j),Ic0,Tc);%uA
    if I_0>=Ic_0
        y(j)=L*T(j)/Rd/d;
    else
        y(j)=T(j)^2/Tc*L/Rd/d;
    end
end

end

