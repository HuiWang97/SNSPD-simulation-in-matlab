function y = rho(I,Ic0,T,Tc,Rd,d)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
y=zeros(size(T));
for j=1:size(T,1)
    Ic_0=Ic(T(j),Ic0,Tc);%uA
    if I<Ic_0
        y(j)=0;
    else
        y(j)=Rd*d;
    end
end
end