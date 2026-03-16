function R = Resistance(T,Tc,I,Ic0,Rd,w,Delta_x)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
y=zeros(size(T));
for j=1:size(T,1)
    Ic_0=Ic(T(j),Ic0,Tc);%uA
    if I<Ic_0
        y(j)=0;
    else
        y(j)=1;
    end
end

T_edge=Tc*sqrt(1-sqrt(I/Ic0));
is_not_sc=find(y~=0);
if isempty(is_not_sc)
    R=0;
else
    x_min=min(is_not_sc);
    x_max=max(is_not_sc);
    if x_min~=1
        d_left=(T(x_min)-T_edge)/(T(x_min)-T(x_min-1));
    else
        d_left=1;
    end
    if x_max~=length(T)
        d_right=(T(x_max)-T_edge)/(T(x_max)-T(x_max+1));
    else
        d_right=length(T);
    end
    R=Rd/w*(sum(y)+d_left+d_right-1)*Delta_x;
end
end