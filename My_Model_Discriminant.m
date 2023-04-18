clear all
clc
close

K_act=210*(10^(-9));
H_IP3=4.0;
K_infty=52*(10^(-6));
K_flux=4.9*(10^(-6));
V_e=10^(-6);
K_e=0.1*(10^(-6));
g=0.51;

K1=K_flux/V_e;
K2=(K_act*g)/V_e;
K3=(g*K_e)/V_e;

c=0:0.001:0.5;

K4=@(c) (((K1/(c^4))*(((K3^2)+(c^2))/((K2^2)+(c^2)))-(1/(c^4)))^(1/4));
p=@(c) ((V_e/((g*K_infty*K4(c))-V_e))^(1/4));
h=@(c) (1/(1+(K4(c)*c)^4));

Tr=@(c) ((2.*K1.*c)./((K2.^2)+(c.^2))*(1./(1+((K4(c).*c).^4)))*((K2.^2)./((K2.^2)+(c.^2)))-((2.*c)./((K3.^2)+(c.^2)))*((K3.^2)./((K3.^2)+(c.^2)))-1);
Det=@(c) ((2.*K1.*c)./((K2.^2)+(c.^2))*(1./(1+((K4(c).*c).^4)))*(((c^2/(K2^2+c^2))-1))+((2.*c)./((K3.^2)+(c.^2)))*(1-((c^2)./((K3.^2)+(c.^2)))) + ((K1*(c^2))/(K2^2 + c^2))*(1/((1+(K4(c)*c)^4)^2))*(4*(K4(c)^4)*(c^3)));
Disc=@(c) ((Tr(c)^2)-(4*Det(c))); %The trace and determinant are both used to help find the discriminant of the Jacobian. 

x0=0.001;
x1=0.4;

options = optimset('Display','iter');

[c1,fval,exitflag,output] = fzero(Disc,x0, options); %The discriminant is solved when equal to zero.
[c2,fval,exitflag,output] = fzero(Disc,x1, options); 

p_left=p(c1) %By substituting the steady states of c into the equation for p, the bifurcation points are found.
p_right=p(c2) %In this case, there is only one bifurcation point.
