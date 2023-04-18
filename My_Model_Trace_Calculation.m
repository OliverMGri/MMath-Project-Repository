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

c=0:0.001:0.5; %The variable c is given a range of values so that it won't become extremely large.

K4=@(c) (((K1/(c^4))*(((K3^2)+(c^2))/((K2^2)+(c^2)))-(1/(c^4)))^(1/4)); %K4 is re-written in terms of c, as noted in the project.
p=@(c) ((V_e/((g*K_infty*K4(c))-V_e))^(1/4)); %p is rearranged as a function of K4.
h=@(c) (1/(1+(K4(c)*c)^4)); %The nullcline of dh/dt.

Tr=@(c) ((2.*K1.*c)./((K2.^2)+(c.^2))*(1./(1+((K4(c).*c).^4)))*((K2.^2)./((K2.^2)+(c.^2)))-((2.*c)./((K3.^2)+(c.^2)))*((K3.^2)./((K3.^2)+(c.^2)))-1); %The trace of the Jacobian is now a function of c.

x0=0.02;
x1=0.30; %These initial conditions are provided to help fzero identify the values of c when the trace equals zero (the steady states).

options = optimset('Display','iter'); %Introducing this option allows the viewer to see the iterations that fzero makes, which makes it easier to see mistakes.

[c1,fval,exitflag,output] = fzero(Tr,x0, options); %The steady states of c are found by solving the trace of the Jacobian when it is equal to zero. This is done for both initial conditions.
[c2,fval,exitflag,output] = fzero(Tr,x1, options); 

c1
c2
p_left=p(c1) %By substituting the steady states of c into the equation for p, the bifurcation points are found.
p_right=p(c2) %The second bifurcation is found here. It is the larger of the two.
