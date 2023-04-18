clear all
clc
close

%ode45 can be used to integrate the ode system.
[t,X]=ode45(@odes,[0 100],[1 1 0.25]);
c = X(:,1); %Variables in terms of t are set.
h = X(:,2);
p = X(:,3);

%Solutions are plotted.
figure 
plot(t,c, 'LineWidth',2);
hold on;
xlabel('$t$','interpreter','latex');
plot(t,h,'-', 'LineWidth',2);
hold on;
plot(t,p,'-', 'LineWidth',2);
l=legend('$c$ : Ca$^{2+}$ concentration','$h$ : Fraction of active IP$_3$Rs', '$p$ : IP$_3$ concentration');
set(l, 'interpreter', 'latex')
%hold off;
set(gca, "FontSize", 16)
exportgraphics(gca,'p_cos_1.png','Resolution',300)

%Now the model is actually introduced.
function Fn = odes(t,x)

K_act=210*(10^(-9));
H_IP3=4.0;
K_infty=52*(10^(-6));
K_flux=4.9*(10^(-6));
V_e=10^(-6);
K_e=0.1*(10^(-6));
g=0.51;

K_inh=K_infty*(x(3)^(H_IP3)/(x(3)^(H_IP3)+1));

K1=K_flux/V_e;
K2=(K_act*g)/V_e;
K3=(g*K_e)/V_e;
K4=V_e/(g*K_inh);

F=K1*(x(2))*((x(1))^2/(K2^2+(x(1))^2))-((x(1))^2/(K3^2+(x(1))^2)); %The derivative of c with respect to t. 
G=(1/(1+(K4*(x(1)))^4))-(x(2)); %The derivative of h with respect to t.
H= 0.1*cos(t); %The derivative of p with respect to t. In this case, it is a cos curve.


Fn=[F; G; H];
end