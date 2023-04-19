clear all
clc
close

list_factory = fieldnames(get(groot,'factory'));
index_interpreter = find(contains(list_factory,'Interpreter'));
for i = 1:length(index_interpreter)
    default_name = strrep(list_factory{index_interpreter(i)},'factory','default');
    set(groot, default_name,'latex');
end %LaTeX form.

veci = 1;

a_var = linspace(0,0.8,50);

cs_1_F=0.25;
cs_2_F=0.05;
v4_1_F=0.29;
v4_2_F=0.18;

roots_1 = [cs_1_F, v4_1_F];
roots_2 = [cs_2_F, v4_2_F];

for a = a_var
    [roots_1, ~, exitflag] = fsolve(@(x)stst_trace(x,a), roots_1);
    Roots_vec_1(veci, :) = roots_1;
    [roots_2, ~, exitflag2] = fsolve(@(x)stst_trace(x,a), roots_2);
    Roots_vec_2(veci, :) = roots_2;
    veci = veci + 1;
end

plot3(a_var, Roots_vec_1(:, 1), Roots_vec_1(:,2))
hold on
plot3(a_var, Roots_vec_2(:, 1), Roots_vec_2(:,2))
xlabel('$\alpha$')
ylabel('c')
zlabel('$v_4$')
%legend('1', '2')
set(gca, "FontSize", 16)
%exportgraphics(gca,'Trace_Steady_States_Graph_fsolve_3D.png','Resolution',300)

function F=stst_trace(x,a)
c=x(1);
v4=x(2);
K_act=210*(10^(-9));
H_IP3=4.0;
K_infty=52*(10^(-6));
K_flux=4.9*(10^(-6));
V_e=10^(-6);
K_e=0.1*(10^(-6));
g=0.51;
betaosc=0.08; %this is tp in the previous model (2)
k_const=1.1;

K1=K_flux/V_e;
K2=(K_act*g)/V_e;
K3=(g*K_e)/V_e;

p = ((v4.*betaosc)./betaosc).*((c + (1 - a).*k_const)./(c + k_const));
K_inh = K_infty.*(p.^(H_IP3)./(p.^(H_IP3)+1));
K4 = V_e./(g.*K_inh);

h = (1./(1 + (K4.*c).^4));
F1 = (K1.*h.*((c.^2)./(K2.^2 + c.^2)) - ((c.^2)./(K3.^2 + c.^2))); %This is the differential equation for c w.r.t. t.

TrJ1 = ((2.*K1.*c)./(K2.^2 + c.^2));
TrJ2 = ((K2.^2)./(K2.^2 + c.^2));
TrJ3 = (1./(1 + (K4.*c).^4));
TrJ4 = ((2.*c)./(K3.^2 + c.^2));
TrJ5 = ((K3.^2)./(K3.^2 + c.^2));
TrJ6 = 1;

TrJ = TrJ1.*TrJ2.*TrJ3 - TrJ4.*TrJ5 - TrJ6;

F = [F1;TrJ];

end