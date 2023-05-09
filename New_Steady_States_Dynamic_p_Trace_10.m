clear all
clc
close

list_factory = fieldnames(get(groot,'factory'));
index_interpreter = find(contains(list_factory,'Interpreter'));
for i = 1:length(index_interpreter)
    default_name = strrep(list_factory{index_interpreter(i)},'factory','default');
    set(groot, default_name,'latex');
end %LaTeX form.

K_act=210*(10^(-9));
H_IP3=4.0;
K_infty=52*(10^(-6));
K_flux=4.9*(10^(-6));
V_e=10^(-6);
K_e=0.1*(10^(-6));
g=0.51;
betaosc=0.08; %this is tp in the previous model (2)
k_const=1.1;
v4 = 0.4*betaosc;
%a = 0.1;

%p =@(c) (v4./betaosc).*((c + (1 - a).*k_const)./(c + k_const));
%K_inh=@(c) K_infty.*(p(c).^(H_IP3)./(p(c).^(H_IP3)+1));

K1=K_flux/V_e;
K2=(K_act*g)/V_e;
K3=(g*K_e)/V_e;
%K4=@(c) V_e./(g.*K_inh(c));

%h =@(c) (1./(1 + (K4(c).*c).^4)); %Derivative set to zero.
%F1 =@(c) (K1.*h(c).*((c.^2)./(K2.^2 + c.^2)) - ((c.^2)./(K3.^2 + c.^2))); %Derivative w.r.t. c.

%TrJ1 =@(c) ((2.*K1.*c)./(K2.^2 + c.^2));
%TrJ2 =@(c) ((K2.^2)./(K2.^2 + c.^2));
%TrJ3 =@(c) (1./(1 + (K4(c).*c).^4));
%TrJ4 =@(c) ((2.*c)./(K3.^2 + c.^2));
%TrJ5 =@(c) ((K3.^2)./(K3.^2 + c.^2));
%TrJ6 =@(c) ones(1);

%F2 =@(c) TrJ1(c).*TrJ2(c).*TrJ3(c) - TrJ4(c).*TrJ5(c) - TrJ6(c);

veci = 1;

a_var = linspace(0,1,50);
%v4_var = linspace(0,0.5,20);

%for v4 = v4_var
cs_1_F=0.5;
cs_2_F=0.01;
for a = a_var
F = @(c) [(K1.*(1./(1 + ((V_e./(g.*(K_infty.*(((v4./betaosc).*((c + (1 - a).*k_const)./(c + k_const))).^(H_IP3)./(((v4./betaosc).*((c + (1 - a).*k_const)./(c + k_const))).^(H_IP3)+1))))).*c).^4)).*((c.^2)./(K2.^2 + c.^2)) - ((c.^2)./(K3.^2 + c.^2)));
    ((2.*K1.*c)./(K2.^2 + c.^2)).*((K2.^2)./(K2.^2 + c.^2)).*(1./(1 + ((V_e./(g.*(K_infty.*(((v4./betaosc).*((c + (1 - a).*k_const)./(c + k_const))).^(H_IP3)./(((v4./betaosc).*((c + (1 - a).*k_const)./(c + k_const))).^(H_IP3)+1))))).*c).^4)) - ((2.*c)./(K3.^2 + c.^2)).*((K3.^2)./(K3.^2 + c.^2)) - ones(1)];

cs_1_F = fsolve(F, cs_1_F);
cs_2_F = fsolve(F, cs_2_F);
cs_1_Fvec(veci) = cs_1_F;
cs_2_Fvec(veci) = cs_2_F;
veci = veci + 1;
end

cs_1_F
cs_2_F
cs_1_Fvec
cs_2_Fvec

plot(a_var, cs_1_Fvec, 'LineWidth', 2)
hold on
%plot(a_var, cs_2_Fvec)
yline(0,'-k')
xlabel('$\alpha$')
ylabel('$c$')
xticks(0:0.2:1)
yticks(-0.2:0.2:1)
xlim([0 1])
ylim([-0.2 1])
set(gca, "FontSize", 24)
exportgraphics(gca,'Trace_Steady_States_Graph_both0_c_a_v404.png','Resolution',300)