clear all
clc
close

list_factory = fieldnames(get(groot,'factory'));
index_interpreter = find(contains(list_factory,'Interpreter'));
for i = 1:length(index_interpreter)
    default_name = strrep(list_factory{index_interpreter(i)},'factory','default');
    set(groot, default_name,'latex');
end %LaTeX form.

a_range = linspace(0,1,301);
c_range = linspace(0,1);
[a_vector, c_vector] = meshgrid(a_range, c_range); %Create a matrix of values for alpha and c.

K_act=210*(10^(-9));
H_IP3=4.0;
K_infty=52*(10^(-6));
K_flux=4.9*(10^(-6));
V_e=10^(-6);
K_e=0.1*(10^(-6));
g=0.51;
betaosc=0.08; %this is tp in the previous model (2)
k_const=1.1;
v4 = 0.5*betaosc;

p = (v4./betaosc).*((c_vector + (1 - a_vector).*k_const)./(c_vector + k_const));
K_inh= K_infty.*(p.^(H_IP3)./(p.^(H_IP3)+1));

K1=K_flux/V_e;
K2=(K_act*g)/V_e;
K3=(g*K_e)/V_e;
K4= V_e./(g.*K_inh);

h = (1./(1 + (K4.*c_vector).^4));
F = (K1.*h.*((c_vector.^2)./(K2.^2 + c_vector.^2)) - ((c_vector.^2)./(K3.^2 + c_vector.^2))); %This is the differential equation for c w.r.t. t.

pcolor(a_vector, c_vector, double(F>0))
colorbar
shading interp
xlabel('$\alpha$')
ylabel('c')

TrJ1 = ((2.*K1.*c_vector)./(K2.^2 + c_vector.^2));
TrJ2 = ((K2.^2)./(K2.^2 + c_vector.^2));
TrJ3 = (1./(1 + (K4.*c_vector).^4));
TrJ4 = ((2.*c_vector)./(K3.^2 + c_vector.^2));
TrJ5 = ((K3.^2)./(K3.^2 + c_vector.^2));
TrJ6 = 1;

TrJ = TrJ1.*TrJ2.*TrJ3 - TrJ4.*TrJ5 - TrJ6;

hold on
pcolor(a_vector, c_vector, double(F>0) + 2*double(TrJ>0))
colorbar
shading interp
xlabel('$\alpha$')
ylabel('c')
set(gca, "FontSize", 16)
exportgraphics(gca,'Trace_Steady_States_Graph_c_alpha.png','Resolution',300)