clear all
clc
close
%addpath(LaTeX_Adapting_Code)

% This script changes all interpreters from tex to latex. 
list_factory = fieldnames(get(groot,'factory'));
index_interpreter = find(contains(list_factory,'Interpreter'));
for i = 1:length(index_interpreter)
    default_name = strrep(list_factory{index_interpreter(i)},'factory','default');
    set(groot, default_name,'latex');
end

K_act=210*(10^(-9));
H_IP3=4.0;
K_infty=52*(10^(-6));
K_flux=4.9*(10^(-6));
V_e=10^(-6);
K_e=0.1*(10^(-6));
g=0.51; %various variables in the differential model of c and h.

betaosc=0.08; %this is tp in the previous model (2)

v4 = 0.3*betaosc; %for simplicity's sake, I just wrote v4 as some multiple of betaosc.

K1=K_flux/V_e;
K2=(K_act*g)/V_e;
K3=(g*K_e)/V_e; %more parameters of the differential model of c and h.

cs_1 = 0.01;
cs_2 = 0.5;

c=0.1; %Setting the acceptible range that we can expect values for c when testing to find the trace.

vi = 1;
alpha_var = linspace(0,1,20);

k_const = 0.01:0.01:3;
a1 = 1 - ((-0.05 + 0.5*k_const)./(k_const));
%a2 = 1 - (betaosc./(v4.*k_const))*0.5.*(0.1 + k_const) + (0.1./k_const);

plot(k_const, a1)
hold on
%plot(k_const, a2)
xlabel('$k_{const}$')
ylabel('$\alpha$')
ylim([0 1])
set(gca, "FontSize", 16)
exportgraphics(gca,'alpha_kconst_ver3.png','Resolution',300)