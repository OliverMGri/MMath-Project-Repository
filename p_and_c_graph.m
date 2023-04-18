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

K1=K_flux/V_e;
K2=(K_act*g)/V_e;
K3=(g*K_e)/V_e;

c=0:0.001:0.5;

K4= (((K1./(c.^4)).*(((K3.^2)+(c.^2))./((K2.^2)+(c.^2)))-(1./(c.^4))).^(1/4)); %K4 is constructed as a function of c.
p= ((V_e./((g.*K_infty.*K4)-V_e)).^(1/4)); %p is constructed as a function of K4 (and therefore it is a function of c).

plot(p, c) %Simply plotting p and c together yields this graph.
xlabel('\(p\)')
ylabel('\(c\)')
xlim([0 0.35]) %The bifurcations that we are interested in must exist within this boundary. The simulations prove this, as there are no constant oscillations if p = 0.1 or 0.35, but there are in-between these values.
set(gca, "FontSize", 16)
exportgraphics(gca,'p_and_c.png','Resolution',300)