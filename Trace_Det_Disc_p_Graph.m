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

K_inh=@(p) (K_infty*(p^(H_IP3)/(p^(H_IP3)+1)));
h_s=@(c_s, p) (1/(1+((V_e/(g*K_inh(p)))*c_s)^4));
fn=@(c_s, p) (K1*h_s(c_s, p)*((c_s)^2/(K2^2+(c_s)^2))-((c_s)^2/(K3^2+(c_s)^2)));

K4=@(c) (((((K1 - 1)*(c^2)) - (K2^2) + (K1*(K3^2)))/((c^6) + ((K2^2)*(c^4))))^(1/4)); %K4 in terms of c.
p=@(c) ((V_e/((g*K_infty*K4(c))-V_e))^(1/4)); %p in terms of k
Found_c = [];
variable_p = [];
j=0.01;
i=1;
while j <= 0.5
    new_p=@(c) p(c) - j;
    Found_c(i) = fsolve(new_p, 0); %An array of c values that are found by solving p = j, where j is a constant between 0 and 0.5.
    variable_p(i) = j;
    j=j+0.005;
    i=i+1;
end

Tr=@(c_s, p) ((2*K1*c_s)/((K2^2)+(c_s^2))*(1/(1+(((V_e/(g*K_inh(p)))*c_s)^4)))*((K2^2)/((K2^2)+(c_s^2)))-((2*c_s)/((K3^2)+(c_s^2)))*((K3^2)/((K3^2)+(c_s^2)))-1);
Det=@(c_s, p) ((2*K1*c_s)/((K2^2)+(c_s^2))*(1/(1+(((V_e/(g*K_inh(p)))*c_s)^4)))*(((c_s^2/(K2^2+c_s^2))-1))+((2*c_s)/((K3^2)+(c_s^2)))*(1-((c_s^2)/((K3^2)+(c_s^2)))) + ((K1*(c_s^2))/(K2^2 + c_s^2))*(1/((1+((V_e/(g*K_inh(p)))*c_s)^4)^2))*(4*((V_e/(g*K_inh(p)))^4)*(c_s^3)));
Disc=@(c_s, p) ((Tr(c_s, p)^2)-(4*Det(c_s, p))); %The trace, determinant and discriminant are formulated as functions of c and p.


for k=1:size(Found_c, 2)
    Tr_values(k) = Tr(Found_c(k), variable_p(k));
    Det_values(k) = Det(Found_c(k), variable_p(k));
    Disc_values(k) = Disc(Found_c(k), variable_p(k)); %The trace, determinant and discriminant are found for the list of c and p values found earlier.
    k=k+1;
end

figure('units','normalized','outerposition',[0 0 1 1]) %A plot of the trace, determinant and discriminant is created. The points where they cross the x-axis are the bifurcation points.
plot(variable_p, Tr_values, 'g', LineStyle='-', LineWidth=2);
hold on;
plot(variable_p, Det_values, 'b', LineStyle='-', LineWidth=2);
hold on;
plot(variable_p, Disc_values, 'r', LineStyle='-', LineWidth=2);
hold on;
plot(variable_p,zeros(size(Found_c,2)), 'k', LineStyle='-', LineWidth=2);
xlabel('$p$', 'Interpreter','latex');
xticks([0 0.1 0.2 0.3 0.4 0.5]);
yticks([-60 -50 -40 -30 -20 -10 0 10 20])
xlim([0 0.5])
%ylim([-10 10])
l=legend('Trace', 'Determinant', 'Discriminant', 'Interpreter','latex');
l.LineWidth=1.5;
l.Location='northeastoutside';
set(gca,'fontsize',16)
set(gca,'linewidth',1.5)
grid on
exportgraphics(gca,'Trace_Det_Disc_p_v2_zoom.png','Resolution',300)