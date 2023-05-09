clear 
clc
close all

list_factory = fieldnames(get(groot,'factory'));
index_interpreter = find(contains(list_factory,'Interpreter'));
for i = 1:length(index_interpreter)
    default_name = strrep(list_factory{index_interpreter(i)},'factory','default');
    set(groot, default_name,'latex');
end

%I would like to test the value of c at numerous values of p.
time_intervals = [0, 100, 200, 300, 400, 500, 600];
values_of_p = [0.15, 0.2, 0.225, 0.25, 0.275, 0.35];  %There are 6 values, as there are 6 time domains.
interval = length(values_of_p);
x0 = [1 1]; %The initial conditions of the ODE.
t = cell(interval,1);
CN = cell(interval, 2);


%This will provide an iterative solution, testing p's values along the timeframe.
for idx = 1 : interval
    p = values_of_p(idx);
    [t{idx},CN{idx}] = ode45(@(t,x) Model(t, x, p), time_intervals(idx:idx+1), x0);
    x0 = CN{idx}(end,:);
  
end

c1 = [CN(:,1)];
c2 = cat(1, c1{:});
c = c2(:,1);
h = c2(:,2);
T = cat(1, t{:});

quantity = zeros(1,length(t));
for i = 1:length(t)
    quantity(i) = size(t{i},1);
end


%finally, we have an array of IP_3 values.
x=zeros(1,length(T));
x(1:quantity(1))= values_of_p(1);
x(quantity(1)+1:quantity(1)+quantity(2))=values_of_p(2);
x(quantity(1)+quantity(2)+1:quantity(1)+quantity(2)+quantity(3))=values_of_p(3);
x(quantity(1)+quantity(2)+quantity(3)+1:quantity(1)+quantity(2)+quantity(3)+quantity(4))=values_of_p(4);
x(quantity(1)+quantity(2)+quantity(3)+quantity(4)+1:quantity(1)+quantity(2)+quantity(3)+quantity(4)+quantity(5))=values_of_p(5);
x(quantity(1)+quantity(2)+quantity(3)+quantity(4)+quantity(5)+1:quantity(1)+quantity(2)+quantity(3)+quantity(4)+quantity(5)+quantity(6))=values_of_p(6);
p1=transpose(x);

%Preparing to plot the solution.
plot(T,c, 'LineWidth',2);
hold on;
xlabel('$t (s)$','interpreter','latex');
plot(T,p1,'-', 'LineWidth',2);
l=legend('$c$ : Ca$^{2+}$ concentration','$p$ : IP$_3$ concentration');
set(l, 'interpreter', 'latex')
set(gca, "FontSize", 14)
exportgraphics(gca,'Step_Function_v2.png','Resolution',300);

%Now we must incorporate the model itself.
function M = Model(t,cn,p)

K_act=210*(10^(-9));
H_IP3=4.0;
K_infty=52*(10^(-6));
K_flux=4.9*(10^(-6));
V_e=10^(-6);
K_e=0.1*(10^(-6));
g=0.51;

K_inh=K_infty*(p^(H_IP3)/(p^(H_IP3)+1));

K1=K_flux/V_e;
K2=(K_act*g)/V_e;
K3=(g*K_e)/V_e;
K4=V_e/(g*K_inh);

Fcn=K1*(cn(2))*((cn(1))^2/(K2^2+(cn(1))^2))-((cn(1))^2/(K3^2+(cn(1))^2));
Gcn=(1/(1+(K4*(cn(1)))^4))-(cn(2));

M=[Fcn; Gcn];
end