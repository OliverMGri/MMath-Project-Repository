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
%v4 = 0.5*betaosc;
a = 0.6;

K1=K_flux/V_e;
K2=(K_act*g)/V_e;
K3=(g*K_e)/V_e;

cs_1_F = 0.01;
cs_2_F = 3;
cs_1_T = 0.01;
cs_2_T = 0.5;

c = 0:0.001:0.5; %Setting the acceptible range that we can expect values for c when testing to find the trace.

veci = 1;

%a_var = linspace(0,1,20);
v4_var = linspace(0.8,0.4,100);

for v4 = v4_var
%for  = a_var
p =@(c) ((v4.*betaosc)./betaosc).*((c + (1 - a).*k_const)./(c + k_const));
K_inh =@(c) K_infty.*(p(c).^(H_IP3)./(p(c).^(H_IP3)+1));
K4 =@(c) V_e./(g.*K_inh(c));

h =@(c) (1./(1 + (K4(c).*c).^4));
F =@(c) (K1.*h(c).*((c.^2)./(K2.^2 + c.^2)) - ((c.^2)./(K3.^2 + c.^2))); %This is the differential equation for c w.r.t. t.

TrJ1 =@(c) ((2.*K1.*c)./(K2.^2 + c.^2));
TrJ2 =@(c) ((K2.^2)./(K2.^2 + c.^2));
TrJ3 =@(c) (1./(1 + (K4(c).*c).^4));
TrJ4 =@(c) ((2.*c)./(K3.^2 + c.^2));
TrJ5 =@(c) ((K3.^2)./(K3.^2 + c.^2));
TrJ6 =@(c) 1;

TrJ =@(c) TrJ1(c).*TrJ2(c).*TrJ3(c) - TrJ4(c).*TrJ5(c) - TrJ6(c);

options = optimset('Display','iter');

[cs_1_F,fval,exitflag,output] = fzero(F, cs_1_F, options); %this should find c when the trace equals zero, which helps us find bifurcation points.
[cs_2_F,fval,exitflag,output] = fzero(F, cs_2_F, options);
[cs_1_T,fval,exitflag,output] = fzero(TrJ, cs_1_T, options); %this should find c when the trace equals zero, which helps us find bifurcation points.
[cs_2_T,fval,exitflag,output] = fzero(TrJ, cs_2_T, options);

cs_1_Fvec(veci) = cs_1_F;
cs_2_Fvec(veci) = cs_2_F;
cs_1_Tvec(veci) = cs_1_T;
cs_2_Tvec(veci) = cs_2_T;
E_F(veci) = exitflag;
veci = veci + 1;
end
%end

cs_1_F
cs_2_F
cs_1_Fvec
cs_2_Fvec

cs_1_T
cs_2_T
cs_1_Tvec
cs_2_Tvec

%plot(v4_var,cs_1_Fvec) %I wanted to plot the values of alpha against one set of bifurcation points.
%hold on
plot(v4_var,cs_2_Fvec) %Choose better initial conditions when between around a=0.1 and 0.35. c=1 or something similar.
hold on
plot(v4_var,cs_1_Tvec) %I wanted to plot the values of alpha against one set of bifurcation points.
hold on
plot(v4_var,cs_2_Tvec)
hold on
yline(0,'-k')
xlabel('$v_4$')
ylabel('$c$')
ylim([-0.1 0.4])
legend({'F', 'Trace 1','Trace 2'}, 'Location', 'northwest')
set(gca, "FontSize", 16)
exportgraphics(gca,'Trace_Steady_States_Graph_c_v4_a06_zoom.png','Resolution',300)
