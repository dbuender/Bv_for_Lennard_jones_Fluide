%%
for i=1:length(T_red)
y=My_over_x(:,i);
temp=300;
hold on
plot(x(1:temp),y(1:temp))

end

%% B_red over T_red
figure(2);
plot(T_red,Bred_over_Tred);
title(['B^*_2 ueber T^*'])
xlabel('$T^* \ /[-]$','Interpreter','Latex');
ylabel('$B_v^* \ /[-]$','Interpreter','Latex');
%% Einfluss sigma
figure();
N_A=6.02214086*10^23;
% LJ Parameter für N_2
epsilon_over_kB=95.05; % [K]
% b_0=2/3 * pi * N_A * sigma^3
b_0=63.78; % [cm^3/mol]
B=(Bred_over_Tred * b_0) ;
T=T_red*epsilon_over_kB;
plot(T,B);
title(['B über T'])
xlabel('$T \ /[K]$','Interpreter','Latex');
ylabel('$B_2 \ /[cm^3/mol]$','Interpreter','Latex');

hold on
% größeres b_0
% b_0=2/3 * pi * N_A * sigma^3
b_0_hb=b_0+20; % [cm^3/mol]
B_hb=(Bred_over_Tred * b_0_hb) ;
plot(T,B_hb);

% kleiners b_0
epsilon_over_kB=epsilon_over_kB-20; % [K]
% b_0=2/3 * pi * N_A * sigma^3
b_0_lb=b_0-20; % [cm^3/mol]
B_lb=(Bred_over_Tred * b_0_lb) ;
plot(T,B_lb);
hold off

legend('\sigma_{N_2}','\sigma > \sigma_{N_2}','\sigma < \sigma_{N_2}')
%% Einfluss epsilion
figure();
N_A=6.02214086*10^23;
% LJ Parameter für N_2
epsilon_over_kB=95.05; % [K]
% b_0=2/3 * pi * N_A * sigma^3
b_0=63.78; % [cm^3/mol]
B=(Bred_over_Tred * b_0) ;
T=T_red*epsilon_over_kB;
plot(T,B);
title(['B über T'])
xlabel('$T \ /[K]$','Interpreter','Latex');
ylabel('$B_2 \ /[cm^3/mol]$','Interpreter','Latex');

hold on
% höheres Epsilon 
epsilon_over_kB_he=epsilon_over_kB+50; % [K]
T_he=T_red*epsilon_over_kB_he;
plot(T_he,B);

% kleineres Epsilon 
epsilon_over_kB_le=epsilon_over_kB-50; % [K]
T_le=T_red*epsilon_over_kB_le;
plot(T_le,B);

legend('\epsilon_{N_2}','\epsilon > \epsilon_{N_2}','\epsilon < \epsilon_{N_2}')
