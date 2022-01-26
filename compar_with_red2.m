% Data for N_2
T_kritt = 126.1; % [K]
p_kritt = 3.394*10^6; % [Pa] Quelle https://www.chemie.de/lexikon/Kritischer_Punkt_%28Thermodynamik%29.html
Rm=8.31446261815324; % [J/(mol K)]

% LJ Parameter für N_2
epsilon_over_kB=95.05; % [K]
% b_0=2/3 * pi * N_A * sigma^3
b_0=63.78*10^-6; % [m^3/mol]

%% z over p_red
figure(3);

% transformation of red to red2 
% red stands for reduzed by LJ Parameters and red2 for reduzed by the
% critikal Point
T_red2=epsilon_over_kB*T_red/T_kritt;
p_red2=2/3*(pi()*Rm*epsilon_over_kB)/(b_0*p_kritt)*p_red;
% The Temparatures at wich z over p_red will plot
T_zplot=[1 1.1 1.2 1.3 1.4 1.5 2];
for i=1:length(T_zplot)
    hold on
    %searching for the nearsest datepoint to tehe requestet Temperature
    l=abs(T_red2-T_zplot(i))==min(abs(T_red2-T_zplot(i)));
    if abs(T_zplot(i)-T_red2(l))>0.1
        warning('nearst Datapoint was quit far away the Graph might be wrong')
    end
    
    plot(p_red2(z_over_pred(l,:)>0.5),z_over_pred(l,z_over_pred(l,:)>0.5),'-');
    
    % write as txt to print in Latex via pgfplots
    table = [p_red2(z_over_pred(l,:)>0.5)' z_over_pred(l,z_over_pred(l,:)>0.5)'];
    %output file
    filename=['z_over_p_red2_at_T_red2_' num2str(T_zplot(i)) '.txt'];
    fid = fopen(filename,'wt'); 
    for ii = 1:size(table,1)
        fprintf(fid,'%g\t',table(ii,:));
        fprintf(fid,'\n');
    end
    fclose(fid);
end
title(['Realgasfaktor (T^* = ' num2str(T_zplot(1)) ' ... '  num2str(T_zplot(length(T_zplot))) ')'])
xlabel('$p_r \ /[-]$','Interpreter','Latex');
ylabel('$z \ /[-]$','Interpreter','Latex');
legend('T_r = 1','T_r = 1.1','T_r = 1.2','T_r = 1.3','T_r = 1.4','T_r = 1.5','T_r = 2')
hold off
