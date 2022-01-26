% to use this scrpit first run main.m

%% B_red over T_red
figure(1);
plot(T_red,Bred_over_Tred);
title(['B^*_2 ueber T^*'])
xlabel('$T^* \ /[-]$','Interpreter','Latex');
ylabel('$B_v^* \ /[-]$','Interpreter','Latex');
%% Influenz of sigma
figure(2)
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

    % write as txt to print in Latex via pgfplots
    table = [T' B'];
    %output file
    fid = fopen('B_over_T.txt','wt'); 
    for ii = 1:size(table,1)
        fprintf(fid,'%g\t',table(ii,:));
        fprintf(fid,'\n');
    end
    fclose(fid);

hold on
% bigger b_0
% b_0=2/3 * pi * N_A * sigma^3
b_0_hb=b_0+20; % [cm^3/mol]
B_hb=(Bred_over_Tred * b_0_hb) ;
plot(T,B_hb);

    % write as txt to print in Latex via pgfplots
    table = [T' B_hb'];
    %output file
    fid = fopen('B_over_T_hb.txt','wt'); 
    for ii = 1:size(table,1)
        fprintf(fid,'%g\t',table(ii,:));
        fprintf(fid,'\n');
    end
    fclose(fid);

% smaller b_0
epsilon_over_kB=epsilon_over_kB-20; % [K]
% b_0=2/3 * pi * N_A * sigma^3
b_0_lb=b_0-20; % [cm^3/mol]
B_lb=(Bred_over_Tred * b_0_lb) ;
plot(T,B_lb);
hold off

    % write as txt to print in Latex via pgfplots
    table = [T' B_lb'];
    %output file
    fid = fopen('B_over_T_lb.txt','wt'); 
    for ii = 1:size(table,1)
        fprintf(fid,'%g\t',table(ii,:));
        fprintf(fid,'\n');
    end
    fclose(fid);

legend('\sigma_{N_2}','\sigma > \sigma_{N_2}','\sigma < \sigma_{N_2}')

%% Influenz epsilion
figure(3)
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
% bigger Epsilon 
epsilon_over_kB_he=epsilon_over_kB+50; % [K]
T_he=T_red*epsilon_over_kB_he;
plot(T_he,B);

    
    % write as txt to print in Latex via pgfplots
    table = [T_he' B'];
    %output file
    fid = fopen('B_over_T_he.txt','wt'); 
    for ii = 1:size(table,1)
        fprintf(fid,'%g\t',table(ii,:));
        fprintf(fid,'\n');
    end
    fclose(fid);

% smaller Epsilon 
epsilon_over_kB_le=epsilon_over_kB-50; % [K]
T_le=T_red*epsilon_over_kB_le;
plot(T_le,B);

hold off
 
    % write as txt to print in Latex via pgfplots
    table = [T_le' B'];
    %output file
    fid = fopen('B_over_T_le.txt','wt'); 
    for ii = 1:size(table,1)
        fprintf(fid,'%g\t',table(ii,:));
        fprintf(fid,'\n');
    end
    fclose(fid);

legend('\epsilon_{N_2}','\epsilon > \epsilon_{N_2}','\epsilon < \epsilon_{N_2}')

