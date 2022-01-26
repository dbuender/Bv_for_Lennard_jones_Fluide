error=0;

% defin the reduced distanc x = sigma/r 
x_diff=0.01;
x_start=x_diff+0;
num_of_kritt_x=1e4; % gives the number of loops, at witch we the programm will stop. Is to avoid an endlessloop if an error occurs 
x=(x_start:x_diff:x_diff*num_of_kritt_x)';

% % defin the reduced Temprature
% T_red_diff=1;
% T_red_start=T_red_diff+0; % so T will not start at zero
% T_red_end=100;
% T_red=T_red_start:T_red_diff:T_red_end;

% defin the reduced Temprature on a logarithmic scale
T_red_start=1;
T_red_end=50;
n_T_red=1000;
c=(log(T_red_end)-log(T_red_start))/n_T_red;
T_red=zeros(1,n_T_red);
T_red(1)=T_red_start;
for n=2:n_T_red
    T_red(n)=T_red(n-1)*exp(c);
end

% % defin the reduced Temprature Manuel
% [T_Boyle_red,error]=find_Boyle_red;
% T_red=[1 2 3 T_Boyle_red 4 5 6];

% define the used variables for the mayer funktion
My_target_diff=1e-10;
My_over_x=zeros(num_of_kritt_x,length(T_red));

% define the used variables for the secondary viralcoefficent
Bred_over_Tred=zeros(1,length(T_red));


%==========================================================================
%            Computaion of B_red by computing Mayers function
%==========================================================================
% Verwendete Gleichung:
% B_red = -3 * intgral von {My-1}*x^2 nach x
% Mayers function My = exp[-4/Tred*(x^-12+x^-6)]


% in the first loop the reduced Temprature will be variied
for i=1:length(T_red)
    
    % the continuous variable for the next loop
    n=1;
    bool_x_reached=0;
    
    My_over_x(n,i)=exp(-4/T_red(i)*(x(n)^-12-x(n)^-6))-1;
    n=2;
    % this Loop should repet till the Mayers funktion tends close enough to
    % zero
   while ~bool_x_reached
       
       My_over_x(n,i)=exp(-4/T_red(i)*(x(n)^-12-x(n)^-6))-1;
       
       % if the the Distanc to zero is small enough
       % wenn der Abstand klein genug ist und eine negative Steigung
       % vorliegt wird abgebrochen
       if My_over_x(n,i)<My_target_diff && My_over_x(n,i)-My_over_x(n-1,i)<0 
           bool_x_reached=1;
%            My_over_x(n+1,i)=-2; % the idea was to use it as cutoff to speed up intigration 
       elseif n==num_of_kritt_x
           bool_x_reached=1;
           error=1;
           warning(['calculations of My was not under set Limit [n=' num2str(n) ',T_red=' num2str(T_red(i))  ']' ])
       end
       
       n=n+1;
       
   end
   
   % calculation of the integral
   Bred_over_Tred(i)=-3*trapz(x,My_over_x(:,i).*x.^2);
   
end
%

%==========================================================================
%            Computaion of z
%==========================================================================


% Variation A
% z = 1/2 + sqrt(1/4 + 2/3 * pi * B_red * p_red / T_red)

% defin the reduced pressure
% p_red_start=0;
% p_red_end=2;
% n_p_red=100;
% p_red_diff=(p_red_end-p_red_start)/n_p_red;
% p_red=p_red_start:p_red_diff:p_red_end;

% define the intervall of the reduced pressure on logarithmic scale
p_red_start=1e-10;
p_red_end=1;
n_p_red=1000;
c=(log(p_red_end)-log(p_red_start))/n_p_red;
p_red=zeros(1,n_p_red);
p_red(1)=p_red_start;
for n=2:n_p_red
    p_red(n)=p_red(n-1)*exp(c);
end

% compute z in a Matrix with T_red constant in each row and p_red in each
% colum
z_over_pred=zeros(length(T_red),length(p_red));
for i=1:length(T_red)
    z_over_pred(i,:) = 1/2 + sqrt(1/4 + 2/3 * pi * Bred_over_Tred(i) .* p_red ./ T_red(i));
end

% as funktion of the Temperture
% z = 1/2 + sqrt(1/4 + 2 * pi * B_red / (3 * T)) p_red

z_over_Tred=zeros(length(T_red),length(p_red));
for i=1:length(T_red)
    z_over_Tred(i,:) = 1/2 + sqrt(1/4 + 2/3 * pi * Bred_over_Tred .* p_red(i) ./ T_red);
end

%==========================================================================
%            Computaion of p_red_kritt
%==========================================================================
% p_red_kritt is the reduced pressure at wich the computen of z over
% Varaiant A is no longer possible
% 
p_red_kritt=zeros(1,length(T_red));
p_red_kritt= -3/(8*pi)*T_red./Bred_over_Tred;

%==========================================================================
%            Computaion of rho_red_kritt
%==========================================================================
% p_red_kritt is the reduced pressure at wich the computen of z over
% Varaiant A is no longer possible
% 
rho_red_kritt=zeros(1,length(T_red));
rho_red_kritt= -3./(4*pi*Bred_over_Tred);

%==========================================================================
%                   Visualisation
%==========================================================================
%% Mayer funktion
figure(1);
for i=1:100:length(T_red)
    y=My_over_x(:,i);
    temp=300;
    hold on
    plot(x(1:temp),y(1:temp))
end
title(['Mayer-Funktion (T^* = ' num2str(T_red_start) ' ... '  num2str(T_red_end) ')'])
xlabel('$x \ /[-]$','Interpreter','Latex');
ylabel('$M_y \ /[-]$','Interpreter','Latex');
hold off

%% B_red over T_red
figure(2);
plot(T_red,Bred_over_Tred);
title(['B^*_2 ueber T^*'])
xlabel('$T^* \ /[-]$','Interpreter','Latex');
ylabel('$B_v^* \ /[-]$','Interpreter','Latex');

% write as txt to print in Latex via pgfplots
table = [T_red' Bred_over_Tred'];
%output file
fid = fopen('B_red_over_T_red.txt','wt'); 
for ii = 1:size(table,1)
    fprintf(fid,'%g\t',table(ii,:));
    fprintf(fid,'\n');
end
fclose(fid);

%% z over p_red
figure(3);
% The Temparatures at wich z over p_red will plot
T_zplot=[1 2 3 3.41 4 5 6];
for i=1:length(T_zplot)
    hold on
    %searching for the nearsest datepoint to tehe requestet Temperature
    l=abs(T_red-T_zplot(i))==min(abs(T_red-T_zplot(i)));
    if abs(T_zplot(i)-T_red(l))>0.1
        warning('nearst Datapoint was quit far away the Graph might be wrong')
    end
    
    plot(p_red(z_over_pred(l,:)>0.5),z_over_pred(l,z_over_pred(l,:)>0.5),'-');
    
    % write as txt to print in Latex via pgfplots
    table = [p_red(z_over_pred(l,:)>0.5)' z_over_pred(l,z_over_pred(l,:)>0.5)'];
    %output file
    filename=['z_over_p_red_at_T_red_' num2str(T_zplot(i)) '.txt'];
    fid = fopen(filename,'wt'); 
    for ii = 1:size(table,1)
        fprintf(fid,'%g\t',table(ii,:));
        fprintf(fid,'\n');
    end
    fclose(fid);
end
title(['Realgasfaktor (T^* = ' num2str(T_zplot(1)) ' ... '  num2str(T_zplot(length(T_zplot))) ')'])
xlabel('$p^* \ /[-]$','Interpreter','Latex');
ylabel('$z \ /[-]$','Interpreter','Latex');
legend('T^* = 1','T^* = 2','T^* = 3','T^*_{Boyle}','T^* = 4','T^* = 5','T^* = 6')
hold off
%% z over T_red
figure(4);
% The Temparatures at wich z over p_red will plot
p_zplot=[0.01 0.1 1];
for i=1:length(p_zplot)
    hold on
    %searching for the nearsest datepoint to tehe requestet pressure
    l=abs(p_red-p_zplot(i))==min(abs(p_red-p_zplot(i)));
    if abs(p_zplot(i)-p_red(l))>0.1
        warning('nearst Datapoint was quit far away the Graph might be wrong')
    end
    
    plot(T_red(z_over_Tred(l,:)>0.5),z_over_Tred(l,z_over_Tred(l,:)>0.5),'-');
    
end
title(['Realgasfaktor (p^* = ' num2str(p_zplot(1)) ' ... '  num2str(p_zplot(length(p_zplot))) ')'])
xlabel('$T^* \ /[-]$','Interpreter','Latex');
ylabel('$z \ /[-]$','Interpreter','Latex');
legend('p^* = 0.01','p^* = 0.1','p^* = 1')
hold off

%% p_red_kritt over T_red
figure(5);
plot(T_red,p_red_kritt);
set(gca, 'YScale', 'log')
title(['p^*_{kritt}'])
xlabel('$T^* \ /[-]$','Interpreter','Latex');
ylabel('$p^* \ /[-]$','Interpreter','Latex');

%% rho_red_kritt over T_red
figure(6);
plot(T_red,rho_red_kritt);
set(gca, 'YScale', 'log')
title(['rho^*_{kritt}'])
xlabel('$T^* \ /[-]$','Interpreter','Latex');
ylabel('$\rho^* \ /[-]$','Interpreter','Latex');

