error=0;

% defin the reduced distanc x = sigma/r 
x_diff=0.01;
x_start=x_diff+0;
num_of_kritt_x=1e4; % gives the number of loops, at witch we the programm will stop. Is to avoid an endlessloop if an error occurs 
x=(x_start:x_diff:x_diff*num_of_kritt_x)';

% defin the reduced Temprature
% T_red_diff=1;
% T_red_start=T_red_diff+0; % so T will not start at zero
% T_red_end=100;
% T_red=T_red_start:T_red_diff:T_red_end;

% defin the reduced Temprature on a logarithmic scale
T_red_start=1;
T_red_end=100;
n_T_red=100;
c=(log(T_red_end)-log(T_red_start))/n_T_red;
T_red=zeros(1,n_T_red);
T_red(1)=T_red_start;
for n=2:n_T_red
    T_red(n)=T_red(n-1)*exp(c);
end


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
       
       %
       % wenn der Abstand klein genug ist und eine negative Steigung
       % vorliegt wird abgebrochen
       if My_over_x(n,i)<My_target_diff && My_over_x(n,i)-My_over_x(n-1,i)<0 
           bool_x_reached=1;
%            My_over_x(n+1,i)=-2; % the idea was to use it as cutoff to speed up intigration 
       elseif n==num_of_kritt_x
           bool_x_reached=1;
           error=1;
           warning(['calculations of My was not under set Limet [n=' num2str(n) ',T_red=' num2str(T_red(i))  ']' ])
       end
       
       n=n+1;
       
   end
   
   % calculation of the integral
   Bred_over_Tred(i)=-3*trapz(x,My_over_x(:,i).*x.^2);
   
end
%

%==========================================================================
%                   Visualisation
%==========================================================================
%% Mayer funktions
figure(1);
for i=1:20:length(T_red)
y=My_over_x(:,i);
temp=300;
hold on
plot(x(1:temp),y(1:temp))
xlabel('$x \ /[-]$','Interpreter','Latex');
ylabel('$M_y \ /[-]$','Interpreter','Latex');
end
hold off

%% B_red over T_red
figure(2);
plot(T_red,Bred_over_Tred);
xlabel('$T^* \ /[-]$','Interpreter','Latex');
ylabel('$B_v^* \ /[-]$','Interpreter','Latex');

