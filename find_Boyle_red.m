function [T_Boyle_red,error]=find_Boyle_red()
% this function will Locate the Boyle Temperature, by using iterative
% approach

error=0;

% defin the reduced distanc x = r/sigma
x_diff=0.01;
x_start=x_diff+0;
num_of_kritt_x=1e4; % gives the number of loops, at witch we the programm will stop. Is to avoid an endlessloop if an error occurs
x=(x_start:x_diff:x_diff*num_of_kritt_x)';

% define the used variables for the mayer funktion
My_target_diff=1e-10;
My_over_x=zeros(num_of_kritt_x,1);

% define the used variables for the secondary viralcoefficent
Bred=-1;
Bred_laststep=-1;
Tdiff=1;
T_Boyle_red=1.5;

while ~(1e-10>abs(Bred))
    if Bred*Bred_laststep<0
        Tdiff=Tdiff/2;
    end
    Bred_laststep=Bred;
    
    if Bred<0
        T_Boyle_red=T_Boyle_red+Tdiff;
    else
        T_Boyle_red=T_Boyle_red-Tdiff;
    end
    
    
    % the continuous variable for the next loop
    n=1;
    bool_x_reached=0;
    % Reset the My_over_x
    My_over_x=zeros(num_of_kritt_x,1);
    
    My_over_x(n)=exp(-4/T_Boyle_red*(x(n)^-12-x(n)^-6))-1;
    n=2;
    % this Loop should repet till the Mayers funktion tends close enough to
    % zero
    while ~bool_x_reached
        
        My_over_x(n)=exp(-4/T_Boyle_red*(x(n)^-12-x(n)^-6))-1;
        
        %
        % wenn der Abstand klein genug ist und eine negative Steigung
        % vorliegt wird abgebrochen
        if My_over_x(n)<My_target_diff && My_over_x(n)-My_over_x(n-1)<0
            bool_x_reached=1;
            %            My_over_x(n+1,i)=-2; % the idea was to use it as cutoff to speed up intigration
        elseif n==num_of_kritt_x
            bool_x_reached=1;
            error=1;
            warning(['calculations of My was not under set Limit [n=' num2str(n) ',T_red=' num2str(T_Boyle_red)  ']' ])
        end
        
        n=n+1;
        
    end
    
    % calculation of the integral
    Bred=-3*trapz(x,My_over_x(:).*x.^2);
end
display(T_Boyle_red)
display(Bred)
end