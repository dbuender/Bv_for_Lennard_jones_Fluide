
%==========================================================================
%            Computaion of specific molar entalphie
%==========================================================================

% The reason for using only every second element, begining with the second
% one, is we use the other to compute the derivative
h_real_red=p_red(2:2:end-1).'*(Bred_over_Tred(2:2:end-1)-T_red(2:2:end-1).*diff(Bred_over_Tred(1:2:end))./diff(T_red(1:2:end)));

figure()
h=surface(T_red(2:2:end-1),p_red(2:2:end-1),h_real_red);
set(h,'LineStyle','none')

%==========================================================================
%            Computaion of specific molar gibsenergie
%==========================================================================

g_real_red=p_red.'*Bred_over_Tred;

figure()
g=surface(T_red,p_red,g_real_red);
set(g,'LineStyle','none')

%==========================================================================
%            Computaion of specific molar Entrophi
%==========================================================================


s_real_red=-p_red(2:2:end-1).'*(diff(Bred_over_Tred(1:2:end))./diff(T_red(1:2:end)));


figure()
s=surface(T_red(2:2:end-1),p_red(2:2:end-1),s_real_red);
set(s,'LineStyle','none')

