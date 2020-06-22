
% Variation B
% This script requiers the main skrkipt to run first, as it funktion as
% an extention.

% Under the asumtion that B_v=B_p*RT we can compute the second viral
% koeffiezent of the pressure explicit viral equation.
% z = 1 + B_red_p * p_red = 1 + B_red_v * 2pi /( T_red * 3) * p_red 
Bred_over_Tred_p=Bred_over_Tred.*2*pi./(T_red.*3);


% compute z in a Matrix with T_red constant in each row and p_red in each
% colum
z_over_pred_VB=zeros(length(T_red),length(p_red));
for i=1:10:length(T_red)
    z_over_pred_VB(i,:) = 1 + Bred_over_Tred_p(i) .* p_red;
end


%% z over p_red
figure(5);
for i=1:10:length(T_red)
    hold on
    plot(p_red,z_over_pred_VB(i,:));
end
title(['Realgasfaktor (T^* = ' num2str(T_red_start) ' ... '  num2str(T_red_end) ')'])
xlabel('$p^* \ /[-]$','Interpreter','Latex');
ylabel('$z \ /[-]$','Interpreter','Latex');
hold off
axis([0 2 0.5 1.2]);