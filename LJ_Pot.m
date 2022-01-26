r_red=0.965:0.01:3;
LJpot=4*1*((1./r_red).^(12)-(1./r_red).^(6));
plot(r_red,LJpot)
axis([0 3 -1 1])

% % write as txt to print in Latex via pgfplots
% table = [r_red' LJpot'];
% %output file
% fid = fopen('LJ_pot.txt','wt'); 
% for ii = 1:size(table,1)
%     fprintf(fid,'%g\t',table(ii,:));
%     fprintf(fid,'\n');
% end
% fclose(fid);

LJforce=24[]