% function newfig(fignum,position)
% FILE D
%
% This function created a new figure, making sure some parameters are set
% correctly and that the figure is completely initialized 
%
% Input Parameters
% ================
% fignum: the number of the figure
% position=[left bottom width height] the location of the figure, in cm.
%       Also used to specify the size of the figure on the paper. 


function newfig(fignum,position)

figure(fignum); delete(fignum); figure(fignum); 
set(gcf,'PaperUnits','centimeters');
set(gcf,'Units','centimeters');
set(gcf,'Position',position);
set(gcf,'PaperSize',[position(3) position(4)]);
set(gcf,'PaperPosition',[0 0 position(3) position(4)]);


