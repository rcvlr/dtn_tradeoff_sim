% function thetaaxes(xaxis,yaxis)
% FILE G
%
% This specialized function formats the x and y axes in case the axes are
% angles in the [-pi,pi] range. 
%
% Input Parameters
% xaxis: if xaxis=1, the x-axis is formatted. 
% yaxis: if yaxis=1, the y-axis is formatted. 

function thetaaxes(xaxis,yaxis)

thetaticks=[-pi -3*pi/4 -pi/2 -pi/4 0 pi/4 pi/2 3*pi/4 pi];
thetaticklabels={'$-\pi$' '$-\frac{3\pi}{4}$' '$-\frac{\pi}{2}$'...
    '$-\frac{\pi}{4}$' '$0$' '$\frac{\pi}{4}$' '$\frac{\pi}{2}$'...
    '$\frac{3\pi}{4}$' '$\pi$'};

if xaxis==1
   set(gca,'XTick',thetaticks); set(gca,'XTickLabel',thetaticklabels);
   set(gca,'XLim',[-pi pi]);
end

if yaxis==1
   set(gca,'YTick',thetaticks); set(gca,'YTickLabel',thetaticklabels);
   set(gca,'YLim',[-pi pi]);
end
    