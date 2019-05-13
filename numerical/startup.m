% startup.m
% FILE B
%
% This script contains default settings, and should be executed at start of
% each routine producing graphical output

% Text
% ATTN LATEX INTERPRETER CAUSES CRASHES
set(groot,'DefaultTextInterpreter','Latex');
set(groot,'DefaultTextFontSize',10);
set(groot,'DefaultTextVerticalAlignment','Top');
set(groot,'DefaultTextHorizontalAlignment','Right');

% Axes Fonts
set(groot,'DefaultAxesFontSize',10);
set(groot,'DefaultAxesFontName','Times');
set(groot,'DefaultAxesTickLabelInterpreter','latex');
   
% Lines
set(groot,'DefaultLineLineWidth',1);
set(groot,'DefaultLineMarkerSize',4);
set(groot,'DefaultAxesLinewidth',1);
set(groot,'DefaultAxesColorOrder',[0 0 0]);

% Page size
set(groot,'DefaultFigurePosition',[10 10 10 10]);
set(groot,'DefaultFigurePaperUnits','centimeters');
set(groot,'DefaultFigureUnits','centimeters');
