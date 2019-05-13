% study1.m
% FILE 31
%
% In this study we run the elliptical FR setting for various combinations
% of the ellipse parameters and compute the Cp versus Vp curves. This is
% the first set of results appearing in the paper
%
% Initialization
clear;
startup;

% FIXED PARAMETERS
% Node parameters
lambda=1;  
v0=1;      
r0=1;      
node_pars=[lambda v0 r0];

% numerical parameters
% Each data point: around 12 minutes and max out on memory!
N=256;      
S=400;      
targetpoints=250; 

verbose=0;  

% environment parameters
potential_type=1;  
cost_type=1;  
density_type=1;  

% VARIABLE PARAMETERS
arange=0.2:0.2:1;
epsilonrange=[0:0.1:0.9 0.99];
I=length(arange);
J=length(epsilonrange);

% WE INITIALIZE THE OUTPUT
Vp_num_1=zeros(I,J);
Cp_num_1=zeros(I,J);
locations_1=zeros(I,J);

% WE RUN SCRIPTS
for i=1:I
    for j=1:J
      
        FR_type=[3 arange(i) epsilonrange(j)];  
        
        % We produce quick update
        st1=num2str(i);
        st2=num2str(length(arange));
        st3=num2str(j);
        st4=num2str(length(epsilonrange));
        disp(['Starting iteration: ' st1  '/' st2 ', ' st3 '/' st4 ' at time']);
        disp(datetime('now','Format','HH:mm:ss'));
        
        % We find B and L
        a=arange(i);
        epsilon=epsilonrange(j);
        beta=a*sqrt(1-epsilon^2);
        area=pi*a*beta;
        B=1.1*a*(1+epsilon);
        L=round(2*B*sqrt(targetpoints/area));
        num_pars=[N B L S];
        
        metricsout=metrics...
         (node_pars,FR_type,potential_type,cost_type,density_type,num_pars,verbose);

        % We produce output
        Vp_num_1(i,j)=metricsout.Vp;
        Cp_num_1(i,j)=metricsout.Cp;   
        locations_1(i,j)=length(metricsout.locations(:,1));
        disp(['Vp and Cp are: ' num2str([metricsout.Vp metricsout.Cp])]);
        disp(['Locations used: ' num2str(locations_1(i,j))]);
        
    end
end

save('saves/study1');




