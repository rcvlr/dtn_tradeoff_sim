% study3.m
% FILE 33
%
% In this study we run the disk FR setting for various combinations
% of the disk radius and Theta_w

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
% Each data point: around 17 minutes and max out on memory!
N=256;
S=400;
targetpoints=250; 

verbose=0;

% environment parameters
potential_type=1;  
cost_type=1; 

% VARIABLE PARAMETERS
krange=4:4:32;
arange=[0.01 0.5:0.5:4];          
I=length(krange);
J=length(arange);

% WE INITIALIZE OUTPUT
Vp_num_3=zeros(I,J);
Cp_num_3=zeros(I,J);
locations_3=zeros(I,J);

% WE RUN SCRIPTS
for i=1:I

    density_type=[2 N krange(i)];  
        
    for j=1:J
      
        FR_type=[1 arange(j)];  
        
        % We produce quick update
        st1=num2str(i);
        st2=num2str(length(krange));
        st3=num2str(j);
        st4=num2str(length(arange));
        disp(['Starting iteration: ' st1  '/' st2 ', ' st3 '/' st4 ' at time']);
        disp(datetime('now','Format','HH:mm:ss'));
        
        % We find B and L
        B=1.1*arange(j);
        L=round(sqrt(4*1.1^2*targetpoints/pi));
        num_pars=[N B L S];
        
        metricsout=metrics(node_pars,FR_type,potential_type,cost_type,density_type,num_pars,verbose);
       
        % We produce output
        Vp_num_3(i,j)=metricsout.Vp;
        Cp_num_3(i,j)=metricsout.Cp;        
        locations_3(i,j)=length(metricsout.locations(:,1));
        disp(['Locations used: ' num2str(locations_3(i,j))]);
        disp(['Vp and Cp are: ' num2str([metricsout.Vp metricsout.Cp])]);
    end
end

save('saves/study3');



