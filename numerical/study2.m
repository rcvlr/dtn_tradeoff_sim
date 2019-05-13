% study2.m
% FILE 32
%
% In this study we run the elliptical FR setting for various combinations
% of the node density lambda and the node bending rate r0

% Initialization
clear;
startup;

% FIXED PARAMETERS
v0=1;        
alpha=1;     
epsilon=0.7; 
FR_type=[3 alpha epsilon];  

% numerical parameters
% Each data point: around 12 minutes and max out on memory!
N=256;      
S=400;     
targetpoints=250; 
beta=alpha*sqrt(1-epsilon^2);
area=pi*alpha*beta;
B=1.1*alpha*(1+epsilon);
L=round(2*B*sqrt(targetpoints/area));
num_pars=[N B L S]; 

potential_type=1;  
cost_type=1; 
density_type=1;  

verbose=0; 

% VARIABLE PARAMETERS
lambdarange=[0.01 0.25 0.5:0.5:5];
r0range=[0.1 1 2 4];
I=length(lambdarange);
J=length(r0range);

% WE INITIALIZE OUTPUT
Vp_num_2=zeros(I,J);
Cp_num_2=zeros(I,J);
locations_2=zeros(I,J);

% WE RUN SCRIPTS
for i=1:I
    
    lambda=lambdarange(i);
       
    for j=1:J
         
        r0=r0range(j); 
        node_pars=[lambda v0 r0];
        
        % We produce quick update
        st1=num2str(i);
        st2=num2str(length(lambdarange));
        st3=num2str(j);
        st4=num2str(length(r0range));
        disp(['Starting iteration: ' st1  '/' st2 ', ' st3 '/' st4 ' at time']);
        disp(datetime('now','Format','HH:mm:ss'));   
        
        metricsout=metrics(node_pars,FR_type,potential_type,cost_type,density_type,num_pars,verbose);

        % We produce output
        Vp_num_2(i,j)=metricsout.Vp;
        Cp_num_2(i,j)=metricsout.Cp;  
        locations_2(i,j)=length(metricsout.locations(:,1));
        disp(['Locations used: ' num2str(length(metricsout.locations(:,1)))]);
        disp(['Vp and Cp are: ' num2str([metricsout.Vp metricsout.Cp])]);
    end
end

save('saves/study2');
