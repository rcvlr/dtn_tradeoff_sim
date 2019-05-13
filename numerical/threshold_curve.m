% function [bbx,bby,tx,ty,der]=threshold_curve(potential_type,FR_type,theta,thetaprime,s)
% FILE 6
%
% This function returns the threshold curve where the potential
% exceeds a threshold, in particular, the following holds:
%
%                U(thetaprime,r)=U(theta,0)
%
% Read the paper for details. 
%
% Input Parameters
% ================
% potential_type: A parameter specifying the type of potential. See the
%      potential routine for information.
% FR_type: A parameter specifying the type of forwarding region. See the
%      boundary routine for information.
% theta: the direction of travel of the node.
% thetaprime: the direction of travel of nodes we are considered for eligibility.
% s: a parametrization of the curve, of length S. it is an increasing
%    sequence of numbers within [0,1].
%
% Output Parameters 
% =================
% bbx,bby(1,S): the cartesian coordinates of the curve for given s(k).
% tx,ty(1,S): the unit vector that is perpenticular to the curve, for given s(k).
% der(1,S): the incremental length covered if s(k) changes by ds is der(k)*ds.

function [bbx,bby,tx,ty,der]=threshold_curve(potential_type,FR_type,theta,thetaprime,s)

% Initialization
S=length(s);
bbx=zeros(1,S); bby=zeros(1,S); tx=zeros(1,S); ty=zeros(1,S); der=zeros(1,S);

if potential_type==1 && abs(thetaprime)<abs(theta)  
    
    for m=1:S
        phi=pi*(2*s(m)-1);    % s is proportional to phi
        r=boundary(FR_type,phi);
        bbx(m)=r*cos(phi);
        bby(m)=r*sin(phi);
        ds=10^(-5);
        phiplus=pi*(2*(s(m)+ds)-1);
        rplus=boundary(FR_type,phiplus);
        bbxplus=rplus*cos(phiplus);
        bbyplus=rplus*sin(phiplus);
        dist=( (bbxplus-bbx(m))^2 + (bbyplus-bby(m))^2 )^0.5;
        der(m)=dist/ds;
        ty(m)=-(bbxplus-bbx(m))/dist;
        tx(m)=(bbyplus-bby(m))/dist;     
    end
    
elseif potential_type==2 % the curve is a line segment on the y axis  
                          % (0<s<0.5) and the boundary of the FR in the 1st 
                          % and 4th quadrants (0.5<s<1).
    lim1=boundary(FR_type,pi/2);
    lim2=-boundary(FR_type,-pi/2);
    for m=1:S
         if s(m)<0.5 % we are in straight line
             bbx(m)=0;
             bby(m)=lim1-2*s(m)*(lim1-lim2);
             tx(m)=-1;
             ty(m)=0;
             der(m)=2*(lim1-lim2);
         else  % we are in the non-straight part
             phi=-pi/2+2*pi*(s(m)-0.5);
             bbx(m)=boundary(FR_type,phi)*cos(phi);
             bby(m)=boundary(FR_type,phi)*sin(phi);
             ds=10^(-5);
             phiplus=-pi/2+2*pi*(s(m)+ds-0.5);
             rplus=boundary(FR_type,phiplus);
             bbxplus=rplus*cos(phiplus);
             bbyplus=rplus*sin(phiplus);
             dist=( (bbxplus-bbx(m))^2 + (bbyplus-bby(m))^2 )^0.5;
             der(m)=dist/ds;
             ty(m)=-(bbxplus-bbx(m))/dist;
             tx(m)=(bbyplus-bby(m))/dist;  
         end         
    end
end










