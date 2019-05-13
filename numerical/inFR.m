% function inFR=inFR(FR_type,x,y)
% FILE 5
%
% This function specifies if a given point belongs to a given FR. 
%
% Input Parameters
% ================
% FR_type: The Forwarding Region used. See the boundary.m routine for a
%     description of possible choices.
% x,y: the cartesian coordinates of the point under investigation.
%
% Output Parameters
% =================
% inFR: This is 1 if the point belongs to the FR and 0 otherwise. 

function inFR=inFR(FR_type,x,y)

phi=atan2(y,x);
if (x^2+y^2)<(boundary(FR_type,phi))^2
     inFR=1;
else, inFR=0;
end
