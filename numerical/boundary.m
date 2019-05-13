% function boundary=boundary(FR_type,phi)
% FILE 1
%
% This function returns the value of b(phi), where b describes 
% the boundary of a FR, for various choices of the FR.
%
% Input Parameters
% ================
% FR_type: The Forwarding Region used. Different choices exist, as follows:
%     FR_type=[1 R]: Forwarding region is a circle of radius R, therefore
%               boundary(phi)=R.
%     FR_type=[2 a b]: Forwarding region is a cardioid with parameters a,b:
%               boundary(phi)=a+b*cos(phi).
%     FR_type=[3 a epsilon]: Forwarding region is an ellipse with parameters a, epsilon:
%               boundary(phi)=a*(1-epsilon^2)/(1-epsilon*cos(phi)).
% phi: The angle phi for which the function b will be calculated
%
% Output Parameters
% =================
% boundary: the value b(phi)
% 
% NB: MUST BE SYMMETRIC WITH RESPECT TO xx' AXIS, OTHERWISE THE CALCULATION
% OF r_D MUST BE MODIFIED

function boundary=boundary(FR_type,phi)
if FR_type(1)==3, boundary=FR_type(2)*(1-(FR_type(3))^2)/(1-FR_type(3)*cos(phi)); % the most-commonly used is placed first
elseif FR_type(1)==2, boundary=FR_type(2)+FR_type(3)*cos(phi);
elseif FR_type(1)==1, boundary=FR_type(2); 
end