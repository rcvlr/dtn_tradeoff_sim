% function cost=cost(cost_type,x,y)
% FILE 3
% 
% This function calculates the cost of transmitting at a receiver with
% cartesian coordinates (with respect to the current holder) (x,y)
%
% Input Parameters
% ================
% cost_type: The type of cost used, as follows (d=sqrt(x^2+y^2)):
%    cost_type=1: cost=d^2
%    cost_type=2: cost=d^4
% x,y: the cartesian coordinates of the receiver
%
% Output Parameters
% =================
% cost: the transmission cost

function cost=cost(cost_type,x,y)

if cost_type==1
    cost=x^2+y^2;
elseif cost_type==2
    cost=(x^2+y^2)^2;
end