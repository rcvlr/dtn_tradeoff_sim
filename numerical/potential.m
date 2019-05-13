% function potential=potential(potential_type,theta,x,y)
% FILE 2
%
% This function returns the potential of a node based on its direction of 
% travel and location, as well as the potential type. 
%
% Input Parameters
% ================
% potential_type: A parameter specifying the type of potential, as follows:
%     potential_type=1:  potential=-|theta|
%     potential_type=2:  potential= progress over cost
% theta: the direction of travel of the node
% x,y: the cartesian coordinates of the node 
%
% Output Parameters
% =================
% potential: the potential of the node

function potential=potential(potential_type, theta,x,y)

if potential_type==1
    potential=-abs(theta);
elseif potential_type==2
    if x==0 && y==0, potential=0;
    else
        potential=x/(x^2+y^2);
    end
end