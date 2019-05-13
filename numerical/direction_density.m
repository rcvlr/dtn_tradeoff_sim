% function direction_density(density_type,theta)
% FILE 4
%
% This function returns the density of the direction with which nodes are
% traveling.
%
% NB: the direction density function must be an EVEN function, otherwise
% errors are introduced in the metrics routine!!!
%
% Input Parameters
% ================
% density_type: A parameter specifying the type of density, as follows:
%     density_type=1:  uniform, i.e., 1/(2pi)
%     density_type=[2 N k]: nonuniform. In particular the density is equal to 
%         N/(16*pi*k) in some directions around the four main directions
%         of the compass and zero in the rest. See the paper for details. This
%         peculiar arrangement is made to minimize discretization errors. 
%         N: must be a multiple of 8, and equal to the length of thetas in
%         the metrics routine (so there is an even number of thetas in each
%         quadrant).
%         k: must be between 1 and N/8. Theta_W=(pi/2)*k/(N/8)=4k*pi/N
% theta: the direction of travel of the node 
%
% Output Parameters
% =================
% direction_density: the direction density in a given direction

function direction_density=direction_density(density_type,theta)

if density_type==1
    direction_density=1/(2*pi);
elseif density_type(1)==2
    N=density_type(2);
    k=density_type(3);   
    if abs(theta-(pi/2)*round(theta/(pi/2)))<=2*pi*k/N
        direction_density=N/(16*pi*k);
    else
        direction_density=0;
    end
else
    direction_density=1/(2*pi);
    disp('ERROR, unknown density type, using uniform density');
end