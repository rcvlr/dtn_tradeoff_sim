% function better_angles=better_angle(N,i)
% FILE 7
% 
% This function provides the indices of angles that are STRICTLY better 
% than a given one i.
% 
% Input Parameters
% ================
% N: total angles
% i: given index

function better_angles=better_angles(N,i)

if i<=floor((N+1)/2)
    better_angles=i+1:(N-i);
else
    better_angles=(N-i+2):i-1;
end
    