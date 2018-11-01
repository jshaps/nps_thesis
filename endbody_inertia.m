function [J] = endbody_intertia(m,r)
%UNTITLED Summary of this function goes here
%   Formulation of inertial tensor of solid spherical end body
i = (2/5)*m*r^2*ones(1,3);

J = diag(i);



end

