clear 
close all
clc

% define gravitational constant and radius of the earth
global mu Re 
mu = 398600; %km^3/s^2
Re = 6378;
R_com = Re + 500;       % orbital radius of com

%-------------
m_drv = 1000;       % mass of debris removal vehicle (drv)
r_drv = 10e-3;      % radius of drv (km)

m_hst = 10000;      % mass of debris (hst)
r_hst = 5e-3;       % radius of hst (km)

%--------- determine mass properties of the system
m_e = (m_drv*m_hst)/(m_drv+m_hst);      % equivalent mass of system (massless tether)

% This parameter will eventually be variable
l_tether = 50;      % initial tether length [km]

% length of tether from hst to COM
l_hst = (m_drv*l_tether)/(m_hst + m_drv);     % [km]
% length of tether from drv to COM
l_drv = l_tether - l_hst;               %[km]

% determine inertial properties
J_hst = endbody_inertia(m_hst,r_hst);
J_drv = endbody_inertia(m_drv,r_drv);

% also will eventually be variable
% Inertia of system in body frame
% y-axis aligned with tether
J_sys = [m_hst*l_hst^2+m_drv*l_drv^2    0   0;...
        0   J_drv(1)+J_hst(1)   0;...
        0   0   m_hst*l_hst^2+m_drv*l_drv^2];
    

t_gam = 15*pi/180;      % tether pitch angle