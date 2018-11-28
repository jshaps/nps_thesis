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
size_drv = [m_drv;r_drv];

m_hst = 10000;      % mass of debris (hst)
r_hst = 5e-3;       % radius of hst (km)
size_hst = [m_hst;r_hst];

%--------- determine mass properties of the system
m_e = (m_drv*m_hst)/(m_drv+m_hst);      % equivalent mass of system (massless tether)

% This parameter will eventually be variable
l_tether = 50;      % initial tether length [km]

% length of tether from hst to COM
l_hst = (m_drv*l_tether)/(m_hst + m_drv);     % [km]
% length of tether from drv to COM
l_drv = l_tether - l_hst;               %[km]

J_sys = dumbbell_inertia(size_hst,size_drv,l_tether);
J_sys_final = dumbbell_inertia(size_hst,size_drv,10);  % find system inertia with 10 km tether

% Define time to Perform Analysis
    Orbits = 5;                              % q=number of orbits
    t0 = 0;                                  % initial time                      [s]
    w0 = sqrt(mu/R_com^3);                    % orbital angular velocity of com    [rad/s]
    n = w0; % note: for circular orbits (AS IN THIS INSTANCE) w0 = n
    T = (2*pi)/w0;                            % orbital period (s)
    tf = (Orbits*T);                    % final time [s]. orbital period*number of orbits
    
% define orital characteristics
COMcoe = getOrbitalElements();
[r0_com,v0_com] = sv_from_coe(COMcoe);        % define initial state vector of system
    
 momentexample_shap(J_sys,tf)    