function [coe] = getOrbitalElements()
% orbital elements
%{
 this example: for COM of tethered system
%
  mu   - gravitational parameter (km^3;s^2)
  coe  - orbital elements [h e RA incl w TA]
         where
             h    = angular momentum (km^2/s)
             e    = eccentricity
             RA   = right ascension of the ascending node (rad)
             incl = inclination of the orbit (rad)
             w    = argument of perigee (rad)
             TA   = true anomaly (rad)
 %}

% packing order for consistency with curtis codes
% h    = coe(1);
% e    = coe(2);
% RA   = coe(3);
% incl = coe(4);
% w    = coe(5);
% TA   = coe(6);

global mu Re

% defined
ra = 500;   % apogee altitude
rp = 500;   % perigee altitude
MA = 23.78*pi/180;         % Mean anomaly at given epoch
incl = 28.5*pi/180;         % inclination [rad]
RA = 80*pi/180;       % right ascension of ascending node [rad]
w = 65*pi/180;             % argument of perigee [rad]

% calculated
a = Re + (ra + rp)/2;       % semi major axis (km)
e = (ra - rp)/(ra + rp);    % eccentricity
p = a*(1-e^2);              % semi-latus rectum
h = sqrt(mu*p);             % magnitude of orbital angular momentum
E = kepler_E(e, MA);         % Eccentric anomaly

%   use eccentric anomaly to calculate true anomaly
TA = 2*atan(sqrt((1+e)/(1-e))*tan(E/2));
TA = unwrap(TA);

% pack data
coe = [h e RA incl w TA];

end

