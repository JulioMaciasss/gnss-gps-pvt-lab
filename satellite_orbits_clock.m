function [satp, satv, dtS] = satellite_orbits_clock(t, Eph, sat)

% Modified by @Omar Garcia Crespillo.
%  - Now it is using the toe.
%  - Limited to GPS processing
% -----------------------------------------------------------
% SYNTAX:
%   [satp, satv] = satellite_orbits(t, Eph, sat, sbas);
%
% INPUT:
%   t = clock-corrected GPS time
%   Eph  = ephemeris matrix
%   sat  = satellite index
%
% OUTPUT:
%   satp = satellite position (X,Y,Z)
%   satv = satellite velocity
%
% DESCRIPTION:
%   Computation of the satellite position (X,Y,Z) and velocity by means
%   of its ephemerides.

%----------------------------------------------------------------------------------------------
%                           goGPS v0.4.3
%
% Copyright (C) 2009-2014 Mirko Reguzzoni, Eugenio Realini
%----------------------------------------------------------------------------------------------
%
%    This program is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this program.  If not, see <http://www.gnu.org/licenses/>.
%----------------------------------------------------------------------------------------------

% Constants
Omegae_dot = 7.2921151467e-5;
F = -4.442807633e-10;

%% get ephemerides
roota     = Eph(4);
ecc       = Eph(6);
omega     = Eph(7);
cuc       = Eph(8);
cus       = Eph(9);
crc       = Eph(10);
crs       = Eph(11);
i0        = Eph(12);
IDOT      = Eph(13);
cic       = Eph(14);
cis       = Eph(15);
Omega0    = Eph(16);
Omega_dot = Eph(17);
toe       = Eph(18);
time_eph  = Eph(32);
% Parameters for clock correction
af2 = Eph(2);
af0 = Eph(19);
af1 = Eph(20);
toc = Eph(21);

%-------------------------------------------------------------------------------
% ALGORITHM FOR THE COMPUTATION OF THE SATELLITE COORDINATES (IS-GPS-200E)
%-------------------------------------------------------------------------------

%eccentric anomaly
[Ek, n] = ecc_anomaly(t, Eph);


A = roota*roota;             %semi-major axis
tk = check_t(t - toe);       %time from the ephemeris reference epoch

fk = atan2(sqrt(1-ecc^2)*sin(Ek), cos(Ek) - ecc);    %true anomaly
phik = fk + omega;                           %argument of latitude

uk = phik                + cuc*cos(2*phik) + cus*sin(2*phik); %corrected argument of latitude
rk = A*(1 - ecc*cos(Ek)) + crc*cos(2*phik) + crs*sin(2*phik); %corrected radial distance
ik = i0 + IDOT*tk        + cic*cos(2*phik) + cis*sin(2*phik); %corrected inclination of the orbital plane

%satellite positions in the orbital plane
x1k = cos(uk)*rk;
y1k = sin(uk)*rk;


%corrected longitude of the ascending node
Omegak = Omega0 + (Omega_dot - Omegae_dot)*tk - Omegae_dot*toe;

%satellite Earth-fixed coordinates (X,Y,Z)
xk = x1k*cos(Omegak) - y1k*cos(ik)*sin(Omegak);
yk = x1k*sin(Omegak) + y1k*cos(ik)*cos(Omegak);
zk = y1k*sin(ik);

% output variables
satp(1,1) = xk;
satp(2,1) = yk;
satp(3,1) = zk;

dt = check_t(t - toc);
%% Clock Correction =====================================
dtR = 0;
%dtR = F * ecc * sqrt(a)*sin(Ek)
dtS = af0+af1*(t-toc)+af2*(t-toc)^2 + dtR; % 1. CALCULAR LA CORRECCION DEL RELOJ AQUI
%% Relativistic correction
% 2. Calcular la correccion relativista despues


%% =======================================================

%-------------------------------------------------------------------------------
% ALGORITHM FOR THE COMPUTATION OF THE SATELLITE VELOCITY (as in Remondi,
% GPS Solutions (2004) 8:181-183 )
%-------------------------------------------------------------------------------
if (nargout > 1)
    Mk_dot = n;
    Ek_dot = Mk_dot/(1-ecc*cos(Ek));
    fk_dot = sin(Ek)*Ek_dot*(1+ecc*cos(fk)) / ((1-cos(Ek)*ecc)*sin(fk));
    phik_dot = fk_dot;
    uk_dot = phik_dot + 2*(cus*cos(2*phik)-cuc*sin(2*phik))*phik_dot;
    rk_dot = A*ecc*sin(Ek)*Ek_dot + 2*(crs*cos(2*phik)-crc*sin(2*phik))*phik_dot;
    ik_dot = IDOT + 2*(cis*cos(2*phik)-cic*sin(2*phik))*phik_dot;
    Omegak_dot = Omega_dot - Omegae_dot;
    x1k_dot = rk_dot*cos(uk) - y1k*uk_dot;
    y1k_dot = rk_dot*sin(uk) + x1k*uk_dot;
    xk_dot = x1k_dot*cos(Omegak) - y1k_dot*cos(ik)*sin(Omegak) + y1k*sin(ik)*sin(Omegak)*ik_dot - yk*Omegak_dot;
    yk_dot = x1k_dot*sin(Omegak) + y1k_dot*cos(ik)*cos(Omegak) - y1k*sin(ik)*ik_dot*cos(Omegak) + xk*Omegak_dot;
    zk_dot = y1k_dot*sin(ik) + y1k*cos(ik)*ik_dot;
    
    satv(1,1) = xk_dot;
    satv(2,1) = yk_dot;
    satv(3,1) = zk_dot;
end
