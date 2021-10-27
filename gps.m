clear all;
clc;
load eph.dat;
load rcvr.dat;
%% GPS Const
c = 299792458.0;    %% speed of light
Wedot = 7.2921151467e-5;   %% Earth rotate rate
mu = 3.986005e+14;         %% Earth universal gravitation constant
F = -4.442807633e-10;      %% Relativistic correction term 

eph = sortrows(eph,2);     %% sort eph and rcvr according to SV Id
rcvr = sortrows(rcvr,2);



SV_position = zeros(size(eph,1),4);
SV_clock_bias = [];
%% Calculate Satellite position
for (i = 1: size(rcvr,1))
    rcvr_tow = rcvr(i,1);            %% receiver time of week            
    pr = rcvr(i,3);                  %% pseudorange
    t = rcvr_tow - pr / c;           %% system time at which signal left SV
    
    svid = eph(i,2);                 %% SV PRN number
    toc = eph(i,3);                  %% reference time of clock parameters
    toe = eph(i,4);                  %% reference time of ephmeris parameters
    af0 = eph(i,5);                  %% clock correction coefficient - clock bias
    af1 = eph(i,6);                  %% clock correction coefficient - drift
    af2 = eph(i,7);                  %% clock correction coefficient - frequency drift
    ura = eph(i,8);                  %% user range accuracy
    e   = eph(i,9);                  %% eccentricity
    sqrta = eph(i,10);               %% sqaure root of semi-major axis 
    dn  = eph(i,11);                 %% mean motion correction
    m0 =  eph(i,12);                 %% mean anomaly at reference time
    w =   eph(i,13);                 %% argument of perigee
    omg0 = eph(i,14);                %% right ascension
    i0 = eph(i,15);                  %% inclination angle at reference time
    odot = eph(i,16);                %% rate of right ascension
    idot = eph(i,17);                %% rate of inclination angle
    cus = eph(i,18);                 %% argument of latitude correction, sine
    cuc = eph(i,19);                 %% argument of latitude correction, cosine
    cis = eph(i,20);                 %% inclination correction, sine
    cic = eph(i,21);                 %% inclination correction, cosine
    crs = eph(i,22);                 %% radius correction, sine
    crc = eph(i,23);                 %% radius correction, cosine
    iod = eph(i,24);                 %% issue of data number
    
    A = sqrta.^2;                    %% semi-major axis
    n0 = sqrt(mu/A.^3);              %% computed mean motion
    tk = t - toe;                    %% Time from ephemeris reference epoch
    if(tk >= 302400)
        tk = tk - 604800;
    elseif(tk <= -302400)
        tk = tk + 604800;
    end;
    n = n0 + dn;                     %% corrected mean motion
    mk = m0 + n*tk;                  %% mean anomaly
    
    syms x;
    eqn = mk == x - e*sin(x);        %% Kepler's Equation for Eccentric anomaly
    Ek = vpasolve(eqn, x);
    Ek = double(Ek);
    
    vk = atan2((sqrt(1-e^2)*sin(Ek)/(1-e*cos(Ek))),((cos(Ek)-e)/(1-e*cos(Ek)))); %% True anomaly
    Phk = vk + w;                    %% Argument of latitude
    
    du = cus*sin(2*Phk) + cuc*cos(2*Phk);    %% Arugment of latitude correction
    dr = crs*sin(2*Phk) + crc*cos(2*Phk);    %% Radius correction
    di = cis*sin(2*Phk) + cic*cos(2*Phk);    %% Inclination correction
    
    uk = Phk + du;                   %% Corrected Argument of latitude
    rk = A * (1 - e * cos(Ek)) + dr; %% Corrected Radius
    ik = i0 + di + idot * tk;        %% Corrected Inclination
    
    xk_prime = rk * cos(uk);         %% Positions in orbital plane
    yk_prime = rk * sin(uk);
    
    Omega = omg0 + (odot - Wedot)* tk - Wedot * toe;  %% Corrected longitude of ascending node
    
    xk = xk_prime * cos(Omega) - yk_prime * cos(ik) * sin(Omega);  %% Earth-fixed coordinates
    yk = xk_prime * sin(Omega) + yk_prime * cos(ik) * cos(Omega);
    zk = yk_prime * sin(ik);
    
    SV_position(i,1) = svid;
    SV_position(i,2) = xk;
    SV_position(i,3) = yk;
    SV_position(i,4) = zk;
    
    %% Satellite clock bias
    delta_tr = F * e * sqrta * sin(Ek);                          %% relativistic correction term
    delta_tsv = af0 + af1*(t - toc) + af2*(t-toc)^2 + delta_tr;  %% satellite clock bias
    SV_clock_bias(i,1) = delta_tsv;
end


x_s = SV_position(:,2);  
y_s = SV_position(:,3);
z_s = SV_position(:,4);

%% Use Least Square to solve receiver position
x_0 = [-2694685.473, -4293642.366, 3857878.924, 0.0]';   %% initial state include position in ECEF coordinate frame and user clock bias
d_x = ones(4,1);    
i = 0; %% iterator
x_record = x_0;
while norm(d_x(1:3)) >= 1e-4
    rho = rcvr(:,3) + c * SV_clock_bias;    %% measurement is corrected pseudorange: pseudorange + satellite clock bias
    rho_ = sqrt((x_s - x_0(1)).^2 + (y_s - x_0(2)).^2 + (z_s - x_0(3)).^2) + c * x_0(4);  %% last predicted pseudorange
    d_rho = rho_ - rho;                 %% minimize the difference between measurement and prediction

    a_x = (x_s-x_0(1))./sqrt((x_s-x_0(1)).^2 + (y_s-x_0(2)).^2 + (z_s-x_0(3)).^2);
    a_y = (y_s-x_0(2))./sqrt((x_s-x_0(1)).^2 + (y_s-x_0(2)).^2 + (z_s-x_0(3)).^2);
    a_z = (z_s-x_0(3))./sqrt((x_s-x_0(1)).^2 + (y_s-x_0(2)).^2 + (z_s-x_0(3)).^2);
    H = [a_x a_y a_z ones(size(rcvr,1),1)]; 
    
    d_x = pinv(H) * d_rho;
    d_x(4) = d_x(4) / (-c);             %% receiver clock term is -c*dt
    x_0 = x_0 + d_x;                    %% update state variable
    x_record(:,i+2) = x_0;
    i = i+1; 

end

%% plot iteration function
titsize = 10;
lgdsize = 10;
fontsize = 16;
makersize =10;

figure 
set(gcf,'position',[0 0 2100*0.5 900*0.8]);

disp(x_0);
it = linspace(0, 3, 4);
plot_x =  x_record(1,:);
plot_y =  x_record(2,:);
plot_z =  x_record(3,:);
plot_t =  x_record(4,:);

subplot(4, 1, 1);
plot(it, plot_x,'*r');
axis tight
ax=gca; ax.FontSize=fontsize;
title(' ','FontSize',titsize);
ylabel ('x (m)');
xlabel ('iteration number');
legend('boxoff')

subplot(4, 1, 2);
plot(it, plot_y,'*r');
axis tight
ax=gca; ax.FontSize=fontsize;
title(' ','FontSize',titsize);
ylabel ('y (m)');
xlabel ('iteration number');
legend('boxoff')

subplot(4, 1, 3);
plot(it, plot_z,'*r');
axis tight
ax=gca; ax.FontSize=fontsize;
title(' ','FontSize',titsize);
ylabel ('z (m)');
xlabel ('iteration number');
legend('boxoff')

subplot(4, 1, 4);
plot(it, plot_t,'*r');
axis tight
ax=gca; ax.FontSize=fontsize;
title(' ','FontSize',titsize);
ylabel ('t (s)');
xlabel ('iteration number');
legend('boxoff')


