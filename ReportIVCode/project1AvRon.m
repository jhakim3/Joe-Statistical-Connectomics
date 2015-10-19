% Av-Ron two-variable simplification of the Hodgkin-Huxley model 
% A minimal biophysical model for an excitable and oscillatory neuron. Biol Cyb 65:487-500 (1991) 
% One neuron
% VNa = +55 mV, VK = -72 mv, VL = -49.4 mv
% resting membrane potential is -60 mV 

function main()
clc; clear all; close all;

global pmbm iext;

%fixed parameters in the model:  Cm =1uF/cm2, mp=3 (for m_inf), wp=4 and s=1.3 (for W)
%i                 1       2     3      4       5      6     7        8       9      10       11
%pmbm[i]          GNa     GK    GL     VNA     VK     VL    gamma    a^w    V1/2^w   a^m     V1/2^w
pmbm  =  [120,    36,    0.3,   55,  -72,   -49.4,   0.2,   0.045,   -55,   0.055,   -33,      0]';

% external current parameters
% amplitude*(t > tstart) -amplitude*(t > tstop)


%8.793916
colorlist={'r','g','b','m'};

q=figure;
time2=[6.1,8.7,8.8,14];
for x=1:4

t2 = time2(x)
%t2=8.793916
%t2=8.793917
%t2=1000

f = @(t) 0*(t >= 0) + 15*(t > 5) -15*(t > 6) + 15*(t > t2) -15*(t > t2+1);
t = linspace(0,25.01,25*100+1); %(0,25.01,25*100+1);

iext=f(t);
 

% set integration time, initial conditions and options for solver
tspan = linspace(0,25,2500);
y0=[-59.4,0.40215];
options = odeset('Jacobian', @mbmodejac, 'Vectorized', 'on', 'InitialStep', 0.01, 'MaxStep', 0.01);

%  use ode45 solver 
[t,y] = ode45(@mbmode, tspan, y0, options);

disp('Done!')

% plot V  vs. t 
%q=figure;
h = plot(t, y(:,1),char(colorlist(x)));
xlabel('time [ms]');
ylabel('V [mV]'); 
ylim([-80 60]);
set(gca,'YTick',-80:20:60);

%title(strcat('Input delay = ',num2str(t2-5),' ms'))
title('Identification of end of refractory period')
hold on
if x==3
plot(t2,55,strcat('*',char(colorlist(x))))    
else
plot(t2,50,strcat('*',char(colorlist(x))))
end

%plot(t,f(t))
hgexport(q,'-clipboard')
end

% add W variable and replot
% 
% figure;
% [ax, h1, h2] = plotyy(t, y(:,1), t, y(:,2));
% axes(ax(1)); ylabel('V [mV]'); ylim([-80 60]);set(gca,'YTick',-80:20:60);
% axes(ax(2)); ylabel('w')
% xlabel('Time, ms.');
% legend('W','V');
% 
end

function ydot = mbmode(t, y)

% returns derivative of the state vector for Av-Ron model equations 
% model parameters in a column vector pmbm
%  index       1    2  3    4    5  6     7     8    9     10     11
%  pmbm[i] = [gNa, gK, gL, VNa, VK, VL, gamma, a^w, V.5^w, a^m,  V.5^w]'
%  iext is a column vector of external curren amplitude in each time step 

%       ydot = mbmode(t,y);     returns dy/dt eval at t,y
%       jac = mbmodejac(t,y)    returns the Jacobian at t,y

global pmbm iext
wp=4;
mp=3;
s=1.3;

%equation (11) and (12) Av-Ron 1991
	ydot = zeros(2,1);
	   ydot(1) = (iext(round((t+0.01)/0.01)) - pmbm(1)*(m_inf(y(1),pmbm)^mp)*(1-y(2)).*(y(1)-pmbm(4))...%I_ext - I_Na
           - pmbm(2)*((y(2)/s)^wp).*(y(1)-pmbm(5)) ... %I_K
           - pmbm(3).*(y(1)-pmbm(6))); %I_L
	ydot(2) = (W_inf(y(1),pmbm)-y(2))./tau(y(1),pmbm);
    
return 
end


function jac = mbmodejac(t,y0)
% returns Jacobian (using numjac) at t,y for Av-Ron eq.

t = 0; y = y0; fac = [];

%computes the Jacobian
[jac,fac]=numjac(@mbmode, t, y, mbmode(t,y), [1e-5;1e-5],fac,0);
end


function min = m_inf(v,pmbm)
% m_inf in eq.(14) 
   min = (1 + exp( -2*pmbm(10).*(v - pmbm(11)) ) )^-1;
return
end


function tau = tau(v,pmbm)
%tau in eq. (15) Av-Ron 1991

   tau = (pmbm(7)*exp(  pmbm(8).*(v - pmbm(9)) ) + ...
          pmbm(7)*exp( -pmbm(8).*(v - pmbm(9))   ) )^-1;
return
end


function win = W_inf(v,pmbm)
%W_inf in eq. (13) Av-Ron 1991
   win = (1 + exp( -2*pmbm(8).*(v - pmbm(9)) ) )^-1;
return
end
