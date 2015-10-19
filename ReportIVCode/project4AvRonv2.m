% Av-Ron two-variable simplification of the Hodgkin-Huxley model
% A minimal biophysical model for an excitable and oscillatory neuron. Biol Cyb 65:487-500 (1991)
% Two neurons + synapses
% VNa = +55 mV, VK = -72 mv, VL = -49.4 mv
% resting membrane potential is -60 mV

function main()
clc; clear all; close all;

global pmbm iext psyn;

%fixed parameters in the model:  Cm =1uF/cm2, mp=3 (for m_inf), wp=4 and s=1.3 (for W)
%iext             1       2     3      4       5      6     7        8       9      10       11
%pmbm[i]          GNa     GK    GL     VNA     VK     VL    gamma    a^w    V1/2^w   a^m     V1/2^w
pmbm  =  [120,    36,    0.3,   55,  -72,   -49.4,   0.2,   0.045,   -55,   0.055,   -33,      0]';
psyn = [1.85, 55, 2, 1.85, 55, 2, 10];

dt=0.01;      % time increment (integration step dt)
t_end =30;  % 1sec
pulse =1;     % pulse width

period = 5;
steps = 20;
tstep = [period/steps: period/steps: period];
t0 = 2;

pulsedt=pulse/dt; % pulse width expressed in dt

i1_list = [];
for i = 1:steps
    f1=zeros(1,t_end/dt);
    f1(t0/dt:t0/dt+pulsedt) = 15;
    f1(t0/dt + period*2/dt + tstep(i)/dt: t0/dt + period*2/dt + tstep(i)/dt + pulsedt) = 15;
    f1=[f1,0];  % add one element at the end (array index cannot start form 0 in Matlab)
    if(i == 1)
        i1_list = f1;
    else
        i1_list = [i1_list;f1];
    end
end
% external current parameters
% amplitude*(t > tstart) -amplitude*(t > tstop)

%f1 = @(t) 0*(t >= 0) + 15*(t > 5) -15*(t > 6) + 15*(t > 18) -15*(t > 19);
f2 = @(t) 0*(t >= 0);% + 15*(t > 5) -15*(t > 6);
%t = linspace(0,25.01,25*100+1); %(0,25.01,25*100+1);
t = [0:dt:t_end];

for i = 1:steps
    iext=[i1_list(i,:) ; f2(t)];
    
    % set integration time, initial conditions and options for solver
    %tspan = linspace(0,25,2500);
    tspan = [0:dt:t_end];
    y0=[-59.4,0.40215,-59.4,0.40215,0,0,0,0];
    options = odeset('Jacobian', @mbmodejac, 'Vectorized', 'on', 'InitialStep', 0.01, 'MaxStep', 0.01);
    
    % now use ode45 solver
    [t1,y] = ode45(@mbmode, tspan, y0, options);
    
    disp('Done!')
    
    % plot V vs. t neuron #1
    figure;
    %subplot(3,1,1);
    plot(t1,y(:,1));
    hold on;
    plot(t1, y(:,3), '--');    
    hold on;
    plot(t1, iext(1,:), 'r');

    
    legend('Neuron 1', 'Neuron 2', 'Ext. Curr.');
    xlabel('time [ms]');
    ylabel('V [mV]');
    ylim([-80 60]);
    set(gca,'YTick',-80:40:60);
    title(strcat('Delay = ',num2str(tstep(i))))
 %   print(['D:\Users\Richard\Documents\2015 Fall\ICM\Project 4 Anderson\' num2str(tstep(i)*10) '.png'], '-dpng');
    % subplot(3,1,2);
    % plot(t, iext(1,:),'r');
    % ylim([0 20]);
    % xlabel('t [ms]');
    % ylabel('Iext');
    % title('Iext vs. t  neuron #1');
    %
    %
    % %plot V vs. t neuron #2
    % subplot(3,1,3);
    % plot(t, y(:,3));
    % xlabel('time [ms]');
    % ylabel('V [mV]');
    % ylim([-80 60]);
    % set(gca,'YTick',-80:40:60);
    % title('V vs. t neuron #2');
    
    % add W variable and replot neuron #2
    %figure;
    %[ax, h4, h5] = plotyy(t, y(:,3), t, y(:,4));
    %axes(ax(1)); ylabel('V [mV]'); ylim([-80 60]);set(gca,'YTick',-80:20:60);
    %axes(ax(2)); ylabel('w')
    %xlabel('Time, ms.');
    %legend('W','V');
    %title('V, W vs. t (neuron #2)');
end

end

function ydot = mbmode(t, y)

% returns derivative of the state vector for Av-Ron model equations
% model parametere are in a column vector pmbm
%  index       1    2  3    4    5  6     7     8    9     10     11
%  pmbm[i] = [gNa, gK, gL, VNa, VK, VL, gamma, a^w, V.5^w, a^m,  V.5^w]'
%  iext is a column vector of external curren amplitude in each time step

%       ydot = mbmode(t,y);     returns dy/dt eval at t,y
%       jac = mbmodejac(t,y)    returns the Jacobian at t,y

global pmbm iext psyn;
wp=4;
mp=3;
s=1.3;

%equation (11) and (12) Av-Ron 1991
ydot = zeros(8,1);
ydot(1) = (iext(1,round((t+0.01)/0.01)) - pmbm(1)*(m_inf(y(1),pmbm)^mp)*(1-y(2)).*(y(1)-pmbm(4))...%I_ext - I_Na
    - pmbm(2)*((y(2)/s)^wp).*(y(1)-pmbm(5)) ... %I_K
    - pmbm(3).*(y(1)-pmbm(6))...%I_L
    - psyn(1).*y(5).*(y(3)-psyn(2))); %I_syn
ydot(2) = (W_inf(y(1),pmbm)-y(2))./tau(y(1),pmbm);
ydot(5) = (-y(5) + y(6))./psyn(3);
ydot(6) = (-y(6) + (y(3) > psyn(7)));

ydot(3) = (iext(2,round((t+0.01)/0.01)) - pmbm(1)*(m_inf(y(3),pmbm)^mp)*(1-y(4)).*(y(3)-pmbm(4))...%I_ext - I_Na
    - pmbm(2)*((y(4)/s)^wp).*(y(3)-pmbm(5)) ... %I_K
    - pmbm(3).*(y(3)-pmbm(6))...%I_L
    - psyn(4).*y(7).*(y(1)-psyn(5))); %I_syn
ydot(4) = (W_inf(y(3),pmbm)-y(4))./tau(y(3),pmbm);
ydot(7) = (-y(7) + y(8))./psyn(6);
ydot(8) = (-y(8) + (y(1) > psyn(7)));

return
end


function jac = mbmodejac(t,y0)
%returns Jacobian (using numjac) at t,y for Av-Ron eq.

t = 0; y = y0; fac = [];

%computes the Jacobian
[jac,fac]=numjac(@mbmode, t, y, mbmode(t,y), [1e-5;1e-5;1e-5;1e-5;1e-5;1e-5;1e-5;1e-5],fac,0);
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
