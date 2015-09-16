function x = est14_Parlikar(MAP,C,b2bP,Period,tau)
% EST02_WK  CO estimator 2: Windkessel 1st-order LTI RC circuit model
%
%   In:   PP  <nx1> vector  --- beat-to-beat pulse pressure
%         HR  <nx1> vector  --- beat-to-beat heart rate
%
%   Out:  X    <nx1> vector  --- estimated CO (uncalibrated)
% 
%   Written by James Sun (xinsun@mit.edu) on Nov 19, 2005.
%alpha = 2; %For triangular pulse shape
    
    %modPP = alpha * (MAP-Pdias);
 
    x = C * (b2bP/Period + MAP/tau);
