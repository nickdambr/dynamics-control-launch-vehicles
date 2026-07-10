function [Abar, Bbar, cbar] = lti_zoh(dt, Vc)
% Exact ZOH discretisation of the LTI GFOLD log-mass system (van Loan / expm).
%   State    xi = [x; y; vx; vy; z],  z = ln(m)
%   Control  w  = [ux; uy; sigma],    u = T/m,  sigma >= ||u||
%   Dynamics xdot=vx, ydot=vy, vxdot=ux, vydot=uy-1, zdot=-Vc*sigma   (LTI)
%
%   INPUT
%     dt - interval length (nondim)
%     Vc - V_ref/c (Tsiolkovsky number)
%   OUTPUT
%     Abar (5x5), Bbar (5x3), cbar (5x1) - discrete-time matrices, CONSTANT
%       across the grid because the system is time-invariant:
%       xi_{k+1} = Abar*xi_k + Bbar*w_k + cbar.
%
% The change of variables z=ln(m), u=T/m linearises the dynamics exactly, so
% the appendix-A ZOH reduces to a single matrix exponential of the augmented
% system    expm([A B c; 0]*dt) = [Abar Bbar cbar; 0 I]
% computed once -- no per-interval integration and no singular mass row.
    arguments
        dt (1,1) double {mustBePositive, mustBeFinite}
        Vc (1,1) double {mustBeFinite}
    end
    A = zeros(5);    A(1,3) = 1;  A(2,4) = 1;             % x<-vx, y<-vy
    B = zeros(5,3);  B(3,1) = 1;  B(4,2) = 1;  B(5,3) = -Vc;   % vx<-ux, vy<-uy, z<--Vc*sigma
    c = [0; 0; 0; -1; 0];                                 % gravity on vy
    E = expm([A, B, c; zeros(4,9)] * dt);
    Abar = E(1:5, 1:5);   Bbar = E(1:5, 6:8);   cbar = E(1:5, 9);
end
