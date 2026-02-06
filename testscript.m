[x1, x2]=solveCoupled(1,1,.1,.1)
function [x1, x2]=solveCoupled(w1,w2,k1,k2)
%solves pair of position-coupled harmonic oscillators with individual frequencies
%w1,w2,and corresponding coupling strengths k1 (osc 1 to osc2) and k2 (osc 2 to osc1).
    A=[0,1,0,0;-w1^2,0,k1,0;0,0,0,1;k2,0,-w2^2,0];
    evolve=@(t) (expm(A*t)*[1,zeros(1,3)]');
    sol=arrayfun(evolve,linspace(0,0.67,1000));%(1:2,:);
    x1=sol(1,:);x2=sol(3,:);
end