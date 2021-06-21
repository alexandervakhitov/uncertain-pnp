% Generated using GBSolver generator Copyright Martin Bujnak,
% Zuzana Kukelova, Tomas Pajdla CTU Prague 2008.
% 
% Please refer to the following paper, when using this code :
%     Kukelova Z., Bujnak M., Pajdla T., Automatic Generator of Minimal Problem Solvers,
%     ECCV 2008, Marseille, France, October 12-18, 2008

function [x y z t] = GB_Solver_3Order_4Variable_Symmetry(a1, b1, c1, d1) 
%loading template
template_name = 'template_pnp';
load(template_name);

%solver parameter
settings.debug      = 0;        % debug mode : very slow
settings.full_QR1   = 0;        % turn on sparse QR
settings.p_inverter = 'all';    % use 'all' to use the 'best' inverter
settings.real       = 1;        % 1: only output real solutions 

%call solver
[sols,stats] = opnp_orig.solver_pfold([a1;b1;c1;d1],settings);

%selecting solutions
x = sols(1,:);
y = sols(2,:);
z = sols(3,:);
t = sols(4,:);
end
