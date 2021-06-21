%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This toolbox illustrates how to use the Uncertainty-aware PnPL
% algorithms described in:
%
%       A. Vakhitov, L. Ferraz, A. Agudo, F. Moreno-Noguer
%       Uncertainty-Aware Camera Pose Estimation from Points and Lines 
%       CVPR 2021
%
% 
% This program is free software: you can redistribute it and/or modify
% it under the terms of the version 3 of the GNU General Public License
% as published by the Free Software Foundation.
% 
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
% General Public License for more details.       
% You should have received a copy of the GNU General Public License
% along with this program. If not, see <http://www.gnu.org/licenses/>.
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This script contains the plots for the 2D+3D noise experiments from 
% the main paper (Figure 3, middle row)
% To reproduce, run the script after main_ordinary_3d_noise_3d.m
%
function plot_results(matpath, save_pref, is_title)

addpath plot_funcs;

if nargin == 0
    matpath = 'results/ordinary3Dresults';
    save_pref = 'ord3d';
    is_title = false;
end
load(matpath);

npts  =  int32(npts);

close all;
yrange= [0 2];

i= 0; w= 375; h= 300; wnarr = 340;

wl = 800;
hl = 110;
dol = 1;


yrange= [min([method_list(:).deleted_med_r]) max([method_list(:).deleted_med_r])];
yrange = [0 3];
plotTitle = '';
xTitle = '';
if is_title
    plotTitle = 'Mean Rotation';
    xTitle = 'Number of Points';
end
figure('color','w','position',[w*i,100,w,h]);%i=i+1;
xdrawgraph(npts,yrange,method_list,'deleted_mean_r',plotTitle,...
    xTitle,'Rotation Err. (deg.)', save_pref, dol, wl, hl);
figure('color','w','position',[w*i,100+h,wnarr,h]);i=i+1;
xdrawgraph(npts,yrange,method_list,'deleted_med_r','',...
    '','', save_pref, dol, wl, hl);
yrange = [0 1.5];

% yrange= [min([method_list(:).deleted_med_t]) max([method_list(:).deleted_med_t])];
% yrange = [0 2];
% yrange = [0 15];
if is_title
    plotTitle = 'Mean Translation';
end
figure('color','w','position',[w*i,100,w,h]);%i=i+1;
xdrawgraph(npts,yrange,method_list,'deleted_mean_t',plotTitle,...
    xTitle,'Translation Err. (%)', save_pref, dol, wl, hl);

figure('color','w','position',[w*i,100+h,wnarr,h]);i=i+1;
% yrange = [0 5];
xdrawgraph(npts,yrange,method_list,'deleted_med_t','',...
    '','', save_pref, dol, wl, hl);


yrange= [0 1.2*max([method_list(:).deleted_med_e])];

figure('color','w','position',[w*i,100,w,h]);%i=i+1;
xdrawgraph(npts,yrange,method_list,'deleted_mean_e','Mean L2 Error',...
    'Number of Points','L2 error');

figure('color','w','position',[w*i,100+h,w,h]);%i=i+1;
xdrawgraph(npts,yrange,method_list,'deleted_med_e','Median L2 Error',...
    'Number of Points','L2 error');

yrange= [0 min(max(1,2*max([method_list(:).pfail])),100)];
figure('color','w','position',[w*i,100+2*h,w,h]);%i=i+1;
xdrawgraph(npts*100,yrange,method_list,'pfail','No solution x method',...
    'Number of Points','% method fails');
i=i+1;

yrange= [0 1.2*max([method_list(:).deleted_med_c])];

figure('color','w','position',[w*i,100,w,h]);%i=i+1;
xdrawgraph(npts,yrange,method_list,'deleted_mean_c','Mean Cost',...
    'Number of Points','Cost (ms)', save_pref, dol, wl, hl);

yrange = [0 15];
figure('color','w','position',[w*i,100+h,w,h]);i=i+1;
xdrawgraph(npts,yrange,method_list,'deleted_med_c','Points',...
    'Number of Points','Med. runtime (ms)', save_pref, dol, wl, hl);