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
% This script contains the plots for the 2D noise experiments 
% from the main paper (Figure 3, top row)
%
% To reproduce, run the script after main_ordinary_3d_noise_2d.m

function plot(matpath, save_pref)
if nargin == 0
    matpath = 'results/ordinary2Dresults';
    save_pref = 'ord2d';
end
addpath plot_funcs;

load(matpath);

npts = int32(npts);

close all;

w= 375; 
h= 300; 
wnarr = 340;

i= 0; 
wl = 700;
hl = 100;
dol = 1;

%Rotation
yrange = [0 1.5];
figure('color','w','position',[w*i,100,w,h]);
xdrawgraph(npts, yrange, method_list, 'deleted_mean_r', 'Mean Rotation',...
           '','Rotation Err. (deg.)', save_pref, dol, wl, hl);
%Do not draw a legend any more
dol = 0;
figure('color','w','position',[w*i,100+h,wnarr,h]);
i=i+1;
xdrawgraph(npts,yrange,method_list,'deleted_med_r','Median Rotation',...
           '','', save_pref, dol, wl, hl);

%Translation
yrange = [0 1.2*max([method_list(:).deleted_mean_t])];
figure('color','w','position',[w*i,100,w,h]);
xdrawgraph(npts,yrange,method_list,'deleted_mean_t','Mean Translation',...
           '','Translation Err. (%)', save_pref, dol, wl, hl);

figure('color','w','position',[w*i,100+h,wnarr ,h]);
i=i+1;
xdrawgraph(npts,yrange,method_list,'deleted_med_t','Median Translation',...
           '','', save_pref, dol, wl, hl);

yrange= [0 1.2*max([method_list(:).deleted_med_c])];
figure('color','w','position',[w*i,100,w,h]);
xdrawgraph(npts,yrange,method_list,'deleted_mean_c','Mean Time',...
    'Number of Points','Time (ms)', save_pref, dol, wl, hl);
