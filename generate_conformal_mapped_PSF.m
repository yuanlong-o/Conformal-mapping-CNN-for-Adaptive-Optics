clc, clear
close all
%% generate conformal mapping grid for a distorted PSF
%  last update: 3/18/2019. YZ


psf = imread('test_PSF.png');
psf = double(psf);
psf = psf - min(psf(:));
psf = psf / max(psf(:));
N_in = size(psf, 1);
%% conformal mapping part, we use free sampling strategies
% input radius sampling number
N_out_r = N_in; 
% input angular sampling number
N_out_a = N_in; 

% upsampling ratio
super_rat = 5;

% sampling start point, since 0 could not be reached
sample_start = 0.1;

% min and max angular range
angle_min = -pi;
angle_max = pi;

% sampling grid 
angular_array = linspace(angle_min + 1e-2, angle_max - 1e-2, (N_out_a + 1) * super_rat);
rad_array = linspace(sample_start, 1, (N_out_r + 1) * super_rat);

% 
mapping_matrix = cell(N_out_r, N_out_r);
mapped_psf = zeros(N_out_r, N_out_r);

X_in = zeros(N_out_r, N_out_a);
Y_in = X_in;
X_out = zeros(N_out_r, N_out_a);
Y_out = X_out;
% for draw


I_out_cut = zeros(N_in);
% tes 1
I_out_cut(end / 2 - 1 : end / 2 + 2, 1 : end / 2) = 1; % test 1
% test 2
% I_out_cut(end / 2 - 1 : end / 2 + 2, :) = 1;
% test 3
% I_out_cut(1 : end / 2 , end / 2 - 1 : end / 2 + 2) = 1;

%% show mapping grid

angular_array_grid = linspace(angle_min + 1e-2, angle_max - 1e-2, (N_out_a + 1));
rad_array_grid = linspace(sample_start, 1, (N_out_r + 1));
for i = 1 : N_out_r + 1
    for j = 1 : N_out_a + 1

        rad_inter = i; % put to 0-1
        angular_inter = j;

        rad =rad_array_grid(rad_inter);
        angular =angular_array_grid(angular_inter);

        psi = rad * cos(angular);
        eta = rad * sin(angular);

        psi_round = round((psi + 1) / 2 * N_in);
        eta_round  = round((eta + 1) / 2 * N_in);

        psi_round  = min(max(psi_round , 1), N_in);
        eta_round  = min(max(eta_round , 1), N_in);

        X_in(i, j) = i;
        Y_in(i, j) = j;
        X_out(i, j) = psi;
        Y_out(i, j) = eta;
    end
end

    
%% show mapped PSF
I_out_cutT = I_out_cut.';
for i = 1 : N_out_r
    for j = 1 : N_out_a
        buf = [];
        for k = 1 : super_rat
            rad_inter = ((i - 1) * super_rat + k); % put to 0-1
            angular_inter = ((j - 1) * super_rat + k);

            rad = rad_array(rad_inter);
            angular = angular_array(angular_inter);

            psi = rad * cos(angular);
            eta = rad * sin(angular);

            psi_round = round((psi + 1) / 2 * N_in);
            eta_round  = round((eta + 1) / 2 * N_in);

            psi_round  = min(max(psi_round , 1), N_in);
            eta_round  = min(max(eta_round , 1), N_in);

            buf = [buf; psi_round , eta_round];          
            % record the mapping matrix
            mapped_psf(i, j) = mapped_psf(i, j) + psf(psi_round, eta_round);

        end
        mapping_matrix(i, j) = {buf};
        mapped_psf(i, j) = mapped_psf(i, j)  / super_rat;

    end
end
    
figure, imagesc(mapped_psf.') % well this is because the imagesc is transposed...
axis equal
axis off
colormap(hot)

%% plot
figure(1), subplot(1, 2, 1), imagesc(psf), axis equal, axis off, colormap(hot), title('distorted PSF')
subplot(1, 2, 2), imagesc(mapped_psf.'), axis equal, axis off, colormap(hot), title('conformal mapped PSF')

    
figure(2);
subplot(1, 2, 1),
plotMesh(conformalSetupAxes(gca, [1, N_out_r], [1, N_out_a]),X_in,Y_in);
axis off
title('traditional grid')

subplot(1, 2, 2),
plotMesh(conformalSetupAxes(gca, [-1, 1], [-1, 1]),X_out,Y_out);
axis off
title('conformal grid')

%% function
function ax = conformalSetupAxes(ax, xData, yData)
% conformalSetupOutputAxes Set up axes in the output/'z' plane.
%
% Supports conformal transformation example, ConformalMappingImageExample.m
% ("Exploring a Conformal Mapping").

% Copyright 2005-2013 The MathWorks, Inc. 

set(ax, 'PlotBoxAspectRatio',[1 1 1],...
		'XLimMode','manual',...
        'YLimMode','manual');
set(ax,'XLim',xData);
set(ax,'YLim',yData);
% set(ax,'Xlabel',text('String','Re(z)'));
% set(ax,'Ylabel',text('String','Im(z)'));
% set(ax,'Title',text('String','Output Plane'));
end

function plotMesh(varargin)
% Plots a mesh on the axes AX with color COLOR, via calls
% to 'LINE'.
%
%  PLOTMESH(AX,X,Y,COLOR) accepts two M-by-N arrays X and Y
%  -- like the output of MESHGRID.
%
%  PLOTMESH(AX,XY,COLOR) accepts a single M-by-N-by-2 array XY
%  -- like the output of TFORMFWD.

%  based on the the different grid index, I would change the color

if nargin == 2
  ax = varargin{1};
  XY = varargin{2};

  X = XY(:,:,1);
  Y = XY(:,:,2);
else
  ax = varargin{1};
  X  = varargin{2};
  Y  = varargin{3};
end
  
color_1 = linspace(0, 1, size(X,1));
color_2 = linspace(0, 1, size(X,2));

for k = 1:size(X,1)
%     line(X(k,:),Y(k,:),'Parent',ax,'Color', [0,  color_1(k), 0]);
    line(X(k,:),Y(k,:),'Parent',ax,'Color', [0, 0, 0], 'Linewidth', 1);
end

for k = 1:size(X,2)
%     line(X(:,k),Y(:,k),'Parent',ax,'Color',[color_2(k), 0.1, 1 - color_2(k)]);
line(X(:,k),Y(:,k),'Parent',ax,'Color',[0, 0, 0], 'Linewidth', 1);
end

end