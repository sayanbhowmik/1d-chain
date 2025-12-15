%% Generate data
clear;
clc

addpath('../../src')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% energies = zeros(201,1);d
% dists = linspace(1,20,201)';
% 
% for i = 1:201
%     dist = dists(i);
%     S = msparc_1d_chain(dist, 32, 0.01, 10, 'GGA_PBE');
%     energies(i) = S.Etotal; 
% end

% save_array = [dists, energies];

% save('1d_32atoms_1:0.1:20dist_GGA_PBE_energies.mat', 'save_array')

energies = zeros(21,1);
dists = linspace(1,20,21)';

for i = 1:21
    dist = dists(i);
    S = msparc_1d_chain(dist, 32, 0.01, 10, 'GGA_PBE');
    energies(i) = S.Etotal; 
end


plot(dists, energies);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% energies = zeros(201,1);
% dists = linspace(1,20,201)';
% 
% for i = 1:201
%     dist = dists(i);
%     S = msparc_1d_chain(dist, 32, 0.01, 10, 'LDA_PW');
%     energies(i) = S.Etotal; 
% end
% 
% save_array = [dists, energies];
% 
% save('1d_32atoms_1:0.1:20dist_LDA_PW_energies.mat', 'save_array')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% energies = zeros(201,1);
% dists = linspace(1,20,201)';
% 
% for i = 1:201
%     dist = dists(i);
%     S = msparc_1d_chain(dist, 32, 0.01, 10, 'None');
%     energies(i) = S.Etotal; 
% end
% 
% save_array = [dists, energies];
% 
% save('1d_32atoms_1:0.1:20dist_NoXC_energies.mat', 'save_array')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Generate beautiful plots
% Define size in inches (e.g., single column width ~3.5 inches)
width = 7;
height = 4.5;
f = figure('Units', 'inches', 'Position', [1 1 width height]);
% f.OuterPosition  = [0,0,18,14];
tiledlayout(1,1, 'Padding', 'compact', 'TileSpacing', 'compact');

plot(dists, energies,'b-','LineWidth', 1.2);

xlabel('Interatomic distance (bohr)', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('$E$ (ha)', 'Interpreter', 'latex', 'FontSize', 14);

box on;               % Enclose the plot in a box
set(gca, ...
    'LineWidth', 1.5, ...       % Thicker axis lines
    'FontSize', 16, ...         % Tick label size
    'FontName', 'Helvetica', ...% Professional font (sans-serif is usually best for plots)
    'TickLabelInterpreter', 'latex'); % LaTeX for tick numbers
set(gca, 'YTickLabel', num2str(energies, '%.4f'));

% xlim([min(dists) max(dists)]);
xlim([6 max(dists)]);

% exportgraphics(f, 'Sigma2.0.pdf', 'Resolution', 300);