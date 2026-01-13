%% Generate data
clear;
% clc

addpath('src')

ang_to_bohr = 1.88973;

r1 = 4.0 * ang_to_bohr;
box = 8.0 * ang_to_bohr;
% box = 128 * ang_to_bohr;

% S = msparc_1d_chain(4., -1, 128, 0.01, 10, 'HSE');

S = msparc_1d_chain(r1, -1, box, 0.01, 10, 'HSE');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Generate beautiful plots

% Example data (replace with your own vectors)
eigvals = S.EigVal;
occ_index = S.occ > 0.75;
occupied = eigvals(occ_index);
occupied_is = linspace(1, length(occupied), length(occupied));
unoccupied = eigvals(~occ_index);
unoccupied_is = linspace(occupied_is(end)+1, occupied_is(end)+1+length(unoccupied), length(unoccupied));

figure;
hold on;
plot(occupied_is, occupied, 'bo', 'MarkerFaceColor','none', 'LineStyle','none'); % blue circles
plot(unoccupied_is, unoccupied, 'r^', 'MarkerFaceColor','none', 'LineStyle','none'); % red triangles
hold off;

xlabel('Index (i)');
ylabel('$\epsilon_i$', 'Interpreter','latex');
axis square
grid on;
legend({'occupied','unoccupied'}, 'Location','northwest', 'FontSize', 14);