clear;
clc

addpath('../src')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% energies = zeros(201,1);
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