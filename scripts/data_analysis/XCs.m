energies = load('scripts/data_gen/1d_32atoms_1:0.1:20dist_GGA_PBE_energies.mat').save_array;

dists = energies(50:end,1);
PBE_energies = energies(50:end,2);

energies = load('scripts/data_gen/1d_32atoms_1:0.1:20dist_LDA_PW_energies.mat').save_array;
LDA_energies = energies(50:end,2);

energies = load('scripts/data_gen/1d_32atoms_1:0.1:20dist_NoXC_energies.mat').save_array;
NoXC_energies = energies(50:end,2);
shift = min(PBE_energies) - min(NoXC_energies);
NoXC_energies = NoXC_energies + shift;

figure
plot(dists, PBE_energies);
hold on;
plot(dists, LDA_energies, '--');
plot(dists, NoXC_energies, ':');
xlabel('Spacing (Bohr');
ylabel('Energy (Ha)');
legend('PBE Energies', 'LDA Energies', 'No XC Energies');