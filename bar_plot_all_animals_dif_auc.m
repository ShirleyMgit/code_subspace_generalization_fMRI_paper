clear all
close all
clc

exp_root = 'C:\Users\User\Documents\fMRI_EXP\';
more_code = fullfile(exp_root, 'cells_data_muller');
data_folder = fullfile(more_code, 'wetransfer_neil-s-data_2023-09-25_1057');
addpath(genpath(more_code));
addpath(genpath(data_folder));

place_cells_rats = [4, 6, 7, 8];
num_place_rats = length(place_cells_rats);
dif_place_cells = -ones(num_place_rats, 1);
v_num_place = -ones(num_place_rats, 1);
for n=1:num_place_rats
    load(['place_cells_animal',num2str(place_cells_rats(n)),'.mat'])

    env1_data = squeeze(cells_array(1, :, :));
    env1_data(isnan(env1_data)) = 0;
    env2_data = squeeze(cells_array(2, :, :));
    env2_data(isnan(env2_data)) = 0;

    cumsum_var_place11 = cal_projection_plot(transpose(env1_data), transpose(env1_data));
    cumsum_var_place12 = cal_projection_plot(transpose(env1_data), transpose(env2_data));

    num_place_cells = length(cumsum_var_place11);
    v_num_place(n) = num_place_cells;
    disp(num_place_cells)
    pc_fraction_place = linspace(0,1,length(cumsum_var_place11)+1);

    place_effect = (sum(cumsum_var_place11) - sum(cumsum_var_place12))/(num_place_cells);
    dif_place_cells(n) = place_effect ;

end

grid_rats = [1, 3, 4];
num_grid_rats = length(grid_rats);
dif_grid_cells = -ones(num_grid_rats,1);
v_num_grid = -ones(num_grid_rats,1);
for n=1:num_grid_rats
    load(['grid_animal',num2str(grid_rats(n)),'.mat'])

    env1_data = squeeze(cells_array(1, :, :));
    env1_data(isnan(env1_data)) = 0;
    env2_data = squeeze(cells_array(2, :, :));
    env2_data(isnan(env2_data)) = 0;

    cumsum_var_grid11 = cal_projection_plot(transpose(env1_data), transpose(env1_data));
    cumsum_var_grid12 = cal_projection_plot(transpose(env1_data), transpose(env2_data));

    num_grid_cells = length(cumsum_var_grid11);
    v_num_grid(n) = num_grid_cells;

    grid_effect = (sum(cumsum_var_grid11) - sum(cumsum_var_grid12))/num_grid_cells;
    dif_grid_cells(n) = grid_effect;

end

bar([0.8,1,1.2, 2.3,2.5,2.7,2.9], [dif_grid_cells; dif_place_cells])
xticklabels([]);

figure(2)
subplot(2,1,1)
plot(v_num_grid, dif_grid_cells,'*')
subplot(2,1,2)
plot(v_num_place, dif_place_cells,'*')

[h, p] = ttest2(dif_grid_cells, dif_place_cells, 'Tail', 'left')