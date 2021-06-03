from os.path import dirname, join, expanduser
from cc3d.CompuCellSetup.CC3DCaller import CC3DCaller

import setup_cells

import winsound
import numpy as np
import pandas as pd
# import matplotlib.pyplot as plt


winsound.Beep(440, 2000)
# CHECKLIST:
# mutate suriving genomes                       ###COMPLETED###
# align genome trajectories                     ###COMPLETED###
# add option for hexagonal mode
# implement network genome description
# investigate smallest mobile unit
# investigate effects of scaling                ### Fitness

# scaling first
# using symmetry in patterns to make stuff controable


def calculate_scores(return_object, func, space, grid):
    starts = return_object['result']['initial data']
    diffs = return_object['result']['final data'] - starts
    # if func == 'total displacement':
    #     return np.sum(np.sum(diffs ** 2, 1) ** 0.5)
    # elif func == 'average displacement':
    #     return np.sum((np.sum(diffs, 0) / diffs.shape[0]) ** 2) ** 0.5
    # elif func == 'vertical displacement':
    #     return diffs.sum(0)[1]
    if func == 'COM':  # == 'COMdist':
        # Calculate average magnitude of COM displacement vectors.
        vec_grid = np.zeros((grid, grid, 2))
        for x in range(0, grid):
            for y in range(0, grid):
                mask = [i for i in range(0, diffs.shape[0]) if x * space < starts[i][0] < (x + 1) * space and
                        y * space < starts[i][1] < (y + 1) * space]
                if len(mask) == 0:
                    vec_grid[x][y][:] = 0
                else:
                    cell_vecs = diffs[mask][:]
                    vec_grid[x][y][:] = np.sum(cell_vecs / len(mask), 0)
                    # Calculate average magnitude of COM displacement vectors.
                    # vec_grid[x][y] = np.sum(np.sum(cell_vecs ** 2, 1) ** 0.5) / len(mask)
                    # Calculate magnitude of average COM displacement vector.
                    # vec_grid[x][y][:] = np.sum(np.sum(cell_vecs / len(mask), 0) ** 2) ** 0.5
        return vec_grid


def measure_fitness(scores_dict):
    measures = np.zeros((len(scores_dict.keys()), 2))
    angles = np.zeros(len(scores_dict.keys()))
    for i in scores_dict.keys():
        vecs = np.zeros((len(scores_dict[i]), 2))
        for j, v in enumerate(scores_dict[i]):
            vecs[j, :] = v
        measures[i, 0] = ((vecs.sum(0) / vecs.shape[0]) ** 2).sum() ** 0.5  # magnitude of average vec
        measures[i, 1] = ((vecs ** 2).sum(1) ** 0.5).sum() / vecs.shape[0]  # average of vec magnitudes
        expected_vec = vecs.sum(0) / vecs.shape[0]
        angles[i] = np.arctan2(expected_vec[1], expected_vec[0])
    return measures, angles


def permute_genome(genome, rotations):
    if rotations == 0:
        return genome
    aligned_genome = np.zeros(genome.shape)
    s = genome.shape[0]
    r = lambda x: [x[2] - 1 - x[1], x[0], x[2]]
    for i in range(0, s):
        for j in range(0, s):
            r_map = [i, j, s]
            for _ in range(0, rotations):
                r_map = r(r_map)
            aligned_genome[r_map[0], r_map[1]] = genome[i, j]
    return aligned_genome


def apply_crossover(parents, parent_angles, phen_dims, phen_mode):
    new_genome = np.zeros((1, parents.shape[1]))
    angle_diff = parent_angles[1] - parent_angles[0]
    symmetry = 0
    multiplier = 0
    if phen_mode == 'car':
        symmetry = 4
        multiplier = 1
    elif 'hex' in phen_mode:
        symmetry = 2
        multiplier = 2
    num_rotations = int(np.mod(np.round(np.mod(angle_diff, 2 * np.pi) / (2 * np.pi / symmetry), 0), symmetry))
    aligned_parent = permute_genome(parents[1, :].reshape(phen_dims), int(multiplier * num_rotations))
    parents[1, :] = aligned_parent.reshape(1, -1)
    for j in range(0, parents.shape[1]):
        new_genome[0, j] = parents[np.random.choice([0, 1], 1)[0], j]
    return new_genome


def apply_mutation(gen, mut_rate, key_dict, skip):
    if mut_rate == 0:
        return gen
    mut_gen = 1 * gen
    for i in range(0, mut_gen.shape[0]):
        if i >= skip:
            for j in range(0, mut_gen.shape[1]):
                if np.random.random() < mut_rate:
                    mut_gen[i, j] = np.random.choice(list((set(key_dict.keys()) | {0}) - {gen[i, j]}))
    return mut_gen


def evolve_next_gen(order, current_gen, survive, offspring, mutate, key_dict, angles, phen_dims, phen_mode):
    population_size = current_gen.shape[0]
    survive_num = int(np.ceil(survive * population_size))
    offspring_num = int(np.floor(offspring * population_size))
    # random_num = population_size - survive_num - offspring_num
    new_gen = np.zeros(current_gen.shape)
    for i in range(0, population_size):
        if i < survive_num:
            new_gen[i, :] = current_gen[order[i], :]
        elif i < survive_num + offspring_num:
            parent_inds = np.random.choice(survive_num, 2, replace=False)
            new_gen[i, :] = apply_crossover(current_gen[order[parent_inds], :],
                                            angles[order[parent_inds]], phen_dims, phen_mode)
            # for j in range(0, current_gen.shape[1]):
            #     ind = np.random.choice([0, 1], 1)
            #     new_gen[i, j] = current_gen[order[parent_inds[ind[0]]]][j]
        else:
            new_gen[i, :] = np.random.choice([0] + list(key_dict.keys()), new_gen.shape[1])
        # mutate_mask = np.sign(np.sign(np.random.random(new_gen.shape) - mutate) + 1)
    new_mut_gen = apply_mutation(new_gen, mutate, key_dict, survive_num)
    return new_mut_gen


def get_genome(seed):
    # if isinstance(seed, str):
    #     gene = np.zeros((1, len(seed)))
    #     for i, num in enumerate(seed):
    #         gene[0, i] = int(num)
    #     return gene
    if isinstance(seed, list):
        genes = np.zeros((len(seed), len(seed[0])))
        for g in range(0, len(seed)):
            for i, num in enumerate(seed[g]):
                genes[g, i] = int(num)
        return genes
    else:
        return np.reshape(seed, (1, -1))


def save_trajectories(traj_list, root, traj_lab, file_name):
    x_data = np.zeros((len(traj_list[0][0]), len(traj_list)))
    y_data = np.zeros((len(traj_list[0][0]), len(traj_list)))
    for t, traj in enumerate(traj_list):
        for s in range(0, len(traj[0])):
            x_data[s, t] = traj[0][s]
            y_data[s, t] = traj[1][s]
    pd.DataFrame(x_data).to_csv(root + traj_lab + 'xData_' + file_name, header=False, index=False)
    pd.DataFrame(y_data).to_csv(root + traj_lab + 'yData_' + file_name, header=False, index=False)


def main():
    # Set up file locations
    root = r'C:/CompuCell3D-py3-64bit/lib/site-packages/ProtocellSimulator/'
    project_name = 'ProtoCellularSim'
    simulation_file_path = join(root, project_name + '/' + project_name + '.cc3d')
    init_cells_file_path = join(root, project_name + '/' + 'init_cell_field.piff')
    variables_file = join(root, project_name + '/' + 'variables.csv')
    root_output_folder = join(root, project_name + '/' + 'Results')
    abbrevs = {'ProtoCellularSim': 'proto'}
    gen_file_prefix = 'gen_'
    # record_file = 'fitness_records.csv'
    # record_file = 'fit_rec_55_COM.csv'
    # rec_file_prefix = 'fit_'

    # Set up parallel processing / multi-threading parameters
    # cpus = 4
    # threads_per_cpu = cpus
    # work_nodes = cpus * threads_per_cpu
    work_nodes = 4  # ############################################################################ 1 / 16
    grid_dim = int(work_nodes ** 0.5)
    phenomes_per_sim = 1  # ###################################################################### 1 / 4
    repeats = int(work_nodes / phenomes_per_sim)

    # Set up simulation parameters
    cell_diam = 10
    actuate_scale = 2
    actuate_period = 100
    num_actuation_cycles = 10  # ################################################################### 10 / 100
    key2name = {3: 'ProtocellPassive', 4: 'ProtocellActiveA', 5: 'ProtocellActiveB'}
    num_active_types = np.array([int('Active' in x) for x in key2name.values()]).sum()
    max_mcs = ((num_active_types + 2) * num_actuation_cycles + 1) * actuate_period  # num_active_types + 2 ???
    phenome_scaling = [1, 1]  # ################################################# [1, 1] / [2, 2] / [3, 3]
    scaling_mode = 'exp'  # ########################################################'tes' / 'exp'
    phenome_mode = 'car'  # ######################################################## 'car' / 'hexH' 'hexV'
    phenome_dims = [5, 5]
    phenome_format = ''
    for d in phenome_dims:
        phenome_format += str(d)
    zone_size = int(10 * np.ceil(np.sum([(cell_diam * phenome_scaling[p] * actuate_scale * phenome_dims[p]) ** 2
                                         for p in range(0, len(phenome_dims))]) ** 0.5 / 10))
    zone_size = 400  # ################################################# 120 / 200
    corner_offset = [int((zone_size - 2 - phenome_dims[p] * cell_diam * phenome_scaling[p]) / 2) + 1
                     for p in range(0, len(phenome_dims))]
    name2dims = {}
    for name in key2name.values():
        name2dims[name] = [cell_diam] * 2
    sim_dim = int(grid_dim * zone_size)

    # Write variables to file
    pd.DataFrame(np.array([['cell diam', cell_diam],
                           ['actuate scale', actuate_scale],
                           ['actuate period', actuate_period],
                           ['num active types', num_active_types],
                           ['max mcs', max_mcs],
                           ['work nodes', 1],  # ###################### work_nodes
                           ['sim dim', sim_dim]]),
                 None, ['Name', 'Value']).to_csv(variables_file, index=False)

    # Set up evolution parameters
    track_flag = False
    evolve_flag = True
    randomise_flag = True

    # rep_seed = ['5555433540354445544405540']

    # rep_seed = ['5550055444554455544455544']  # good, one cell is sub-optimal
    # rep_seed = ['4444444444544455545555555']
    # rep_seed = ['4000444044540455545555555']  # <<< crafted1

    # rep_seed = ['5350305554554405444054044']  # <<< first success
    # traj_label = 'trial1_'
    # rep_seed = ['0333550354405444554430044']  # <<< second attempt
    # rep_seed = ['5355535455540454404440444']  # <<< third attempt
    # traj_label = 'trial3_'
    # rep_seed = ['5350055444554445544455543']  # <<< trial 4
    # traj_label = 'trial4_'
    # rep_seed = ['4444444444544455545555555']
    # traj_label = 'custom_'
    rep_seed = ['4444444444544455555555555']
    traj_label = 'custom2r_'
    # rep_seed = ['4455544455444554445544555']
    # traj_label = 'custom2u_'
    # rep_seed = ['5555555555544454444444444']
    # traj_label = 'custom2l_'
    # rep_seed = ['5554455444554445544455544']
    # traj_label = 'custom2d_'
    # rep_seed = ['4000050000000000000000000']  # <<< two cell example
    # traj_label = 'twoCell_'

    # rep_seed = ['0350500554303403334050044']  # <<< combination
    # rep_seed = '0404035350430040340534003'  # <<<< random
    init_gen = 0
    number_of_generations = 20  # ########################################################### 10 / 50
    fitness_function = 'COM'
    survival_prop = 0.2  # 0.2
    offspring_prop = 0.5  # 0.4
    mutation_rate = 0.01
    align_genomes = False
    pop_size = 8  # ############################################## 32

    # Begin simulations
    record_file = phenome_mode + phenome_format + '_'
    if phenome_scaling[0] == 1 and phenome_scaling[1] == 1:
        record_file += 'sca'
    else:
        record_file += scaling_mode
    record_file += str(phenome_scaling[0]) + str(phenome_scaling[1]) + '_act' + str(num_actuation_cycles)
    if evolve_flag:
        record_file = 'ev_fit_' + record_file + ('_rep' + str(repeats) +
                                                 '_sur' + str(int(100 * survival_prop)) +
                                                 '_cro' + str(int(100 * offspring_prop)) +
                                                 '_mut' + str(int(100 * mutation_rate)) +
                                                 '_rot' + str(int(align_genomes)))
    if init_gen == 0:
        with open(join(root, record_file + '.csv'), 'w') as f_out:
            f_out.write('code,gen,fit\n')  # _MA,fit_AM
    trajectories = []
    for gen_num, sim in enumerate([simulation_file_path] * number_of_generations):
        if evolve_flag:
            gen_fname = gen_file_prefix + str(phenome_format) + '_' + str(gen_num + init_gen) + '.csv'
            generation = pd.read_csv(join(root_output_folder, gen_fname), header=None).to_numpy()
            if gen_num == 0 and randomise_flag:
                generation = np.zeros((pop_size, int(np.prod(phenome_dims))))
                generation = evolve_next_gen(np.arange(0, generation.shape[0], 1), generation, 0, 0,
                                             mutation_rate, key2name, np.zeros(generation.shape[0]), phenome_dims, phenome_mode)
        else:
            generation = get_genome(rep_seed)
        # fitnesses = [0 for _ in range(0, generation.shape[0])]
        fit_vecs = {}
        for g in range(0, generation.shape[0]):
            fit_vecs[g] = []
        for i in range(0, int(generation.shape[0] / phenomes_per_sim)):
            grid_arrangement = np.array([int(i * phenomes_per_sim + np.floor(j / repeats))
                                         for j in range(0, work_nodes)]).reshape((grid_dim, grid_dim))
            setup_cells.setup_init_cells(init_cells_file_path, generation, grid_arrangement, tuple(phenome_dims), zone_size,
                             corner_offset, key2name, name2dims, phenome_scaling, phenome_mode, scaling_mode)
            # key_matrix = np.reshape(generation[[i]][:], tuple(phenome_dims))  # [0], dims[1]))
            # if not key_matrix.sum() == 0:
            # place_cells(init_cells_file_path, key_matrix, key2name, name2dims)
            cc3d_caller = CC3DCaller(cc3d_sim_fname=sim,
                                     screenshot_output_frequency=0,
                                     output_dir=join(root_output_folder, abbrevs[project_name]),
                                     result_identifier_tag=i)
            ret_value = cc3d_caller.run()
            if track_flag:
                trajectories.append([ret_value['result']['trackX'], ret_value['result']['trackY']])
            fit_vec_grid = calculate_scores(ret_value, fitness_function, zone_size, grid_dim)
            # for index in set([a for a in grid_arrangement.reshape(-1)]):
            for x in range(0, fit_vec_grid.shape[0]):
                for y in range(0, fit_vec_grid.shape[1]):
                    fit_vecs[grid_arrangement[x][y]].append(fit_vec_grid[x][y])
            # fitnesses[i] = calculate_scores(ret_value, fitness_function, zone_size, grid_dim)
        fit_measures, fit_angles = measure_fitness(fit_vecs)
        # fitnesses = fit_measures.sum(1)
        fitnesses = fit_measures[:, 0] / num_actuation_cycles
        if evolve_flag:
            if not align_genomes:
                fit_angles = np.zeros(generation.shape[0])
            ordered_genomes = np.argsort(-1 * fitnesses)
            next_gen = evolve_next_gen(ordered_genomes, generation, survival_prop, offspring_prop,
                                       mutation_rate, key2name, fit_angles, phenome_dims, phenome_mode)
            new_gen_fname = gen_file_prefix + str(phenome_format) + '_' + str(gen_num + init_gen + 1) + '.csv'
            pd.DataFrame(next_gen).to_csv(join(root_output_folder, new_gen_fname),
                                          header=False, index=False)
        with open(join(root, record_file + '.csv'), 'a') as f_out:
            for i in range(0, generation.shape[0]):
                gene_code = ''
                for j in range(0, generation.shape[1]):
                    gene_code += str(int(generation[i, j]))
                f_out.write(gene_code + ',' + str(gen_num) + ',' +
                            str(fitnesses[i]) + '\n')  # + ',' + str(fit_measures[i, 1])
    if track_flag:
        save_trajectories(trajectories, root, traj_label, record_file + '.csv')
        # plot_trajectories(trajectories, zone_size)
    winsound.Beep(440, 2000)


if __name__ == '__main__':
    main()
