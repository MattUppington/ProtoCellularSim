from cc3d.CompuCellSetup.CC3DCaller import CC3DCaller

import setup_cells

import os
import json
import winsound
import numpy as np
import pandas as pd


def setup_init_cells(file_path, gen, arrangement, phen_dim, zone, offsets, key_name, name_dims,
                     phenome_scaling, phenome_mode, scaling_mode):
    index = 0
    with open(file_path, 'w') as f_out:
        index = write_hollow_box_to_piff(f_out, index, 'Wall', [0, 0], [zone, zone])
        for x in range(0, arrangement.shape[0]):
            for y in range(0, arrangement.shape[1]):
                # index = hollow_box_2d(f_out, index, 'Wall', [x * zone, y * zone], [zone, zone])
                key_matrix = np.reshape(gen[[arrangement[x][y]], :], phen_dim)  # ][
                start = [x * zone + offsets[0], y * zone + offsets[1]]
                index = place_cells(f_out, index, key_matrix, key_name, name_dims, start,
                                    phenome_scaling, phenome_mode, scaling_mode)


def write_hollow_box_to_piff(file, ind, pixel_type, corner, dims):
    new_ind = ind
    shifts = np.array([[0, 0, 0, 0],
                       [0, 0, 1, 1],
                       [1, 1, 0, 0],
                       [1, 1, 1, 1]])
    edges = np.array([[0, 1, 0, 0],
                      [0, 0, 0, 1],
                      [1, 1, 0, 1],
                      [0, 1, 1, 1]])
    edges[:, [0, 1]] *= dims[0] - 2
    edges[:, [2, 3]] *= dims[1] - 2
    for side in range(0, 4):
        line = str(new_ind) + ' ' + pixel_type
        for elem in range(0, 4):
            axis = int(np.floor(elem / 2))
            line += ' ' + str(corner[axis] + shifts[side, elem] + edges[side, elem])
        line += ' 0 0\n'
        file.write(line)
        new_ind += 1
    return new_ind


def place_cells(file, ind, key_mat, key_name_dict, name_dims_dict, offset, phenome_scaling,
                phenome_mode, scaling_mode):
    # with open(fname, 'w') as f_out
    #     index = 0
    #     offset = 2
    shunt = None
    if phenome_mode == 'car':
        shunt = [0, 0]
    elif phenome_mode == 'hexH':
        shunt = [1, 0]
    elif phenome_mode == 'hexV':
        shunt = [0, 1]
    expand = [1, 1]
    tesselate = [1, 1]
    if scaling_mode == 'exp':
        expand = phenome_scaling
    elif scaling_mode == 'tes':
        tesselate = [key_mat.shape[0], key_mat.shape[1]]
    new_ind = ind
    for x in np.arange(0, key_mat.shape[0]):
        for y in np.arange(0, key_mat.shape[1]):
            key = key_mat[x][y]
            if not key == 0:
                cell_type = key_name_dict[key]
                for p1 in range(0, phenome_scaling[0]):
                    xx = expand[0] * x + p1 * tesselate[0]  # # # # # # # # # # # # # # # # # # # # # # #
                    for p2 in range(0, phenome_scaling[1]):
                        yy = expand[1] * y + p2 * tesselate[1]
                        x_coords = name_dims_dict[cell_type][0] * (xx * np.ones(2) + np.array([0, 1])) + offset[0]
                        y_coords = name_dims_dict[cell_type][1] * (yy * np.ones(2) + np.array([0, 1])) + offset[1]
                        x_coords += shunt[0] * np.mod(y, 2) * int(phenome_scaling[0] *
                                                               name_dims_dict[cell_type][0] / 2)
                        y_coords += shunt[1] * np.mod(x, 2) * int(phenome_scaling[1] *
                                                               name_dims_dict[cell_type][1] / 2)
                        line = (str(new_ind) + ' ' + cell_type +
                                ' ' + str(int(x_coords[0])) + ' ' + str(int(x_coords[1])) +
                                ' ' + str(int(y_coords[0])) + ' ' + str(int(y_coords[1])) + ' 0 0\n')
                        file.write(line)
                        new_ind += 1
    return new_ind


def calculate_zone_vectors(return_object, zone, grid):
    starts = return_object['result']['initial data']
    diffs = return_object['result']['final data'] - starts
    vec_grid = np.zeros((grid[0], grid[1], 2))
    for x in range(0, grid[0]):
        for y in range(0, grid[1]):
            mask = [i for i in range(0, diffs.shape[0]) if x * zone[0] < starts[i][0] < (x + 1) * zone[0] and
                    y * zone[1] < starts[i][1] < (y + 1) * zone[1]]
            if len(mask) == 0:
                vec_grid[x][y][:] = 0
            else:
                cell_vectors = diffs[mask][:]
                vec_grid[x][y][:] = np.sum(cell_vectors / len(mask), 0)
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


def evolve_next_gen(order, current_gen, survive, offspring, mutate, key_dict, angles, phen_dims, config):
    # population_size = current_gen.shape[0]
    survive_num = int(np.ceil(survive * config['population size']))
    offspring_num = int(np.floor(offspring * config['population size']))
    new_gen = np.zeros(current_gen.shape)
    for i in range(0, config['population size']):
        if i < survive_num:
            new_gen[i, :] = current_gen[order[i], :]
        elif i < survive_num + offspring_num:
            parent_inds = np.random.choice(survive_num, 2, replace=False)
            new_gen[i, :] = apply_crossover(current_gen[order[parent_inds], :],
                                            angles[order[parent_inds]], phen_dims)
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


def main(config):
    dep_config = setup_cells.get_dep_configs(config)
    root = os.getcwd()
    simulation_file_path = os.path.join(os.path.join(root, config['project name']), config['project name'] + '.cc3d')
    sim_output_directory = os.path.join(os.path.join(root, config['project name']), config['sim output subdir'])
    init_cells_file_path = os.path.join(os.path.join(root, config['project name']), 'init_cell_field.piff')
    generations_folder = os.path.join(root, config['generation prefix']['folder'])

    root_output_folder = os.path.join(root, 'Results')

    # Set up parallel processing / multi-threading parameters
    # cpus = 4
    # threads_per_cpu = cpus
    # work_nodes = cpus * threads_per_cpu
    work_nodes = 1  # ############################################################################ 1 / 16
    # grid_dim = [4, 4]
    phenomes_per_sim = 1  # ###################################################################### 1 / 4
    repeats = int(work_nodes / phenomes_per_sim)

    # Set up simulation parameters
    num_actuation_cycles = 50  # ################################################################### 10 / 100
    key2name = {3: 'ProtocellPassive', 4: 'ProtocellActiveA', 5: 'ProtocellActiveB'}
    phenome_dims = [5, 5]


    # Write variables to file
    # pd.DataFrame(np.array([['cell diam', cell_diam],
    #                        ['actuate scale', actuate_scale],
    #                        ['actuate period', actuate_period],
    #                        ['num active types', num_active_types],
    #                        ['max mcs', max_mcs],
    #                        ['work nodes', 16],  # ###################### work_nodes
    #                        ['sim dim', sim_dim]]),
    #              None, ['Name', 'Value']).to_csv(variables_file, index=False)

    # Set up evolution parameters
    track_flag = True
    evolve_flag = False
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
    survival_prop = 0.2  # 0.2
    offspring_prop = 0.5  # 0.4
    mutation_rate = 0.01
    align_genomes = False

    phenotype_dimension_string = (str(config['phenotype dimensions'][0]).zfill(2) +
                                  str(config['phenotype dimensions'][1]).zfill(2) + '_')
    record_file = config['phenotype layout'] + phenotype_dimension_string
    if config['phenotype scale'][0] == 1 and config['phenotype scale'][1] == 1:
        record_file += 'sca'
    else:
        record_file += config['scaling mode']
    record_file += (str(config['phenotype scale'][0]).zfill(2) + str(config['phenotype scale'][1]).zfill(2) +
                    '_cyc' + str(config['num cycles']))
    if evolve_flag:
        record_file = 'ev_fit_' + record_file + ('_rep' + str(repeats) +
                                                 '_sur' + str(int(100 * survival_prop)) +
                                                 '_cro' + str(int(100 * offspring_prop)) +
                                                 '_mut' + str(int(100 * mutation_rate)) +
                                                 '_rot' + str(int(align_genomes)))
    if config['initial generation'] == 0:
        with open(os.path.join(root, record_file + '.csv'), 'w') as f_out:
            f_out.write('code,gen,fit\n')  # _MA,fit_AM
    trajectories = []
    # Begin simulations
    winsound.Beep(440, 2000)
    for gen_num, sim in enumerate([simulation_file_path] * config['num generations']):
        if evolve_flag:
            gen_filename = (config['generation prefix']['file'] + phenotype_dimension_string +
                            str(gen_num + config['initial generation']) + '.csv')
            generation = pd.read_csv(os.path.join(generations_folder, gen_filename), header=None).to_numpy()
            if gen_num == 0 and randomise_flag:
                generation = np.zeros((config['population size'], int(np.prod(config['phenotype dimensions']))))
                generation = evolve_next_gen(np.arange(0, generation.shape[0], 1), generation, 0, 0,
                                             mutation_rate, key2name, np.zeros(generation.shape[0]), phenome_dims, config)
        else:
            generation = get_genome(rep_seed)
        # fitnesses = [0 for _ in range(0, generation.shape[0])]
        end_vectors = {}
        for g in range(0, generation.shape[0]):
            end_vectors[g] = []
        for i in range(0, int(generation.shape[0] / phenomes_per_sim)):
            zone_genome_index = np.array([int(i * phenomes_per_sim + np.floor(j / repeats))
                                         for j in range(0, dep_config['num zones'])]).reshape(config['grid dimensions'])
            # setup_init_cells(init_cells_file_path, generation, zone_genome_index, tuple(phenome_dims),
            #                  dep_config['zone size'], dep_config['corner offset'], key2name, name2dims,
            #                  phenome_scaling, phenome_mode, scaling_mode)
            setup_cells.initialise_piff(init_cells_file_path, generation, zone_genome_index, config)
            cc3d_caller = CC3DCaller(cc3d_sim_fname=sim,
                                     screenshot_output_frequency=0,
                                     output_dir=sim_output_directory,
                                     result_identifier_tag=i)
            ret_value = cc3d_caller.run()
            if track_flag:
                trajectories.append([ret_value['result']['trackX'], ret_value['result']['trackY']])
            zone_end_vectors = calculate_zone_vectors(ret_value, dep_config['zone size'], config['grid dimensions'])
            # for index in set([a for a in zone_genome_index.reshape(-1)]):
            for x in range(0, zone_end_vectors.shape[0]):
                for y in range(0, zone_end_vectors.shape[1]):
                    end_vectors[zone_genome_index[x][y]].append(zone_end_vectors[x][y])
            # fitnesses[i] = calculate_scores(ret_value, fitness_function, zone_size, grid_dim)
        fit_measures, fit_angles = measure_fitness(end_vectors)
        # fitnesses = fit_measures.sum(1)
        fitnesses = fit_measures[:, 0] / config['num cycles']
        if evolve_flag:
            if not align_genomes:
                fit_angles = np.zeros(generation.shape[0])
            ordered_genomes = np.argsort(-1 * fitnesses)
            next_gen = evolve_next_gen(ordered_genomes, generation, survival_prop, offspring_prop,
                                       mutation_rate, key2name, fit_angles, phenome_dims, config)
            new_gen_filename = (config['generation file prefix'] + phenotype_dimension_string +
                                str(gen_num + config['initial generation'] + 1) + '.csv')
            pd.DataFrame(next_gen).to_csv(os.path.join(root_output_folder, new_gen_filename),
                                          header=False, index=False)
        with open(os.path.join(root, record_file + '.csv'), 'a') as f_out:
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
    with open('config.json') as config_file:
        raw_config = json.load(config_file)
    main(raw_config)
