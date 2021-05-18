from cc3d.CompuCellSetup.CC3DCaller import CC3DCaller

import setup_cells

import os
import json
import winsound
import numpy as np
import pandas as pd


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


def evolve_next_gen(order, current_gen, survive, offspring, mutate, angles, config):
    population_size = current_gen.shape[0]
    survive_num = int(np.ceil(survive * population_size))
    offspring_num = int(np.floor(offspring * population_size))
    new_gen = np.zeros(current_gen.shape)
    for i in range(0, population_size):
        if i < survive_num:
            new_gen[i, :] = current_gen[order[i], :]
        elif i < survive_num + offspring_num and survive_num > 1:
            parent_indices = np.random.choice(survive_num, 2, replace=False)
            new_gen[i, :] = apply_crossover(current_gen[order[parent_indices], :],
                                            angles[order[parent_indices]], config)
        else:
            new_gen[i, :] = np.random.choice([0] + list(config['key name signal'].keys()), new_gen.shape[1])
        # mutate_mask = np.sign(np.sign(np.random.random(new_gen.shape) - mutate) + 1)
    new_mut_gen = apply_mutation(new_gen, mutate, survive_num, config)
    return new_mut_gen


def apply_crossover(parents, parent_angles, config):
    new_genome = np.zeros((1, parents.shape[1]))
    angle_diff = parent_angles[1] - parent_angles[0]
    symmetry = 0
    multiplier = 0
    if config['phenotype layout'] == 'car':
        symmetry = 4
        multiplier = 1
    elif 'hex' in config['phenotype layout']:
        symmetry = 2
        multiplier = 2
    num_rotations = int(np.mod(np.round(np.mod(angle_diff, 2 * np.pi) / (2 * np.pi / symmetry), 0), symmetry))
    aligned_parent = permute_genome(parents[1, :].reshape(config['phenotype dimensions']),
                                    int(multiplier * num_rotations))
    parents[1, :] = aligned_parent.reshape(1, -1)
    for j in range(0, parents.shape[1]):
        new_genome[0, j] = parents[np.random.choice([0, 1], 1)[0], j]
    return new_genome


def permute_genome(genome, rotations):
    if rotations == 0:
        return genome
    aligned_genome = np.zeros(genome.shape)
    s = genome.shape[0]

    def rot(x):
        return [x[2] - 1 - x[1], x[0], x[2]]
    # r = lambda x: [x[2] - 1 - x[1], x[0], x[2]]
    for i in range(0, s):
        for j in range(0, s):
            r_map = [i, j, s]
            for _ in range(0, rotations):
                r_map = rot(r_map)
            aligned_genome[r_map[0], r_map[1]] = genome[i, j]
    return aligned_genome


def apply_mutation(gen, mut_rate, skip, config):
    if mut_rate == 0:
        return gen
    mut_gen = 1 * gen
    for i in range(0, mut_gen.shape[0]):
        if i >= skip:
            for j in range(0, mut_gen.shape[1]):
                if np.random.random() < mut_rate:
                    mut_gen[i, j] = np.random.choice(list((set(config['key name signal'].keys()) | {0}) - {gen[i, j]}))
    return mut_gen


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


def measure_fitness(end_vector_dict_list):
    measures = np.zeros((len(end_vector_dict_list.keys()), 1))  # 2
    angles = np.zeros(len(end_vector_dict_list.keys()))
    for i in end_vector_dict_list.keys():
        vectors = np.zeros((len(end_vector_dict_list[i]), 2))
        for j, v in enumerate(end_vector_dict_list[i]):
            vectors[j, :] = v
        measures[i, 0] = ((vectors.sum(0) / vectors.shape[0]) ** 2).sum() ** 0.5  # magnitude of average vec
        # measures[i, 1] = ((vectors ** 2).sum(1) ** 0.5).sum() / vectors.shape[0]  # average of vec magnitudes
        expected_vec = vectors.sum(0) / vectors.shape[0]
        angles[i] = np.arctan2(expected_vec[1], expected_vec[0])
    return measures, angles


def save_trajectories(trajectory_list, root, trajectory_lab, file_name):
    x_data = np.zeros((len(trajectory_list[0][0]), len(trajectory_list)))
    y_data = np.zeros((len(trajectory_list[0][0]), len(trajectory_list)))
    for t, traj in enumerate(trajectory_list):
        for s in range(0, len(traj[0])):
            x_data[s, t] = traj[0][s]
            y_data[s, t] = traj[1][s]
    pd.DataFrame(x_data).to_csv(root + trajectory_lab + 'xData_' + file_name, header=False, index=False)
    pd.DataFrame(y_data).to_csv(root + trajectory_lab + 'yData_' + file_name, header=False, index=False)


def main(config):
    dep_config = setup_cells.get_dep_configs(config)
    root = os.getcwd()
    simulation_file_path = os.path.join(os.path.join(root, config['project name']), config['project name'] + '.cc3d')
    sim_dump_directory = os.path.join(os.path.join(root, config['project name']), config['sim dump subdir'])
    init_cells_file_path = os.path.join(os.path.join(root, config['project name']), 'init_cell_field.piff')
    generations_folder = os.path.join(root, config['generation prefix']['folder'])

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
    survival_prop = 0.2  # 0.2
    offspring_prop = 0.5  # 0.4
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
        record_file = 'ev_fit_' + record_file + ('_rep' + str(config['repeats']) +
                                                 '_sur' + str(int(100 * config['survive proportion'])) +
                                                 '_cro' + str(int(100 * config['offspring proportion'])) +
                                                 '_mut' + str(int(100 * config['mutation rate'])) +
                                                 '_rot' + str(int(align_genomes)))
    if config['initial generation'] == 0:
        with open(os.path.join(os.path.join(root, config['generation prefix']['folder']),
                               record_file + '.csv'), 'w') as f_out:
            f_out.write('code,gen,fit\n')  # _MA,fit_AM
    trajectories = []
    # Begin simulations
    winsound.Beep(440, 2000)
    for gen_num, sim in enumerate([simulation_file_path] * config['num generations']):
        if evolve_flag:
            gen_filename = (config['generation prefix']['file'] + phenotype_dimension_string +
                            str(gen_num + config['initial generation']) + '.csv')
            if gen_num == 0 and randomise_flag:
                generation = np.zeros((dep_config['population size'], int(np.prod(config['phenotype dimensions']))))
                generation = evolve_next_gen(np.arange(0, generation.shape[0], 1), generation, 0, 0,
                                             config['mutation rate'], np.zeros(generation.shape[0]), config)
            else:
                generation = pd.read_csv(os.path.join(generations_folder, gen_filename), header=None).to_numpy()
        else:
            generation = get_genome(rep_seed)
        end_vectors = {}
        for g in range(0, generation.shape[0]):
            end_vectors[g] = []
        for i in range(0, int(generation.shape[0] / config['phenotypes per simulation'])):
            zone_genome_index = np.array([int(i * config['phenotypes per simulation'] +
                                              np.floor(j / config['repeats'])) for j in
                                          range(0, dep_config['num zones'])]).reshape(dep_config['grid dimensions'])
            setup_cells.initialise_piff(init_cells_file_path, generation, zone_genome_index, config)
            cc3d_caller = CC3DCaller(cc3d_sim_fname=sim,
                                     screenshot_output_frequency=0,
                                     output_dir=sim_dump_directory,
                                     result_identifier_tag=i)
            ret_value = cc3d_caller.run()
            if track_flag:
                trajectories.append([ret_value['result']['trackX'], ret_value['result']['trackY']])
            zone_end_vectors = calculate_zone_vectors(ret_value, dep_config['zone size'], dep_config['grid dimensions'])
            for x in range(0, zone_end_vectors.shape[0]):
                for y in range(0, zone_end_vectors.shape[1]):
                    end_vectors[zone_genome_index[x][y]].append(zone_end_vectors[x][y])
        fitness, predicted_directions = measure_fitness(end_vectors)
        fitness_per_cycle = fitness[:, 0] / config['num cycles']
        if evolve_flag:
            if not align_genomes:
                predicted_directions = np.zeros(generation.shape[0])
            ordered_genomes = np.argsort(-1 * fitness_per_cycle)
            next_gen = evolve_next_gen(ordered_genomes, generation, survival_prop, offspring_prop,
                                       config['mutation rate'], predicted_directions, config)
            new_gen_filename = (config['generation prefix']['file'] + phenotype_dimension_string +
                                str(gen_num + config['initial generation'] + 1) + '.csv')
            pd.DataFrame(next_gen).to_csv(os.path.join(generations_folder, new_gen_filename),
                                          header=False, index=False)
        with open(os.path.join(os.path.join(root, config['generation prefix']['folder']),
                               record_file + '.csv'), 'a') as f_out:
            for i in range(0, generation.shape[0]):
                gene_code = ''
                for j in range(0, generation.shape[1]):
                    gene_code += str(int(generation[i, j]))
                f_out.write(gene_code + ',' + str(gen_num) + ',' +
                            str(fitness_per_cycle[i]) + '\n')  # + ',' + str(fitness[i, 1])
    if track_flag:
        save_trajectories(trajectories, root, traj_label, record_file + '.csv')
    winsound.Beep(440, 2000)


if __name__ == '__main__':
    with open('config.json') as config_file:
        raw_config = json.load(config_file)
    main(raw_config)
