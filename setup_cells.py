import numpy as np


def get_dep_configs(config):
    zone_size = [config['cell diameter'] * config['actuate scale'] * config['phenotype dimensions'][i] *
                 config['phenotype scale'][i] * config['zone phenotype ratio']
                 for i in range(0, len(config['phenotype dimensions']))]
    corner_offset = [int((zone_size[i] - 2 - config['phenotype dimensions'][i] *
                          config['cell diameter'] * config['phenotype scale'][i]) / 2) + 1
                     for i in range(0, len(config['phenotype dimensions']))]
    num_zones = config['grid dimensions'][0] * config['grid dimensions'][1]
    sim_dims = [int(zone_size[i] * config['grid dimensions'][i]) for i in range(0, len(zone_size))]
    cycle_dur = config['actuate period'] * max([len(config['key name signal'][k]['signal'])
                                                for k in config['key name signal'].keys()])
    max_mcs = int(config['settle delay'] + cycle_dur * config['num cycles'])
    return {'zone size': zone_size,
            'corner offset': corner_offset,
            'num zones': num_zones,
            'sim dimensions': sim_dims,
            'cycle duration': cycle_dur,
            'max mcs': max_mcs}


def initialise_piff(piff_file_path, gen, arrangement, config):
    dep_config = get_dep_configs(config)
    index = 0
    with open(piff_file_path, 'w') as piff_file:
        # index = write_hollow_box_to_piff(piff_file, index, 'Wall', [0, 0], dep_config['zone size'])
        for x in range(0, arrangement.shape[0]):
            for y in range(0, arrangement.shape[1]):
                box_start = [x * config['zone size'][0], y * config['zone size'][1]]
                index = write_hollow_box_to_piff(piff_file, index, 'Wall', box_start, dep_config['zone size'])
                key_matrix = np.reshape(gen[[arrangement[x][y]], :], config['phenotype dimensions'])  # ][
                start = [x * dep_config['zone size'][0] + dep_config['corner offset'][0],
                         y * dep_config['zone size'][1] + dep_config['corner offset'][1]]
                index = place_cells(piff_file, index, key_matrix, start, config)


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


def place_cells(file, ind, key_mat, offset, config):
    shunt = None
    if config['phenotype layout'] == 'car':
        shunt = [0, 0]
    elif config['phenotype layout'] == 'hexH':
        shunt = [1, 0]
    elif config['phenotype layout'] == 'hexV':
        shunt = [0, 1]
    expand = [1, 1]
    tessellate = [1, 1]
    if config['scaling mode'] == 'exp':
        expand = config['phenotype scale']
    elif config['scaling mode'] == 'tes':
        tessellate = [key_mat.shape[0], key_mat.shape[1]]
    name_dims_dict = {}
    for key in config['key name signal'].keys():
        name_dims_dict[config['key name signal'][key]['name']] = [config['cell diameter']] * 2
    new_ind = ind
    for x in np.arange(0, key_mat.shape[0]):
        for y in np.arange(0, key_mat.shape[1]):
            key = key_mat[x][y]
            if not key == 0:
                cell_type = config['key name signal'][key]['name']
                for p1 in range(0, config['phenotype_scale'][0]):
                    xx = expand[0] * x + p1 * tessellate[0]  # # # # # # # # # # # # # # # # # # # # # # #
                    for p2 in range(0, config['phenotype_scale'][1]):
                        yy = expand[1] * y + p2 * tessellate[1]
                        x_coordinates = name_dims_dict[cell_type][0] * (xx * np.ones(2) + np.array([0, 1])) + offset[0]
                        y_coordinates = name_dims_dict[cell_type][1] * (yy * np.ones(2) + np.array([0, 1])) + offset[1]
                        x_coordinates += shunt[0] * np.mod(y, 2) * int(config['phenotype scale'][0] *
                                                                       name_dims_dict[cell_type][0] / 2)
                        y_coordinates += shunt[1] * np.mod(x, 2) * int(config['phenotype scale'][1] *
                                                                       name_dims_dict[cell_type][1] / 2)
                        line = (str(new_ind) + ' ' + cell_type +
                                ' ' + str(int(x_coordinates[0])) + ' ' + str(int(x_coordinates[1])) +
                                ' ' + str(int(y_coordinates[0])) + ' ' + str(int(y_coordinates[1])) + ' 0 0\n')
                        file.write(line)
                        new_ind += 1
    return new_ind
