import numpy as np


def setup_init_cells(file_path, gen, arrangement, phen_dim, zone, offsets, key_name, name_dims,
                     phenome_scaling, phenome_mode, scaling_mode):
    index = 0
    with open(file_path, 'w') as f_out:
        for x in range(0, arrangement.shape[0]):
            for y in range(0, arrangement.shape[1]):
                index = hollow_box_2d(f_out, index, 'Wall', [x * zone, y * zone], [zone, zone])
                key_matrix = np.reshape(gen[[arrangement[x][y]], :], phen_dim)  # ][
                start = [x * zone + offsets[0], y * zone + offsets[1]]
                index = place_cells(f_out, index, key_matrix, key_name, name_dims, start,
                                    phenome_scaling, phenome_mode, scaling_mode)


def hollow_box_2d(file, ind, pixel_type, corner, dims):
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
                    xx = expand[0] * x + p1 * tesselate[0]  # # # # # # # ##  # # # # # # # # # # #  # # #
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