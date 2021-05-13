import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import math
import string


def convert_colour(hex_string):
    hex = [str(n) for n in range(0, 10)] + [l for l in string.ascii_lowercase[:6]]
    hex2dec = {}
    for x, h in enumerate(hex):
        hex2dec[h] = x
    rgb_col = [0, 0, 0]
    for c in range(0, len(rgb_col)):
        rgb_col[c] = (hex2dec[hex_string[2 * c + 1]] + 16 * hex2dec[hex_string[2 * c]]) / 256
    return rgb_col


def movement(a):
    r_min = 1
    r_max = 1.25992104989
    r = (1 - a) * r_min + a * r_max
    return (r_max - r_min) / 2 + r * math.log((r + r_min) / (r + r_max), math.e)


def actuation_cycles(e_col, l_col, p_col):
    plot_settings = {'Passive': {'sym': 's-', 'col': convert_colour(p_col)},
                     'Late': {'sym': 'o-', 'col': convert_colour(l_col)},
                     'Early': {'sym': '^-', 'col': convert_colour(e_col)}}
    labels = ['Early', 'Late', 'Passive']
    mcs = [0, 50, 100, 150, 200, 250, 300, 350, 400]
    a = 1.25992104989
    b = (a + 1) / 2
    radius_p = [1, 1, 1, 1, 1, 1, 1, 1, 1]
    radius_1 = [1, b, a, a, a, b, 1, 1, 1]
    radius_2 = [1, 1, 1, b, a, a, a, b, 1]

    fig, axs = plt.subplots(2, 1)
    axs[0].plot(mcs, radius_p, plot_settings['Passive']['sym'],
                label='Passive', color=plot_settings['Passive']['col'])
    axs[0].plot(mcs, radius_2, plot_settings['Late']['sym'],
                label='Late', color=plot_settings['Late']['col'])
    axs[0].plot(mcs, radius_1, plot_settings['Early']['sym'],
                label='Early', color=plot_settings['Early']['col'])
    axs[0].set(xlabel='MCS', ylabel='Radius')
    axs[0].set_title('Variation of radii over actuation cycles')
    axs[0].legend()


    motion = [0]
    motion.append(motion[-1] + movement(0))
    motion.append(motion[-1] - movement(1))
    motion.append(motion[-1] - movement(1))
    motion.append(motion[-1] + movement(0))
    axs[1].plot([0, 100, 200, 300, 400], motion, 'k^-')
    axs[1].set(xlabel='MCS', ylabel='$x_{COM}$')
    # axs[1].set_title('Graph showing motion of centre of mass for two agent system')

    # for ax in axs.flat:
    #     ax.label_outer()
    plt.show()


def show_fitness_per_gen(data_list, survive_prop):
    colours = {0: convert_colour('929e00'), 1: convert_colour('ff86e8'), 2: convert_colour('3c75a2'), 3: 'm'}
    for i, data in enumerate(data_list):
        pop_size = len(data[data['gen'] == data['gen'].unique()[0]])
        survive_num = int(np.ceil(survive_prop * pop_size))
        g_nums = sorted(list(data['gen'].unique()) * survive_num)
        best_fits = []
        averages = []
        errors = np.zeros((2, len(data['gen'].unique())))
        for g in data['gen'].unique():
            generation = data[data['gen'] == g]
            fitnesses = [fit for fit in generation['fit'].sort_values(0, False)[0:survive_num]]
            best_fits += fitnesses
            averages.append(sum(fitnesses) / len(fitnesses))
            errors[0, g] = averages[-1] - np.min(fitnesses)
            errors[1, g] = np.max(fitnesses) - averages[-1]
        # plt.plot(g_nums, best_fits, colours[i] + '.', linewidth=0.5, label='Trial ' + str(i+1))
        plt.errorbar(np.arange(i / 5, len(averages) + i / 5, 1), averages, color=colours[i], yerr=errors,
                     linewidth=2, elinewidth=1, label='Example trial ' + str(i+1)) #
        # for b in range(0, survive_num):
        #     data['best'] = data['gen'].diff(periods=1)
        # data['best'] = data['gen'].diff(periods=survive_num)
        # data.iloc[0:survive_num, data.columns.get_loc('best')] = 1
        # plot_values = [data['fit'][data['gen'] == i & data['best'] == 1] for i in data['gen'].unique()]
        # plt.plot(data['gen'][data['best'] == 1], data['fit'][data['best'] == 1], colours[i])
        # plt.boxplot(plot_values)
    plt.xlabel('Generation')
    plt.ylabel('Fitness')
    # plt.title('Distributions of best fitness scores per generation')
    plt.axis([0, 50, 0, 1.2])
    plt.legend(loc='lower right')
    plt.show()


def my_basic_edit_dist(string1, string2):
    dist = None
    if len(string1) == len(string2):
        dist = 0
        for i in range(0, len(string1)):
            if not string1[i] == string2[i]:
                dist += 1
    return dist


def show_diversity(data):
    x_values = []
    y_values = []
    for gen_num in data['gen'].unique():
        print('Calculating Diversity: generation = ' + str(gen_num))
        codes = data['code'][data['gen'] == gen_num].to_numpy()
        for i in range(0, len(codes) - 1):
            for j in range(i + 1, len(codes)):
                x_values.append(gen_num)
                y_values.append(my_basic_edit_dist(codes[i], codes[j]))
    heatmap, xedges, yedges = np.histogram2d(x_values, y_values, bins=(11, 25))
    # plt.plot(x_values, y_values, 'kx')
    extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]
    plt.imshow(heatmap.T, extent=extent, origin='lower')
    plt.show()


def compare_fitness(data):
    codes = data['code'].unique()
    fig, axs = plt.subplots(len(codes))
    fig.suptitle('Vertically stacked subplots')
    limits = [np.floor(data['fit'].min()), np.ceil(data['fit'].max())]
    my_bins = np.arange(limits[0] - 0.5, limits[1] + 0.5, 1)
    for c in range(0, len(codes)):
        axs[c].hist(data['fit'][data['code'] == codes[c]], bins=my_bins)
    # for x in data_df['code'].unique():
        # blah = data_df['gen'][data_df['code'] == x]
    plt.show()


# cumulative number of appearances of genes across generations
# distribution of inter edit distances
# implement evolution switch
# remove gen_num and i index from root output directory to eliminate redudant directories in Results folder
# include evolution flag in generation and fitness file names for differentiation / to avoid premature overwriting

def state_abundance(data_list, survive_prop, colour_list):
    counter_list = []
    fig, axs = plt.subplots(1, 1)
    colours = {'4': colour_list[0], '5': colour_list[1], '3': colour_list[2], '0': colour_list[3]}
    labels = {'0': 'No cell', '3': 'Passive cell', '4': 'Active cell (Early)', '5': 'Active cell (Late)'}
    symbols = {'0': 'x-', '3': 's-', '4': '^-', '5': 'o-'}
    for i, data in enumerate(data_list):
        pop_size = len(data[data['gen'] == data['gen'].unique()[0]])
        counter_list.append({'0': np.zeros(data['gen'].max() + 1),
                             '3': np.zeros(data['gen'].max() + 1),
                             '4': np.zeros(data['gen'].max() + 1),
                             '5': np.zeros(data['gen'].max() + 1)})
        for gen_num in data['gen'].unique():
            print('data set ' + str(i + 1) + ' generation ' + str(gen_num))
            # codes = data['code'][data['gen'] == gen_num].to_numpy()
            generation = data[data['gen'] == gen_num]
            survive_num = int(np.ceil(survive_prop * pop_size))
            # counter_list[i] = {'0': 0, '3': 0, '4': 0, '5': 0}
            for genome in generation['code'][0:survive_num]:
                for index in genome:
                    counter_list[i][index][gen_num] += 1
        for ind in counter_list[i].keys():
            axs.plot(data['gen'].unique(), counter_list[i][ind] / (25 * survive_num), symbols[ind], color=colours[ind], label=labels[ind])
    axs.set(xlabel='Generation', ylabel='Proportion of Genome')
    axs.set_title('Average Abundance of Cell States per Generation')
    axs.axis([0, 50, 0, 0.5])
    axs.legend()
    plt.show()


def neighbour_distribution(data_list, phenome_dims, phenome_mode):
    norm_x_coords = np.tile(np.expand_dims(np.arange(0, phenome_dims[0]), 1), [1, phenome_dims[1]])
    norm_y_coords = np.tile(np.expand_dims(np.arange(0, phenome_dims[1]), 0), [phenome_dims[0], 1])
    shunt = None
    if phenome_mode == 'car':
        shunt = [0, 0]
    elif phenome_mode == 'hexH':
        shunt = [0.5, 0]
    elif phenome_mode == 'hexV':
        shunt = [0, 0.5]
    x_coords = norm_x_coords + shunt[0] * np.mod(norm_y_coords, 2)
    y_coords = norm_y_coords + shunt[1] * np.mod(norm_x_coords, 2)

    inds = {0: 0, 3: 1, 4: 2, 5: 3}
    meta_neighbour_mat = np.zeros((len(inds), len(inds)))

    hist_list = {}
    for i in inds.keys():
        for j in inds.keys():
            hist_list[str(i) + str(j)] = []
    final_gen = None
    for data in data_list:
        final_gen = data['code'][data['gen'] == data['gen'].max()]
        offspring_num = int((0.2 + 0.4) * len(final_gen))
        for genome in final_gen[0: offspring_num]:
            phenome = np.zeros((1, len(genome)))
            for slot, index in enumerate(genome):
                phenome[0, slot] = int(index)
            phenome = phenome.reshape(phenome_dims)
            neigh_mat = np.zeros((len(inds), len(inds)))
            counts = np.zeros((len(inds), len(inds)))
            for i1 in range(phenome.shape[0]):
                for j1 in range(phenome.shape[1]):
                    for i2 in range(phenome.shape[0]):
                        for j2 in range(phenome.shape[1]):
                            dist = p_dist(np.array([x_coords[i1, j1], y_coords[i1, j1]]),
                                          np.array([x_coords[i2, j2], y_coords[i2, j2]]), 2)
                            if dist != 0:
                                neigh_mat[inds[phenome[i1, j1]], inds[phenome[i2, j2]]] += dist
                                counts[inds[phenome[i1, j1]], inds[phenome[i2, j2]]] += 1
                                hist_list[str(int(phenome[i1, j1])) + str(int(phenome[i2, j2]))].append(dist)
            meta_neighbour_mat += neigh_mat / counts
    for i in range(0, len(inds)):
        for j in range(i + 1, len(inds)):
            meta_neighbour_mat[i, j] = 0
    labels = ['Empty', 'Passive', 'Early', 'Late']
    plot_heatmap(meta_neighbour_mat / len(final_gen) / len(data_list), labels)
    # fig2, axs = plt.subplots(len(inds), len(inds))
    # for i in inds.keys():
    #     for j in inds.keys():
    #        axs[inds[i], inds[j]].hist(hist_list[str(i) + str(j)], bins=np.arange(0, 8, 0.2))
    plt.show()


def p_dist(x1, x2, p):
    return np.sum((x2 - x1) ** p) ** (1/p)


def plot_heatmap(matrix, labels):
    fig1, ax = plt.subplots()
    ax.imshow(matrix)
    ax.set_xticks(np.arange(0, matrix.shape[0]))
    ax.set_yticks(np.arange(0, matrix.shape[1]))

    ax.set_xticklabels(labels)
    ax.set_yticklabels(labels)
    plt.setp(ax.get_xticklabels(), rotation=45, ha='right', rotation_mode='anchor')
    for j in range(0, matrix.shape[1]):
        for i in range(j, matrix.shape[0]):
            _ = ax.text(j, i, np.round(matrix[i, j], 1), ha='center', va='center')


def plot_trajectories(x_data, y_data, window_sizeX, window_sizeY, vec_nums, text_pos, scale_bar):
    end_vecs = np.concatenate((np.expand_dims(x_data[-1, :] - x_data[0, :], 0),
                               np.expand_dims(y_data[-1, :] - y_data[0, :], 0)), 0)
    fitness = np.sum((end_vecs.mean(1) / 50) ** 2) ** 0.5
    av_dist = np.mean(np.sum(end_vecs ** 2, 0) ** 0.5)
    print('distances')
    print(np.sum(end_vecs ** 2, 0) ** 0.5)
    print('average distance covered = ' + str(av_dist))
    angles = np.arctan2(end_vecs[1, :], end_vecs[0, :]) * 180 / np.pi
    angles2 = np.mod(angles, 360)
    std_angle = np.min([np.std(angles), np.std(angles2)])
    print('standard deviation of directions = ' + str(std_angle))
    print('angles: ')
    print(angles)
    print(angles2)
    print('average direction:')
    print(np.mean(angles))
    print(np.mean(angles2))
    plot_settings = {'traj': {'sym': '-', 'col': [1, 0, 1]},
                     'vec': {'sym': '--', 'col': [0.5, 0.5, 0.5]},
                     'start': {'sym': 'd', 'col': [0, 0, 0]},
                     'final': {'sym': '*', 'col': [0, 0, 0]}}
    labels = ['Trajectories', 'Start positions', 'End positions', 'End vectors']
    fig1, ax1 = plt.subplots()
    fig2, ax2 = plt.subplots()
    for c in range(0, x_data.shape[1]):
        ax1.plot(x_data[:, c], y_data[:, c], plot_settings['traj']['sym'])  # , label=labels[0])  #color=plot_settings['traj']['col']
        if c in vec_nums:
            ax2.plot(x_data[:, c], y_data[:, c], plot_settings['traj']['sym']) #, label=labels[0])  #color=plot_settings['traj']['col']
        # else:
        ax2.plot(x_data[[0, -1], c], y_data[[0, -1], c], plot_settings['vec']['sym'],
                 color=plot_settings['vec']['col'], label=labels[3])
        ax1.plot(x_data[0, c], y_data[0, c], plot_settings['start']['sym'],
                 color=plot_settings['start']['col'], label=labels[1])
        ax1.plot(x_data[-1, c], y_data[-1, c], plot_settings['final']['sym'],
                 color=plot_settings['final']['col'], label=labels[2])
        # ax2.plot(x_data[[0, -1], c], y_data[[0, -1], c], plot_settings['vec']['sym'],
        #          color=plot_settings['vec']['col'], label=labels[3])
        ax2.plot(x_data[0, c], y_data[0, c], plot_settings['start']['sym'],
                 color=plot_settings['start']['col'], label=labels[1])
        ax2.plot(x_data[-1, c], y_data[-1, c], plot_settings['final']['sym'],
                 color=plot_settings['final']['col'], label=labels[2])
        labels = ['', '', '', '']
        # plt.plot(traj[0], traj[1], colours[4] + '-', label='Trial ' + str(i))
        # plt.plot(traj[0][0], traj[1][0], colours[4] + 'd')
        # plt.plot(traj[0][-1], traj[1][-1], colours[4] + '*')
    # ax1.set_title('Example trajectories of two-cell network')
    ax2.text(text_pos[0], text_pos[1],
             '$D = $' + str(np.round(av_dist, 1)) + '$\mu m$ \n $\sigma = $' + str(np.round(std_angle, 1)) + '$^\circ$',
             fontsize=16)
    ax2.plot([scale_bar[0], scale_bar[0] + 10], [scale_bar[1], scale_bar[1]],
             color=[0, 0, 0], linewidth=5, label='Cell Diameter')
    ax1.axis([window_sizeX[0], window_sizeX[1], window_sizeY[0], window_sizeY[1]])
    ax1.legend()
    # ax2.set_title('Sampled End Vectors of Cell Network Trajectories')
    ax2.axis([window_sizeX[0], window_sizeX[1], window_sizeY[0], window_sizeY[1]])
    ax2.legend()
    plt.show()


def summarise_trajectories(root, file_name, scaling_mode):
    upper_scale_lim = 3
    magnitude_mat = np.zeros((upper_scale_lim, upper_scale_lim))
    direction_mat = np.zeros((upper_scale_lim, upper_scale_lim))
    mag_dev_mat = np.zeros((upper_scale_lim, upper_scale_lim))
    dir_dev_mat = np.zeros((upper_scale_lim, upper_scale_lim))
    fig, axs = plt.subplots(1, 2)
    for i in range(1, upper_scale_lim + 1):
        for j in range(1, upper_scale_lim + 1):
            if i == 1 and j == 1:
                temp_scale_mode = 'sca'
            else:
                temp_scale_mode = scaling_mode
            xData = pd.read_csv(root + 'xData_' + file_name + temp_scale_mode +
                                str(i) + str(j) + '.csv', header=None).to_numpy()
            yData = pd.read_csv(root + 'yData_' + file_name + temp_scale_mode +
                                str(i) + str(j) + '.csv', header=None).to_numpy()
            end_vectors = np.concatenate((np.expand_dims(xData[-1, :] - xData[0, :], 1),
                                          np.expand_dims(yData[-1, :] - yData[0, :], 1)), 1)
            average_vec = end_vectors.sum(0, keepdims=True) / end_vectors.shape[0]
            magnitude_mat[i - 1, j - 1] = np.sum(average_vec ** 2) ** 0.5
            mag_dev_mat[i - 1, j - 1] = np.var(np.sum(end_vectors ** 2, 1) ** 0.5) ** 0.5
            average_angle = np.arctan2(average_vec[0, 1], average_vec[0, 0])
            direction_mat[i - 1, j - 1] = np.mod(average_angle, 2 * np.pi)
            directions = np.arctan2(end_vectors[:, 1], end_vectors[:, 0])
            angle_diffs_sym = directions - average_angle
            angle_diffs_pos = np.mod(directions, 2 * np.pi) - np.mod(average_angle, 2 * np.pi)
            angle_diffs = np.zeros(directions.shape)
            for a in range(angle_diffs.shape[0]):
                if np.abs(angle_diffs_pos[a]) <= np.abs(angle_diffs_sym[a]):
                    angle_diffs[a] = angle_diffs_pos[a]
                else:
                    angle_diffs[a] = angle_diffs_sym[a]
            dir_dev_mat[i - 1, j - 1] = np.var(angle_diffs) ** 0.5
            axs[0].plot(direction_mat[i - 1, j - 1], magnitude_mat[i - 1, j - 1], 'x', color=[i / 3, 0, j / 3])
            axs[1].plot(direction_mat[i - 1, j - 1], magnitude_mat[i - 1, j - 1], 'x', color=[i / 3, 0, j / 3])
            # ax.errorbar(direction_mat[i - 1, j - 1], magnitude_mat[i - 1, j - 1], 'kx', yerr=mag_dev_mat[i - 1, j - 1])
            axs[0].plot(direction_mat[i - 1, j - 1] * np.array([1, 1]), magnitude_mat[i - 1, j - 1] +
                        np.array([1, -1]) * mag_dev_mat[i - 1, j - 1], '--', color=[i / 3, 0, j / 3])
            axs[1].plot(direction_mat[i - 1, j - 1] + np.array([1, -1]) * dir_dev_mat[i - 1, j - 1],
                        magnitude_mat[i - 1, j - 1] * np.array([1, 1]), '--', color=[i / 3, 0, j / 3])
            axs[0].annotate('(' + str(i) + ', ' + str(j) + ')',
                            xy=(direction_mat[i - 1, j - 1], magnitude_mat[i - 1, j - 1]))
            axs[1].annotate('(' + str(i) + ', ' + str(j) + ')',
                            xy=(direction_mat[i - 1, j - 1], magnitude_mat[i - 1, j - 1]))
            # ax.text(magnitude_mat[i - 1, j - 1], direction_mat[i - 1, j - 1], 3, ha='center', va='center')
    plt.show()


            # angle_diffs_raw = directions - np.arctan2(average_vec[0, 1], average_vec[0, 0])
            # angle_diffs = angle_diffs_raw.copy()
            # for a in range(0, angle_diffs.shape[0]):
            #     if np.pi < angle_diffs_raw[a]:
            #         angle_diffs[a] = angle_diffs_raw[a] - 2 * np.pi
            #     elif angle_diffs_raw[a] <= -np.pi:
            #         angle_diffs[a] = 2 * np.pi + angle_diffs_raw[a]
            # end_magnitudes = np.sum(end_vectors ** 2, 1) ** 0.5
            # average_arc = angle_diffs * end_magnitudes / sum(end_magnitudes)




def main():
    early_colour = 'e8d268'
    late_colour = '8738d2'
    passive_colour = '27cbbc'
    nocell_colour = '222222'
    # actuation_cycles(early_colour, late_colour, passive_colour)

    root = 'C:/CompuCell3D-py3-64bit/lib/site-packages/MySimulations/SavedConfigurations/'

    data_df1 = pd.read_csv(root + 'oldSimulations/ev_fit_rec_55_COM_success.csv')
    data_df2 = pd.read_csv(root + 'oldSimulations/ev_fit_rec_55_COM_attempt2.csv')
    data_df3 = pd.read_csv(root + 'oldSimulations/ev_fit_rec_55_COM_thirdTry.csv')
    show_fitness_per_gen([data_df1, data_df2, data_df3], 0.2)  # , data_df3, data_df4], 0.4)
    # state_abundance([data_df3], 0.2, [convert_colour(early_colour),
    #                                   convert_colour(late_colour),
    #                                   convert_colour(passive_colour),
    #                                   convert_colour(nocell_colour)])  # data_df1, data_df2,

    # neighbour_distribution([data_df1], [5, 5], 'car')  # , data_df2, data_df3

    # xData = pd.read_csv(root + 'trajectories/xData_twoCell_car55_sca11.csv', header=None).to_numpy()
    # yData = pd.read_csv(root + 'trajectories/yData_twoCell_car55_sca11.csv', header=None).to_numpy()
    # xData = pd.read_csv(root + 'trajectories/xData_trial1_car55_sca11.csv', header=None).to_numpy()
    # yData = pd.read_csv(root + 'trajectories/yData_trial1_car55_sca11.csv', header=None).to_numpy()
    # xData = pd.read_csv(root + 'trajectories/xData_custom2_car55_sca11aaa.csv', header=None).to_numpy()
    # yData = pd.read_csv(root + 'trajectories/yData_custom2_car55_sca11aaa.csv', header=None).to_numpy()

    traj_label = ['twoCell', 'trial1', 'trial4', 'trial3', 'custom2r', 'custom2u', 'custom2l', 'custom2d']
    lab_num = 3
    scale = ['', 'sca11', 'exp22', 'exp33']
    scale_num = 1
    view = {'twoCell': {'center': [0, -5], 'range': 9, 'vecs': [2, 5, 7],
                        'text': [2.5, -12], 'scale': [-2.5, -12.5]},
            'trial1': {'center': [-10, -12], 'range': 18, 'vecs': [2, 3, 7],
                       'text': [-26, -10], 'scale': [-26, -11]},
            'trial4': {'center': [0, -23], 'range': 28, 'vecs': [1, 2, 7],
                       'text': [-24, -12], 'scale': [-22, -3]},
            'trial3': {'center': [-27, 0], 'range': 33, 'vecs': [6, 7, 9],
                       'text': [-55, -25], 'scale': [-55, -28]},
            'custom2r': {'center': [48, 0], 'range': 50, 'vecs': [],
                         'text': [0, 0], 'scale': [0, 0]},
            'custom2u': {'center': [0, 0], 'range': 10, 'vecs': [],
                         'text': [0, 0], 'scale': [0, 0]},
            'custom2l': {'center': [0, 0], 'range': 10, 'vecs': [],
                         'text': [0, 0], 'scale': [0, 0]},
            'custom2d': {'center': [0, 0], 'range': 10, 'vecs': [],
                         'text': [0, 0], 'scale': [0, 0]}}

    xData = pd.read_csv((root + 'trajectories/' + traj_label[lab_num] + '_xData_car55_' +
                         scale[scale_num] + '_act50.csv'),
                        header=None).to_numpy()
    yData = pd.read_csv((root + 'trajectories/' + traj_label[lab_num] + '_yData_car55_' +
                         scale[scale_num] + '_act50.csv'),
                        header=None).to_numpy()
    # xData = pd.read_csv(root + 'trajectories/trial4_xData_car55_sca11_act50.csv', header=None).to_numpy()
    # yData = pd.read_csv(root + 'trajectories/trial4_yData_car55_sca11_act50.csv', header=None).to_numpy()
    # center = [-20, 0]
    center = view[traj_label[lab_num]]['center']
    # range = 25
    range = view[traj_label[lab_num]]['range']
    traj_vec = view[traj_label[lab_num]]['vecs']
    text_pos = view[traj_label[lab_num]]['text']
    scale = view[traj_label[lab_num]]['scale']
    # plot_trajectories(xData[:, 0:10] - xData[0, 0], yData[:, 0:10] - yData[0, 0],
    #                   [center[0] - range, center[0] + range],
    #                   [center[1] - range, center[1] + range], traj_vec, text_pos, scale)

    # summarise_trajectories(root + 'trajectories/', 'trial3_car55_', 'exp')



    # data_df1 = pd.read_csv(root + 'ev_fit_rep4_phe55_sca11_sur40_cro50_mut1_rot0aaa.csv')
    # data_df2 = pd.read_csv(root + 'ev_fit_rep4_phe55_sca11_sur40_cro50_mut5_rot0bbb.csv')
    # data_df3 = pd.read_csv(root + 'ev_fit_rep4_phe55_sca11_sur40_cro50_mut5_rot0ccc.csv')
    # data_df4 = pd.read_csv(root + 'ev_fit_rep4_phe55_sca11_sur40_cro50_mut1_rot1aaa.csv')


    # show_diversity(data_df1)
    # compare_fitness(data_df2)
    # need to add separate (histogram) plots for repeat / non-evolved data


if __name__ == '__main__':
    main()
