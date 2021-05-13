from cc3d.core.PySteppables import *
from cc3d import CompuCellSetup
from pathlib import Path

import numpy as np
import pandas as pd
import random


time_mcs_ratio = 1


class FileReaderSteppable(SteppableBasePy):
    def __init__(self, frequency=1):
        SteppableBasePy.__init__(self, frequency)

    def start(self):
        global_data_file = r'C:/CompuCell3D-py3-64bit/lib/site-packages/MySimulations/ProtoCellSim/variables.csv'
        df = pd.read_csv(global_data_file)
        # self.set_max_mcs(int(df['Value'][df['Name'] == 'max mcs']))
        self.shared_steppable_vars['cell diam'] = float(df['Value'][df['Name'] == 'cell diam'])
        self.shared_steppable_vars['actuate period'] = float(df['Value'][df['Name'] == 'actuate period'])
        self.shared_steppable_vars['actuate scale'] = float(df['Value'][df['Name'] == 'actuate scale'])
        self.shared_steppable_vars['num active types'] = float(df['Value'][df['Name'] == 'num active types'])
        # self.shared_steppable_vars['light frequency'] = float(df['Value'][df['Name'] == 'light frequency'])


class ConstraintInitializerSteppable(SteppableBasePy):
    def __init__(self, frequency=1):
        SteppableBasePy.__init__(self, frequency)
        self.my_cell_radius = None

    def cell_setup(self, cell):
        cell.targetVolume = self.shared_steppable_vars['neutral volume']
        cell.lambdaVolume = 20.0
        self.lengthConstraintPlugin.setLengthConstraintData(cell, 20, 2 * self.my_cell_radius,
                                                            2 * self.my_cell_radius)
        cell.targetSurface = self.shared_steppable_vars['neutral surface']
        cell.lambdaSurface = 2.0

    def start(self):
        self.my_cell_radius = self.shared_steppable_vars['cell diam'] / 2
        self.shared_steppable_vars['neutral volume'] = int(np.floor(np.pi * self.my_cell_radius ** 2))
        self.shared_steppable_vars['neutral surface'] = 8 * self.my_cell_radius
        for cell in self.cell_list_by_type(self.PROTOCELLPASSIVE):  # , self.PROTOCELLACTIVE1, self.PROTOCELLACTIVE2):
            # cell.targetVolume = self.shared_steppable_vars['neutral volume']
            # cell.lambdaVolume = 20.0
            # self.lengthConstraintPlugin.setLengthConstraintData(cell, 20, 2 * self.my_cell_radius,
            #                                                     2 * self.my_cell_radius)
            # cell.targetSurface = self.shared_steppable_vars['neutral surface']
            # cell.lambdaSurface = 2.0
            self.cell_setup(cell)
            cell.dict['Phenotype'] = 'Passive'
            cell.type = self.PROTOCELL
        for cell in self.cell_list_by_type(self.PROTOCELLACTIVEA):
            # if cell.type == self.PROTOCELLACTIVE1:
            self.cell_setup(cell)
            cell.dict['Phenotype'] = 'Active1'
            cell.type = self.PROTOCELL
        for cell in self.cell_list_by_type(self.PROTOCELLACTIVEB):
            # elif cell.type == self.PROTOCELLACTIVE2:
            self.cell_setup(cell)
            cell.dict['Phenotype'] = 'Active2'
            cell.type = self.PROTOCELL


class CellMotilitySteppable(SteppableBasePy):
    def __init__(self, frequency=50):
        SteppableBasePy.__init__(self, frequency)
        self.my_max_force = 12.0 / time_mcs_ratio

    def start(self):
        ext_pot_alg = self.get_xml_element('Ext Pot Alg')
        if ext_pot_alg.cdata == 'CenterOfMassBased':
            self.my_max_force *= self.shared_steppable_vars['neutral volume']

    def step(self, mcs):
        for cell in self.cell_list_by_type(self.PROTOCELL):
            angle = 2 * np.pi * random.uniform(0, 1)
            magnitude = self.my_max_force * random.uniform(0, 1)
            cell.lambdaVecX = magnitude * np.cos(angle)
            cell.lambdaVecY = magnitude * np.sin(angle)


class ActuationSteppable(SteppableBasePy):
    def __init__(self, frequency=10):
        SteppableBasePy.__init__(self, frequency)
        self.my_period = None
        self.my_scale = None
        self.my_rate_per_vol = None
        self.my_max_vol = None
        self.my_num_act = None

    def start(self):
        self.my_period = self.shared_steppable_vars['actuate period'] * time_mcs_ratio
        self.my_scale = self.shared_steppable_vars['actuate scale']
        self.my_rate_per_vol = (self.my_scale - 1) / (self.my_period / self.frequency)
        self.my_max_vol = self.shared_steppable_vars['neutral volume'] * self.my_scale
        self.my_num_act = self.shared_steppable_vars['num active types']

    def step(self, mcs):

        for cell in self.cell_list_by_type(self.PROTOCELL):
            # print('--------------------------ACTUATION STEPPABLE-----------------------------')
            if 'Active' in cell.dict['Phenotype']:
                active_type = int(cell.dict['Phenotype'][-1])
                act_stage = np.mod(np.floor(mcs / self.my_period), self.my_num_act + 2)
                switch = (active_type <= act_stage <= active_type + self.my_num_act - 1)
                coeff = np.sign(int(switch) - 0.5)
                # illumination = green_field[cell.xCOM, cell.yCOM, 0]
                # threshold = 5
                # coeff = np.sign(int(illumination >= threshold) - 0.5)
                new_vol = (cell.targetVolume + coeff * self.my_rate_per_vol *
                           self.shared_steppable_vars['neutral volume'])
                new_vol_bounded = saturate(new_vol, self.shared_steppable_vars['neutral volume'], self.my_max_vol)
                cell.targetVolume = new_vol_bounded
                new_diam = vol2diam(new_vol_bounded)
                lambda_len = self.lengthConstraintPlugin.getLambdaLength(cell)
                self.lengthConstraintPlugin.setLengthConstraintData(cell, lambda_len, new_diam, new_diam)
                cell.targetSurface = 4 * new_diam
                cell.lambdaVolume = 10
                cell.lambdaSurface = 10


class LinkUpdaterSteppable(SteppableBasePy):
    def __init__(self, frequency=10):
        SteppableBasePy.__init__(self, frequency)
        self.FPPplugin = None

    def start(self):
        self.FPPplugin = CompuCell.getFocalPointPlasticityPlugin()

    def step(self, mcs):
        # print('--------------------------LINK UPDATER STEPPABLE-----------------------------')
        checked_off = []
        # red_field = self.field.RedLight
        for cell in self.cell_list_by_type(self.PROTOCELL):
            if cell.id not in checked_off:
                # checked_off.append(cell.id)
                # photocleave = (red_field[cell.xCOM, cell.yCOM, cell.zCOM] == 1)
                length1 = self.lengthConstraintPlugin.getTargetLength(cell)
                for fppd in FocalPointPlasticityDataList(self.FPPplugin, cell):
                    if fppd.neighborAddress.type == self.PROTOCELL:  # and fppd.neighborAddress.id not in checked_off:
                        # if photocleave:
                        #     self.FPPplugin.deleteFocalPointPlasticityLink(cell, fppd.neighborAddress)
                        #     continue
                        length2 = self.lengthConstraintPlugin.getTargetLength(fppd.neighborAddress)
                        dist = (length1 + length2) / 2 + 1  # + 0
                        lam_link = fppd.lambdaDistance
                        self.FPPplugin.setFocalPointPlasticityParameters(cell, fppd.neighborAddress,
                                                                         lam_link, dist, 100)
                        # checked_off.append(fppd.neighborAddress.id)


class SecretionSteppable(SteppableBasePy):
    def __init__(self, frequency=10):
        SteppableBasePy.__init__(self, frequency)
        self.FPPplugin = None

    def start(self):
        self.FPPplugin = CompuCell.getFocalPointPlasticityPlugin()

    def step(self, mcs):
        crosshair = self.get_field_secretor('CameraCrosshair')
        for cell in self.cell_list_by_type(self.PROTOCELL):
            num_links = len(FocalPointPlasticityDataList(self.FPPplugin, cell))
            crosshair.secreteOutsideCellAtBoundary(cell, num_links)
            crosshair.secreteInsideCellAtBoundary(cell, num_links)


class LightControlSteppable(SteppableBasePy):
    def __init__(self, frequency=10):
        SteppableBasePy.__init__(self, frequency)
        self.my_pixel_dims = [8, 8]
        self.my_threshold = 5
        self.my_light_switch = 0

    def step(self, mcs):
        red_field = self.field.RedLight
        green_field = self.field.GreenLight
        # crosshair_field = self.field.CameraCrosshair
        red_field[:, :, 0] = 0
        self.my_light_switch = (np.mod(mcs, 2 * self.shared_steppable_vars['light frequency']) /
                                (2 * self.shared_steppable_vars['light frequency'] - 1) >= 0.5)
        green_field[:, :, 0] = int(10 * self.my_light_switch)
        # for cell in self.cell_list_by_type(self.PROTOCELL):
        #     center = [cell.xCOM, cell.yCOM, 0]
        #     if crosshair_field[center[0], center[1], center[2]] > self.my_threshold:
        #         coord = [int(np.floor(center[0] / self.my_pixel_dims[0])),
        #                  int(np.floor(center[1] / self.my_pixel_dims[1])), 0]
        #         red_field[self.my_pixel_dims[0] * coord[0]:self.my_pixel_dims[0] * (coord[0] + 1),
        #                   self.my_pixel_dims[1] * coord[1]:self.my_pixel_dims[1] * (coord[1] + 1), 0] = 10


class RecordSteppable(SteppableBasePy):
    def __init__(self, frequency=100):
        SteppableBasePy.__init__(self, frequency)
        self.my_trackX = []
        self.my_trackY = []

    def step(self, mcs):
        # print('--------------------------RECORD STEPPABLE-----------------------------')
        start_time = self.shared_steppable_vars['actuate period']
        if start_time <= mcs < start_time + self.frequency:
            for cell in self.cell_list_by_type(self.PROTOCELL):
                cell.dict['initial data'] = np.array([cell.xCOM, cell.yCOM, cell.zCOM])
        if np.mod(mcs, 4 * self.shared_steppable_vars['actuate period']) == 0:
            x_coord = 0
            y_coord = 0
            for cell in self.cell_list_by_type(self.PROTOCELL):
                x_coord += cell.xCOM
                y_coord += cell.yCOM
            self.my_trackX.append(x_coord / len(self.cell_list_by_type(self.PROTOCELL)))
            self.my_trackY.append(y_coord / len(self.cell_list_by_type(self.PROTOCELL)))

                #cell.lambdaVecY = -10.0
        # for i, cell in enumerate(self.cell_list_by_type(self.PROTOCELL)):
        #     if i == 1:
        #         self.my_track2.append(cell.yCOM)
        #         self.my_track3.append(cell.targetVolume)
        #     elif i == 0:
        #         self.my_track1.append(cell.yCOM)

        # elif 2500 <= mcs < 2500 + self.frequency:
        #     for cell in self.cell_list_by_type(self.PROTOCELL):
        #         cell.dict['actuated_data'] = np.array([cell.xCOM, cell.yCOM, cell.zCOM])
        # elif 4000 <= mcs < 4000 + self.frequency:
        #     for cell in self.cell_list_by_type(self.PROTOCELL):
        #         cell.dict['final_data'] = np.array([cell.xCOM, cell.yCOM, cell.zCOM])

    def finish(self):
        for cell in self.cell_list_by_type(self.PROTOCELL):
            cell.dict['final data'] = np.array([cell.xCOM, cell.yCOM, cell.zCOM])
        outputs = ['initial data', 'final data']
        pg = CompuCellSetup.persistent_globals
        return_dict = {}
        for label in outputs:
            return_dict[label] = np.zeros((len(self.cell_list_by_type(self.PROTOCELL)), 2))
            for i, cell in enumerate(self.cell_list_by_type(self.PROTOCELL)):
                return_dict[label][i, :] = cell.dict[label][:2]
        return_dict['trackX'] = self.my_trackX
        return_dict['trackY'] = self.my_trackY
        # return_dict['track3'] = self.my_track3
        # print(self.my_track1[-1])
        # print(self.my_track2[-1])
        # print(self.my_track3)
        pg.return_object = return_dict

        # output_dir = self.output_dir
        # for d in range(0, len(outputs
        #     pg.return_object[outputs[d]] =
        #     output_path = Path(output_dir).joinpath('cell_' + outputs[d] + '.csv')
        #         with open(output_path, 'w') as fout:
        #             for cell in self.cell_list_by_type(self.PROTOCELL):
        #                 fout.write(str(cell.id).zfill(3) + ' {} {}\n'.format(cell.dict[outputs[d]][0],
        #                                                                      cell.dict[outputs[d]][1]))


class PhotocleaveSteppable(SteppableBasePy):
    def __init__(self, frequency=10):
        SteppableBasePy.__init__(self, frequency)
        decay_period = 60
        growth_period = 150
        self.my_max_density = 3.0
        self.my_decay_rate = self.my_max_density / (decay_period / frequency)
        self.my_growth_rate = self.my_max_density / (growth_period / frequency)

    def step(self, mcs):
        red_field = self.field.RedLight
        for cell in self.cell_list_by_type(self.PROTOCELL):
            current_linker_density = self.adhesionFlexPlugin.getAdhesionMoleculeDensity(cell, 'Linker')
            if red_field[cell.xCOM, cell.yCOM, cell.zCOM] == 1:
                new_density = current_linker_density - self.my_decay_rate
            else:
                new_density = current_linker_density + self.my_growth_rate
            new_density = np.array([new_density, 0]).max()
            new_density = np.array([new_density, self.my_max_density]).min()
            self.adhesionFlexPlugin.setAdhesionMoleculeDensity(cell, 'Linker', new_density)


class CollectionSteppable(SteppableBasePy):
    def __init__(self, frequency=10):
        SteppableBasePy.__init__(self, frequency)
        self.my_lattice_dims = None

    def start(self):
        lattice_dims = np.array([[0, 0, 0]]).T
        for x, y, z in self.every_pixel():
            new_coord = np.array([[x, y, z]]).T
            lattice_dims = np.max(np.concatenate((lattice_dims, new_coord), 1), 1, keepdims=True)
        self.my_lattice_dims = list(lattice_dims[:, 0] + 1)
        attr_field = self.field.ATTR
        attr_field[int(self.my_lattice_dims[0] / 2) - 10:int(self.my_lattice_dims[0] / 2) + 10,
                   int(self.my_lattice_dims[1] / 2) - 10:int(self.my_lattice_dims[1] / 2) + 10, 0] = 5
        # for x, y, z in self.every_pixel():
        #     attr_field[x, y, z] = np.abs(x - self.my_lattice_dims[0] / 2) + np.abs(y - self.my_lattice_dims[1] / 2)


def vol2diam(vol):
    return 2 * (vol / np.pi)**0.5


# def radius2volume(rad):
#     return np.pi * rad**2


def saturate(val, min_bound, max_bound):
    val1 = np.array([val, min_bound]).max()
    return np.array([val1, max_bound]).min()
