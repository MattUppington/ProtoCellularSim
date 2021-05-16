# from cc3d import CompuCellSetup
from ProtoCellularSimSteppables import *

from MySimulations.setup_cells import get_dep_configs

import os
import json


def configure_simulation():
    config_filename = os.path.join(os.getcwd(), 'config.json')
    with open(config_filename) as config_file:
        config = json.load(config_file)
    dep_config = get_dep_configs(config)
    cell_diameter = config['cell diameter']
    sim_dim = dep_config['sim dimensions']
    max_mcs = dep_config['max mcs']
    work_nodes = config['work nodes']

    from cc3d.core.XMLUtils import ElementCC3D
    compucell3d_element = ElementCC3D("CompuCell3D", {"Revision": "20210123", "Version": "4.2.4"})
    meta_data_element = compucell3d_element.ElementCC3D("Metadata")
    meta_data_element.ElementCC3D("NumberOfProcessors", {}, str(work_nodes))
    meta_data_element.ElementCC3D("DebugOutputFrequency", {}, "10")

    potts_element = compucell3d_element.ElementCC3D("Potts")
    potts_element.ElementCC3D("Dimensions", {"x": str(sim_dim[0]), "y": str(sim_dim[1]), "z": "1"})
    potts_element.ElementCC3D("Steps", {}, str(max_mcs))
    potts_element.ElementCC3D("Temperature", {}, "3.0")
    potts_element.ElementCC3D("NeighborOrder", {}, "1")
    potts_element.ElementCC3D("Boundary_x", {}, "Periodic")
    potts_element.ElementCC3D("Boundary_y", {}, "Periodic")

    plugin_element = compucell3d_element.ElementCC3D("Plugin", {"Name": "CellType"})
    plugin_element.ElementCC3D("CellType", {"TypeId": "0", "TypeName": "Medium"})
    plugin_element.ElementCC3D("CellType", {"Freeze": "", "TypeId": "1", "TypeName": "Wall"})
    plugin_element.ElementCC3D("CellType", {"TypeId": "2", "TypeName": "Protocell"})
    plugin_element.ElementCC3D("CellType", {"TypeId": "3", "TypeName": "ProtocellPassive"})
    plugin_element.ElementCC3D("CellType", {"TypeId": "4", "TypeName": "ProtocellActiveA"})
    plugin_element.ElementCC3D("CellType", {"TypeId": "5", "TypeName": "ProtocellActiveB"})

    compucell3d_element.ElementCC3D("Plugin", {"Name": "Volume"})

    compucell3d_element.ElementCC3D("Plugin", {"Name": "Surface"})

    plugin_element_1 = compucell3d_element.ElementCC3D("Plugin", {"Name": "ExternalPotential"})
    plugin_element_1.ElementCC3D("Algorithm", {"id": "Ext Pot Alg"}, "PixelBased")

    compucell3d_element.ElementCC3D("Plugin", {"Name": "CenterOfMass"})

    compucell3d_element.ElementCC3D("Plugin", {"Name": "NeighborTracker"})

    compucell3d_element.ElementCC3D("Plugin", {"Name": "PixelTracker"})

    plugin_element_2 = compucell3d_element.ElementCC3D("Plugin", {"Name": "Contact"})
    plugin_element_2.ElementCC3D("Energy", {"Type1": "Medium", "Type2": "Medium"}, "0.0")
    plugin_element_2.ElementCC3D("Energy", {"Type1": "Medium", "Type2": "Wall"}, "0.0")
    plugin_element_2.ElementCC3D("Energy", {"Type1": "Medium", "Type2": "Protocell"}, "50.0")  # 20
    plugin_element_2.ElementCC3D("Energy", {"Type1": "Wall", "Type2": "Wall"}, "10.0")
    plugin_element_2.ElementCC3D("Energy", {"Type1": "Wall", "Type2": "Protocell"}, "30.0")
    plugin_element_2.ElementCC3D("Energy", {"Type1": "Protocell", "Type2": "Protocell"}, "200.0")  # 50
    plugin_element_2.ElementCC3D("NeighborOrder", {}, "1")  # 2

    plugin_element_3 = compucell3d_element.ElementCC3D("Plugin", {"Name": "FocalPointPlasticity"})
    plugin_element_3.ElementCC3D("Local")
    parameters_element = plugin_element_3.ElementCC3D("Parameters", {"Type1": "Protocell", "Type2": "Protocell"})
    parameters_element.ElementCC3D("Lambda", {}, "500")  # 100
    parameters_element.ElementCC3D("ActivationEnergy", {"id": "link energy"}, "-0")
    parameters_element.ElementCC3D("TargetDistance", {"id": "protocell diameter"}, str(cell_diameter))
    parameters_element.ElementCC3D("MaxDistance", {}, "100")
    parameters_element.ElementCC3D("MaxNumberOfJunctions", {}, "10")
    plugin_element_3.ElementCC3D("NeighborOrder", {}, "1")

    compucell3d_element.ElementCC3D("Plugin", {"Name": "LengthConstraint"})

    plugin_element_4 = compucell3d_element.ElementCC3D("Plugin", {"Name": "Connectivity"})
    plugin_element_4.ElementCC3D("Penalty", {}, "100000")

    compucell3d_element.ElementCC3D("Plugin", {"Name": "Secretion"})

    steppable_element = compucell3d_element.ElementCC3D("Steppable", {"Type": "PIFInitializer"})
    steppable_element.ElementCC3D("PIFName", {}, "init_cell_field.piff")

    CompuCellSetup.setSimulationXMLDescription(compucell3d_element)


CompuCellSetup.register_steppable(steppable=FileReaderSteppable())
CompuCellSetup.register_steppable(steppable=ConstraintInitializerSteppable())
CompuCellSetup.register_steppable(steppable=ActuationSteppable())
CompuCellSetup.register_steppable(steppable=LinkUpdaterSteppable())
CompuCellSetup.register_steppable(steppable=RecordSteppable())

# CompuCellSetup.register_steppable(steppable=SecretionSteppable())
# CompuCellSetup.register_steppable(steppable=LightControlSteppable())
# CompuCellSetup.register_steppable(steppable=CellMotilitySteppable())
# CompuCellSetup.register_steppable(steppable=PhotoCleaveSteppable())
# CompuCellSetup.register_steppable(steppable=CollectionSteppable())

# <XMLScript Type="XMLScript">Simulation/ProtoCellSim.xml</XMLScript>
configure_simulation()
CompuCellSetup.run()
