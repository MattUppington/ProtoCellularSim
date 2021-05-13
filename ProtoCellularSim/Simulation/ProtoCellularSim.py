from cc3d import CompuCellSetup
from ProtoCellularSimSteppables import *


def configure_simulation():
    import pandas as pd
    variables_file_path = r'C:/CompuCell3D-py3-64bit/lib/site-packages/MySimulations/ProtoCellSim/variables.csv'
    variables = pd.read_csv(variables_file_path)
    cell_diam = variables['Value'][variables['Name'] == 'cell diam'].values[0]
    # actuate_scale = variables['Value'][variables['Name'] == 'actuate scale']
    # actuate_period = None
    # num_active_types = variables['Value'][variables['Name'] == 'num active types']
    max_mcs = variables['Value'][variables['Name'] == 'max mcs'].values[0]
    work_nodes = variables['Value'][variables['Name'] == 'work nodes'].values[0]
    sim_dim = variables['Value'][variables['Name'] == 'sim dim'].values[0]

    from cc3d.core.XMLUtils import ElementCC3D
    CompuCell3DElmnt = ElementCC3D("CompuCell3D", {"Revision": "20200518", "Version": "4.2.1"})
    MetadataElmnt = CompuCell3DElmnt.ElementCC3D("Metadata")
    MetadataElmnt.ElementCC3D("NumberOfProcessors", {}, str(work_nodes))  # str(work_nodes)
    MetadataElmnt.ElementCC3D("DebugOutputFrequency", {}, "10")

    PottsElmnt = CompuCell3DElmnt.ElementCC3D("Potts")
    PottsElmnt.ElementCC3D("Dimensions", {"x": str(sim_dim), "y": str(sim_dim), "z": "1"})
    PottsElmnt.ElementCC3D("Steps", {}, str(max_mcs))
    PottsElmnt.ElementCC3D("Temperature", {}, "3.0")
    PottsElmnt.ElementCC3D("NeighborOrder", {}, "1")
    PottsElmnt.ElementCC3D("Boundary_x", {}, "Periodic")
    PottsElmnt.ElementCC3D("Boundary_y", {}, "Periodic")

    PluginElmnt = CompuCell3DElmnt.ElementCC3D("Plugin", {"Name": "CellType"})
    PluginElmnt.ElementCC3D("CellType", {"TypeId": "0", "TypeName": "Medium"})
    PluginElmnt.ElementCC3D("CellType", {"Freeze": "", "TypeId": "1", "TypeName": "Wall"})
    PluginElmnt.ElementCC3D("CellType", {"TypeId": "2", "TypeName": "Protocell"})
    PluginElmnt.ElementCC3D("CellType", {"TypeId": "3", "TypeName": "ProtocellPassive"})
    PluginElmnt.ElementCC3D("CellType", {"TypeId": "4", "TypeName": "ProtocellActiveA"})
    PluginElmnt.ElementCC3D("CellType", {"TypeId": "5", "TypeName": "ProtocellActiveB"})

    CompuCell3DElmnt.ElementCC3D("Plugin", {"Name": "Volume"})

    CompuCell3DElmnt.ElementCC3D("Plugin", {"Name": "Surface"})

    PluginElmnt_1 = CompuCell3DElmnt.ElementCC3D("Plugin", {"Name": "ExternalPotential"})
    PluginElmnt_1.ElementCC3D("Algorithm", {"id": "Ext Pot Alg"}, "PixelBased")

    CompuCell3DElmnt.ElementCC3D("Plugin", {"Name": "CenterOfMass"})

    CompuCell3DElmnt.ElementCC3D("Plugin", {"Name": "NeighborTracker"})

    CompuCell3DElmnt.ElementCC3D("Plugin", {"Name": "PixelTracker"})

    PluginElmnt_2 = CompuCell3DElmnt.ElementCC3D("Plugin", {"Name": "Contact"})
    PluginElmnt_2.ElementCC3D("Energy", {"Type1": "Medium", "Type2": "Medium"}, "0.0")
    PluginElmnt_2.ElementCC3D("Energy", {"Type1": "Medium", "Type2": "Wall"}, "0.0")
    PluginElmnt_2.ElementCC3D("Energy", {"Type1": "Medium", "Type2": "Protocell"}, "50.0") #20
    PluginElmnt_2.ElementCC3D("Energy", {"Type1": "Wall", "Type2": "Wall"}, "10.0")
    PluginElmnt_2.ElementCC3D("Energy", {"Type1": "Wall", "Type2": "Protocell"}, "30.0")
    PluginElmnt_2.ElementCC3D("Energy", {"Type1": "Protocell", "Type2": "Protocell"}, "200.0") #50
    PluginElmnt_2.ElementCC3D("NeighborOrder", {}, "1")  # 2

    PluginElmnt_3 = CompuCell3DElmnt.ElementCC3D("Plugin", {"Name": "FocalPointPlasticity"})
    PluginElmnt_3.ElementCC3D("Local")
    ParametersElmnt = PluginElmnt_3.ElementCC3D("Parameters", {"Type1": "Protocell", "Type2": "Protocell"})
    ParametersElmnt.ElementCC3D("Lambda", {}, "500") #100
    ParametersElmnt.ElementCC3D("ActivationEnergy", {"id": "link energy"}, "-0")
    ParametersElmnt.ElementCC3D("TargetDistance", {"id": "protocell diameter"}, str(cell_diam))
    ParametersElmnt.ElementCC3D("MaxDistance", {}, "100")
    ParametersElmnt.ElementCC3D("MaxNumberOfJunctions", {}, "10")
    PluginElmnt_3.ElementCC3D("NeighborOrder", {}, "1")

    CompuCell3DElmnt.ElementCC3D("Plugin", {"Name": "LengthConstraint"})

    PluginElmnt_4 = CompuCell3DElmnt.ElementCC3D("Plugin", {"Name": "Connectivity"})
    PluginElmnt_4.ElementCC3D("Penalty", {}, "100000")

    CompuCell3DElmnt.ElementCC3D("Plugin", {"Name": "Secretion"})

    SteppableElmnt = CompuCell3DElmnt.ElementCC3D("Steppable", {"Type": "PIFInitializer"})
    SteppableElmnt.ElementCC3D("PIFName", {}, "init_cell_field.piff")

    CompuCellSetup.setSimulationXMLDescription(CompuCell3DElmnt)


CompuCellSetup.register_steppable(steppable=FileReaderSteppable())
CompuCellSetup.register_steppable(steppable=ConstraintInitializerSteppable())
CompuCellSetup.register_steppable(steppable=ActuationSteppable())
CompuCellSetup.register_steppable(steppable=LinkUpdaterSteppable())
CompuCellSetup.register_steppable(steppable=RecordSteppable())

# CompuCellSetup.register_steppable(steppable=SecretionSteppable())
# CompuCellSetup.register_steppable(steppable=LightControlSteppable())
# CompuCellSetup.register_steppable(steppable=CellMotilitySteppable())
# CompuCellSetup.register_steppable(steppable=PhotocleaveSteppable())
# CompuCellSetup.register_steppable(steppable=CollectionSteppable())

# <XMLScript Type="XMLScript">Simulation/ProtoCellSim.xml</XMLScript>
configure_simulation()
CompuCellSetup.run()
