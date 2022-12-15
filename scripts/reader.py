import numpy as np
import pandas as pd

def get_truth_particles(tree, evN):
    particles=pd.DataFrame()

    particle_px = tree.array("particle_px")[evN]/1e3
    particle_py = tree.array("particle_py")[evN]/1e3
    particle_pz = tree.array("particle_pz")[evN]/1e3
    particle_e = tree.array("particle_e")[evN]/1e3

    particle_phi=np.arctan2(particle_py,particle_px)
    
    particle_pt = np.sqrt(np.square(particle_px)+np.square(particle_py))
    particle_theta=np.arctan2(particle_pt,particle_pz)
    particle_eta=-np.log(np.tan(particle_theta/2.))

    
    particle_pdgid = tree.array("particle_pdgid")[evN]
    particle_pdgid_name=np.asarray([get_pdgid_name(pdgid) for pdgid in particle_pdgid])
    particle_class=np.asarray([get_class(pdgid) for pdgid in particle_pdgid])
    particle_class_name=np.asarray([get_class_name[cl] for cl in particle_class])
    #particle_color=np.asarray([get_truth_particle_color[cl] for cl in particle_class])
    
    particles['pT'] = particle_pt
    particles['Phi'] = particle_phi
    particles['Eta'] = particle_eta
    #particles['Class'] = particle_class
    particles['Class'] = particle_class_name
    particles['Particle'] = particle_pdgid_name
    #particles['Color'] = particle_color

    return particles

def get_calorimeter_info(tree, evN):
    calorimeter=pd.DataFrame()

    #detector
    cell_eta = tree.array("cell_eta")[evN]
    cell_phi = map_phi_in_pi(tree.array("cell_phi")[evN])
    cell_e = tree.array("cell_e")[evN]/1e3
    cell_layer = tree.array("cell_layer")[evN]
    #cell_eT=cell_e/np.cosh(cell_eta)[evN]
    cell_topo_idx = tree.array("cell_topo_idx")[evN]
    cell_parent_idx = tree.array("cell_parent_idx")[evN]

    calorimeter['Eta']   = cell_eta
    calorimeter['Phi']   = cell_phi
    calorimeter['En']    = cell_e
    calorimeter['Topo']  = cell_topo_idx
    calorimeter['Parent']  = cell_parent_idx
    calorimeter['Layer'] = cell_layer

    return calorimeter

def get_tracks_info(tree, evN):
    tracks=pd.DataFrame()

    track_qoverp = tree.array("track_qoverp")[evN]
    track_phi = tree.array("track_phi")[evN]
    track_theta = tree.array("track_theta")[evN]
    track_eta=-np.log(np.tan(track_theta/2.))
    track_pt=abs(1/track_qoverp)*np.sin(track_theta)/1e3 #to GeV

    track_not_reg=tree.array("track_not_reg")[evN]    
    condition_reg=np.where(track_not_reg==0)


    tracks["pT"] = track_pt[condition_reg]
    tracks["Eta0"] = track_eta[condition_reg]
    tracks["Eta1"] = tree.array("track_eta_layer_0")[evN][condition_reg]
    tracks["Eta2"] = tree.array("track_eta_layer_1")[evN][condition_reg]
    tracks["Eta3"] = tree.array("track_eta_layer_2")[evN][condition_reg]
    tracks["Eta4"] = tree.array("track_eta_layer_3")[evN][condition_reg]
    tracks["Eta5"] = tree.array("track_eta_layer_4")[evN][condition_reg]
    tracks["Eta6"] = tree.array("track_eta_layer_5")[evN][condition_reg]
    tracks["Phi0"] = map_phi_in_pi(track_phi[condition_reg])
    tracks["Phi1"] = map_phi_in_pi(tree.array("track_phi_layer_0")[evN][condition_reg])
    tracks["Phi2"] = map_phi_in_pi(tree.array("track_phi_layer_1")[evN][condition_reg])
    tracks["Phi3"] = map_phi_in_pi(tree.array("track_phi_layer_2")[evN][condition_reg])
    tracks["Phi4"] = map_phi_in_pi(tree.array("track_phi_layer_3")[evN][condition_reg])
    tracks["Phi5"] = map_phi_in_pi(tree.array("track_phi_layer_4")[evN][condition_reg])
    tracks["Phi6"] = map_phi_in_pi(tree.array("track_phi_layer_5")[evN][condition_reg])

    return tracks



def get_jet_info(tree, evN, collection):
    jets = pd.DataFrame()

    jets["pT"] = tree.array(collection+"_jet_pt")[evN]/1e3 #in GeV
    jets["Eta"] = tree.array(collection+"_jet_eta")[evN]
    jets["Phi"] = map_phi_in_pi(tree.array(collection+"_jet_phi")[evN])
    constituents = tree.array(collection+"_jet_constituents_list")[evN]
    n = np.array( [len(c) for c in constituents])

    jets["n"] = n

    return jets



def get_pdgid_name(element):
    if element in [3222,-3222, 3112,-3112]:
        return "Sigma"
    elif element in [211,-211]:
        return "Pion"
    elif element in [321, -321]:
        return "Kaon"
    elif element in [2212,-2212]:
        return "Proton"
    elif element in [3312, -3312]:
        return "Csi"
    elif element in [3334, -3334]:
        return "Omega"
    elif element in [130, 310, 311]:
        return "Kaon"
    elif element in [2112, -2112]:
        return "Neutron"
    elif element in [3322,-3322]:
        return "Csi"
    elif element in [3122, -3122]:
        return "Lambda"
    elif element in [22]:
        return "Photon"
    elif element in [11,-11]:
        return "Electron"
    elif element in [13, -13]:
        return "Muon"
    else:
        return "Unknown"

def get_class(element):
    list_Gamma = [22]
    list_ChHad = [3222 ,-3222 ,211   ,-211  ,321   ,-321  ,2212  ,-2212 ,3112  ,-3112 ,3312  ,-3312 ,3334  ,-3334 ]
    list_NuHad = [130  ,310   ,311   ,2112  ,-2112 ,-3322 ,3322  ,3122  ,-3122]
    list_El    = [11   ,-11]
    list_Mu    = [13   ,-13]
    map_class=[list_Gamma, list_NuHad, list_ChHad, list_El, list_Mu]
    
    particle_class=-1
    for i,m in enumerate(map_class):
        if element in m:
            particle_class=i
            break
    if  particle_class==-1: return None

    return particle_class

get_class_name ={
    0:"Photon",
    1:"Neutral Hadron",
    2:"Charged Hadron",
    3:"El",
    4:"Mu"
}
get_truth_particle_color = {
    "Photon":"#AC9700",
    "Neutral Hadron":"#5D5D5D",
    "Charged Hadron":"#B60000",
    "Electron":"#007307",
    "Muon":"#6600A4"
}

calo_map ={
    0:"ECAL 1",
    1:"ECAL 2",
    2:"ECAL 3",
    3:"HCAL 4",
    4:"HCAL 5",
    5:"HCAL 6",
}
def set_range(config):
    
    eta_max = config["eta max"]
    phi_max = config["phi max"]
    if phi_max=="PI": phi_max=np.pi

    eta_min = config["eta min"]
    phi_min = config["phi min"]
    if phi_min=="-PI": phi_min=-np.pi
    

    return max(-3,eta_min), min(3,eta_max), max(-np.pi,phi_min), min(np.pi,phi_max) 




def map_phi_in_pi(arr):
    return (arr+np.pi)%(2*np.pi)-np.pi