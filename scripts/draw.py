import numpy as np
import pandas as pd
import plotly.graph_objects as go
from reader import *

def draw_calo(figure, config, lay, calorimeter):
    col_idx, row_idx = (lay)//3+1,(lay)%3+1
    cbarlocs_y = [.164,.5,.836][::-1]
    cbarlocs_x = [.45, 1.01]
        
    cbar_max = max(calorimeter["En"])
    if lay<3: cbar_max=max(calorimeter[calorimeter["Layer"]<3]["En"])
    else: cbar_max = max(calorimeter[calorimeter["Layer"]>2]["En"])

    calo_colo = "Blues"
    if(lay<3): calo_colo="Greens"

    eta_bins=np.linspace(-3,3,config["granularity"][lay])
    phi_bins=np.linspace(-np.pi,np.pi,config["granularity"][lay])
    layer, eta_edges, phi_edges = np.histogram2d(   calorimeter[calorimeter["Layer"]==lay]["Eta"], 
                                                            calorimeter[calorimeter["Layer"]==lay]["Phi"], 
                                                            weights=calorimeter[calorimeter["Layer"]==lay]["En"],
                                                            bins=(eta_bins,phi_bins), 
                                                            range=[[-3,3], [-np.pi, np.pi]])
    eta_centers=(eta_edges[:-1]+eta_edges[1:])/2
    phi_centers=(phi_edges[:-1]+phi_edges[1:])/2
    layer[np.where(layer<1e-6)]=1e-6

    if config["zscale"]=="log": 
        z_txt="log( E [GeV] )"
        layerT = layer.T
        layerT[np.where(layerT<=1e-6)]=1e-6
        layerT=np.log10(layerT)
    elif config["zscale"]=="lin": 
        z_txt="GeV"
        layerT = layer.T

    figure.add_trace(
        go.Heatmap(
            x = eta_centers,
            y = phi_centers,
            z = layerT,
            hovertemplate=
                    "<b>Cell Layer "+str(lay+1)+"</b><br><br>" +
                    "Eta: %{x:.3f}<br>" +
                    "Phi: %{y:.3f}<br>" +
                    "Energy: %{z:.3f} "+z_txt+"<br>" +
                    "<extra></extra>",
            type = 'heatmap',
            colorscale = calo_colo,
            colorbar=dict( tickfont_size=30, 
                           len=1/3,
                            x=cbarlocs_x[col_idx-1],  
                            y=cbarlocs_y[row_idx-1],
                            title_side="right",  
                            title=z_txt),         
            zmin=0,
            #zmax=cbar_max,
            name=str(calo_map[lay]),
        ),
        row=row_idx,
        col=col_idx,        
    )
    
def draw_tracks(figure, lay, tracks):
    col_idx, row_idx = (lay)//3+1,(lay)%3+1

    delta_phi=tracks["Phi"+str(lay+1)]-tracks["Phi0"]
    sym=np.array([ "star-triangle-down" for i in range(len(delta_phi))])
    sym[np.where(( (delta_phi>0) & (abs(delta_phi)<np.pi) ) | ( (delta_phi<0) & (abs(delta_phi)>np.pi) ) )]="star-triangle-up"

    figure.add_trace( #Extrapolation
        go.Scatter(
            x=tracks["Eta"+str(lay+1)],
            y=tracks["Phi"+str(lay+1)],
            text=tracks["pT"],
            mode='markers',
            name="Track Extrap. Lay",
            hovertemplate=
                "<b>Track Extrapolation Layer "+str(lay+1)+"</b><br><br>" +
                "Eta: %{x:.3f}<br>" +
                "Phi: %{y:.3f}<br>" +
                "pT: %{text:.3f} GeV<br>" +
                "<extra></extra>",
            marker=dict(
                size=10,
                symbol=sym,
                color="black",
                sizeref=1,
                line=dict(
                    color="black",
                    width=0
                ),  
            ),
        ),    
        row=row_idx,
        col=col_idx,
    )

    figure.add_trace( #PV
        go.Scatter(
            x=tracks["Eta0"],
            y=tracks["Phi0"],
            text=tracks["pT"],
            mode='markers', 
            name="Track PV",
            hovertemplate=
                "<b>Track PV</b><br><br>" +
                "Eta: %{x:.3f}<br>" +
                "Phi: %{y:.3f}<br>" +
                "pT: %{text:.3f} GeV<br>" +
                "<extra></extra>",
            marker=dict(
                size=6,
                symbol="x",
                color="black",
                sizeref=1,
                line=dict(
                    color="black",
                    width=0
                ),  
            ),
        ),    
        row=row_idx,
        col=col_idx,
    )

    for trk_idx in tracks.T:             #Loop over tracks
        for lay_idx in range(lay+1):  #Loop over previous layers
            deltaPhi_lay = abs(tracks.T[trk_idx]["Eta"+str(lay_idx+1)] - tracks.T[trk_idx]["Eta"+str(lay_idx)])

            eta_i, eta_f = tracks["Eta"+str(lay_idx)][trk_idx], tracks["Eta"+str(lay_idx+1)][trk_idx]
            phi_i, phi_f = tracks["Phi"+str(lay_idx)][trk_idx], tracks["Phi"+str(lay_idx+1)][trk_idx]
            
            if(deltaPhi_lay > np.pi):
                eta_m = ( eta_i * (np.pi - abs(phi_f)) + eta_f * (np.pi - abs(phi_i)) ) / (2 * np.pi - abs(phi_f) - abs(phi_i))
                phi_m = np.pi*np.sign(tracks["Phi"+str(lay_idx)][trk_idx])

                figure.add_shape(
                    type="line",
                    x0=eta_i,
                    y0=phi_i,
                    x1=eta_m,
                    y1=phi_m,
                    line=dict(
                        color="black",
                        width=1,
                    ),
                    row=row_idx,
                    col=col_idx,
                )
                figure.add_shape(
                    type="line",
                    x0=eta_m,
                    y0=phi_m,
                    x1=eta_f,
                    y1=phi_f,
                    line=dict(
                        color="black",
                        width=1,
                    ),
                    row=row_idx,
                    col=col_idx,
                )
            else:
                figure.add_shape(
                    type="line",
                    x0=eta_i,
                    y0=phi_i,
                    x1=eta_f,
                    y1=phi_f,
                    line=dict(
                        color="black",
                        width=1,
                    ),
                    row=row_idx,
                    col=col_idx,
                )

def draw_particles(figure, config, lay, particles):
    col_idx, row_idx = (lay)//3+1,(lay)%3+1
    
    classes = config["draw particles"]

    for cl in classes:
        if len(particles[particles["Class"]==cl])>0:
            figure.add_trace( 
                go.Scatter(
                    x=particles[particles["Class"]==cl]["Eta"], 
                    y=particles[particles["Class"]==cl]["Phi"],
                    name=cl,# | "+str(calo_map[lay]),
                    text=particles[particles["Class"]==cl]["Particle"],
                    marker_size=particles[particles["Class"]==cl]["pT"],
                    hovertemplate=
                    "<b>%{text}</b><br><br>" +
                    "Eta: %{x:.3f}<br>" +
                    "Phi: %{y:.3f}<br>" +
                    "pT: %{marker.size:.3f} GeV<br>" +
                    "<extra></extra>",
                    mode='markers',
                    marker=dict(
                        symbol="circle",
                        color=get_truth_particle_color[cl],
                        size=100*particles[particles["Class"]==cl]["pT"],
                        #line=dict(
                        #    color=color,
                        #    width=5,
                        #),
                        sizemode="area",
                        sizeref=0.005
                    ),
                    opacity=.9,
                #marker={'sizemode':'area',
                #            'sizeref':.05}
                    #color='Color',
                    #size='Pt [GeV]', hover_data=['Particle']
                ),
            row=row_idx,
            col=col_idx,

            )

def draw_truthjets(figure, config, lay, jets):
    col_idx, row_idx = (lay)//3+1,(lay)%3+1    

    r0 = config["R cone"]
    truthjet_color="#A34A00"
    

    if not config["only cone"]:
        figure.add_trace( 
            go.Scatter(
                x=jets['Eta'],
                y=jets['Phi'],
                name="Truth Jet",# | "+str(calo_map[lay]),
                marker_size=jets["pT"],
                hovertemplate=
                "<b>"+jet_name+"</b><br><br>" +
                "Eta: %{x:.3f}<br>" +
                "Phi: %{y:.3f}<br>" +
                "pT: %{marker.size:.3f} GeV<br>" +
                "<extra></extra>",
                mode='markers',
                marker=dict(
                    color=truthjet_color, 
                    size=jets["pT"],
                    #line=dict(
                    #    color=color,
                    #    width=0
                    #),
                    sizemode="area",
                    sizeref=0.003
                ),
                opacity=0.9
            ),
            row=row_idx,
            col=col_idx,
        )

    for index, row in jets.iterrows():
        figure.add_shape(
            type="circle",            
            x0=row["Eta"]-r0, y0=row["Phi"]-r0,
            x1=row["Eta"]+r0, y1=row["Phi"]+r0,
            fillcolor=None,
            opacity=1,     
            line=dict(
                color=truthjet_color,
                width=4,
                dash="solid",
            ),   
            row=row_idx,
            col=col_idx,
        )





def draw_jets(figure, config, lay, jets, collection):
    col_idx, row_idx = (lay)//3+1,(lay)%3+1    

    r0 = config["R cone"]
    jet_color = config["draw "+collection+"jets"]["color"]
    jet_name = config["draw "+collection+"jets"]["name"]
    if not config["only cone"]:
        figure.add_trace( 
            go.Scatter(
                x=jets['Eta'],
                y=jets['Phi'],
                name=jet_name,
                marker_size=jets["pT"],
                hovertemplate=
                "<b>"+jet_name+"</b><br><br>" +
                "Eta: %{x:.3f}<br>" +
                "Phi: %{y:.3f}<br>" +
                "pT: %{marker.size:.3f} GeV<br>" +
                "<extra></extra>",
                mode='markers',
                marker=dict(
                    color=jet_color, 
                    size=jets["pT"],
                    #line=dict(
                    #    color=color,
                    #    width=0
                    #),
                    sizemode="area",
                    sizeref=0.01
                ),
                opacity=0.9
            ),
            row=row_idx,
            col=col_idx,
        )

    for index, row in jets.iterrows():
        figure.add_shape(
            type="circle",            
            x0=row["Eta"]-r0, y0=row["Phi"]-r0,
            x1=row["Eta"]+r0, y1=row["Phi"]+r0,
            fillcolor=None,
            opacity=1,     
            line=dict(
                color=jet_color,
                width=4,
                dash=config["draw "+collection+"jets"]["line"],
            ),   
            row=row_idx,
            col=col_idx,
        )