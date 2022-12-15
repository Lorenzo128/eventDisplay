#import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
#import ROOT
import uproot3 as ur
import json
import os
import sys

#import plotly
#import chart_studio.plotly as py
#import plotly.graph_objects as go
#import plotly.express as px
from plotly.subplots import make_subplots

from reader import *
from draw import *

if __name__ == "__main__":

    print("===========================================================================")
    print("============================== EVENT DISPLAY ==============================")
    print("===========================================================================\n")


    config_path = sys.argv[1]

    with open(config_path, 'r') as fp:
        config = json.load(fp)

    if not os.path.exists('./htmls'):
        os.mkdir('./htmls')


    

    #SET RANGE DISPLAYED
    eta_min, eta_max, phi_min, phi_max = set_range(config)
    EventTree=ur.open(config["path_to_truth"])["Low_Tree"]
    EventNumber = config["number"]

    
    #if config["draw particles"]:
        #print("Drawing Truth Infos")
        
        #for cl in config["draw particles"]:
        #    exec("truth_{0} = truth_particles[truth_particles[\"Class Name\"] ==\"{0}\"]".format(cl))

    if config["draw calo"]:
        granularity = config["granularity"]
        calorimeter=get_calorimeter_info(EventTree ,EventNumber)
    if config["draw tracks"]:
        tracks=get_tracks_info(EventTree ,EventNumber)
    if config["draw particles"]:
        particles=get_truth_particles(EventTree ,EventNumber)
    if config["draw truthjets"]:
        truth_jets=get_jet_info(EventTree ,EventNumber, "true")
    if config["draw pflowjets"]:
        pflow_jets=get_jet_info(EventTree ,EventNumber, "pflow")
    
    #f=10
    if config["draw calo"]:
        fig = make_subplots(rows=3, cols=2, shared_xaxes=True, shared_yaxes=True, horizontal_spacing=0.117, vertical_spacing=0.02)
        for lay in range(6):
            #if lay == 0: continue
            draw_calo(fig, config, lay, calorimeter)
            if config["draw tracks"]:
                if(len(tracks)<1): continue
                draw_tracks(fig, lay, tracks)
            if config["draw particles"]:
                draw_particles(fig, config, lay, particles)
            if config["draw truthjets"]:
                draw_jets(fig, config, lay, truth_jets,"truth")
            if config["draw pflowjets"]:
                draw_jets(fig, config, lay, pflow_jets,"pflow")



        fig.update_layout(height=2500, width=2050)
        fig['layout']['yaxis']['title']='φ'
        fig['layout']['yaxis3']['title']='φ'
        fig['layout']['yaxis5']['title']='φ'
        fig['layout']['xaxis5']['title']='η'
        fig['layout']['xaxis6']['title']='η'

    else:
        fig = make_subplots(rows=1, cols=1, shared_xaxes=True, shared_yaxes=True)
        if config["draw particles"]:
                draw_particles(fig, config, 0, particles)
        if config["draw particles"]:
            draw_particles(fig, config, 0, particles)
        if config["draw truthjets"]:
            draw_jets(fig, config, 0, truth_jets,"truth")
        if config["draw pflowjets"]:
            draw_jets(fig, config, 0, pflow_jets,"pflow")
        fig.update_layout(height=900, width=1000)
        fig['layout']['yaxis']['title']='φ'
        fig['layout']['xaxis']['title']='η'

    fig.layout.plot_bgcolor = '#fff'    
    fig.layout.paper_bgcolor = '#fff'
    #fig.update_layout(legend= {'itemsizing': 'constant'})

    fig.update_xaxes(range=[eta_min, eta_max],title_font=dict(family="Times New Roman", size=60), tickfont_size=50,dtick = (eta_max-eta_min)/5)
    fig.update_yaxes(range=[phi_min, phi_max],title_font=dict(family="Times New Roman", size=60), tickfont_size=50,dtick = (phi_max-phi_min)/5)
    fig.update_xaxes(showgrid=True, gridcolor="grey", linewidth=2, zeroline=True, zerolinewidth=4, zerolinecolor='grey', showline=True, mirror=True, linecolor='black', )
    fig.update_yaxes(showgrid=True, gridcolor="grey", linewidth=2, zeroline=True, zerolinewidth=4, zerolinecolor='grey', showline=True, mirror=True, linecolor='black', )

    fig.update_layout(
        showlegend=False
    )



    if not config["name fig"]: 
        name = "figure"
    else:
        name = config["name fig"]

    if config["save html"]: 
        fig.write_html(name+".html")

    if config["save image"]: 
        fig.write_image(name+'.pdf', scale=8)
        fig.write_image(name+'.png', scale=8)
