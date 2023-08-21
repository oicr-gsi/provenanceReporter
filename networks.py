# -*- coding: utf-8 -*-
"""
Created on Wed Jun 14 16:29:36 2023

@author: rjovelin
"""


import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import io
import base64




def get_node_labels(block_workflows, workflow_names):
    '''
    (dict, dict) -> dict
    
    Returns a dictionary with node labels (workflows) for each block (sample pair)
    and bmpp parent workflows
    
    Parameters
    ----------
    - block_workflows (dict): Dictionary of workflow run ids organized by sample pair and bmpp parent workflows
    - workflow_names (dict): Dictionary with workflow name for each workflow run id in project
    '''
    
    labels = {}
    for block in block_workflows:
        labels[block] = {}
        for bmpp in block_workflows[block]:
            labels[block][bmpp] = []
            for workflow in block_workflows[block][bmpp]:
                workflow_name = workflow_names[workflow][0]
                workflow_name = workflow_name.split('_')[0]
                if workflow_name.lower() == 'varianteffectpredictor':
                    workflow_name = 'VEP'
                elif workflow_name.lower() == 'bammergepreprocessing':
                    workflow_name = 'bmpp'
                labels[block][bmpp].append(workflow_name)
       
    return labels


def make_adjacency_matrix(block_workflows, parent_workflows):
    '''
    (dict, dict) -> dict    
    
    Returns a dictionary of matrices showing relationships among workflows for each
    sample and bmpp parent workflows
    
    Parameters
    ----------
    - block_workflows (dict): Dictionary of workflow run ids organized by sample pair and bmpp parent workflows
    - parent_workflows (dict): Dictionary with parent-children workflow relationships
    '''
          
    matrix = {}
    
    for block in block_workflows:
        matrix[block] = {}
        for bmpp in block_workflows[block]:
            M = []
            for i in block_workflows[block][bmpp]:
                m = []
                for j in block_workflows[block][bmpp]:
                    if i not in parent_workflows and j not in parent_workflows:
                        m.append(0)
                    elif i == j:
                        m.append(0)
                    elif i in parent_workflows:
                        if j in parent_workflows[i]:
                            m.append(1)
                        else:
                            m.append(0)
                    elif i not in parent_workflows:
                        if i in parent_workflows[j]:
                            m.append(1)
                        else:
                            m.append(0)
                M.append(m)
            matrix[block][bmpp] = M

    return matrix




def convert_figure_to_base64(figure):
    '''
    (matplotlib.figure.Figure) -> str
    
    Returns a base64 string representation of figure
    
    Parameters
    ----------
    - figure (matplotlib.figure.Figure): Figure showing workflow relationships
    '''
    
    my_stringIObytes = io.BytesIO()
    figure.savefig(my_stringIObytes, format='png')
    
    my_stringIObytes.seek(0)
    my_base64_jpgData = base64.b64encode(my_stringIObytes.read()).decode('utf-8')
      
    return my_base64_jpgData



def show_graph(adjacency_matrix, mylabels):
    '''
    (dict, dict) -> matplotlib.figure.Figure
    
    Returns a matplotlib figure of workflow relationships
        
    Parameters
    ----------
    - matrix (dict): Dictionary of matrices of workflow relationships for each
                     sample pair (block) and bmpp parent workflows
    - labels (dict): Dictionary with node labels (workflows) for each block (sample pair)
                     and bmpp parent workflows
    '''
    



    figure = plt.figure()
    #figure.set_size_inches(2.5, 2)

    # add a plot to figure (N row, N column, plot N)
    ax = figure.add_subplot(1, 1, 1)



    
    rows, cols = np.where(adjacency_matrix == 1)
    edges = zip(rows.tolist(), cols.tolist())
    gr = nx.Graph()
    
    gr.add_edges_from(edges)
    
    #nx.draw(gr, node_size=500, labels=mylabels, with_labels=True)
    
    nodes = list(gr)
    N = {}
    for i in nodes:
        N[i] = mylabels[i]
    
    nx.draw(gr, node_size=1200, node_color='#ffe066', font_size = 14, with_labels=True, labels=N, linewidths=2)
    

    # write title
    #ax.set_title(title, size = 14)
        
    # # add space between axis and tick labels
    # ax.yaxis.labelpad = 18
    # ax.xaxis.labelpad = 18
    
    # # do not show lines around figure  
    # ax.spines["top"].set_visible(False)    
    # ax.spines["bottom"].set_visible(False)    
    # ax.spines["right"].set_visible(False)    
    # ax.spines["left"].set_visible(False)  
    
    # # do not show ticks
    # plt.tick_params(axis='both', which='both', bottom=False, top=False,
    #                 right=False, left=False, labelleft=False, labelbottom=False,
    #                 labelright=False, colors = 'black', labelsize = 12, direction = 'out')  
    
    # # set up same network layout for all drawings
    # Pos = nx.spring_layout(G)
    # # draw edges    
    # nx.draw_networkx_edges(gr, pos=1, width=0.7, edge_color='grey', style='solid',
    #                         alpha=0.4, ax=ax, arrows=False, node_size=5,
    #                         nodelist=AllNodes, node_shape='o')
    
    #nx.draw_networkx_edges(gr, pos=nx.spring_layout(gr), width=0.7, edge_color='grey', style='solid',
    #                        alpha=0.4, ax=ax, arrows=False, node_size=5,
    #                        node_shape='o')
    
    # draw all nodes, color according to degree
    # nodelist = sorted(degree.keys())
    # node_color = [degree[i] for i in nodelist]
    
   
    # nodes = nx.draw_networkx_nodes(G, pos=Pos, with_labels=False, node_size=5,
    #                                node_color=node_color, node_shape='o', alpha=0.3,
    #                                linewidths=0, edgecolors='grey', ax=None,
    #                                nodelist=nodelist, cmap=cmap)
    # nodes.set_clim(min(node_color), max(node_color)+1) 

    # # add discrete color bar for node degree
    # divider = make_axes_locatable(ax)
    # cax = divider.append_axes("bottom", size="5%", pad=0.05)
    

    # save figure    
    plt.tight_layout()
    #figure.savefig(Outputfile, bbox_inches = 'tight')
    plt.close()

    return figure




def plot_workflow_network(matrix, labels):
    '''
    (dict, dict) -> dict
    
    Returns a dictionary with plots of workflow relationships for each sample pair and 
    bmpp parent workflows
    
    Parameters
    ----------
    - matrix (dict): Dictionary of matrices of workflow relationships for each
                     sample pair (block) and bmpp parent workflows
    - labels (dict): Dictionary with node labels (workflows) for each block (sample pair)
                     and bmpp parent workflows
    '''
        
        
    F = {}
    # convert to numpy 2-D array
    for block in matrix:
        F[block] = {}
        for bmpp in matrix[block]:
            matrix[block][bmpp] = np.array(matrix[block][bmpp])
            figure = show_graph(matrix[block][bmpp], mylabels=labels[block][bmpp])   
            # convert to base64
            F[block][bmpp] = convert_figure_to_base64(figure)
   
    return F


