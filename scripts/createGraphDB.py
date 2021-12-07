import networkx
import pandas as pd
from itertools import combinations
import numpy as np
import math
import json
from networkx import algorithms as netalg
def save_networks(net,species):
    networkx.readwrite.graphml.write_graphml(net,f'./databases/combined_graphdbs/{species}/complex_tf_target_{species}.xml')
    networkx.readwrite.gml.write_gml(net,f'./databases/combined_graphdbs/{species}/complex_tf_target_{species}.gml')
    jsondata=networkx.readwrite.json_graph.node_link_data(net)
    with open(f'./databases/combined_graphdbs/{species}/complex_tf_target_{species}.json', 'w') as outfile:
        json.dump(jsondata, outfile)

def get_network_info(network):
    interaction_types=[]
    node_types=[]
    for n in network.nodes.keys():
        node_type='protein'
        if network.nodes[n].get('is_TF'):node_type='transcription factor'
        elif network.nodes[n].get('is_CTF'):node_type='transcription cofactor'
        else:pass
        node_types.append(node_type)
    for i1,i2 in network.edges.keys():
        interaction_types.append(network[i1][i2]['interaction'])
    return node_types,interaction_types
def calculate_graph_metrics(graph):
            dc=np.median(list(netalg.centrality.degree_centrality(graph).values()))
            tv=netalg.cluster.transitivity(graph)
            avgc=netalg.cluster.average_clustering(graph)
            avgcc=np.median(list(netalg.centrality.closeness_centrality(graph).values()))
            ac=netalg.assortativity.degree_assortativity_coefficient(graph)
            degree=round(graph.number_of_edges()/graph.number_of_nodes(),2)
            #sworld=netalg.smallworld.omega(graph,niter=10,nrand=3)
            return round(dc,2),round(tv,2),round(avgc,2),round(avgcc,2),round(ac,2),degree
def add_corum(network,complexes,complex_names):
    idx=0
    for c in complexes:
        c=str.split(c,';')
        if len(c)==1:
            network.add_edge(c[0],c[0],complex=complex_names[idx],interaction='complex')
        else:
            for prot_pairs in combinations(c,2):
                if network.has_edge(prot_pairs[0],prot_pairs[1]):
                    network[prot_pairs[0]][prot_pairs[1]].update({'complex':network[prot_pairs[0]][prot_pairs[1]]['complex']+';'+complex_names[idx]})
                    network[prot_pairs[1]][prot_pairs[0]].update({'complex':network[prot_pairs[1]][prot_pairs[0]]['complex']+';'+complex_names[idx]})
                else:
                    network.add_edge(prot_pairs[0],prot_pairs[1],complex=complex_names[idx],interaction='complex')
                    network.add_edge(prot_pairs[1],prot_pairs[0],complex=complex_names[idx],interaction='complex')
        idx+=1
    return network
def add_complexome(network,complexes,complex_names):
    idx=0
    for c in complexes:
        c=str.split(c,'|')
        for i in range(len(c)):
            c[i]=c[i][:-3]
        if len(c)==1:
            network.add_edge(c[0],c[0],complex=complex_names[idx],interaction='complex')
        else:
            for prot_pairs in combinations(c,2):
                if network.has_edge(prot_pairs[0],prot_pairs[1]):
                    network[prot_pairs[0]][prot_pairs[1]].update({'complex':network[prot_pairs[0]][prot_pairs[1]]['complex']+';'+complex_names[idx]})
                    network[prot_pairs[1]][prot_pairs[0]].update({'complex':network[prot_pairs[1]][prot_pairs[0]]['complex']+';'+complex_names[idx]})
                else:
                    network.add_edge(prot_pairs[0],prot_pairs[1],complex=complex_names[idx],interaction='complex')
                    network.add_edge(prot_pairs[1],prot_pairs[0],complex=complex_names[idx],interaction='complex')

        idx+=1
    return network
def add_trrust(network,tf,target,attribute):
    for t,t2,a in zip(tf,target,attribute):
        network.add_edge(t,t2,interaction='regulation',regulation_type=a)
        network.nodes[t]['is_TF']=True
    return network
def add_yeastract(network,tf,target,target_gname):
    for t,t2,t3 in zip(tf,target,target_gname):
        if not isinstance(t2,str):
            t2 = t3
        network.add_edge(t,t2,interaction='regulation')
        network.nodes[t]['is_TF']=True
    return network

def add_hftarget(network,tf,target,attribute):
    for t,t2,a in zip(tf,target,attribute):
        network.add_edge(t,t2,interaction='regulation',tissue_type=a)
        network.nodes[t]['is_TF']=True
    return network
def add_animal_tfdb3_info(network,tf,attribute,ctf=False):
    for t,a in zip(tf,attribute):
        if t in network.nodes:
            if ctf:
                network.nodes[t]['is_CTF']=True
            else:
                network.nodes[t]['is_TF']=True
            network.nodes[t]['family']=a
    return network



def createDB(species='human'):
    if species=='human':
        net=networkx.DiGraph()
        corum=pd.read_csv('./databases/corum/complex.txt',delimiter='\t')
        corum=corum[corum['Organism']=='Human']
        complexes_names=corum['ComplexName'].values
        complexes=corum['subunits(UniProt IDs)'].values
        net=add_corum(net,complexes,complexes_names)
        #print(net.number_of_edges()/2,net.number_of_nodes())
        complexome=pd.read_csv('./databases/complexome/homo_sapiens_complex.tsv',delimiter='\t')
        complexes_names=complexome['Recommended name'].values
        complexes=complexome['Identifiers (and stoichiometry) of molecules in complex'].values
        net=add_complexome(net,complexes,complexes_names)
        #print(net.number_of_edges()/2,net.number_of_nodes())
        trrust=pd.read_csv('./databases/trrust/homo_sapiens_tf_genetarget_processed.csv',delimiter=',').dropna()
        net=add_trrust(net, trrust['Uniprot TF'].values,trrust['Uniprot Target'].values,trrust['Interaction'].values)
        #print(net.number_of_edges()/2,net.number_of_nodes())
        tftarget=pd.read_csv('./databases/tftarget/homo_sapiens_tf_genetarget_processed.csv',delimiter=',').dropna()
        net=add_hftarget(net, tftarget['Uniprot TF'].values,tftarget['Uniprot Target'].values,tftarget['Tissue'].values)
        #print(net.number_of_edges()/2,net.number_of_nodes())
        animaltfdb3=pd.read_csv('./databases/animal_tfdb3/homo_sapiens_tf_processed.csv',delimiter=',').dropna()
        net=add_animal_tfdb3_info(net,animaltfdb3['Uniprot'].values,animaltfdb3['Family'].values)
        animaltfdb3=pd.read_csv('./databases/animal_tfdb3/homo_sapiens_ctf_processed.csv',delimiter=',').dropna()
        net=add_animal_tfdb3_info(net,animaltfdb3['Uniprot'].values,animaltfdb3['Family'].values,ctf=True)
        #networkx.readwrite.graphml.write_graphml(net,'./results/graphml')
        #networkx.readwrite.gml.write_gml(net,'./results/gml')
    elif species=='mouse':
        net=networkx.DiGraph()
        corum=pd.read_csv('./databases/corum/complex.txt',delimiter='\t')
        corum=corum[corum['Organism']=='Mouse']
        complexes_names=corum['ComplexName'].values
        complexes=corum['subunits(UniProt IDs)'].values
        net=add_corum(net,complexes,complexes_names)
        #print(net.number_of_edges()/2,net.number_of_nodes())
        complexome=pd.read_csv('./databases/complexome/mouse_complex.tsv',delimiter='\t')
        complexes_names=complexome['Recommended name'].values
        complexes=complexome['Identifiers (and stoichiometry) of molecules in complex'].values
        net=add_complexome(net,complexes,complexes_names)
        #print(net.number_of_edges()/2,net.number_of_nodes())
        trrust=pd.read_csv('./databases/trrust/mouse_tf_genetarget_processed.csv',delimiter=',').dropna()
        net=add_trrust(net, trrust['Uniprot TF'].values,trrust['Uniprot Target'].values,trrust['Interaction'].values)
        #print(net.number_of_edges()/2,net.number_of_nodes())
        animaltfdb3=pd.read_csv('./databases/animal_tfdb3/mouse_tf_processed.csv',delimiter=',').dropna()
        net=add_animal_tfdb3_info(net,animaltfdb3['Uniprot'].values,animaltfdb3['Family'].values)
        animaltfdb3=pd.read_csv('./databases/animal_tfdb3/mouse_ctf_processed.csv',delimiter=',').dropna()
        net=add_animal_tfdb3_info(net,animaltfdb3['Uniprot'].values,animaltfdb3['Family'].values,ctf=True)

    elif species=='yeast':
        net=networkx.DiGraph()
        complexome=pd.read_csv('./databases/complexome/yeast_complex.tsv',delimiter='\t')
        complexes_names=complexome['Recommended name'].values
        complexes=complexome['Identifiers (and stoichiometry) of molecules in complex'].values
        net=add_complexome(net,complexes,complexes_names)
        #print(net.number_of_edges()/2,net.number_of_nodes())
        yt=pd.read_csv('./databases/yeastract/yeast_tf_genetarget_processed.csv',delimiter=',')
        net=add_yeastract(net, yt['Uniprot TF'].values,yt['Uniprot Target'].values,yt['Target'].values)
        #print(net.number_of_edges()/2,net.number_of_nodes())
    return net
    



if __name__ == '__main__':
    all_species=['yeast','mouse','human']
    for species in all_species:
        print(species.upper())
        network=createDB(species)
        nodes,edges=get_network_info(network)
        node_types,node_counts=np.unique(nodes,return_counts=True)
        edge_types,edge_counts=np.unique(edges,return_counts=True)
        try:
            complex_edge_index=list(edge_types).index('complex')
            edge_counts[complex_edge_index]=int(edge_counts[complex_edge_index]/2)
        except:pass
        print(node_types,node_counts)
        print(edge_types,edge_counts)
        print(f'Nodes: {np.sum(node_counts)}\t Edges: {np.sum(edge_counts)}')
        print(calculate_graph_metrics(network))
        #save_networks(network,species)