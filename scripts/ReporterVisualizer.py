import plotly.figure_factory as ff
import plotly.graph_objects as go
import numpy as np
from scripts.utils import is_equal,get_folder_name
from datetime import datetime
import os
import yaml
import logging
import pandas as pd
import requests

logging.getLogger("requests").setLevel(logging.WARNING)
logging.getLogger("urllib3").setLevel(logging.WARNING)

class ReporterVisualizer():
    def __init__(self,compared_complex=None,settings=None,input1=None,input2=None):
        self.compared_complex = compared_complex
        self.input1 = input1
        self.input2 = input2
        self.input1_name=settings['INPUT1']['NAME']
        self.input2_name=settings['INPUT2']['NAME']
        self.settings=settings
        self.complexes={}
        self.use_log=settings['REPORTER VISUALIZER SETTINGS']['USE LOG']
        self.theme=settings['REPORTER VISUALIZER SETTINGS']['THEME'].lower()
        self.savefigs=settings['REPORTER VISUALIZER SETTINGS']['SAVE FIGURES']
        self.figure_counter=1
        if self.theme=='dark':
            self.bgcolor='#121212'
            self.fontcolor='white'
            self.plotlytheme='plotly_dark'
        elif self.theme=='light':
            self.bgcolor='#FFFFFF'
            self.fontcolor='black'
            self.plotlytheme='simple_white'
        else:
            self.bgcolor='#121212'
            self.fontcolor='white'
            self.plotlytheme='plotly_dark'
        self.netdb=None
        if is_equal(settings['DATA TYPE'],'protein'):
            self.protein_only=True
        elif is_equal(settings['DATA TYPE'],'complex'):
            self.protein_only=False
        else:
            logging.critical('Unsupported data type')

        if self.protein_only:
            self.base_elem='protein'
        else:
            self.base_elem='complex'
    def write(self,figure,height='90%'):
        if self.savefigs:
            figure.write_image(f'{self.directory}/figures/{str(self.figure_counter)}_{figure.layout.title["text"]}.png',scale=4)
            figure.write_image(f'{self.directory}/figures/{str(self.figure_counter)}_{figure.layout.title["text"]}.svg')
            self.figure_counter += 1
        figure.update_layout(autosize=True,width=None)
        self.dashboard.write(figure.to_html(full_html=False, include_plotlyjs='cdn',default_height=height))
    def write_string(self,string):
        string='<div>'+string+'</div>'
        self.dashboard.write(string)

    def add_background_color(self):

        string='''<style>
        body {
        background-color: '''+self.bgcolor+''';
        }
        .center {
        text-align: center;
        color: '''+self.fontcolor+''';
        }
        .center2 {
        text-align: center;
        color: '''+self.fontcolor+''';
        font-size: large;
        }
        .main{
        color: '''+self.fontcolor+''';
        font-weight: bold;
        margin-left: 2%;
        margin-right: 2%;
        }
        .subTF{
        color: DeepSkyBlue;
        font-weight: bold;
        margin-left: 2%;
        margin-right: 2%;
        }
        .subCTF{
        color: BlueViolet;
        font-weight: bold;
        margin-left: 2%;
        margin-right: 2%;
        }
        .sub{
            color:'''+self.fontcolor+''';
            font-weight: bold;
            margin-left: 2%;
            margin-right: 2%;
        }
        .sub2{
            color:'''+self.fontcolor+''';
            margin-left: 2%;
            margin-right: 2%;
        }
        </style>'''
        self.dashboard.write(string)
    def run_all(self):
        from plotly.subplots import make_subplots
        import networkx
        directory=get_folder_name(self.settings)
        self.directory=directory
        
        f = open(f'{directory}/settings.yaml', 'w+')
        yaml.dump(self.settings, f, allow_unicode=True)
        f.close()
        if not self.protein_only:
            self.compared_complex.network_plot.write_html(f'{directory}/complex_network.html')
        self.dashboard = open(f'{directory}/dashboard.html','w')
        species=self.settings['SPECIES'].lower()
        net=networkx.readwrite.gml.read_gml(f'./databases/combined_graphdbs/{species}/complex_tf_target_{species}.gml')
        self.netdb=net
        self.disease_gene_db=pd.read_csv('./databases/disease_gene/all_gene_disease_associations.tsv',delimiter='\t')
        self.drug_gene_db=pd.read_csv('./databases/drug_gene/drug_gene.csv',delimiter=',')
        self.add_background_color()
        self.add_complex_info()
        self.print_statistics()
        fig = make_subplots(rows=1, cols=2)
        self.plot_distributions(fig)
        self.plot_violin(fig)
        proteins=self.plot_abundance_difference()
        self.write_GO()
        proteins2=proteins
        while not self.create_target_subnetwork(proteins2):
            proteins2=proteins2[:-1]
        idx=1
        for p in proteins:
            self.create_target_subnetwork([p],idx=idx)
            idx+=1
            
        proteins2=proteins
        while not self.create_target_subnetwork(proteins2,centered=True):
            proteins2=proteins2[:-1]
        idx=1
        for p in proteins:
            self.create_target_subnetwork([p],idx=idx,centered=True)
            idx+=1       
        self.plot_graphs()
        self.dashboard.close()
    def create_target_subnetwork(self,proteins,idx=-1,centered=False):
        
        if not self.protein_only:
            proteins=proteins[:5]
            proteins_new=[]
            for p in proteins:
                tp=str.split(p,',')
                for p2 in tp:
                    proteins_new.append(p2)
            proteins=proteins_new
        else:
            proteins=proteins[:15]
        net=self.netdb
        nodes_to_include=list(proteins)
        def filter_subnetwork(network,node1,node2):
            exclude=True
            edge=network.get_edge_data(node1,node2)
            # if edge['interaction']=='regulation' and not self.settings['REPORTER VISUALIZER SETTINGS']['REGULATION']:
            #     if 'complex' in edge:pass
            #     else:return True

            # if edge['interaction']=='complex' and not self.settings['REPORTER VISUALIZER SETTINGS']['COMPLEX']:
            #     if 'regulation_type' in edge or 'tissue_type' in edge:pass
            #     else:return True

            if self.settings['SPECIES'].lower()!='human':return False
            if 'tissue_type' in edge:
                                    
                    tissue=str.split(edge['tissue_type'],',')
                    for t in self.settings['REPORTER VISUALIZER SETTINGS']['TISSUE']:
                        #print(t)
                        if t in tissue:
                            exclude=False
                            break
            elif 'complex' in edge:exclude=False
            return exclude
            

        def filter_network(nodes,degrees=5,predecessors=True):
            
                
            if degrees==0:
                return
            else:
                temp_nodes=[]
                for n in nodes:
                    try:
                        neighbors=net.successors(n)
                        for n2 in neighbors:
    
                            if filter_subnetwork(net,n,n2):continue
                            temp_nodes.append(n2)
                            nodes_to_include.append(n2)
                        if predecessors:
                            neighbors=net.predecessors(n)
                            for n2 in neighbors:
                                if filter_subnetwork(net,n,n2):continue
                                temp_nodes.append(n2)
                                nodes_to_include.append(n2)
                    except:pass
                if not temp_nodes:
                    return
                filter_network(temp_nodes,degrees-1)


        filter_network(proteins,1)
        nodes_to_include=np.unique(nodes_to_include)
        #if idx!=-1:
            #print(len(nodes_to_include))
        if len(nodes_to_include)>150 and len(proteins)>1:
            #print(len(nodes_to_include),len(proteins))
            return False
        #print(len(nodes_to_include),len(proteins))
        net=net.subgraph(nodes_to_include)
        from pyvis.network import Network
        
        nt = Network(height='750px', width='100%', bgcolor=self.bgcolor, font_color=self.fontcolor,directed=True,heading="Interactions (blue->complex, red->gene target) of proteins with significant abundance change (light green)")
        #nt.from_nx(net)
        for n in net.nodes:
            if n in proteins:
                color='rgb(20,250,40)'
            else:
                color='rgb(10,122,70)'
            nt.add_node(n,color=color)
        for e in net.edges:
            
            if centered and (e[0] not in proteins and e[1] not in proteins):
                continue
            if net[e[0]][e[1]]['interaction']=='complex' and self.settings['REPORTER VISUALIZER SETTINGS']['COMPLEX']:
                color='rgb(20,60,250)'
                label=net[e[0]][e[1]]['complex']
                nt.add_edge(e[0],e[1],color=color,title=label)
            elif net[e[0]][e[1]]['interaction']=='regulation' and self.settings['REPORTER VISUALIZER SETTINGS']['REGULATION']:
                color='rgb(240,50,50)'
                try:label=net[e[0]][e[1]]['regulation_type']
                except:label=''
                nt.add_edge(e[0],e[1],color=color,title=label)
        add=''
        if centered:add='_centered'
        if idx==-1:
            nt.write_html(f'{self.directory}/complex_target_subnetwork{add}.html')
        else:
            nt.write_html(f'{self.directory}/complex_target_subnetworks/no{str(idx)}{add}.html')
        return True
    def add_complex_info(self):
        corum=pd.read_csv('./databases/corum/complex.txt',delimiter='\t')
        corum=corum[corum['Organism']==self.settings['SPECIES'].lower().capitalize()]
        complexes_names=corum['ComplexName'].values
        complexes=corum['subunits(UniProt IDs)'].values
        complexes_new=[]
        for c in complexes:
            c=str.split(c,';')
            complexes_new.append(c)

        complexome=pd.read_csv(f'./databases/complexome/{self.settings["SPECIES"].lower()}_complex.tsv',delimiter='\t')
        complexes_names=np.hstack((complexes_names,complexome['Recommended name'].values))
        complexes=complexome['Identifiers (and stoichiometry) of molecules in complex'].values
        for c in complexes:
            c=str.split(c,'|')
            temp_c=[]
            for i in range(len(c)):
                c[i]=c[i][:-3]
                temp_c.append(c[i])
            complexes_new.append(temp_c)

        self.complexes['names']=complexes_names
        self.complexes['proteins']=complexes_new

    def get_best_matching_real_complex(self,complex):
        complex=str.split(complex,',')
        def complex_similarity_score(complexproteins1,complexproteins2):
            complex1=set(complexproteins1)
            complex2=set(complexproteins2)
            intersect=len(complex1.intersection(complex2))
            total_size=max(len(complex1),len(complex2))
            return (intersect/total_size)*100
        match_percent=[]
        
        for c in self.complexes['proteins']:
            match_percent.append(complex_similarity_score(c,complex))
        if np.max(match_percent)>=50:
            best_match=self.complexes['names'][np.argmax(match_percent)]
            match_b=round(np.max(match_percent),2)
        else:
            best_match='no match'
            match_b=''
        return best_match,match_b
    def get_uniprot_info(self,protein):
        upinfo={'name':'N/A','info':'N/A','disease':'N/A','gofunction':'N/I','genename':'N/A'}
        try:
            r = requests.get(f'https://www.uniprot.org/uniprot/{protein}.xml')
            data=r.text
            import xml.etree.ElementTree as ET
            root = ET.fromstring(data)
            root=root[0]
            for child in root:
                if 'protein' in child.tag:
                    upinfo['name']=child[0][0].text
                if 'gene' in child.tag:
                    upinfo['genename']=child[0].text
                if 'comment' in child.tag:
                    attrib=child.attrib
                    
                    if attrib['type']=='function':
                        upinfo['info']=child[0].text
                    elif attrib['type']=='disease':
                        upinfo['disease']=child[0].text
                        #print(child[0].attrib)
                        #print(child[0]['name'])
        except:
            pass
        return upinfo
    def write_GO(self):
        idx=1
        string='<div><h3 class="main">Descriptions of proteins/complexes with largest change in abundance (blue-> transcription factor, violet-> transcription cofactor)</h3>'
        for proteincomplex in self.most_interesting:
            best_match,match_percent=self.get_best_matching_real_complex(proteincomplex)
            if self.protein_only:
                bmstring=''
            else:
                bmstring=f'|| Best matching real complex: {best_match} | Match: {match_percent}'
            string+=f'<hr class="main"><h4 class="main">No.{str(idx)} {bmstring}</h4>'
            for p in proteincomplex.split(','):
                upinfo=self.get_uniprot_info(p)
                
                disease_gene_subset=self.disease_gene_db[self.disease_gene_db['geneSymbol']==upinfo['genename'].upper()]
                if len(disease_gene_subset)>0:
                    disease_gene_subset=disease_gene_subset.sort_values('score',axis=0,ascending=False)
                    upinfo['disease']=disease_gene_subset['diseaseName'].iloc[0]
                    
                dinfo='N/A'
                drug_gene_subset=self.drug_gene_db[self.drug_gene_db['GeneName']==upinfo['genename'].upper()]
                if len(drug_gene_subset)>0:
                    dinfo=' | '.join(list(drug_gene_subset['Drug']))
                pclass="sub"
                try:
                    
                    if self.netdb.nodes[p].get('is_TF'):
                        pclass="subTF"
                    elif self.netdb.nodes[p].get('is_CTF'):
                        pclass="subCTF"
                except:pass
                string+=f'''<p class="{pclass}">{p} {upinfo["name"]}</p><p class="sub2">
                Info: {upinfo["info"]} 
                <br>Disease: {upinfo["disease"]} 
                <br>Approved/Patented Drugs: {dinfo}
                <br>GO Function: {upinfo["gofunction"]}</p>'''
            idx+=1
        string+='</div>'
        self.dashboard.write(string)
    def plot_graphs(self):
        string1='<div class="center2"><a class="center2" href="./complex_target_subnetwork.html">Protein target subnetwork</a></div>'
        string2='<div class="center2"><a class="center2" href="./complex_network.html">Complex network comparison</a></div>'
        stringp='<p></p>'
        self.dashboard.write(stringp)
        self.dashboard.write(string1)
        self.dashboard.write(stringp)
        if not self.protein_only:
            self.dashboard.write(string2)
    def print_statistics(self):
        def calc_stat1(string):
            return round((self.compared_complex.graph1[string]/self.compared_complex.graph2[string])*100-100,1)
        def calc_stat2(string):
            return round((self.compared_complex.input1_statistics[string]/self.compared_complex.input2_statistics[string])*100-100,1)
        ps='<h3 class="center">'
        pe='</h3>'
        ps1='<h1 class="center">'
        pe1='</h2>'
        if self.protein_only:
            n1=self.compared_complex.protein_only_numbers1
            n2=self.compared_complex.protein_only_numbers2
        else:
            n1=len(self.compared_complex.complexproteins1)
            n2=len(self.compared_complex.complexproteins2)

        string=f'{ps1}Dashboard for comparison {str.replace(self.directory,"./results/","")}{pe1}'
        string=string+f'''{ps}<b>Number of filtered {self.input1_name}</b> {self.base_elem}(e)s: {n1}{pe}
        {ps}<b>Number of filtered {self.input2_name}</b>  {self.base_elem}(e)s: {n2}{pe}
        {ps}<b>Species: {self.settings['SPECIES'].lower().capitalize()}</b>{pe}
        {ps}<b>Filters: {' & '.join(self.settings['FILTER SETTINGS']['FILTER'])}</b>{pe}
        '''
        

        if not self.protein_only:
            fig=go.Figure()
            fig0=go.Figure()
            fig0 = go.Figure(data=[go.Table(header=dict(values=['Metrics','Median','Mean']),
                 cells=dict(values=[
                     ['Complex similarity','Assymetric complex similarity (treated in untreated)','Assymetric complex similarity (untreated in treated)','Protein number difference'],
                     [round(np.median(np.max(self.compared_complex.complexproteins_similarities,0)),1),round(np.median(np.max(self.compared_complex.assymetric_complex_similarities,0)),1),round(np.median(np.max(self.compared_complex.assymetric_complex_similarities2,0)),1),round(np.median(self.compared_complex.protein_differences),1)],
                      [round(np.mean(np.max(self.compared_complex.complexproteins_similarities,0)),1),round(np.mean(np.max(self.compared_complex.assymetric_complex_similarities,0)),1),round(np.mean(np.max(self.compared_complex.assymetric_complex_similarities2,0)),1),round(np.mean(self.compared_complex.protein_differences),1) ],
                      ]))
                     ])
            fig0.update_layout(template=self.plotlytheme,height=400,title='Complex similarity metrics between the two conditions')
            fig = go.Figure(data=[go.Table(header=dict(values=['Metrics',self.input1_name, self.input2_name, 'Difference (%)']),
                 cells=dict(values=[
                     ['Degree centrality','Transitivity', 'Average clustering', 'Closeness centrality','Assortativity', 'Similarity'],
                     [self.compared_complex.graph1["degree_centrality"], self.compared_complex.graph1["transitivity"], self.compared_complex.graph1["average_clustering"], self.compared_complex.graph1["closeness_centrality"],self.compared_complex.graph1["assortavity"],self.compared_complex.graph1["similarity"]],
                      [self.compared_complex.graph2["degree_centrality"], self.compared_complex.graph2["transitivity"], self.compared_complex.graph2["average_clustering"], self.compared_complex.graph2["closeness_centrality"],self.compared_complex.graph2["assortavity"],self.compared_complex.graph2["similarity"]],
                      [calc_stat1('degree_centrality'), calc_stat1('transitivity'), calc_stat1('average_clustering'), calc_stat1('closeness_centrality'),calc_stat1('assortavity'),0],
                      ]))
                     ])
            fig.update_layout(template=self.plotlytheme,height=400,title='Network comparison')
        fig2 = go.Figure(data=[go.Table(header=dict(values=['Metrics',self.input1_name, self.input2_name, 'Difference (%)']),
                 cells=dict(values=[
                     ['Max','Min', 'Mean', 'Median','STD', 'IQR','Kurtosis','Skewness'],
                      [self.compared_complex.input1_statistics["max"], self.compared_complex.input1_statistics["min"], self.compared_complex.input1_statistics["mean"], self.compared_complex.input1_statistics["median"],self.compared_complex.input1_statistics["standard_deviation"],self.compared_complex.input1_statistics["iqr"],self.compared_complex.input1_statistics["kurtosis"],self.compared_complex.input1_statistics["skewness"]],
                      [self.compared_complex.input2_statistics["max"], self.compared_complex.input2_statistics["min"], self.compared_complex.input2_statistics["mean"], self.compared_complex.input2_statistics["median"],self.compared_complex.input2_statistics["standard_deviation"],self.compared_complex.input2_statistics["iqr"],self.compared_complex.input2_statistics["kurtosis"],self.compared_complex.input2_statistics["skewness"]],
                      [calc_stat2('max'), calc_stat2('min'), calc_stat2('mean'), calc_stat2('median'),calc_stat2('standard_deviation'),calc_stat2('iqr'),calc_stat2('kurtosis'),calc_stat2('skewness')],
                      ]))
                     ])
        fig2.update_layout(template=self.plotlytheme,height=400,title='Abundance')
        self.write_string(string)
        if not self.protein_only:
            self.write(fig0,height='42%')
            self.write(fig,height='42%')
        self.write(fig2,height='42%')

    def plot_distributions(self,fig):
        if self.protein_only:
            data=[self.compared_complex.protein_abundance1,self.compared_complex.protein_abundance2]
        else:
            data=[self.compared_complex.complexabundance1,self.compared_complex.complexabundance2]
        if self.use_log:
            if not self.protein_only:
                a=np.log(data[0]+1)
                b=np.log(data[1]+1)
            else:
                a=np.log(data[0])
                b=np.log(data[1])      
            a[a==-np.inf]=0
            b[b==-np.inf]=0
            data=[a,b]

        figdist = ff.create_distplot(data, [self.input1_name,self.input2_name])
        figdist.update_layout(template=self.plotlytheme,title=f'Distribution of {self.base_elem} log abundances')
        self.write(figdist)
    def plot_violin(self,fig):
        if self.protein_only:
            data=[self.compared_complex.protein_abundance1,self.compared_complex.protein_abundance2]
        else:
            data=[self.compared_complex.complexabundance1,self.compared_complex.complexabundance2]
        if self.use_log:
            if not self.protein_only:
                a=np.log(data[0]+1)
                b=np.log(data[1]+1)
            else:
                a=np.log(data[0])
                b=np.log(data[1])   
            a[a==-np.inf]=0
            b[b==-np.inf]=0
            data=[a,b]

        figv = go.Figure()
        figv.add_trace(go.Violin(y=data[0],
                            x0=f'{self.input1_name}',
                            name=f'{self.input1_name}',
                            box_visible=True,
                            meanline_visible=True))
        figv.add_trace(go.Violin(y=data[1],
                            x0=f'{self.input2_name}',
                            name=f'{self.input2_name}',
                            box_visible=True,
                            meanline_visible=True))

        figv.update_layout(yaxis_zeroline=False,template=self.plotlytheme,title=f'Distribution of {self.base_elem} log abundances')
        self.write(figv)
    def plot_abundance_difference(self,limit=15):
        if not self.protein_only:
            ordered_proteincomplex_list=[]
            ordered_proteincomplex_list_full=[]
            for pc in self.compared_complex.protein_only_all_proteins_list_ordered:
                tempname=','.join(pc)
                ordered_proteincomplex_list_full.append(tempname)
                tempname=tempname[:40]
                ordered_proteincomplex_list.append(tempname)
                if len(ordered_proteincomplex_list)==limit:
                    break
        else:
            ordered_proteincomplex_list=self.compared_complex.protein_only_all_proteins_list_ordered[:limit]
            ordered_proteincomplex_list_full=self.compared_complex.protein_only_all_proteins_list_ordered[:limit]
        fig = go.Figure([go.Bar(x=ordered_proteincomplex_list, 
        y=self.compared_complex.protein_only_diff_abundance_ordered[:limit])])
        fig.update_layout(template=self.plotlytheme,title=f'Largest {self.base_elem} abundance differences')
        self.write(fig)
        self.most_interesting=ordered_proteincomplex_list_full
        return self.most_interesting
        #plt.bar(x=np.arange(len(proteins_tf)),height=abundance_tf,tick_label=proteins_tf)
#plt.ylabel('Abundance difference for TFs (treated-untreated)')


