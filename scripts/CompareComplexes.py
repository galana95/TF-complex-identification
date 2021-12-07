import numpy as np
import scipy.stats as stats
from scripts.ProteinComplexSimulation import ProteinComplexes
from scripts.utils import is_equal
import logging
#from memory_profiler import profile

class CompareProteinComplexes():
    def __init__(self,input1=None,input2=None,settings=None):
        self.theme=settings['REPORTER VISUALIZER SETTINGS']['THEME'].lower()
        if self.theme=='dark':
            self.bgcolor='#121212'
            self.fontcolor='white'
        elif self.theme=='light':
            self.bgcolor='#FFFFFF'
            self.fontcolor='black'
        else:
            self.bgcolor='#121212'
            self.fontcolor='white'
        if is_equal(settings['DATA TYPE'],'protein'):
            self.protein_data_only=True
        elif is_equal(settings['DATA TYPE'],'complex'):
            self.protein_data_only=False
        else:
            logging.critical('Unsupported data type')
        if self.protein_data_only:
            input1_prot_names=list(input1.iloc[:,0].values)
            input2_prot_names=list(input2.iloc[:,0].values)
            try:
                input1_abundances=list(input1['Simulation abundance'].values)
                input2_abundances=list(input2['Simulation abundance'].values)
            except:
                input1_abundances=list(input1.iloc[:,1].values)
                input2_abundances=list(input2.iloc[:,1].values)
            self.add_protein_only_data(input1_prot_names,input2_prot_names,input1_abundances,input2_abundances)
            logging.info(f'N1:{len(input1_prot_names)}, N2:{len(input2_prot_names)}')
        else:
            cp1=ProteinComplexes(input1)
            cp2=ProteinComplexes(input2)
            self.complexproteins1 =cp1.simulation_data['List of proteins'].values
            self.complexproteins2 =cp2.simulation_data['List of proteins'].values
            self.complexabundance1 =cp1.simulation_data['Average abundance'].values
            self.complexabundance2 =cp2.simulation_data['Average abundance'].values
            num_comps=len(self.complexproteins2)*len(self.complexproteins1)
            logging.info(f'N1:{len(self.complexproteins1)}, N2:{len(self.complexproteins2)}')
            if num_comps>1000*1000:
                logging.info('Lot of comparisons, expect slower performance')
            if num_comps>10000*1000:
                logging.info('Lot of comparisons, expect very slow performance')
            if num_comps>10000*10000:
                logging.warning('Lot of comparisons, expect extremely slow performance with considerable memory usage')
    def run_comparisons(self):
        if self.protein_data_only:
            self.get_protein_only_numbers()
            self.calculate_basic_statistics(self.protein_abundance1,self.protein_abundance2)
            self.compare_protein_abundance()
        else:
            self.split_protein_complexes()
            self.calculate_all_similarities()
            n1,n2=self.get_protein_numbers_in_each_complex()
            self.calculate_basic_statistics(n1,n2)
            self.calculate_network_based_similarities()
            self.compare_complex_abundance()
            self.create_visualize_graph()


    def add_protein_only_data(self,input1,input2,abundance1,abundance2,settings=None):
        self.proteins1=input1
        self.proteins2=input2
        self.protein_abundance1=abundance1
        self.protein_abundance2=abundance2
    def complex_similarity_score(self,complexproteins1,complexproteins2):
        complex1=set(complexproteins1)
        complex2=set(complexproteins2)
        intersect=len(complex1.intersection(complex2))
        total_size=max(len(complex1),len(complex2))
        return (intersect/total_size)*100
    def assymetric_complex_similarity_score(self,complexproteins1,complexproteins2):
        complex1=set(complexproteins1)
        complex2=set(complexproteins2)
        not_in_complex2=len(complex1.difference(complex2))
        not_in=100-(not_in_complex2/len(complex1))*100
        return not_in
    def calc_protein_differences(self,complexproteins1,complexproteins2):
        complex1=set(complexproteins1)
        complex2=set(complexproteins2)
        intersect=len(complex1.intersection(complex2))
        total_size=max(len(complex1),len(complex2))
        return total_size-intersect
    def split_protein_complexes(self):
        for i in range(len(self.complexproteins1)):
            self.complexproteins1[i]=self.complexproteins1[i].split(",")
        for i in range(len(self.complexproteins2)):
            self.complexproteins2[i]=self.complexproteins2[i].split(",")
    
    def calculate_all_similarities(self,all=True):
        complex_similarities=np.zeros((len(self.complexproteins1),len(self.complexproteins2)))
        if all:
            assymetric_complex_similarities=np.zeros((len(self.complexproteins1),len(self.complexproteins2)))
            assymetric_complex_similarities2=np.zeros((len(self.complexproteins1),len(self.complexproteins2)))
        protein_differences=np.zeros((len(self.complexproteins1),len(self.complexproteins2)))
        for i in range(len(self.complexproteins1)):
            for j in range(len(self.complexproteins2)):
                if all:
                    assymetric_complex_similarities[i,j]=self.assymetric_complex_similarity_score(self.complexproteins1[i],self.complexproteins2[j])
                    assymetric_complex_similarities2[i,j]=self.assymetric_complex_similarity_score(self.complexproteins2[j],self.complexproteins1[i])
                complex_similarities[i,j]=self.complex_similarity_score(self.complexproteins1[i],self.complexproteins2[j])
                protein_differences[i,j]=self.calc_protein_differences(self.complexproteins1[i],self.complexproteins2[j])
        if all:
            self.assymetric_complex_similarities=assymetric_complex_similarities
            self.assymetric_complex_similarities2=assymetric_complex_similarities2
        else:
            self.assymetric_complex_similarities=np.nan
            assymetric_complex_similarities2=np.nan
        self.complexproteins_similarities=complex_similarities
        self.protein_differences=protein_differences
    def get_protein_numbers_in_each_complex(self):

        complex1_protein_nums=[]
        complex2_protein_nums=[]
        for i in range(len(self.complexproteins1)):
            complex1_protein_nums.append(len(self.complexproteins1[i]))
        for i in range(len(self.complexproteins2)):
            complex2_protein_nums.append(len(self.complexproteins2[i]))
        return complex1_protein_nums,complex2_protein_nums
    def get_protein_only_numbers(self):
        self.protein_only_numbers1=len(set(self.proteins1))
        self.protein_only_numbers2=len(set(self.proteins2))
        return len(set(self.proteins1)),len(set(self.proteins2))
    def calculate_basic_statistics(self,complex1_protein_nums,complex2_protein_nums):
        self.input1_statistics={'min':round(np.min(complex1_protein_nums),1),'max':round(np.max(complex1_protein_nums),1),'mean':round(np.mean(complex1_protein_nums),1),
        'median':round(np.median(complex1_protein_nums),1),'iqr':round(stats.iqr(complex1_protein_nums),1),'standard_deviation':round(np.std(complex1_protein_nums),1),
        'skewness':round(stats.skew(complex1_protein_nums),1),'kurtosis':round(stats.kurtosis(complex1_protein_nums),1)}
        self.input2_statistics={'min':round(np.min(complex2_protein_nums),1),'max':round(np.max(complex2_protein_nums),1),'mean':round(np.mean(complex2_protein_nums),1),
        'median':round(np.median(complex2_protein_nums),1),'iqr':round(stats.iqr(complex2_protein_nums),1),'standard_deviation':round(np.std(complex2_protein_nums),1),
        'skewness':round(stats.skew(complex2_protein_nums),1),'kurtosis':round(stats.kurtosis(complex2_protein_nums),1)}

    
    def calculate_network_based_similarities(self):
        import networkx as nx
        from itertools import combinations
        from networkx import algorithms as netalg
        def calculate_graph_metrics(graph):
            dc=np.mean(list(netalg.centrality.degree_centrality(graph).values()))
            tv=netalg.cluster.transitivity(graph)
            avgc=netalg.cluster.average_clustering(graph)
            avgcc=np.mean(list(netalg.centrality.closeness_centrality(graph).values()))
            ac=netalg.assortativity.degree_assortativity_coefficient(graph)
            #sworld=netalg.smallworld.omega(graph,niter=10,nrand=3)
            return round(dc,2),round(tv,2),round(avgc,2),round(avgcc,2),round(ac,2)
        G1 = nx.Graph()
        G2 = nx.Graph()
        #all_proteins=[]
        # for complex in self.complexproteins2:
        #     for protein in complex:
        #         if protein not in all_proteins:
        #             protein.append(all_proteins)
        # for protein in all_proteins:
        #     net.add_node(protein, label=protein)
        #     cpglens.append(len(cpg))

        for complex in self.complexproteins1:
            if len(complex)<2:continue
            for prot_pairs in combinations(complex,2):
                G1.add_edge(prot_pairs[0],prot_pairs[1], weight=1)
        for complex in self.complexproteins2:
            if len(complex)<2:continue
            for prot_pairs in combinations(complex,2):
                G2.add_edge(prot_pairs[0],prot_pairs[1], weight=1)
        gm2=calculate_graph_metrics(G2)
        gm1=calculate_graph_metrics(G1)
        #similarity=netalg.similarity.graph_edit_distance(G1,G2,timeout=30)
        #try:similarity=round(similarity,2)
        #except:pass
        similarity=np.nan
        self.graph1={'graph':G1,'similarity':similarity,'degree_centrality':gm1[0],'transitivity':gm1[1],'average_clustering':gm1[2],
        'closeness_centrality':gm1[3],'assortavity':gm1[4]}
        self.graph2={'graph':G2,'similarity':similarity,'degree_centrality':gm2[0],'transitivity':gm2[1],'average_clustering':gm2[2],
        'closeness_centrality':gm2[3],'assortavity':gm2[4]}

        return G1,G2
    
    def create_visualize_graph(self,keep_only_proteins=50):
        from pyvis.network import Network
        from itertools import combinations
        import random
        net = Network(height='750px', width='100%',  bgcolor=self.bgcolor, font_color=self.fontcolor,heading="Proteins of complexes with significant abundance differences. Red & Blue -> proteins present only in one condition, Green -> proteins present in both.")
    
        divider=len(self.complexproteins2)+len(self.complexproteins1)
        if len(self.complexproteins1)>keep_only_proteins:
            divider=keep_only_proteins
            #idxscp1=random.sample(range(0, len(self.complexproteins1)-1), keep_only_proteins)
            args=np.flip(np.argsort(self.complexabundance1))
            args=args[:keep_only_proteins]
            temp_complexproteins1=self.complexproteins1[args]
            temp_complexabundance1=self.complexabundance1[args]
        else:
            temp_complexproteins1=self.complexproteins1
            temp_complexabundance1=self.complexabundance1
        if len(self.complexproteins2)>keep_only_proteins:
            #idxscp2=random.sample(range(0, len(self.complexproteins2)-1), keep_only_proteins)
            args=np.flip(np.argsort(self.complexabundance2))
            args=args[:keep_only_proteins]
            temp_complexproteins2=self.complexproteins2[args]
            temp_complexabundance2=self.complexabundance2[args]
        else:
            temp_complexproteins2=self.complexproteins2
            temp_complexabundance2=self.complexabundance2
        all_proteins=[]
        all_abundances={}
        all_numbers={}
        cp1_proteins=[]
        cp2_proteins = []
        idx=0
        for complex in temp_complexproteins2:
            for protein in complex:
                if temp_complexabundance2[idx]<1:
                    abund=0
                else:
                    abund=np.log2(temp_complexabundance2[idx])
                try:
                    all_abundances[protein]+=abund
                    all_numbers[protein]+=1
                except:
                    all_abundances[protein]=abund
                    all_numbers[protein]=1
                if protein not in cp2_proteins:
                    cp2_proteins.append(protein)
                if protein not in all_proteins:
                    all_proteins.append(protein)
            idx+=1
        idx=0
        for complex in temp_complexproteins1:
            for protein in complex:
                if temp_complexabundance1[idx]<1:
                    abund=0
                else:
                    abund=np.log2(temp_complexabundance1[idx])
                try:
                    all_abundances[protein]+=abund
                    all_numbers[protein]+=1
                except:
                    all_abundances[protein]=abund
                    all_numbers[protein]=1
                if protein not in cp1_proteins:
                    cp1_proteins.append(protein)
                if protein not in all_proteins:
                    all_proteins.append(protein)
            idx+=1

        color=''
        for protein in all_proteins:
            if protein in cp1_proteins and protein in cp2_proteins:
                color='rgb(5,255,50)'
            elif protein in cp2_proteins:
                color='rgb(255,50,50)'
            elif protein in cp1_proteins:
                color='rgb(10,60,255)'
            size=all_abundances[protein]
            abundvals=np.asarray(list(all_abundances.values()))
            abundvals=abundvals[abundvals>0]
            if size==0:size=np.min(abundvals)/2
            net.add_node(protein, label=protein,color=color,size=((size/all_numbers[protein])+size)/2+10)

            #cpglens.append(len(cpg))

        red = list(np.random.random(size=keep_only_proteins*2) * 256)
        green = list(np.random.random(size=keep_only_proteins*2) * 256) 
        blue = list(np.random.random(size=keep_only_proteins*2) * 256) 
        idx=0
        for complex in temp_complexproteins1:
            if len(complex)<2:continue
            for prot_pairs in combinations(complex,2):
                    net.add_edge(prot_pairs[0],prot_pairs[1],color=f'rgb({red[idx]},{green[idx]},{blue[idx]})',width=4)
            idx+=1
        #idx=0
        for complex in temp_complexproteins2:
            if len(complex)<2:continue
            for prot_pairs in combinations(complex,2):
                    net.add_edge(prot_pairs[0],prot_pairs[1],color=f'rgb({red[idx]},{green[idx]},{blue[idx]})',width=4)
            idx+=1
        self.network_plot=net
    def compare_protein_abundance(self):
        proteins1=self.proteins1
        proteins2=self.proteins2
        abund1=self.protein_abundance1
        abund2=self.protein_abundance2
        abundance1=[]
        abundance2=[]
        all_proteins_list=[]
        for i in range(len(proteins1)):
            if proteins1[i] not in all_proteins_list:
                all_proteins_list.append(proteins1[i])
                abundance1.append(abund1[i])
                try:
                    idx=proteins2.index(proteins1[i])
                    abundance2.append(abund2[idx])
                except:
                    abundance2.append(0)
        for i in range(len(proteins2)):
            if proteins2[i] not in all_proteins_list:
                all_proteins_list.append(proteins2[i])
                abundance2.append(abund2[i])
                try:
                    idx=proteins1.index(proteins2[i])
                    abundance1.append(abund1[idx])
                except:
                    abundance1.append(0)
        diff_abundance=np.abs(np.asarray(abundance1)-np.asarray(abundance2))
        args=np.flip(np.argsort(diff_abundance))
        all_proteins_list_ordered=np.asarray(all_proteins_list)[args]
        diff_abundance_ordered=diff_abundance[args]
        diff_abundance_ordered2=np.asarray(abundance1)-np.asarray(abundance2)
        diff_abundance_ordered2=diff_abundance_ordered2[args]
        self.protein_only_all_proteins_list_ordered=all_proteins_list_ordered
        self.protein_only_diff_abundance_ordered=diff_abundance_ordered2
        return all_proteins_list_ordered,diff_abundance_ordered2
    def is_in_complex_list(self,complex,complex_list):
        complex=set(complex)
        for complex2 in complex_list:
            complex2=set(complex2)
            if complex==complex2:
                return True
        return False
    def find_in_complex_list(self,complex,complex_list):
        complex=set(complex)
        idx=0
        for complex2 in complex_list:
            complex2=set(complex2)
            if complex==complex2:
                return int(idx)
            idx+=1
        return np.nan
    
    def compare_complex_abundance(self):
        complex1=self.complexproteins1
        complex2=self.complexproteins2
        abund1=self.complexabundance1
        abund2=self.complexabundance2

        abundance1=[]
        abundance2=[]
        all_proteins_list=[]
        for i in range(len(complex1)):
            if not self.is_in_complex_list(complex1[i],all_proteins_list):
                all_proteins_list.append(complex1[i])
                abundance1.append(abund1[i])
                try:
                    idx=self.find_in_complex_list(complex1[i],complex2)
                    abundance2.append(abund2[idx])
                except:
                    abundance2.append(0)
        for i in range(len(complex2)):
            if not self.is_in_complex_list(complex2[i],all_proteins_list):
                all_proteins_list.append(complex2[i])
                abundance2.append(abund2[i])
                try:
                    idx=self.find_in_complex_list(complex2[i],complex1)
                    abundance1.append(abund1[idx])
                except:
                    abundance1.append(0)

        diff_abundance=np.abs(np.asarray(abundance1)-np.asarray(abundance2))
        args=np.flip(np.argsort(diff_abundance))
        all_proteins_list_ordered=np.asarray(all_proteins_list)[args]
        diff_abundance_ordered=diff_abundance[args]
        diff_abundance_ordered2=np.asarray(abundance1)-np.asarray(abundance2)
        #print(max(diff_abundance_ordered2),min(diff_abundance_ordered2),np.mean(diff_abundance_ordered2))
        diff_abundance_ordered2=diff_abundance_ordered2[args]
        self.protein_only_all_proteins_list_ordered=all_proteins_list_ordered
        self.protein_only_diff_abundance_ordered=diff_abundance_ordered2

        return all_proteins_list_ordered,diff_abundance_ordered2
        

