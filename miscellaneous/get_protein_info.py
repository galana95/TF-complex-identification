import numpy as np
import networkx
import requests
import pickle
import json
species='human'
net=networkx.readwrite.gml.read_gml(f'./databases/combined_graphdbs/{species}/complex_tf_target_{species}.gml')
PROTEINS={}
idx=0
for protein in net.nodes:
    if idx%1000 == 0:
        print(idx)
    idx+=1
    try:
        r = requests.get(f'https://www.uniprot.org/uniprot/{protein}.xml')
        data=r.text
        import xml.etree.ElementTree as ET
        root = ET.fromstring(data)
        root=root[0]
        #root = tree.getroot()
        #root=root._children[0]
        valid_features=np.asarray(["domain","region of interest","modified residue","cross-link","compositionally biased region","zinc finger region"],dtype=str)
        PROTEIN=[]
        for child in root:
            if 'sequence' in child.tag:
                attrib['sequence']=list(child.text)
                PROTEIN.append(attrib)
            if 'feature' in child.tag:
                attrib=child.attrib
                if attrib['type'] not in valid_features:continue
                subchild=child[0]
                position_start=subchild[0].attrib['position']
                try:
                    position_end=subchild[1].attrib['position']
                except:
                    position_end=position_start
                attrib['pos_start']=position_start
                attrib['pos_end']=position_end
                PROTEIN.append(attrib)
        PROTEINS[protein]=PROTEIN
    except:pass
        #print(attrib)
        #break
with open('./protein_data_human.p', 'wb') as handle:
    pickle.dump(PROTEINS, handle)
with open('./protein_data_human.json', 'w') as fp:
    json.dump(PROTEINS, fp)
exit()
AMINO_ACID_CONVERSION={''}
feature_type_counts={"domain":[],"region of interest":[],"modified residue":[],"cross-link":[],"sequence":[]}
for protein in PROTEINS.keys():
    temp_prot=PROTEINS[protein]
    for feature in temp_prot:
        if feature['type'] in valid_features:
            feature_type_counts[feature['type']].append(feature['description'])
        try:
            seq=feature['sequence']
            for s in seq:
                feature_type_counts['sequence'].append(s)
        except:pass
for type in feature_type_counts.keys():
    a,b=np.unique(feature_type_counts[type],return_counts=True)
    print(type,a,b)