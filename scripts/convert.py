import pandas as pd
import numpy as np
def convert_gene_to_uniprot(genes):
    import urllib.parse
    import urllib.request

    url = 'https://www.uniprot.org/uploadlists/'

    params = {
    'from': 'GENENAME',
    'to': 'ID',
    'format': 'tab',
    'query': 'AATF ABL1'
    }

    data = urllib.parse.urlencode(params)
    data = data.encode('utf-8')
    req = urllib.request.Request(url, data)
    with urllib.request.urlopen(req) as f:
        response = f.read()
    response1 = response.text()
    return response


def drug_gene_to_csv(name):
    file = open(name, 'r')
    Lines = file.readlines()
    

    # Strips the newline character
    data=[]
    temp_data=[]
    temp_gene=""
    for line in Lines:
        try:
            line_data=line.split('\t')
            #print(line_data)
            if line_data[1]=='GENENAME':

                temp_gene=line_data[2]
            if line_data[1]=='DRUGINFO' and (line_data[4]=='Approved\n' or line_data[4]=='Patented\n'):
                temp_data.append(temp_gene)
                temp_data.append(line_data[3])
                temp_data.append(line_data[4].replace('\n',''))
                data.append(np.asarray(temp_data))
                temp_data=[]
        except:pass
    data=np.asarray(data)
    print(data.shape)
    print(data[:10])
    df=pd.DataFrame(np.asarray(data),columns=['GeneName','Drug','DrugStatus'])
    df.to_csv('./databases/drug_gene/drug_gene.csv')
        

if __name__ == '__main__':
    import networkx
    import graphistry
    graphistry.register(api=3, protocol="https", server="hub.graphistry.com", username="sz.janni", password="Jedike73172833096")    
    net=networkx.readwrite.gml.read_gml(f'./databases/combined_graphdbs/human/complex_tf_target_human.gml')
    graphistry.bind(source='src', destination='dst', node='nodeid').plot(net)
    exit()
    data=pd.read_csv('./data/Ab_24_prediction.csv',delimiter=';')
    data2=pd.read_csv('./data/yeasgene_to_uniprot.tsv',delimiter='\t')
    for d in range(data.shape[0]):
        data_temp=data['List of proteins'].iloc[d]
        data_temp=str.split(data_temp,',')
        for d2 in range(len(data_temp)):
            data_temp2=data2[data2.iloc[:,0]==data_temp[d2]]
            try:
                data_temp2=str(data_temp2['Entry'].values[0])
            except:continue
            data_temp[d2]=data_temp2
            
        data_temp=','.join(data_temp)
        
        data['List of proteins'].iloc[d]=data_temp
    data.to_csv('./data/Ab_24_prediction_up.csv',index=False)
    exit()
    
    data=pd.read_csv('./data/Ab_18.csv',delimiter='\t',header=None)
    data2=pd.read_csv('./data/yeasgene_to_uniprot.tsv',delimiter='\t')
    for d in range(data.shape[0]):

        data_temp=data2[data2.iloc[:,0]==data.iloc[d,0]]
        try:data.iloc[d,0]=str(data_temp['Entry'].values[0])
        except:pass
    data.to_csv('./data/Ab_18_up.csv',index=False)
    exit()
    
    drug_gene_to_csv('./databases/drug_gene/P1-01-TTD_target_download.txt')
    exit()

    data=pd.read_csv('./databases/animal_tfdb3/homo_sapiens_ctf.txt',delimiter='\t')
    data2=pd.read_csv('./databases/animal_tfdb3/homo_sapiens_ctf_convert.tab',delimiter='\t')
    family=data['Family'].values
    data=data['Symbol'].values
    

    uniprot=[]

    for d in data:
        data_temp=data2[data2.iloc[:,0]==d]
        try:
            uniprot.append(data_temp.iloc[0,1])
        except:uniprot.append(np.nan)
    uniprot=np.asarray(uniprot)
    df=pd.DataFrame(np.vstack((data,uniprot,family)).T,columns=['Gene name','Uniprot','Family'])
    df.to_csv('./databases/animal_tfdb3/homo_sapiens_ctf_processed.csv')


    exit()
    data=pd.read_csv('./databases/yeastract/yeast_tf_gene_target.tsv',delimiter=';')
    data2=pd.read_csv('./databases/yeastract/yeast_TF.tab',delimiter='\t')
    data3=pd.read_csv('./databases/yeastract/yeast_Gene.tab',delimiter='\t')

    uniprot1=[]
    uniprot2=[]
    regulation=[]

    data_1=data.iloc[:,0]
    data_2=data.iloc[:,1]
    for d,d2 in zip(data_1,data_2):
        data_temp_1=data2[data2.iloc[:,0]==d]
        data_temp_2=data3[data3.iloc[:,0]==d2]
        if len(data_temp_1)>0:
            uniprot1.append(data_temp_1.iloc[0,1])
        else:
            uniprot1.append(np.nan)
        if len(data_temp_2)>0:
            uniprot2.append(data_temp_2.iloc[0,1])
        else:
            uniprot2.append(np.nan)
        #try:
        #    uniprot.append(data_temp.iloc[0,1])
        #except:uniprot.append(np.nan)
    uniprot1=np.asarray(uniprot1)
    uniprot2=np.asarray(uniprot2)
    df=pd.DataFrame(np.hstack((data.values[:,:2],uniprot1.reshape(-1,1),uniprot2.reshape(-1,1))),columns=['TF','Target','Uniprot TF','Uniprot Target'])
    df.to_csv('./databases/yeastract/yeast_tf_genetarget_processed.csv')
    exit()

    
    # data=pd.read_csv('./databases/tftarget/homo_sapiens_tf_genetarget.txt',delimiter='\t')
    # data2=pd.read_csv('./databases/tftarget/tftarget_TF.tab',delimiter='\t')
    # data3=pd.read_csv('./databases/tftarget/tftarget_gene.tab',delimiter='\t')

    # uniprot1=[]
    # uniprot2=[]
    # regulation=[]

    # data_1=data.iloc[:,0]
    # data_2=data.iloc[:,1]
    # for d,d2 in zip(data_1,data_2):
    #     data_temp_1=data2[data2.iloc[:,0]==d]
    #     data_temp_2=data3[data3.iloc[:,0]==d2]
    #     if len(data_temp_1)>0:
    #         uniprot1.append(data_temp_1.iloc[0,1])
    #     else:
    #         uniprot1.append(np.nan)
    #     if len(data_temp_2)>0:
    #         uniprot2.append(data_temp_2.iloc[0,1])
    #     else:
    #         uniprot2.append(np.nan)
    #     #try:
    #     #    uniprot.append(data_temp.iloc[0,1])
    #     #except:uniprot.append(np.nan)
    # uniprot1=np.asarray(uniprot1)
    # uniprot2=np.asarray(uniprot2)
    # df=pd.DataFrame(np.hstack((data.values,uniprot1.reshape(-1,1),uniprot2.reshape(-1,1))),columns=['TF','Target','Tissue','Uniprot TF','Uniprot Target'])
    # df.to_csv('./databases/tftarget/homo_sapiens_tf_genetarget_processed.csv')
    # exit()


    data=pd.read_csv('./databases/trrust/mouse_tf_genetarget.tsv',delimiter='\t')
    data2=pd.read_csv('./databases/trrust/trrust_mouse_TF.tab',delimiter='\t')
    data3=pd.read_csv('./databases/trrust/trrust_mouse_Gene.tab',delimiter='\t')

    uniprot1=[]
    uniprot2=[]
    regulation=[]

    data_1=data.iloc[:,0]
    data_2=data.iloc[:,1]
    for d,d2 in zip(data_1,data_2):
        data_temp_1=data2[data2.iloc[:,0]==d]
        data_temp_2=data3[data3.iloc[:,0]==d2]
        if len(data_temp_1)>0:
            uniprot1.append(data_temp_1.iloc[0,1])
        else:
            uniprot1.append(np.nan)
        if len(data_temp_2)>0:
            uniprot2.append(data_temp_2.iloc[0,1])
        else:
            uniprot2.append(np.nan)
        #try:
        #    uniprot.append(data_temp.iloc[0,1])
        #except:uniprot.append(np.nan)
    uniprot1=np.asarray(uniprot1)
    uniprot2=np.asarray(uniprot2)
    df=pd.DataFrame(np.hstack((data.values[:,:3],uniprot1.reshape(-1,1),uniprot2.reshape(-1,1))),columns=['TF','Target','Interaction','Uniprot TF','Uniprot Target'])
    df.to_csv('./databases/trrust/mouse_tf_genetarget_processed.csv')
    exit()

    data=pd.read_csv('./databases/animal_tfdb3/mouse_tf.txt',delimiter='\t')
    data2=pd.read_csv('./databases/animal_tfdb3/animal_tfdb3_mouse.tab',delimiter='\t')
    family=data['Family'].values
    data=data['Symbol'].values
    

    uniprot=[]

    for d in data:
        data_temp=data2[data2.iloc[:,0]==d]
        try:
            uniprot.append(data_temp.iloc[0,1])
        except:uniprot.append(np.nan)
    uniprot=np.asarray(uniprot)
    df=pd.DataFrame(np.vstack((data,uniprot,family)).T,columns=['Gene name','Uniprot','Family'])
    df.to_csv('./databases/animal_tfdb3/mouse_tf_processed.csv')
    exit()


    data=pd.read_csv('./databases/human_tfnames/homo_sapiens_tf.txt',delimiter=';',header=None)
    data2=pd.read_csv('./databases/human_tfnames/human_tf_names.tab',delimiter='\t')

    uniprot=[]
    dataval=[]
    for d in data.values:
        dataval.append(d[0])
        data_temp=data2[data2.iloc[:,0]==d[0]]
        try:
            uniprot.append(data_temp.iloc[0,1])
        except:uniprot.append(np.nan)
    uniprot=np.asarray(uniprot)
    df=pd.DataFrame(np.vstack((dataval,uniprot)).T,columns=['Gene name','Uniprot'])
    df.to_csv('./databases/human_tfnames/homo_sapiens_tf_processed.csv')
    exit()


    data=pd.read_csv('./databases/yeasttract/yeast_tf.csv',delimiter=';')
    data2=pd.read_csv('./databases/yeasttract/yeastract_TF.tab',delimiter='\t')
    data=np.unique(data['Transcription Factor'])
    uniprot=[]
    for d in data:
        d=d[:-1]
        data_temp=data2[data2.iloc[:,0]==d]
        try:
            uniprot.append(data_temp.iloc[0,1])
        except:uniprot.append(np.nan)
    uniprot=np.asarray(uniprot)
    df=pd.DataFrame(np.vstack((data,uniprot)).T,columns=['Gene name','Uniprot'])
    df.to_csv('./databases/yeasttract/yeast_tf_processed.csv')
