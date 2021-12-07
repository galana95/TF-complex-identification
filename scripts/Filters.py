import pandas as pd
import numpy as np
from scripts.utils import is_equal
class Filter():
    def __init__(self,settings):
        self.filter_names = settings['FILTER SETTINGS']['FILTER']
        self.tf_db=settings['FILTER SETTINGS']['TRANSCRIPTION FACTOR DB']
        self.complex_db=settings['FILTER SETTINGS']['COMPLEX DB']
        self.genetarget_db=settings['FILTER SETTINGS']['TARGET DB']
        self.db_join_operator=settings['FILTER SETTINGS']['OPERATOR']
        self.tfcomplex=settings['FILTER SETTINGS']['TFCOMPLEX']
    def run_filters(self,data):
        for f in self.filter_names:
            temp_filter=self.get_filter(f)
            mask=temp_filter(data)
            if f=='partial_tf_complex':
                mask=mask[1]
            elif f=='tf_complex':
                mask=mask[0]

            data=data.iloc[np.asarray(mask).astype(bool),:]

        return data
    def get_filter(self,filt):
        if filt=='tf_complex' or filt=='partial_tf_complex':
            return self.filter_transcription_factor_complexes
        elif filt=='tf':
            return self.filter_transcription_factor_proteins
        elif filt=='minprotnum2':
            return self.minprotnum2
        
    def minprotnum2(self,pred_comp):
        minprotabove=[]
        pred_complexes=pred_comp['List of proteins']
        for proteins in pred_complexes:
            proteins=proteins.split(',')
            complex_length=len(proteins)
            if complex_length>=2:
                minprotabove.append(1)
            else:
                minprotabove.append(0)
        return minprotabove

    def load_tfs(self):
        all_tfs=[]
        #if is_equal(self.db_join_operator,'union'):
        for db in self.tf_db:
            if is_equal(db,'jaspar'):
                jaspar=pd.read_csv('./databases/jaspar/tf.csv')
                jaspar2=jaspar.values
                for j in jaspar2:
                    all_tfs.append(j[0])
            elif is_equal(db,'human_tfnames'):
                htf=pd.read_csv('./databases/human_tfnames/homo_sapiens_tf_processed.csv')
                htf=htf['Uniprot'].values
                for tf in htf:
                    all_tfs.append(tf)
            elif is_equal(db,'yeastract'):
                yt=pd.read_csv('./databases/yeastract/yeast_tf_processed.csv')
                yt=yt['Uniprot'].values
                for y in yt:
                    all_tfs.append(y)
            elif is_equal(db,'animal_tfdb3'):
                atf=pd.read_csv('./databases/animal_tfdb3/homo_sapiens_tf_processed.csv')
                atf=atf['Uniprot'].values
                for tf in atf:
                    all_tfs.append(tf)
                atf=pd.read_csv('./databases/animal_tfdb3/mouse_tf_processed.csv')
                atf=atf['Uniprot'].values
                for tf in atf:
                    all_tfs.append(tf)
                
        return np.unique(all_tfs)
    def load_atfs(self):
        all_tfs=[]
        #if is_equal(self.db_join_operator,'union'):
        for db in self.tf_db:
            if is_equal(db,'jaspar'):
                jaspar=pd.read_csv('./databases/jaspar/tf.csv')
                jaspar2=jaspar.values
                for j in jaspar2:
                    all_tfs.append(j[0])
            elif is_equal(db,'human_tfnames'):
                htf=pd.read_csv('./databases/human_tfnames/homo_sapiens_tf_processed.csv')
                htf=htf['Uniprot'].values
                for tf in htf:
                    all_tfs.append(tf)
            elif is_equal(db,'yeastract'):
                yt=pd.read_csv('./databases/yeastract/yeast_tf_genetarget_processed.csv')
                yt=yt['Uniprot TF'].values
                for y in yt:
                    all_tfs.append(y)
            elif is_equal(db,'animal_tfdb3'):
                atf=pd.read_csv('./databases/animal_tfdb3/homo_sapiens_tf_processed.csv')
                atf=atf['Uniprot'].values
                for tf in atf:
                    all_tfs.append(tf)
                atf=pd.read_csv('./databases/animal_tfdb3/mouse_tf_processed.csv')
                atf=atf['Uniprot'].values
                for tf in atf:
                    all_tfs.append(tf)
                atf=pd.read_csv('./databases/animal_tfdb3/homo_sapiens_ctf_processed.csv')
                atf=atf['Uniprot'].values
                for tf in atf:
                    all_tfs.append(tf)
                atf=pd.read_csv('./databases/animal_tfdb3/mouse_ctf_processed.csv')
                atf=atf['Uniprot'].values
                for tf in atf:
                    all_tfs.append(tf)
                
        return np.unique(all_tfs)
    def load_ctfs(self):
        all_tfs=[]
        #if is_equal(self.db_join_operator,'union'):
        for db in self.tf_db:
            if is_equal(db,'animal_tfdb3'):
                atf=pd.read_csv('./databases/animal_tfdb3/homo_sapiens_ctf_processed.csv')
                atf=atf['Uniprot'].values
                for tf in atf:
                    all_tfs.append(tf)
                atf=pd.read_csv('./databases/animal_tfdb3/mouse_ctf_processed.csv')
                atf=atf['Uniprot'].values
                for tf in atf:
                    all_tfs.append(tf)
                
        return np.unique(all_tfs)

    def filter_transcription_factor_proteins(self,proteins):
        try:
            proteins=proteins['Protein name']
        except:
            proteins=proteins.iloc[:,0]
        if self.tfcomplex.lower()=='ctf':
            tf_data=self.load_ctfs()
        elif self.tfcomplex.lower()=='tf':
            tf_data=self.load_tfs()
        elif self.tfcomplex.lower()=='atf':
            tf_data=self.load_atfs()
        is_tf=[]
        for i in range(len(proteins)):
            if proteins[i] in tf_data:
                is_tf.append(1)
            else:
                is_tf.append(0)
        return is_tf
    def filter_transcription_factor_complexes(self,pred_comp):
        complex_thresholds={1:100,2:100,3:66,4:75,5:60,6:50,7:42,8:37,9:33,10:30,11:27,12:25}
        pred_complexes=pred_comp['List of proteins']

        if self.tfcomplex.lower()=='ctf':
            tf_data=self.load_ctfs()
        elif self.tfcomplex.lower()=='tf':
            tf_data=self.load_tfs()
        elif self.tfcomplex.lower()=='atf':
            tf_data=self.load_atfs()
        is_transcription_factor=[]
        is_partial_transcription_factor=[]
        #pred_complexes=self.simulation_data.iloc[:,2]
        idx=0
        all_proteins=[]
        for proteins in pred_complexes:
            proteins=proteins.split(',')
            complex_length=len(proteins)
            if complex_length>12:
                threshold=25
            else:
                threshold=complex_thresholds[complex_length]
            is_complex_num=0
            for prot in proteins:
                all_proteins.append(prot)
                if prot in tf_data:
                    
                    is_complex_num+=1
            complex_percent=(is_complex_num/complex_length)*100
            if complex_percent>=threshold:
                is_transcription_factor.append(1)
            else:
                is_transcription_factor.append(0)
            if is_complex_num>0:
                is_partial_transcription_factor.append(1)
            else:
                is_partial_transcription_factor.append(0)
            idx+=1
        #self.simulation_data['is_partial_transcription_factor']=is_partial_transcription_factor
        #self.simulation_data['is_transcription_factor']=is_transcription_factor
        #preds=preds[preds['is_partial_transcription_factor']==1]
        #preds=preds['List of proteins']
        return is_transcription_factor,is_partial_transcription_factor
    def filter_proteins(self,is_transcription_factor=None):
        try:
            complexes=self.simulation_data[is_transcription_factor==1]
        except:pass
        complexes=complexes['List of proteins']
        proteins=[]
        for i in range(len(complexes)):
            proteins.append(str.split(complexes.iloc[i],','))
        return proteins