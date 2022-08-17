import pandas as pd
import numpy as np
import cobra
from scipy.stats import norm
import math


def convert2modifiedAbsZ(df):
    zsDf={}
    scaleFactor1=1/norm.ppf(3/4) #1.482602218505602
    scaleFactor2=(math.pi/2)**(1/2) #1.2533141373155001
    for i in df.index:
        data=df.loc[i]
        temp=data.dropna()
        median = np.median(temp)
        MAD = np.median([np.abs(v - median) for v in temp])
        MeanAD = np.mean([np.abs(v - median) for v in temp])
        if MAD ==0:
            absZscore=abs((temp - median)/scaleFactor2*MeanAD)
        else:
            absZscore= abs((temp-median)/(scaleFactor1*MAD))
        zsDf[i]=absZscore.to_dict()
    zsDf=pd.DataFrame(zsDf)
    zsDf=zsDf.transpose()
    return zsDf


def runStep4(averdf, pvdf,threshold=3.5):
    unusedSubsys1=['Unassigned','Transport, mitochondrial', 'Transport, peroxisomal', 'Transport, extracellular', 'Exchange/demand reaction', 'Transport, lysosomal', 'Transport, nuclear','Transport, endoplasmic reticular', 'Transport, golgi apparatus']
    unusedSubsys2=['Unassigned','Exchange/demand reaction']
    
    prepro={}
    
    for meta in pvdf.index:
        data=pvdf.loc[meta]
        data=data.dropna()
        sigsData=data[data<=0.05]
        
        temp=averdf[sigsData.index].loc[meta]
        prepro[meta]=temp.to_dict()
        
    prepro=pd.DataFrame.from_dict(prepro)
    prepro=prepro.transpose()
    
    zdf=convert2modifiedAbsZ(prepro)
    
    record=[]
    ESSENTIALS=['MNXM142','MNXM199','MNXM231','MNXM140','MNXM78','MNXM61','MNXM94','MNXM97', 'MNXM134']
    for meta in zdf.index:
        data=zdf.loc[meta]
        data=data.dropna()
        if len(data)<=2:
            highZsData=data
        else:
            highZsData=data[data>threshold]
        
        highZsData=highZsData.to_dict()
        metabol, path=meta.split('^')
        
        if not metabol[:-2] in ESSENTIALS:
            if path in unusedSubsys1:
                continue
        else:
            if path in unusedSubsys2:
                continue
        
        for G in highZsData:
            record.append([G, metabol, path])
            
    record=pd.DataFrame(record, columns=['gene','metabolite', 'pathway'])
    
    return record


def standardizedCompare(pro, pro_pv, threshold=3.5):
    predictedMGPs=runStep4(pro, pro_pv, threshold)
    model=cobra.io.read_sbml_model('./mgpPrediction/data/Recon2M.2_MNX_Ensembl_Transcript.xml')
    
    nameDictM={}
    nameDictN={}
    for i in model.metabolites:
        nameDictM[i.id]=i.name
        nameDictN[i.name]=i.id
    currencyMetabolites=['CO2', 'H2O', 'NADH', 'FADH2', 'CoA', 'H(+)', 'NH3', 'NAD(+)', 'CoA', 'NADH',
                         'FAD', 'O2', 'NADP(+)', 'H2O2', 'NADPH', 'ATP', 'ADP', 'GDP', 'GTP','dCMP',
                         'dCTP','dUDP','uracil','CTP','UTP','Na(+)','phosphate','guanosine','thymine','dUMP',
                         'cytosine','dUTP', 'GMP','K(+)','adenine','CMP', 'adenosine','carbonic acid','diphosphate',
                         'sulfate','CDP','dGDP', 'dTMP', 'dAMP','uridine','thiamine','dCDP', 'dTTP',
                         'guanine','UMP','dTDP','thymidine','UDP','sulfite','AMP', 'dADP', 'dGTP', 'dATP', 'dCTP','Zn(2+)',
                        'dCMP', 'cytidine']    
    drops=[]
    for idx in predictedMGPs.index:
        test=predictedMGPs['metabolite'][idx]
        if nameDictM[test] in currencyMetabolites:
            drops.append(idx)
    
    predictedMGPs=predictedMGPs.drop(drops)

    return predictedMGPs

