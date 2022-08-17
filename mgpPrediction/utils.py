import pandas as pd
import cobra
import numpy as np

def parsing(data):
    check={}
    for idx in data.index:
        gene, meta, path=list(data.loc[idx])
        pair=(gene, meta[:-2], path)
        if not pair in check:
            check[pair]=[]
        check[pair].append([path])
    temp={}
    for pair in check:
        temp[pair]=['&'.join(x) for x in check[pair]]
        temp[pair]='&'.join(temp[pair])
    return temp
    
def pairing(data):
    model=cobra.io.read_sbml_model('./mgpPrediction/data/Recon2M.2_MNX_Ensembl_Transcript.xml')
    nameDict={x.id[:-2]:x.name for x in model.metabolites}
    record=[]
    for idx in data:
        record.append([idx[0], idx[1],idx[2]])
    result=pd.DataFrame(record, columns=['Target gene', 'Target metabolite MetaNetX ID','Target pathway'])
    result['Target metabolite']=[nameDict[x] for x in result['Target metabolite MetaNetX ID']]
    result=result[['Target gene','Target metabolite', 'Target metabolite MetaNetX ID', 'Target pathway']]
    return result

def createOrganizedOutput(mgps):

    parsed=parsing(mgps)
    result=pairing(parsed)
    
    cancerdf=result

    cancerdf=cancerdf.reset_index(drop=True)
    return cancerdf
        
