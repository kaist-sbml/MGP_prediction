import pandas as pd
import numpy as np


def quantileNormalize(data):
    df = data.copy()
    temp = {}
    
    for c in df:
        temp[c] = df[c].sort_values(na_position='first').values
    df2 = pd.DataFrame(temp)
    rank = df2.median(axis = 1).tolist()

    for c2 in df:
        t = df[c2].rank( pct=True, method='max' ).values
        df[c2] = [ np.nanpercentile( rank, i*100 ) if ~np.isnan(i) else np.nan for i in t ]
    return df
    
def fluxSum(df):
    parsed={}
    for idx in df.index:
        meta=idx.split('^')[0]
        if not meta in parsed:
            parsed[meta]=[]
        parsed[meta].append(idx)
    
    result={}
    for meta in parsed:
        called=df.loc[parsed[meta]].sum().to_dict()
        result[meta]=called
    result=pd.DataFrame(result)
    
    return result.transpose()

def Filter(fluxD, threshold):
    drops=[]
    for i in fluxD.index:
        data=fluxD.loc[i].values
        if abs(max(data)-min(data))<threshold:
            drops.append(i)
    return fluxD.drop(drops)


def FilterQuantileNorm(df, threshold=1e-11):
    fluxSumRaw=fluxSum(df)
    fluxSumRaw=Filter(fluxSumRaw, threshold)
    
    
    fluxSumQuan=quantileNormalize(fluxSumRaw)
    coef=fluxSumQuan/fluxSumRaw
    coef=coef.fillna(0)
    
    temp={}
    for idx in df.index:
        meta=idx.split('^')[0]
        if not meta in coef.index:
            continue
        called=df.loc[idx]*coef.loc[meta]
        temp[idx]=called
    result=pd.DataFrame(temp)
    result=result.transpose()
    result=result.fillna(0)
    return result
