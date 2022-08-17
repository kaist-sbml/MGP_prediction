from mgpPrediction import contributionFluxSum
import pandas as pd
import numpy as np
from mgpPrediction import pathwayConservenormalize
from scipy import stats
import time
import numpy as np


def sumPathWay(df):
    result={}
    sd={}
    for idx in df.index:
        temp=idx.split('^')
        meta=temp[0]
        path=temp[1]
        if len(temp)!=2:
            print("Check the id of metabolites in GEM!")
            raise
        
        if not meta in sd:
            sd[meta]=[]
        sd[meta].append(idx)
    for meta in sd:
        result[meta]=df.loc[sd[meta]].sum().to_dict()
    
    result=pd.DataFrame(result)
    result=result.transpose()
    
    return result

    
def statTest(df, mut):
    result={}
    
    for mg in mut.index:
        temp={}
        posi=[x for x in mut if mut[x][mg]==1]
        nega=set(mut.columns)-set(posi)
        nega=list(nega)
        
        posidf=df[posi]
        negadf=df[nega]
        
        for meta in posidf.index:
            d1=posidf.loc[meta].values
            d2=negadf.loc[meta].values
            
            pv=stats.ranksums(d1, d2)[-1]
            temp[meta]=pv
        result[mg]=temp
    result=pd.DataFrame(result)
    return result

def pathWayAnalysis(fluxdf, metabolites, sigGenes, mutdf):
    refDf=fluxdf.loc[[x for x in fluxdf.index if x.startswith(metabolites+'^')]]
    resultAvg={}
    resultPv={}
    for path in refDf.index:
        temp={}
        temp_pv={}
        for gene in sigGenes:
            posi=[x for x in mutdf if mutdf[x][gene]==1]
            nega=set(mutdf)-set(posi)
            nega=list(nega)

            posi_avg=refDf[posi].loc[path].sum()/len(posi)
            temp[gene]=posi_avg

            posi_d=refDf[posi].loc[path].values
            nega_d=refDf[nega].loc[path].values
            pv=stats.ranksums(posi_d, nega_d)[-1]
            temp_pv[gene]=pv
        resultAvg[path]=temp
        resultPv[path]=temp_pv
    
    return resultAvg, resultPv

def makeItDf(dictdata):
    df=pd.DataFrame.from_dict(dictdata)
    df=df.transpose()
    return df


def runAnalysis(FLUXPATH, MUTPATH):
    st=time.time()
    fluxSumCalculator=contributionFluxSum.ContributionFluxSum(FLUXPATH)
    proR=fluxSumCalculator.calculateCFluxSum('producingRxns')
    print('Step1 finished!', 'Elapsed time:%.3f seconds'%(time.time()-st))
    proR=proR.round(12)
    
    fluxSum2=sumPathWay(proR)
    fluxSum2=pathwayConservenormalize.Filter(fluxSum2, threshold=1e-11)
    
    mutdf=pd.read_csv(MUTPATH, index_col=0)
    netSamples=set(fluxSum2.columns)&set(mutdf.columns)
    
    mutdf=mutdf[list(netSamples)]
    fluxSum2_m=fluxSum2[list(netSamples)]
    fluxSum2=pathwayConservenormalize.quantileNormalize(fluxSum2_m)
    ## keep 0 values when it was 0 before normalization
    isNotZero=(fluxSum2_m!=0)*1
    fluxSum2=fluxSum2*isNotZero
    
    pro=statTest(fluxSum2, mutdf)
    print('Step2 finished!', 'Elapsed time:%.3f seconds'%(time.time()-st))

    producingDirection_m=proR
    
    sigGenes={}
    temp=pro.transpose()
    for meta in temp:
        temp2=list(temp[temp[meta]<0.05].index)
        if len(temp2)>0:
            sigGenes[meta]=temp2
    
    producingDirection=pathwayConservenormalize.FilterQuantileNorm(producingDirection_m)
    ## keep 0 values when it was 0 before normalization
    isNotZero=(producingDirection_m!=0)*1
    producingDirection=producingDirection*isNotZero
    
    producingDirection=producingDirection.round(14)
    
    pro_FINALResult={}
    pro_FINALResult_pv={}
    for meta in sigGenes:
        r1, r2 = pathWayAnalysis(producingDirection, meta, sigGenes[meta], mutdf)
        pro_FINALResult.update(r1)
        pro_FINALResult_pv.update(r2)    
    
    pro_pv_result=makeItDf(pro_FINALResult_pv)
    
    pro_FINALResult=makeItDf(pro_FINALResult)
    print('Step3 finished!', 'Elapsed time:%.3f seconds'%(time.time()-st))
    
    return pro_FINALResult, pro_pv_result
    
    
    
