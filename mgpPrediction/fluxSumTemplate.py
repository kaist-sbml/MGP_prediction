import cobra

def run():
    result={}
    model=cobra.io.read_sbml_model('./mgpPrediction/data/Recon2M.2_MNX_Ensembl_Transcript.xml')
    for m in model.metabolites:
        temp={'pre':[], 'next':[]}
        for r in m.reactions:
            pro=[x.id for x in r.products]
            if m.id in pro:
                temp['pre'].append(r.id)
            else:
                temp['next'].append(r.id)
        result[m.id]=temp
    Template=result
    return Template
