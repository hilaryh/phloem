from cobra import io
import numpy as np
from cobra.flux_analysis.variability import *
from cobra.flux_analysis.parsimonious import pfba
import pandas as pd
from math import nan
import re

def mainCR(minstd,modelpath=None,cplxs=0):
    transcriptome_data='Kim_et_al_2021_transcriptome_latest.csv'
    gene_mapping      ='Model_gene_rxn_mappings_Arabidopsis.csv'
    vfva_sol=None
    # If no preexisting files, generate from stem model
    # stem_modelm,vfva_sol=loadModel(modelpath)
    # Dm, Dstdm = loadTranscriptome(stem_modelm,transcriptome_data,gene_mapping,minstd)
    modelpath='combo_stem.xml'
    stem_modelm=io.read_sbml_model(modelpath)
    stem_modelm.reactions.aPPI_c_CC_SE_l.lower_bound=-1000
    # # Halve SE and pSE maintenance costs...
    # for rxn in stem_modelm.reactions:
    #     if 'SE' in rxn.id and ('ATPase_tx' in rxn.id or 'turnover' in rxn.id):
    #         rxn.lower_bound = 0.5*rxn.lower_bound

    if not vfva_sol:
        vfva_sol=flux_variability_analysis(stem_modelm)
    
    Ddict, Dstddict = loadTranscriptome(stem_modelm,transcriptome_data,gene_mapping,minstd,cplxs)
    stem_modelm,rxn_weighting,rxnitweighting,regrexdifference,coefficients,regrexconstraint,Ddiffp,Ddiffn, Dcons,Dconsdp,Dconsdn = addRegrexCons(stem_modelm,Ddict,Dstddict,gene_mapping,minstd)
    saveRegrexResults(stem_modelm,regrexdifference,Ddict,rxnitweighting,vfva_sol,minstd,modelpath)

def saveRegrexResults(stem_modelm,regrexdifference,Ddict,rxn_weighting,vfva_sol,minstd,modelpath):
    ub={}
    lb={}
    for rxn in stem_modelm.reactions:
        ub[rxn.id]=rxn.upper_bound
        lb[rxn.id]=rxn.lower_bound

    stem_modelm.objective='diel_biomass'
    stem_modelm.reactions.diel_biomass.upper_bound=1000
    # stem_modelm.reactions.diel_biomass.lower_bound=temp['diel_biomass']
    solmax = pfba(stem_modelm)
    print(solmax['diel_biomass'])

    stem_modelm.reactions.PSII_RXN_p_CC_l.upper_bound=0
    stem_modelm.objective='diel_biomass'
    nopsiisolmax = pfba(stem_modelm)
    stem_modelm.reactions.PSII_RXN_p_CC_l.upper_bound=ub['PSII_RXN_p_CC_l']
    regrexobj=stem_modelm.problem.Objective(regrexdifference,
        direction='min')
    stem_modelm.objective=regrexobj
    stem_modelm.reactions.diel_biomass.lower_bound=float(solmax['diel_biomass']*0.68)
    rsol=pfba(stem_modelm)
    
    fva_sol2 = flux_variability_analysis(stem_modelm)
    
    stem_modelm.reactions.diel_biomass.lower_bound=float(solmax['diel_biomass']*0.68)
    # # Reverse H+PPase
    # stem_modelm.reactions.PROTON_PPI_ec_CC_l.upper_bound=-0.1
    # revppisol=pfba(stem_modelm)
    # stem_modelm.reactions.PROTON_PPI_ec_CC_l.upper_bound=ub['PROTON_PPI_ec_CC_l']

    # # Limit mitochondrial bits
    # exclrxns=pd.read_csv('excludefromSE.csv')
    # exclrxnlist=list(exclrxns['0'])
    # for rxn in [x for x in stem_modelm.reactions if '_SE_' in x.id or '_pSE_' in x.id]:
    #     if any([x in rxn.id for x in exclrxnlist]):
    #         rxn.lower_bound=0
    #         rxn.upper_bound=0
    # mitosol=pfba(stem_modelm)
    # mitofva_sol = flux_variability_analysis(stem_modelm)

    data={'Reaction ID':[x.id for x in stem_modelm.reactions],
        # 'cell':['_'.join(x.id.split('_')[-2:]) if ('_l' in x.id or '_d' in x.id) else '' for x in stem_modelm.reactions],
        'Transcript abundance':[Ddict[x.id] if x.id in Ddict.keys() else '' for x in stem_modelm.reactions],
        'Subsystem':[x.notes['SUBSYSTEM']  if 'SUBSYSTEM' in str(x.notes) else '' for x in stem_modelm.reactions],
        'Stoichiometry':[x.reaction for x in stem_modelm.reactions],
        'Lower bound':[x.lower_bound for x in stem_modelm.reactions],
        'Upper bound':[x.upper_bound for x in stem_modelm.reactions],
        'RegrEx':[rsol[x.id] for x in stem_modelm.reactions],
        'FVA_min':[fva_sol2.minimum[x.id] for x in stem_modelm.reactions],
        'FVA_max':[fva_sol2.maximum[x.id] for x in stem_modelm.reactions],
        'Parsimonious':[solmax[x.id] for x in stem_modelm.reactions],
        'FVA_min (parsimonious)':[vfva_sol.minimum[x.id] for x in stem_modelm.reactions],
        'FVA_max (parsimonious)':[vfva_sol.maximum[x.id] for x in stem_modelm.reactions]}
        # 'RegrEx (AVP1 reversed)':[revppisol[x.id] for x in stem_modelm.reactions],
        # 'FVA_min (w/o PSII)':[nopsiifva_sol2.minimum[x.id] for x in stem_modelm.reactions],
        # 'FVA_max (w/o PSII)':[nopsiifva_sol2.maximum[x.id] for x in stem_modelm.reactions],
        # 'RegrEx (w/o PSII)':[nopsiirsol[x.id] for x in stem_modelm.reactions],
        # 'Parsimonious (w/o PSII)':[nopsiisolmax[x.id] for x in stem_modelm.reactions],
        # 'Restr mitochondria':[mitosol[x.id] for x in stem_modelm.reactions],
        # 'FVA_min (mito)':[mitofva_sol.minimum[x.id] for x in stem_modelm.reactions],
        # 'FVA_max (mito)':[mitofva_sol.maximum[x.id] for x in stem_modelm.reactions]}
    df=pd.DataFrame(data)
    from datetime import datetime, date, time
    now = datetime.now().strftime('%Y_%m_%d')
    if modelpath:
        df.to_csv('RegrExFVAsol_'+modelpath.split('.')[0]+str(minstd)+'.csv')
    else:
        df.to_csv('RegrExFVAsol_'+str(minstd)+'.csv')



def addRegrexCons(stem_modelm,Ddict,Dstddict,gene_mapping,minstd):
    # rxns = [rxn.id for rxn in stem_modelm.reactions]
    Drxns = [rxn for rxn in list(Ddict.keys()) if Ddict[rxn]>-1 and '_l' in rxn]
    Ddiffp={}
    Ddiffn={}
    Dcons={}
    Dconsdp={}
    Dconsdn={}
    mapFrame=pd.read_csv(gene_mapping)
    rxn_ids=list(mapFrame['rxn_id'])
    rxnweighting={}
    rxnitweighting={}


    for rxnit in Drxns:
        if not 'MC_l' in rxnit:
            map_no=[rxn_ids[x] for x in range(len(rxn_ids)) if rxn_ids[x] in rxnit]
            if len(map_no)>1:
                if '_'.join(rxnit.split('_')[:-2]) in map_no:
                    map_no=['_'.join(rxnit.split('_')[:-3])]
                else:
                    print(rxnit,map_no)
            rxn=map_no[0]
            sameDrxns = [x.id for x in stem_modelm.reactions if rxn in x.id and any([y in x.id for y in ['_CC_l','_SE_l','_pSE_l']])]
            if any([abs(len(x)-len(sameDrxns[0]))>1 for x in sameDrxns]):
                print('Possible problem - these reactions are lumped together:')
                print(sameDrxns)
        else:
            map_no=['_'.join(rxnit.split('_')[:-3])]
            rxn=map_no[0]
            sameDrxns = [x.id for x in stem_modelm.reactions if rxn =='_'.join(x.id.split('_')[:-3]) and 'MC_l' in x.id]
            if any([abs(len(x)-len(sameDrxns[0]))>1 for x in sameDrxns]):
                print('Possible problem - these reactions are lumped together:')
                print('root = ',rxnit)
                print(sameDrxns)
            rxn=rxnit
            # sameDrxns = [rxn]
        if rxn not in Dcons.keys():
            Ddiffp[rxn] = stem_modelm.problem.Variable('Ddiffp'+rxn)
            Ddiffn[rxn] = stem_modelm.problem.Variable('Ddiffn'+rxn)
            Dconsdp[rxn] = stem_modelm.problem.Constraint(Ddiffp[rxn],lb=0,ub=1000)
            Dconsdn[rxn] = stem_modelm.problem.Constraint(Ddiffn[rxn],lb=0,ub=1000)
            stem_modelm.add_cons_vars([Ddiffp[rxn], Ddiffn[rxn],Dconsdp[rxn],Dconsdn[rxn]])
            Dcoeffs = {}
            for ii in sameDrxns:
                Dcoeffs[stem_modelm.reactions.get_by_id(ii).forward_variable]=1
                Dcoeffs[stem_modelm.reactions.get_by_id(ii).reverse_variable]=1
            Dcoeffs[Ddiffp[rxn]]=-1
            Dcoeffs[Ddiffn[rxn]]=1
            Dcons[rxn] = stem_modelm.problem.Constraint(0,
                lb=Ddict[rxnit],
                ub=Ddict[rxnit])
            stem_modelm.add_cons_vars([Dcons[rxn]])
            stem_modelm.solver.update()
            Dcons[rxn].set_linear_coefficients(coefficients=Dcoeffs)
            if Dstddict[rxnit]==0 or Dstddict[rxnit]<minstd:
                Dstddict[rxnit]=minstd
            rxnweighting[rxn]=1/Dstddict[rxnit]
            rxnitweighting[rxnit]=rxnweighting[rxn]
        else:
            rxnitweighting[rxnit]=rxnweighting[rxn]

    coefficients = {}
    regrexdifference = stem_modelm.problem.Variable('regrexdifference')
    stem_modelm.add_cons_vars(regrexdifference)
    coefficients[regrexdifference]=-1
    for rxn in list(Ddiffp.keys()):
        coefficients[Ddiffp[rxn]]=rxnweighting[rxn]
        coefficients[Ddiffn[rxn]]=rxnweighting[rxn]
    regrexconstraint = stem_modelm.problem.Constraint(0, lb=0, ub=0)
    stem_modelm.add_cons_vars(regrexconstraint)
    stem_modelm.solver.update()
    regrexconstraint.set_linear_coefficients(coefficients=coefficients)
    return stem_modelm,rxnweighting,rxnitweighting,regrexdifference,coefficients,regrexconstraint,Ddiffp,Ddiffn, Dcons,Dconsdp,Dconsdn
        


def loadPrep(stem_modelm,Dfile='matlabD.csv'):
    # rxns=list(stoic_matm.columns)
    rxns = [rxn.id for rxn in stem_modelm.reactions]
    df=pd.read_csv(Dfile)
    rxnm=list(df['rxn'])
    Dm=list(df['D'])
    Dstdm=list(df['Dstd'])
    if set(rxns).difference(set(rxnm)):
        print('Problem - csv may be for different model.')

    Ddict = {}
    Dstddict = {}
    for it,rxn in enumerate(rxnm):
        Ddict[rxn] = Dm[it]
        Dstddict[rxn] = Dstdm[it]

    return Ddict, Dstddict


def loadModel(modelpath):
    # Import model
    stem_model=io.read_sbml_model(modelpath)
    ## Trim constrained reactions
    for rxn in stem_model.reactions:
        if (rxn.upper_bound==0.0) & (rxn.lower_bound==0.0):
            stem_model.remove_reactions([rxn.id])
    stem_model.reactions.diel_biomass.lower_bound=0
    
    # Solve model
    sol=stem_model.optimize()

    # Remove unnecessary reactions
    stem_model.reactions.Phloem_output_tx_d.lower_bound = float(0.68*sol['diel_biomass'])
    stem_model.reactions.Phloem_output_tx_d.upper_bound = float(0.68*sol['diel_biomass'])
    fva_sol = flux_variability_analysis(stem_model)
    nofluxrxns = []
    for rxn in stem_model.reactions:
        if fva_sol.maximum[rxn.id]==0 and fva_sol.minimum[rxn.id]==0:
            nofluxrxns+=[rxn.id]
    stem_model.remove_reactions(nofluxrxns)
    stem_model.reactions.Phloem_output_tx_d.upper_bound = 1000

    return stem_model,fva_sol

def loadTranscriptome(stem_model,transcriptome_data,gene_mapping,minstd,cplxs):
    rxns = [rxn.id for rxn in stem_model.reactions]
    transcriptFrame = pd.read_csv(transcriptome_data)    
    gene_ids    = list(transcriptFrame['Gene ID'])
    gene_names    = list(transcriptFrame['Gene name'])
    wt_transcripts={}
    wt_transcripts[0] = [np.log2(x) if (x!=0 and np.log2(x)>0) else 0 for x in list(transcriptFrame['FPKM P1 CC'])]
    wt_transcripts[1] = [np.log2(x) if (x!=0 and np.log2(x)>0) else 0 for x in list(transcriptFrame['FPKM P2 CC'])] 
    transcriptsMC = {}
    transcriptsMC[0] = [20*np.log2(x) if (x!=0 and np.log2(x)>0) else 0 for x in list(transcriptFrame['FPKM P1'])]
    transcriptsMC[1] = [20*np.log2(x) if (x!=0 and np.log2(x)>0) else 0 for x in list(transcriptFrame['FPKM P2'])]

    mapFrame=pd.read_csv(gene_mapping)
    rxn_ids=list(mapFrame['rxn_id'])
    if '4' in gene_mapping:
        if not cplxs:
            maps=list(mapFrame['All related genes'])
        else:
            maps=list(mapFrame['Genes not in complex'])
            complex_maps = list(mapFrame['Genes in complex2'])
    else:    
        maps=list(mapFrame['Gene ID'])
    altmaps = list(mapFrame['GeneSynonyms'])

    D=[]
    Dstd=[]
    for rxn in rxns:
        Ds=[]
        # Check if reaction is light phase, phloem
        if any([x in rxn for x in ['_CC_l','_SE_l','_pSE_l']]):        
            # Find rxn root in list of core model reactions
            map_no=[x for x in range(len(rxn_ids)) if rxn_ids[x] in rxn]
            # If the reaction is a derivative of a core model reaction, find the associated genes and add the associated transcript abundance of each gene to the list of reaction data
            if map_no:
    #             print(map_no)
    #             print(rxn)
                genes=maps[map_no[0]]
                altgenenames=altmaps[map_no[0]]
                if cplxs:
                    complex_temp = complex_maps[map_no[0]]
                    if str(complex_temp)!='nan':
                        complexes=complex_temp.split(';')
                    else:
                        complexes=''
                if str(genes)!='nan' or str(altgenenames)!='nan':
                    if str(genes)!='nan':
                        gene_list=re.split(',',genes)
                    else:
                        gene_list=[]
                    if str(altgenenames)!='nan':
                        altgene_list=re.split(',',altgenenames)
                    else:
                        altgene_list=[]
                    rxn_locs = [x for x in range(len(gene_ids)) if str(gene_ids[x]) in gene_list]
                    olenrxn_locs=len(rxn_locs)
                    altrxn_locs = [x for x in range(len(gene_names)) if str(gene_names[x]) in altgene_list]
                    rxn_locs = list(set(rxn_locs+altrxn_locs))
                    # if olenrxn_locs!=len(rxn_locs):
                        # print('Potentially different transcript abundance: ',rxn)
                    if rxn_locs:            
                        for loc in rxn_locs:
                            Ds_temp=[wt_transcripts[1][loc]]+[wt_transcripts[0][loc]]
                            Ds=Ds+Ds_temp
                if cplxs:
                    cplx_locs = [[x for x in range(len(gene_ids)) if str(gene_ids[x]) in y] for y in complexes.split(',')]
                    if cplx_locs:
                        for cpx in cplx_locs:
                            Ds_temp = [min([wt_transcripts[1][loc] for loc in cpx])]+[min([wt_transcripts[0][loc] for loc in cpx])]
                            Ds=Ds+Ds_temp

        # Repeat for light MC reactions but use MC transcript data
        if 'MC_l' in rxn:        
            map_no=[x for x in range(len(rxn_ids)) if rxn_ids[x] in rxn]
            if map_no:
                genes=maps[map_no[0]]
                altgenenames=altmaps[map_no[0]]
                if cplxs:
                    complex_temp = complex_maps[map_no[0]]
                    if str(complex_temp)!='nan':
                        complexes=complex_temp.split(';')
                if str(genes)!='nan':
                    gene_list=re.split(',',genes)
                    rxn_locs = [x for x in range(len(gene_ids)) if str(gene_ids[x]) in gene_list]
                    altrxn_locs = [x for x in range(len(gene_names)) if str(gene_names[x]) in altgene_list]
                    rxn_locs = list(set(rxn_locs+altrxn_locs))
                    if rxn_locs:            
                        for loc in rxn_locs:
                            Ds_temp=[transcriptsMC[1][loc]]+[transcriptsMC[0][loc]]
                            Ds=Ds+Ds_temp
                if cplxs:
                    cplx_locs = [[x for x in range(len(gene_ids)) if str(gene_ids[x]) in y] for y in complexes.split(',')]
                    if cplx_locs:
                        for cpx in cplx_locs:
                            Ds_temp = [min([wt_transcripts[1][loc] for loc in cpx])]+[min([wt_transcripts[0][loc] for loc in cpx])]
                            Ds=Ds+Ds_temp
        # If there was no data on a reaction, return a 'nan' otherwise, take the maximum of the abundance data found
        if Ds==[]:
            D.append(-1000)
            Dstd.append(1)
        else:
            D.append(np.max(Ds))
            Dstd.append(min(np.std(Ds),minstd))
    # D[0]=max(D)
    if len(D)!=len(rxns):
        print('Error: D is the wrong length')
    if len(Dstd)!=len(rxns):
        print('Error: Dstd is the wrong length')

    Ddict = {}
    Dstddict = {}
    for it,rxn in enumerate(rxns):
        Ddict[rxn] = D[it]/max(D)
        Dstddict[rxn] = Dstd[it]/max(D)
        if Ddict[rxn]<0:
            Ddict[rxn]=nan
    return Ddict, Dstddict
