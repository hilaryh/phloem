import numpy as np
from datetime import datetime
import cobra
import re
from scipy.optimize import curve_fit
from cobra import flux_analysis
from collections import OrderedDict
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
from matplotlib.colors import LinearSegmentedColormap
from scipy import stats
from sklearn.datasets import *
import math
from graphviz import Digraph
def netMetaboliteStoich(cobra_model,rxnlist):
  '''
  Function to identify net synthesis/consumption of metabolites in a given set
  of reactions
  args: 1) a solved cobra model 2) a list of reactions from which net metabolite
  stoichiometry should be calculated
  output: a dictionay file with metabolite ids as key and net soichiometry as
  value
  '''
  netMet = dict()
  for rxn in rxnlist:
    rxn = cobra_model.reactions.get_by_id(rxn)
    if round(rxn.flux,6)==0:
      print(rxn.id+" flux is 0.")
      netMet = dict()
      break
    for met in rxn.metabolites:
      if netMet.keys().__contains__(met.id):
        netMet[met.id]=netMet[met.id]+((rxn.flux/abs(rxn.flux))*rxn.metabolites.get(met))
      else:
        netMet[met.id]=((rxn.flux/abs(rxn.flux))*rxn.metabolites.get(met))
    tempList = list()
    for k in netMet.keys():
        if netMet[k]==0:
            tempList.append(k)
    for k in tempList:
        del netMet[k]
  return netMet


#Function to print out all reactions generating/consuming a metabolite of inter-
#est
#args: 1) a solved cobra model 2) metabolite ID 3) ID of alternate charged state
#(use "" if none) 4) output file (use "" if no output file is required)
#output: none
def writeMetabSummary(cobra_model, met, Amet, outfile):
  met=cobra_model.metabolites.get_by_id(met)
  if not Amet == "":
    Amet=cobra_model.metabolites.get_by_id(Amet)
  if not outfile=="":
    fout = open(outfile,"w")
    fout.write("rxn ID\treaction\tmetabolite flux\n")
  for rxn in met.reactions:
    sto=rxn.metabolites.get(met)
    if Amet=="":
      sto1=0
    else:
      sto1=rxn.metabolites.get(Amet)
    if outfile=="":
      print(rxn.id+"\t"+rxn.reaction+"\t"+str(rxn.x*(sto+sto1)))
    else:
      fout.write(rxn.id+"\t"+rxn.reaction+"\t"+str(rxn.x*(sto+sto1))+"\n")
  if not outfile=="":
    fout.close()


#Function to calculate night time carbon conversion efficiency in diel model
#args: 1) a solved cobra model 2) day-night accumulation tag (default = "diel-
#-Transfer") 3) reaction ID for night time output (default = phloem_output_tx-
#-2), 4) reaction ID representing CO2 respired
#output: carbon conversion efficiency
def predictCCE(C3_model,accumulation_tag="dielTransfer",output="Phloem_output_tx2",CO2rxn = "CO2_tx2"):
  

  for met in C3_model.metabolites:
    if not met.formula:
      met.formula=""

  Cin = 0
  Cout = 0

  for rxn in C3_model.reactions.query(accumulation_tag):
    if round(rxn.x,5)>0:
      for met in rxn.products:
        #print met
        if met.formula.__contains__("C"):
          #print str(Cin)+"---"+met.id
          #print str(rxn.x)+"\t"+str(rxn.metabolites.get(met))+"\t"+str(int(re.split(r"[C,H]",met.formula)[1]))
          Cin = Cin + (rxn.x * rxn.metabolites.get(met) * int(re.split(r"[C,H]",met.formula)[1]))

  for rxn in C3_model.reactions.query(accumulation_tag):
    if round(rxn.x,5)<0:
      for met in rxn.reactants:
        if met.formula.__contains__("C"):
          #print str(Cout)+"---"+met.id
          #print str(rxn.x)+"\t"+str(rxn.metabolites.get(met))+"\t"+str(int(re.split(r"[C,H]",met.formula)[1]))
          Cout = Cout + (rxn.x * rxn.metabolites.get(met) * int(re.split(r"[C,H]",met.formula)[1]))

  rxn = C3_model.reactions.get_by_id(output)
  for met in rxn.reactants:
    if met.formula.__contains__("C"):
      #print str(Cout)+"---"+met.id
      #print str(rxn.x)+"\t"+str(-1 * rxn.metabolites.get(met))+"\t"+str(int(re.split(r"[C,H]",met.formula)[1]))
      Cout = Cout + (rxn.x * -1 * rxn.metabolites.get(met) * int(re.split(r"[C,H]",met.formula)[1]))

  Cout = Cout + (-1*C3_model.reactions.get_by_id(CO2rxn).x)

  if(not round(Cin,5) == round(Cout,5)):
    print("Error, Cin = "+str(Cin)+" and Cout = "+str(Cout))
    return 0
  else:
    print("Cin = "+str(Cin)+"\tCO2 ="+str(C3_model.reactions.get_by_id(CO2rxn).x)+"\t"+str(1 + ((C3_model.reactions.get_by_id(CO2rxn).x)/Cin)))
    return 1 + ((C3_model.reactions.get_by_id("CO2_tx2").x)/Cin)



####################################################################
#This function generates a tab seperated file that can be used with#
#Cytoscape to visualize metabolic flux                             #
#                                                                  #
#inputs: 1) a cobra model with feasible solution 2) the name of the#
#output file  3)  the number of cells in the model (eg: 2 for diel #
#C3 and 4 for diel C4)                                             #
#                                                                  #
####################################################################

def generateFluxMap(cobra_model,solution, outfile,phases = 2):
    #solution = cobra.flux_analysis.parsimonious.optimize_minimal_flux(cobra_model)
    #solution = cobra.flux_analysis.parsimonious.pfba(cobra_model)          #If the previous line returns error comment it out and uncomment this line instead

    #open output file for writing
    f = open(outfile,"w");

    #use rxnSet to identify reaction that have already been processed
    rxnSet = set()

    mult=set()
    #Looping through all reactions in the model
    for rxn in cobra_model.reactions:
        #Get the ID
        RXN=rxn.id
        #declare a boolean variable multFlag to keep track of whether this reaction is present in multiple models
        multFlag=False

        #check if the reaction has already been processed before and if yes skip this run in the loop
        if(rxnSet.__contains__(RXN)):
            continue
        if rxn.id.__contains__("EX") or rxn.id.__contains__("Transfer"):
            multFlag=False
        #check if the reaction ends with one or two i.e it is present more than once in the model
        elif(["1","2","3","4","5","6","7","8","9"].__contains__(rxn.id[len(rxn.id)-1])):
            #change the id to without the suffix 1-9 and declare it as a reaction which has multiple instances
            RXN = rxn.id[0:len(rxn.id)-1]
            multFlag=True
        elif rxn.id[len(rxn.id)-2:] == "10":
            #change the id to without the suffix 10 and declare it as a reaction which has multiple instances
            RXN = rxn.id[0:len(rxn.id)-2]
            multFlag=True

        #if metabolite has multiple instances
        values = dict()
        status1 = dict()
        status2 = dict()
        if(multFlag):
            tempvalue = list()
            temp1 = list()
            temp2 = list()
            mult.add(RXN)
            #add the reaction we are about to process to the reactions processed list
            for i in range(1,phases+1):
                rxnSet.add(RXN+str(i))
                if(round(float(solution.fluxes.get(RXN+str(i)))*10000000) == 0):
                    tempvalue.append(0)
                    temp1.append("none")
                    temp2.append("none")
                elif(float(solution.fluxes.get(RXN+str(i)))*10000 > 0):
                    tempvalue.append(solution.fluxes.get(RXN+str(i)))
                    temp1.append("produced")
                    temp2.append("consumed")
                elif(float(solution.fluxes.get(RXN+str(i)))*10000 < 0):
                    tempvalue.append(solution.fluxes.get(RXN+str(i)))
                    temp1.append("consumed")
                    temp2.append("produced")
            values[RXN] = tempvalue
            status1[RXN] = temp1
            status2[RXN] = temp2

            #select 1 reaction so that we can identify the reactants and products which can be then used to generate the edge shared_name
            rxn=cobra_model.reactions.get_by_id(RXN+"1")

            for reac in rxn.reactants:
                REAC=reac.id
                if(REAC.__contains__("1")):
                    if(REAC.rindex("1")==len(REAC)-1):
                        REAC=REAC[0:len(REAC)-1]
                    f.write("R_"+RXN+" (reaction-reactant) M_"+REAC)
                    for i in range(1,phases+1):
                        f.write("\t"+str(values[RXN][i-1])+"\t"+str(status2[RXN][i-1]))
                    f.write("\n")
                if(RXN.__contains__("biomass")):
                    f.write("R_"+RXN+" (reaction-product)) M_"+REAC)
                    for i in range(1,phases+1):
                        f.write("\t"+str(values[RXN][i-1])+"\t"+str(status1[RXN][i-1]))
                    f.write("\n")
            for prod in rxn.products:
                PROD=prod.id
                if(PROD.__contains__("1")):
                    if(PROD.rindex("1")==len(PROD)-1):
                        PROD=PROD[0:len(PROD)-1]
                f.write("R_"+RXN+" (reaction-product) M_"+PROD)
                for i in range(1,phases+1):
                    f.write("\t"+str(values[RXN][i-1])+"\t"+str(status1[RXN][i-1]))
                f.write("\n")
            if(RXN.__contains__("biomass")):
                f.write("R_"+RXN+" (reaction-reactant) M_"+REAC)
                for i in range(1,phases+1):
                    f.write("\t"+str(values[RXN][i-1])+"\t"+str(status2[RXN][i-1]))
                f.write("\n")
        else:
            #add the reaction we are about to process to the reactions processed list
            rxnSet.add(RXN)
            if(round(float(solution.fluxes.get(rxn.id))*10000000) == 0):
                value = 0;
                status1= "none";
                status0= "none";
            elif(solution.fluxes.get(rxn.id)*10000 > 0):
                value = solution.fluxes.get(rxn.id)*1000;
                status1= "produced";
                status0= "consumed";
            elif(solution.fluxes.get(rxn.id)*10000 < 0):
                value = solution.fluxes.get(rxn.id)*1000;
                status1= "consumed";
                status0= "produced";

            for reac in rxn.reactants:
                REAC=reac.id
                if(REAC.__contains__("1")):
                    if(REAC.rindex("1")==len(REAC)-1):# or (met.id.rindex("2")==len(rxn.id)-1):
                        REAC=REAC[0:len(REAC)-1]
                f.write("R_%s (reaction-reactant) M_%s\t%s\t%s\t0\tnone\n" % (RXN,REAC,value,status0));
            for prod in rxn.products:
                PROD=prod.id
                if(PROD.__contains__("1")):
                    if(PROD.rindex("1")==len(PROD)-1):# or (met.id.rindex("2")==len(rxn.id)-1):
                        PROD=PROD[0:len(PROD)-1]
                f.write("R_%s (reaction-product) M_%s\t%s\t%s\t0\tnone\n" % (RXN,PROD,value,status1));

    f.close();



####################################################
# This function estimates Rubisco carboxylase flux #
# at which the net CO2 uptake rate is equal to the #
# user defined value                               #
####################################################

def estimateVcFromNetCO2(model,netCO2uptake,Vc_ID="RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p1",CO2in_ID="CO2_tx1",verbose=False):

    
    # Initally constraint Vc flux to half current Vc_ID
    sol = flux_analysis.parsimonious.optimize_minimal_flux(model)

    if sol.fluxes[Vc_ID]==0:
        print('Warning no rubisco activity?')
        co2_vals={}
        co2_vals[sol.fluxes[Vc_ID]]=(netCO2uptake - sol.fluxes[CO2in_ID])/netCO2uptake
        fluxranges=np.logspace(-2,2,6)
        
        for val in fluxranges:
            model.reactions.get_by_id(Vc_ID).upper_bound = val
            model.reactions.get_by_id(Vc_ID).lower_bound = val
            try:
                sol= flux_analysis.parsimonious.optimize_minimal_flux(model)
                co2_vals[sol.fluxes[Vc_ID]]=(netCO2uptake - sol.fluxes[CO2in_ID])/netCO2uptake
            except:
                print('Stuck at ',val)
        popt, pcov = curve_fit(Hyp_func, list(co2_vals.keys()), list(co2_vals.values()),maxfev=10000)
        prev=(-popt[1])/popt[0]
        # model_feasible=0
        # count=0
        # while not model_feasible:
        #     try:
        #         popt, pcov = curve_fit(Hyp_func, list(co2_vals.keys())[count:], list(co2_vals.values()[count:]),maxfev=10000)
        #         prev=(-popt[1])/popt[0]
        #         model.reactions.get_by_id(Vc_ID).lower_bound = prev
        #         model.reactions.get_by_id(Vc_ID).upper_bound = prev
        #         sol = flux_analysis.parsimonious.optimize_minimal_flux(model)
        #         model_feasible=1
        #         count+=1
        #     except:
        #         now = prev + ((netCO2uptake -sol.fluxes[CO2in_ID]))/4
        #         if now<0:
        #             now=0
        #         print(now)
        #         model.reactions.get_by_id(Vc_ID).lower_bound = now
        #         model.reactions.get_by_id(Vc_ID).upper_bound = now
        #         model_feasible=1
        #         sol = flux_analysis.parsimonious.optimize_minimal_flux(model)
        #         co2_vals[sol.fluxes[Vc_ID]]=(netCO2uptake - sol.fluxes[CO2in_ID])/netCO2uptake
        return prev
    else:
        model.reactions.get_by_id(Vc_ID).lower_bound = sol.fluxes[Vc_ID]/2
        model.reactions.get_by_id(Vc_ID).upper_bound = sol.fluxes[Vc_ID]/2

        #perform pFBA
        # If it doesn't work, just start with original Vc_ID and more iterations
        try:
            sol = flux_analysis.parsimonious.optimize_minimal_flux(model)
            max_it=4
        except:
            max_it=6
            

        #set loop counter
        i=0

        #Use a while loop to decrease/increase Vc flux until net CO2 rate is similar to given value (or loop counter hits 10)
        prev = sol.fluxes[Vc_ID]
        co2_vals={}
        while((netCO2uptake - sol.fluxes[CO2in_ID])/netCO2uptake > 0.001 and i<max_it):
            i=i+1
            prev = sol.fluxes[Vc_ID]
            # Increment in Vc flux is set by given netCo2 uptake - model predicted CO2 uptake rate in previous pFBA run
            now = prev + ((netCO2uptake -sol.fluxes[CO2in_ID]))
            model.reactions.get_by_id(Vc_ID).lower_bound = now
            model.reactions.get_by_id(Vc_ID).upper_bound = now
            sol = flux_analysis.parsimonious.optimize_minimal_flux(model)
            co2_vals[sol.fluxes[Vc_ID]]=(netCO2uptake - sol.fluxes[CO2in_ID])/netCO2uptake
            if verbose:
                print("----"+str(i)+"----")
                print("Vc flux ="+str(now))
                print("net CO2 uptake ="+str(sol.fluxes[CO2in_ID]))
                print("Target CO2 uptake ="+str(netCO2uptake))
                
        if (-(netCO2uptake - sol.fluxes[CO2in_ID])/netCO2uptake) > 0.001:
            print('Warning '+CO2in_ID+' too large:', sol.fluxes[CO2in_ID])
            while(-(netCO2uptake - sol.fluxes[CO2in_ID])/netCO2uptake > 0.001 and i<10):
                i=i+1
                prev = sol.fluxes[Vc_ID]
                # Increment in Vc flux is set by given netCo2 uptake - model predicted CO2 uptake rate in previous pFBA run
                now = prev + ((netCO2uptake -sol.fluxes[CO2in_ID]))/4
                if now<0:
                    now=0
                print(i, now)
                model.reactions.get_by_id(Vc_ID).lower_bound = now
                model.reactions.get_by_id(Vc_ID).upper_bound = now
                sol = flux_analysis.parsimonious.optimize_minimal_flux(model)
                co2_vals[sol.fluxes[Vc_ID]]=(netCO2uptake - sol.fluxes[CO2in_ID])/netCO2uptake
                if verbose:
                    print("----"+str(i)+"----")
                    print("Vc flux ="+str(now))
                    print("net CO2 uptake ="+str(sol.fluxes[CO2in_ID]))
                    print("Target CO2 uptake ="+str(netCO2uptake))
        print(list(co2_vals.keys()), list(co2_vals.values()))
        popt, pcov = curve_fit(Hyp_func, list(co2_vals.keys()), list(co2_vals.values()),maxfev=10000)
        prev=(-popt[1])/popt[0]
        model_feasible=0
        count=0
        while not model_feasible:
            try:
                popt, pcov = curve_fit(Hyp_func, list(co2_vals.keys())[count:], list(co2_vals.values()[count:]),maxfev=10000)
                prev=(-popt[1])/popt[0]
                model.reactions.get_by_id(Vc_ID).lower_bound = prev
                model.reactions.get_by_id(Vc_ID).upper_bound = prev
                sol = flux_analysis.parsimonious.optimize_minimal_flux(model)
                model_feasible=1
                count+=1
            except:
                now = prev + ((netCO2uptake -sol.fluxes[CO2in_ID]))/4
                if now<0:
                    now=0
                print(now)
                model.reactions.get_by_id(Vc_ID).lower_bound = now
                model.reactions.get_by_id(Vc_ID).upper_bound = now
                model_feasible=1
                sol = flux_analysis.parsimonious.optimize_minimal_flux(model)
                co2_vals[sol.fluxes[Vc_ID]]=(netCO2uptake - sol.fluxes[CO2in_ID])/netCO2uptake
        return prev

def Hyp_func(x, a,b):
    return a*x+b
######################################################
# This function estimates biomass/phloem output flux #
# at which the net CO2 uptake rate is equal to the   #
# user defined value                                 #
####################################################
def estimateOutputFromNetCO2(model,netCO2uptake,Output_ID="diel_biomass",Vc_ID="RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p1",CO2in_ID="CO2_tx1",verbose=False):
    # Initally constraint Vc flux to net CO2 uptake rate
    print(model.reactions.get_by_id("Photon_tx_l").upper_bound)
    print(model.reactions.get_by_id("ATPase_tx_MC_l").upper_bound)
    model.reactions.get_by_id(Vc_ID).lower_bound = netCO2uptake
    model.reactions.get_by_id(Vc_ID).upper_bound = netCO2uptake
    print(model.reactions.get_by_id(Vc_ID).upper_bound)
    #perform pFBA
    sol = flux_analysis.parsimonious.pfba(model)

    #unconstrain Vc
    model.reactions.get_by_id(Vc_ID).lower_bound = 0
    model.reactions.get_by_id(Vc_ID).upper_bound = 1000

    #set loop counter
    i=0

    #Use a while loop to increase Vc flux until net CO2 rate is similar to given value (or loop counter hits 10)
    while((netCO2uptake - sol.fluxes[CO2in_ID])/netCO2uptake > 0.001):# and i<10):
        i=i+1
        prev = sol.fluxes[Output_ID]
        # Increment in Vc flux is set by given netCo2 uptake - model predicted CO2 uptake rate in previous pFBA run
        now = prev + (prev*((netCO2uptake - sol.fluxes[CO2in_ID])/netCO2uptake))
        model.reactions.get_by_id(Output_ID).lower_bound = now
        model.reactions.get_by_id(Output_ID).upper_bound = now

        flux_analysis.parsimonious.pfba(model)
        if verbose:
            print("----"+str(i)+"----")
            print("Vc flux ="+str(model.reactions.get_by_id(Vc_ID).x))
            print("net CO2 uptake ="+str(sol.fluxes[CO2in_ID]))
            print("Target CO2 uptake ="+str(netCO2uptake))
            print("Before:"+str(prev))
            print("After:"+str(now))
    return prev

######################################################
# This function estimates biomass/phloem output flux #
# at which the Vc flux is equal to the user defined  #
# value                                              #
######################################################
def estimateOutputFromVc(model,Vc,Output_ID="diel_biomass",Vc_ID="RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p1",verbose=False):
    # Initally constraint Vc flux to net CO2 uptake rate
    model.reactions.get_by_id(Vc_ID).lower_bound = Vc/2
    model.reactions.get_by_id(Vc_ID).upper_bound = Vc/2
    if verbose:
        print(model.reactions.get_by_id("Photon_tx1").upper_bound)
        print(model.reactions.get_by_id("ATPase_tx1").upper_bound)
        print(model.reactions.get_by_id(Vc_ID).upper_bound)
    #perform pFBA
    sol = flux_analysis.parsimonious.pfba(model)

    #unconstrain Vc
    model.reactions.get_by_id(Vc_ID).lower_bound = 0
    model.reactions.get_by_id(Vc_ID).upper_bound = 1000

    #set loop counter
    i=0

    #Use a while loop to increase Vc flux until net CO2 rate is similar to given value (or loop counter hits 10)
    while((Vc - sol.fluxes[Vc_ID])/Vc > 0.001):# and i<10):
        i=i+1
        prev = sol.fluxes[Output_ID]
        # Increment in Vc flux is set by given netCo2 uptake - model predicted CO2 uptake rate in previous pFBA run
        now = prev + (prev*((Vc - model.reactions.get_by_id(Vc_ID).x)/Vc))
        model.reactions.get_by_id(Output_ID).lower_bound = now
        model.reactions.get_by_id(Output_ID).upper_bound = now

        sol = flux_analysis.parsimonious.pfba(model)
        if verbose:
            print("----"+str(i)+"----")
            print("Vc flux ="+str(sol.fluxes[Vc_ID]))
            print("Target Vc flux ="+str(Vc))
            print("Before:"+str(prev))
            print("After:"+str(now))
    return prev


#####################################################################
# This function prints information on CO2 generating fluxes and can #
# be used to study dark respiration                                 #
# inputs: 1) an FBA solved model, 2) a list of fluxs to ignore, 3)  #
# the tag used to mark night time metabolites ex: "_night", 4)      #
# a list of solution objects                                        #
#####################################################################
def printDarkRespirationFluxes(model,rxn2avoid=["CO2_tx1","CO2_ec1","CO2_mc1","CO2_pc1","GCVMULTI_RXN_m1"],night_tag = "2",custom_solutions=[]):
  for met in model.metabolites.query("CARBON_DIOXIDE"):
    if met.id.__contains__(night_tag):
      continue
    for rxn in met.reactions:
      if rxn2avoid.__contains__(rxn.id):
        continue
      if len(custom_solutions) == 0:
        if (rxn.metabolites.get(met) > 0 and rxn.x > 0) or (rxn.metabolites.get(met) < 0 and rxn.x < 0):
          print(rxn.id+"\t"+rxn.reaction+"\t"+str(rxn.x*rxn.metabolites.get(met)))
        else:
          if (rxn.metabolites.get(met) > 0 and rxn.upper_bound > 0) or (rxn.metabolites.get(met) < 0 and rxn.lower_bound < 0):
            tot=0
            for sol in custom_solutions:
              tot=tot+abs(sol.x_dict.get(rxn.id))
              if tot==0:
                continue

            print(rxn.id+"\t"+rxn.reaction,)
            for sol in custom_solutions:
              print("\t"+str(rxn.metabolites.get(met)*sol.x_dict.get(rxn.id)),)
            print("")

#####################################################################
# This function generates ATP budgets for a given flux distribution #
# inputs: 1) an FBA model, 2) a dictionary object with reaction ids #
# as keys and reaction fluxes as values, 3) name of output file (op-#
# -tional), 4) Option to show plots, 5) If choosing to show plot, c-#
# -hoose wether to use percentage or absolute values in the plot. 6)#
# Provide a day or night indicator tag to specify day or night ATP  #
# summary 7) a destination file to save plot to 8) a dictionary to  #
# specify colour for fluxes in plot                                 #
#####################################################################

def generateATPbudget(model,solution,outfile="",show_plot=True,percentage=False,day_or_night_tag="_MC_l",save_plot_to="temp.png",pltrange=[],colourDict={},keep_same=[],outdict=None):
  if outfile!="":
    fout = open(outfile,"w")
  ATPdict = dict()
  total = 0
  for p in ("c","p","m","x"):
    met=model.metabolites.get_by_id("ATP_"+p+day_or_night_tag)
    met1=model.metabolites.get_by_id("aATP_"+p+day_or_night_tag)
    for rxn in met.reactions:
      if rxn.id.__contains__("ATP_AMP_mc") or rxn.id.__contains__("ATP_ADP_mc") or rxn.id.__contains__("ATP_pc") or rxn.id.__contains__("AMP_ATP_xc") or rxn.id.__contains__("ATP_ADP_Pi_pc"):
        continue
      sto=rxn.metabolites.get(met)
      sto1=rxn.metabolites.get(met1)
      if not sto1:
        sto1=0
      if outfile!="":
        fout.write(rxn.id+"\t"+rxn.reaction+"\t"+str(solution[rxn.id]*(sto+sto1))+"\t"+met.compartment+"\n")
      ATPdict[rxn.id]=solution[rxn.id]*(sto+sto1)
      if solution[rxn.id]*(sto+sto1) > 0:
        total = total + (solution[rxn.id]*(sto+sto1))
  if outfile!="":
    fout.close()

  tempDict = dict()
  for rxn in ATPdict.keys():
    tempDict[rxn]=abs(ATPdict[rxn])

  #sort ATPdict by values
  sorted_by_value = sorted(tempDict.items(), key= lambda x:x[1],reverse=True)

  ATPdict2 = OrderedDict()
  ATPdict2["Others-pos"]=0
  ATPdict2["Others-neg"]=0
  baseline = dict()
  pos_base=0
  neg_base=0
  i=0
  other_thresh=[sorted_by_value[x][1] for x in range(17,1,-1) if len([y for y in sorted_by_value if y[1]>sorted_by_value[x][1]])<=17]
  other_thresh=max(other_thresh+[5e-5])
  for TEMP in sorted_by_value:
    rxn = TEMP[0]
    if ATPdict[rxn]>0:
      if ATPdict[rxn] < other_thresh:
        if percentage:
          ATPdict2["Others-pos"]=ATPdict2["Others-pos"]+float(ATPdict[rxn]*100)/total
        else:
          ATPdict2["Others-pos"]=ATPdict2["Others-pos"]+ATPdict[rxn]
        continue
      base = pos_base
      if percentage:
        ATPdict2[rxn]=float(ATPdict[rxn]*100)/total
        pos_base = pos_base + float(ATPdict[rxn]*100)/total
      else:
        pos_base = pos_base + ATPdict[rxn]
        ATPdict2[rxn]=ATPdict[rxn]
    else:
      if abs(ATPdict[rxn]) < other_thresh:
        if percentage:
          ATPdict2["Others-neg"]=ATPdict2["Others-neg"]+float(ATPdict[rxn]*100)/total
        else:
          ATPdict2["Others-neg"]=ATPdict2["Others-neg"]+ATPdict[rxn]
        continue
      base = neg_base
      if percentage:
        ATPdict2[rxn]=float(ATPdict[rxn]*100)/total
        neg_base = neg_base + float(ATPdict[rxn]*100)/total
      else:
        neg_base = neg_base + ATPdict[rxn]
        ATPdict2[rxn]=ATPdict[rxn]
    i=i+1
    baseline[rxn]=base
  baseline["Others-pos"]=pos_base
  baseline["Others-neg"]=neg_base
  # Bring common reactions to the front
  tempDict=ATPdict2.copy()
  for rxn in keep_same:
    for key in ATPdict2.keys():
      if rxn in key:
        tempDict.move_to_end(key,last=False)
  ATPdict2=tempDict
  if show_plot:
    plt.rcParams.update({'font.size': 10}) #sets a global fontsize
    plt.rcParams['xtick.major.size'] = 5 # adjusts tick line length and width
    plt.rcParams['xtick.major.width'] = 1
    plt.rcParams['ytick.major.size'] = 5
    plt.rcParams['ytick.major.width'] = 1
    plt.rcParams['axes.linewidth']=2 # makes axes line thicker
    thisplot = plt.figure(figsize=(3,4))
    NUM_COLORS=len(ATPdict2.keys())
    if NUM_COLORS<=20:
        cm = plt.get_cmap('tab20')
    else:
        cm = plt.get_cmap('gist_ncar')
    for i,rxn in enumerate(ATPdict2.keys()):
      if colourDict.keys().__contains__(rxn):
        plt.bar(1,ATPdict2[rxn],width=0.1,bottom=baseline[rxn],label=rxn,color=colourDict[rxn])
      else:
        plt.bar(1,ATPdict2[rxn],width=0.1,bottom=baseline[rxn],label=rxn,color=cm(1.*i/NUM_COLORS))
    plt.xlim(0.8,1.2)
    if percentage:
      plt.ylabel("ATP produced/consumed (%)")
    else:
      plt.ylabel("ATP produced/consumed (in moles)")
    handles, labels = plt.gca().get_legend_handles_labels()
    labels2=list(set(labels)-set(["Others-neg","Others-pos"]))+list(["Others-neg","Others-pos"])
    handles2=[handles[labels.index(i)] for i in labels2]
    lgd=plt.legend(handles2,labels2,bbox_to_anchor=(1,1))
    plt.axhline(0,linestyle="--",color="black")
    plt.tight_layout
    if pltrange:
        plt.ylim(-pltrange,pltrange)
    plt.savefig(save_plot_to, bbox_extra_artists=(lgd,), bbox_inches='tight',format="png")
    plt.show()
    if outdict:
        return thisplot, ATPdict
    else:
        return thisplot


#####################################################################
# This function generates ATP budgets for a given flux distribution #
# inputs: 1) an FBA model, 2) a dictionary object with reaction ids #
# as keys and reaction fluxes as values, 3) name of output file (op-#
# -tional), 4) Option to show plots, 5) If choosing to show plot, c-#
# -hoose wether to use percentage or absolute values in the plot 6) #
# Provide a day or night indicator tag to specify day or night NAD(-#
# -P)H summary 7) a destination file to save plot to 8) a dictionary#
# to specify colour for fluxes in plot                              #
#####################################################################
def generateNADHNADPHbudget(model,solution,outfile="",show_plot=True,percentage=False,day_or_night_tag="1",save_plot_to="temp",colourDict={}):
    if outfile!="":
        fout = open(outfile,"w")
    Reddict = dict()
    total = 0
    for red in ["NADPH","NADH"]:
        for p in ("c","p","m","x"):
            if len(model.metabolites.query(red+"_"+p+day_or_night_tag))==0:
                continue
            met=model.metabolites.get_by_id(red+"_"+p+day_or_night_tag)
            for rxn in met.reactions:
                sto=rxn.metabolites.get(met)
                sto1=0#rxn.metabolites.get(met1)
                if outfile!="":
                    fout.write(rxn.id+"\t"+rxn.reaction+"\t"+str(solution[rxn.id]*(sto+sto1))+"\t"+met.compartment+"\n")
                Reddict[rxn.id]=solution[rxn.id]*(sto+sto1)
                if solution[rxn.id]*(sto+sto1) > 0:
                    total = total + (solution[rxn.id]*(sto+sto1))
    if outfile!="":
        fout.close()

    tempDict = dict()
    for rxn in Reddict.keys():
      tempDict[rxn]=abs(Reddict[rxn])

    #sort by values
    sorted_by_value = sorted(tempDict.items(), key= lambda x:x[1],reverse=True)



    Reddict2 = dict()
    Reddict2["Others-pos"]=0
    Reddict2["Others-neg"]=0
    baseline = dict()
    pos_base=0
    neg_base=0
    i=0
    other_thresh=[sorted_by_value[x][1] for x in range(17,1,-1) if len([y for y in sorted_by_value if y[1]>sorted_by_value[x][1]])<=17]
    other_thresh=max(other_thresh+[5e-5])
    for TEMP in sorted_by_value:
        rxn = TEMP[0]
        if Reddict[rxn]>0:
            if Reddict[rxn] < other_thresh:
                if percentage:
                    Reddict2["Others-pos"]=Reddict2["Others-pos"]+float(Reddict[rxn]*100)/total
                else:
                    Reddict2["Others-pos"]=Reddict2["Others-pos"]+Reddict[rxn]
                continue
            base = pos_base
            if percentage:
                Reddict2[rxn]=float(Reddict[rxn]*100)/total
                pos_base = pos_base + float(Reddict[rxn]*100)/total
            else:
                pos_base = pos_base + Reddict[rxn]
                Reddict2[rxn]=Reddict[rxn]
        else:
            if abs(Reddict[rxn]) < other_thresh:
                if percentage:
                    Reddict2["Others-neg"]=Reddict2["Others-neg"]+float(Reddict[rxn]*100)/total
                else:
                    Reddict2["Others-neg"]=Reddict2["Others-neg"]+Reddict[rxn]
                continue
            base = neg_base
            if percentage:
                Reddict2[rxn]=float(Reddict[rxn]*100)/total
                neg_base = neg_base + float(Reddict[rxn]*100)/total
            else:
                neg_base = neg_base + Reddict[rxn]
                Reddict2[rxn]=Reddict[rxn]
        i=i+1
        baseline[rxn]=base
    baseline["Others-pos"]=pos_base
    baseline["Others-neg"]=neg_base

    if show_plot:
        plt.rcParams.update({'font.size': 10}) #sets a global fontsize
        plt.rcParams['xtick.major.size'] = 5 # adjusts tick line length and width
        plt.rcParams['xtick.major.width'] = 1
        plt.rcParams['ytick.major.size'] = 5
        plt.rcParams['ytick.major.width'] = 1
        plt.rcParams['axes.linewidth']=2 # makes axes line thicker
        thisplot = plt.figure(figsize=(3,4))
        NUM_COLORS=len(Reddict2.keys())
        if NUM_COLORS<=20:
            cm = plt.get_cmap('tab20')
        else:
            cm = plt.get_cmap('gist_ncar')
        for i,rxn in enumerate(Reddict2.keys()):
            if colourDict.keys().__contains__(rxn):
              plt.bar(1,Reddict2[rxn],width=0.1,bottom=baseline[rxn],label=rxn,color=colourDict[rxn])
            else:
              plt.bar(1,Reddict2[rxn],width=0.1,bottom=baseline[rxn],label=rxn,color=cm(1.*i/NUM_COLORS))
        plt.xlim(0.8,1.2)
        if percentage:
            plt.ylabel("NAD(P)H produced/consumed (%)")
        else:
            plt.ylabel("NAD(P)H produced/consumed (in moles)")
        handles, labels = plt.gca().get_legend_handles_labels()
        labels2=list(set(labels)-set(["Others-neg","Others-pos"]))+list(["Others-neg","Others-pos"])
        handles2=[handles[labels.index(i)] for i in labels2]
        lgd=plt.legend(handles2,labels2,bbox_to_anchor=(1,1))
        plt.axhline(0,linestyle="--",color="black")
        plt.tight_layout
        plt.savefig(save_plot_to, bbox_extra_artists=(lgd,), bbox_inches='tight',format="png")
#         plt.savefig(save_plot_to, bbox_inches='tight',format="png")
        plt.show()
        return thisplot

def generateuniformbudget(model,solution,metin='ATP',altmets=['aATP'],cells=['_MC_','_CC_','_SE_','_pSE_'],lightphase=['l'],diflims=[],outfile="",show_plot=True,percentage=False,save_plot_to=None,colourDict={},outdict=None,numindict=14,names=False,fixedy=1.05,ylims=None,thresh=5e-5,coloradj=4,percell='',filetype='png',comps=("c","p","m","x")):
  addnames(model)
  # Things to add:
  # Percentages instead of flux values
  rxnexceptions = ['PROTON_PPI_ec','turnover','ATPase_tx','_ec_']
  if len(comps)==1:
      rxnexceptions = rxnexceptions+['_pc','_mc_','_xc_','_ec_','_vc_']
  cellmult = [1,20,4,10/3]
  # Collect ATP dict for all cells
  metdicts={}
  total={}
  for cell in cells:
    for phase in lightphase:
      celltag=cell+phase
      metdicts[celltag]={}
      total[celltag] = 0
      for p in comps:
        thismet=model.metabolites.get_by_id(metin+"_"+p+celltag)
        thisalts=[]
        for altmet in altmets:
          if any([altmet+"_"+p+celltag==met.id for met in model.metabolites]):
            thisalts+=[model.metabolites.get_by_id(altmet+"_"+p+celltag)]
        for rxn in thismet.reactions:
          if (len(rxn.id.split('_')[-3])>1 and not any([y in rxn.id for y in rxnexceptions])) or solution[rxn.id]==0:
            continue
          sto=rxn.metabolites.get(thismet)
          altsto=0
          for altmet in thisalts:
            if altmet in rxn.metabolites:
              altsto+=rxn.metabolites.get(altmet)
          metdicts[celltag][rxn.id]=solution[rxn.id]*(sto+altsto)
          if percell:
              if 'MC' in celltag:
                cellmultiplier = cellmult[0]
              elif 'CC' in celltag:
                cellmultiplier = cellmult[1]
              elif 'pSE' in celltag:
                cellmultiplier = cellmult[3]
              elif 'SE' in celltag:
                cellmultiplier = cellmult[2]
              else:
                  print('Couldnt recognise cell: ',celltag)
              metdicts[celltag][rxn.id]=metdicts[celltag][rxn.id]*cellmultiplier
          if solution[rxn.id]*(sto+altsto) > 0:
            total[celltag] += (solution[rxn.id]*(sto+altsto))
        altrxns = list(set([list(x.reactions)[ii] for x in thisalts for ii in range(len(x.reactions))]).difference(set(thismet.reactions)))
        for rxn in altrxns:
          if ((len(rxn.id.split('_')[-3])>1 and not any([y in rxn.id for y in rxnexceptions])) or solution[rxn.id]==0):
            continue
          altsto=0
          for altmet in thisalts+[thismet]:
            if altmet in rxn.metabolites:
                altsto+=rxn.metabolites.get(altmet)
          metdicts[celltag][rxn.id]=solution[rxn.id]*(sto+altsto)
          if percell:
            if 'MC' in celltag:
                cellmultiplier = cellmult[0]
            elif 'CC' in celltag:
                cellmultiplier = cellmult[1]
            elif 'pSE' in celltag:
                cellmultiplier = cellmult[3]
            elif 'SE' in celltag:
                cellmultiplier = cellmult[2]
            else:
                print('Couldnt recognise cell: ',celltag)
            metdicts[celltag][rxn.id]=metdicts[celltag][rxn.id]*cellmultiplier
          if solution[rxn.id]*(altsto) > 0:
            total[celltag] += (solution[rxn.id]*(sto+altsto))
    
  if names:
    namedict={}
    for rxn in set([x for cell in list(metdicts.keys())for x in list(metdicts[cell].keys())]):
      tempname=model.reactions.get_by_id(rxn).name.replace('--','-').replace(',',' ')
      if ':' in tempname:
        tempname=tempname.split(':')[1]
      if tempname=='NIL' or not tempname:
        tempname=rxn
      tempname=tempname[0].upper()+tempname[1:]
      if not any([y in rxn for y in rxnexceptions]):
        tempname=tempname+' ('+rxn.split('_')[-3]+')'
      namedict['_'.join(rxn.split('_')[:-2])]=tempname+' '*(62-len(tempname))
  # Work out threshold for 'other' reactions
  # Want no more than 17 different reactions per cell.
  not17={}
  for cell in list(metdicts.keys()):
    not17[cell]=sorted(metdicts[cell].items(),key=lambda x: abs(x[1]),reverse=True)[numindict:]
  combnot17 = []
  allrxns = []
  rxnstosort=list(set([x.replace(cells[1],cells[0]) for x in list(metdicts[list(metdicts.keys())[1]].keys())]+list(metdicts[list(metdicts.keys())[0]].keys())))
#   sortedrxnlist=sorted(rxnstosort,key=lambda x:sum([metdicts[y][x.replace('_'.join(x.split('_')[-2:]),'_'.join(y.split('_')[-2:]))] if x.replace('_'.join(x.split('_')[-2:]),'_'.join(y.split('_')[-2:])) in metdicts[y].keys() else 0 for y in list(metdicts.keys())]))
  sortedrxnlist=sorted(rxnstosort,key=lambda x:sum([metdicts[y][x.replace('_'.join(x.split('_')[-2:]),'_'.join(y.split('_')[-2:]))] if x.replace('_'.join(x.split('_')[-2:]),'_'.join(y.split('_')[-2:])) in metdicts[y].keys() else 0 for y in ['_CC_l']]))
#   print(sortedrxnlist)
#   print([metdicts['_CC_l'][x.replace('MC','CC')] if x.replace('MC','CC') in metdicts['_CC_l'].keys() else 0 for x in sortedrxnlist])
#   print([metdicts['_MC_l'][x] if x in metdicts['_MC_l'].keys() else 0 for x in sortedrxnlist])
  for cell in list(not17.keys()):
    combnot17+=['_'.join(x.split('_')[:-2]) for x in [y[0] for y in not17[cell]]]
    allrxns +=['_'.join(x.split('_')[:-2]) for x in sortedrxnlist]
  toprxns = set(allrxns).difference(set(combnot17))
  toprxns = sorted(list(toprxns),key=lambda x:sortedrxnlist.index(x+list(metdicts.keys())[0]))
  if 'MC_l' in list(metdicts.keys()):
    midpoints = [rxn for rxn in toprxns if (rxn+'_CC_l' not in metdicts['_CC_l'].keys() or metdicts['_CC_l'][rxn+'_CC_l']==0) and metdicts['_MC_l'][rxn+'_MC_l']>0]
  else:
    midpoints = [rxn for rxn in toprxns if (rxn+'_CC_l' not in metdicts['_CC_l'].keys() or metdicts['_CC_l'][rxn+'_CC_l']==0)]
#   print([rxn for rxn in toprxns if metdicts['_MC_l'][rxn+'_MC_l']>0])
#   print([rxn+': '+str(metdicts['_MC_l'][rxn+'_MC_l']) for rxn in toprxns])
  if not midpoints:
    #   print(2)
      midpoints = [rxn for rxn in toprxns if rxn+'_CC_l' in metdicts['_CC_l'].keys() and metdicts['_CC_l'][rxn+'_CC_l']>0]
#   print(midpoints)
  midpoint = max(0,toprxns.index(midpoints[0])-1)
#   print(toprxns[midpoint])
#   toprxns = sorted(list(toprxns),key=lambda x:metdicts[celltag][x])
  alldict=metdicts.copy()
  # Replace all reactions with low flux with 'other', make baseline dict
  for cell in list(metdicts.keys()):
    cellrxns=list(metdicts[cell].keys())
    metdicts[cell]['Others-pos']=0
    metdicts[cell]['Others-neg']=0
    for rxn in cellrxns:
      if not any([x in rxn for x in toprxns]):
        # print(rxn)
        rxnval=metdicts[cell].pop(rxn)
        if rxnval>0:
          metdicts[cell]['Others-pos']+=rxnval
        elif rxnval<0:
          # print(metdicts[cell].keys())
          metdicts[cell]['Others-neg']+=rxnval

  metdicts2={}
  for cell in list(metdicts.keys()):
    metdicts2[cell]={}
    for rxn in toprxns:
      if any([rxn+cell==y for y in list(metdicts[cell].keys())]):
        metdicts2[cell][rxn]=metdicts[cell][rxn+cell]
      else:
        metdicts2[cell][rxn]=0
    for rxn in ['Others-pos','Others-neg']:
      metdicts2[cell][rxn]=metdicts[cell][rxn]
  baseline={}
  maxneg=0
  maxpos=0
  for cell in list(metdicts2.keys()):
    baseline[cell]={}
    posbase=0
    negbase=0
    for rxn in toprxns[midpoint::-1]:
      if metdicts2[cell][rxn]>thresh:
        base=posbase
        posbase+=metdicts2[cell][rxn]
      elif metdicts2[cell][rxn]<-thresh:
        base=negbase
        negbase+=metdicts2[cell][rxn]
      else:
        base=posbase
      baseline[cell][rxn]=base
    #   print(rxn,metdicts2[cell][rxn],baseline[cell][rxn],posbase,negbase)
    for rxn in toprxns[(midpoint+1):]:
      if metdicts2[cell][rxn]>thresh:
        base=posbase
        posbase+=metdicts2[cell][rxn]
      elif metdicts2[cell][rxn]<-thresh:
        base=negbase
        negbase+=metdicts2[cell][rxn]
      else:
        base=posbase
      baseline[cell][rxn]=base
    baseline[cell]["Others-pos"]=posbase
    baseline[cell]["Others-neg"]=negbase
    # print(rxn,metdicts2[cell][rxn],baseline[cell][rxn],posbase,negbase)
    if posbase>maxpos and cell not in diflims:
        maxpos=posbase+metdicts[cell]['Others-pos']
    if maxneg>negbase and cell not in diflims:
        maxneg=negbase+metdicts[cell]['Others-neg']
  if show_plot:
    plt.rcParams.update({'font.size': 16}) #sets a global fontsize
    plt.rcParams['xtick.major.size'] = 5 # adjusts tick line length and width
    plt.rcParams['xtick.major.width'] = 1
    plt.rcParams['ytick.major.size'] = 5
    plt.rcParams['ytick.major.width'] = 1
    plt.rcParams['axes.linewidth']=2 # makes axes line thicker
    mpf=max((midpoint-coloradj)/(len(toprxns)+2),0.05)
    # print(mpf) 
    cdict = {'blue':   [[0.0,  174.0/255, 174.0/255],
                        [mpf-0.05,  236/255, 236/255],
                        [mpf+0.05,  220/255, 220/255],
                        [1.0,  84/255, 84/255]],
                'red':    [[0.0,  110.0/255, 110.0/255],
                        [mpf-0.05, 218/255, 218/255],
                        [mpf+0.05,  220/255, 220/255],
                        [1.0,  48/255, 48/255]],
                'green':  [[0.0,  100.0/255, 100.0/255],
                        [mpf-0.05,  218/255, 218/255],
                        [mpf+0.05,  243/255, 243/255],
                        [1.0,  162/255, 162/255]]}
    newcmp = LinearSegmentedColormap('testCmap', segmentdata=cdict, N=256)
    cm=newcmp
    for cell in list(metdicts.keys()):
      thisplot = plt.figure(figsize=(3,4))
    #   NUM_COLORS=len(metdicts2[cell].keys())
      NUM_COLORS=len(toprxns)+2
      for i,rxn in enumerate((['Others-neg']+toprxns+['Others-pos'])[-1::-1]):
        if abs(metdicts2[cell][rxn])>thresh:
          if colourDict.keys().__contains__(rxn):
            plt.bar(1,metdicts2[cell][rxn],width=0.1,bottom=baseline[cell][rxn],label=rxn,color=colourDict[rxn])
          else:
            mp=0 #midpoint-2
            plt.bar(1,metdicts2[cell][rxn],width=0.1,bottom=baseline[cell][rxn],label=rxn,color=cm(i/NUM_COLORS))
            # plt.bar(1,metdicts2[cell][rxn],width=0.1,bottom=baseline[cell][rxn],label=rxn,color=cm((0.5*i+0.5*NUM_COLORS+mp)/(NUM_COLORS-mp)-(0.5*NUM_COLORS+mp)/(NUM_COLORS-mp)))
      plt.xlim(0.8,1.2)
      plt.tick_params('x',labelbottom=0,bottom=0)
      plt.box(False)
      if percentage:
        plt.ylabel(metin+" produced/consumed (%)")
      else:
        plt.ylabel(metin+" produced/consumed \n(µmol/m$^2$/s)")
      handles, labels = plt.gca().get_legend_handles_labels()
      if names:
        lgd=plt.legend(handles,[namedict[x] if x in list(namedict.keys()) else x for x in labels],bbox_to_anchor=(1.05,1),loc='upper left',frameon=False)
      else:
        lgd=plt.legend(handles,labels,bbox_to_anchor=(1.05,1),loc='upper left',frameon=False)
      plt.axhline(0,linestyle="--",color="black")
      plt.title(cell[1:])
      if fixedy and cell not in diflims:
        plt.ylim(fixedy*maxneg,fixedy*maxpos)
      if ylims:
        plt.ylim(ylims[0],ylims[1])
      plt.tight_layout
      if save_plot_to:
        if not outfile:
          thisoutfile='Figures/budget'+metin+cell+percell+'.'+filetype
        else:
          thisoutfile='Figures/'+outfile+'_budget'+metin+cell+percell+'.'+filetype
        plt.savefig(thisoutfile, bbox_extra_artists=(lgd,), bbox_inches='tight',format=filetype)
      plt.show()
    if outdict:
        return thisplot, metdicts2,alldict
    else:
        return thisplot
  return metdicts2

def addnames(stem_model):
    for rxn in stem_model.reactions:
        if 'PROTON_PPI' in rxn.id:
            rxn.name = 'H$^+$-PP$_i$ase'
        if 'PROTON_ATPase_c' in rxn.id:
            rxn.name = 'H$^+$-ATPase'
        if 'ATPPHOSPHORIBOSYLTRANS' in rxn.id:
            rxn.name = 'ATP-phosphoribosyl transferase'
        if 'UDPKIN' in rxn.id:
            rxn.name = 'UDP kinase'
        if 'GDPKIN' in rxn.id:
            rxn.name = 'GDP kinase'
        if '6PFRUCTPHOS_RXN' in rxn.id:
            rxn.name = '6-phosphofructokinase'
        if 'Mitochondrial_ATP_Sy' in rxn.id:
            rxn.name = 'Mitochondrial ATP synthase'
        if 'Plastidial_ATP_Sy' in rxn.id:
            rxn.name = 'Plastidial ATP synthase'
        if 'Protein_turnover' in rxn.id:
            rxn.name = 'Protein turnover'
        if 'RNA_turnover' in rxn.id:
            rxn.name = 'RNA turnover'
        if 'ATPase_tx' in rxn.id:
            rxn.name='Cell maintenance'
        if 'PYRROLINECARBDEHYDROG_RXN_NADP' in rxn.id:
            rxn.name = 'Pyrroline5-carboxylate reductase'
        if 'RXN_8001' in rxn.id:
            rxn.name = 'Histidinol dehydrogenase'
        if 'ATPPHOSPHORIBOSYLTRANS' in rxn.id:
            rxn.name = 'ATP phosphoribosyltransferase'
        if 'PYRUVDEH_RX' in rxn.id:
            rxn.name = 'Pyruvate dehydrogenase'
        if '2OXOGLUTARATEDEH_RXN' in rxn.id:
            rxn.name = '2-oxoglutarate dehydrogenase'
        if 'HYDROXYPYRUVATE_REDUCTASE_RXN_NAD' in rxn.id:
            rxn.name = 'Hydroxypyruvate reductase (NAD+)'
        if 'MALATE_DEHYDROGENASE_NADP_RXN' in rxn.id:
            rxn.name = 'Malate dehydrogenase (NADP+)'
        if 'RXN66_3' in rxn.id:
            rxn.name = 'Acetaldehyde dehydrogenase'
        if 'SUCROSE_ec' in rxn.id:
            rxn.name = 'Sucrose transport'
        if 'X5P_Pi_pc' in rxn.id:
            rxn.name = 'Plastidial Xylulose-5-phosphate-Pi shuttle'
        if 'Pi_PROTON_mc' in rxn.id:
            rxn.name = 'Mitochondrial Pi-H$^+$ shuttle'
        if 'H_pc' in rxn.id:
            rxn.name = 'Plastidial H$^+$ channel'
        if 'GAP_Pi_pc' in rxn.id:
            rxn.name = 'Plastidial GAP-Pi shuttle'
    return
#####################################################################
# This function converts flux to grams based on the formula of prod-#
#-ucts or reactants(if there are no products or if the user chooses #
# to consider only the reactants)				    #
# inputs: 1) a cobra model, 2) a reaction ID, 3) flux value (option-#
#-al) and 4) a boolean trigger (optional): False to use products for#
# calculation and True to use the reactants instead 		    #
#####################################################################
def convertFluxMol2Grams(model,rxnID,flux=0,noProd=True):
    rxn = model.reactions.get_by_id(rxnID)
    mass = 0
    if flux == 0:
        flux = rxn.flux
    if len(rxn.products)!= 0 and not noProd:
        for met in rxn.products:
            print(met.id)
            print(met.formula)
            if met.formula == "" or met.formula == "NA":
                noProd = True
            if not noProd:
                mass = mass + (met.formula_weight*rxn.metabolites.get(met))
    if noProd or len(rxn.products)==0:
        mass = 0
        for met in rxn.reactants:
            #print(met.id)
            #print(met.formula)
            if met.formula == "" or met.formula == "NA":
                print("PROBLEM: no formula for "+met.id)
                return
            else:
                mass = mass + abs(met.formula_weight*rxn.metabolites.get(met))
    mass = mass*flux
    return mass

def getInchiHcount(met):
    '''
    This function counts the number of protons when metabolite is uncharged
    Input: 1) metabolite with INCHI in sbml notes
    Output: Number of protons (int)
    Author: Sanu Shameer (sanushameer@gmail.com)
    '''
    temp = 0
    temp = met.notes.get("INCHI").split("/")[1].split("H")[1].split("O")[0].split("N")[0].split("P")[0].split("S")[0]
    return int(temp)

def getFormulaFromInchiAndCharge(met):
    '''
    This function generates metabolite formula based on INCHI and charge
    Input: 1) cobra metabolite object
    Output: metabolite formula
    Author: Sanu Shameer (sanushameer@gmail.com)
    '''
    temp = ""
    for elem in ["C","H","N","O","P","S"]:
        if elem == "H" and elem in met.formula:
            number = ""
            if met.elements.get(elem)>1:
                number = getInchiHcount(met)+met.charge
                if int(number) == number:
                    number = str(int(number))
                else:
                    number = str(number)
            temp = temp+elem+number
        elif elem in met.formula:
            number = ""
            if met.elements.get(elem)>1:
                number = str(met.elements.get(elem))
            temp = temp+elem+number
    return temp

def updateFormulaFromInchiAndCharge(met):
    '''
    This function generates updated metabolite formula based on INCHI and charge
    Input: 1) cobra metabolite object
    Output:
    Author: Sanu Shameer (sanushameer@gmail.com)
    '''
    temp = getFormulaFromInchiAndCharge(met)
    met.formula = temp
    return

def findAlternateMetStates(cobra_model):
    '''
    This function checks for "a" and "b" prefixes to find alternate charge state
    -s of metabolites and returns a python dictionary with the alternate charge
    state as the key and their default charge state as the value
    Input: 1) a cobra model
    Output: 1) a python dictionary with a cobra metabolite as the key and value
    '''

    fractionMets=dict()
    for met in cobra_model.metabolites:
        prefix=""
        a=re.search("^a{1,3}",met.id)
        anion=""
        if a:
            anion=a.group(0)
            prefix=anion
        b=re.search("^b{1,3}",met.id)
        basic=""
        if b:
            basic=b.group(0)
            prefix=basic
        if (not prefix == ""):
            try:
                mainMet = cobra_model.metabolites.get_by_id(met.id.replace(prefix,"",1))
                fractionMets[met]=mainMet
            except:
                continue
    return fractionMets

###function to take the fluxes of reactions at different solutions and to find correlations between them
#takes as an input the model to find the reactions names and a dictionary with FBA/pFBA solutions
#from Sanu
#returns:
    #rxnvalues - a dictionary with reactions as keys and the its fluxes as values
    #clustered - a set with the the reactions that cluster together or are single as frozensets
    #grouping - a dictionary with 1 to number of groupings (cluster of more than 1 reaction) as keys
    # and the reactions that cluster together as values

def pearsons_correlation(cobra_model, solutiondict):
    keyList = list(solutiondict.keys())
    rxnvalues=dict()

    for rxn in cobra_model.reactions:
        values=list()
        for i in keyList:
            values.append(solutiondict[i][rxn.id])
        rxnvalues[rxn.id]=values

    PCCdict=dict()
    i=0;
    for rxn1 in cobra_model.reactions:
        print("Processing "+str((i+1))+" of "+str(len(cobra_model.reactions)))
        i+=1
        PCC=dict()
        for rxn2 in cobra_model.reactions:
            PCC[rxn2.id] = stats.pearsonr(rxnvalues[rxn1.id],rxnvalues[rxn2.id])
        PCCdict[rxn1.id]=PCC

    usedSet=set()
    clustered=set()

    for rxn1 in PCCdict.keys():
        if rxn1 not in usedSet:
            usedSet.add(rxn1)
            tempset=set()
            tempset.add(rxn1)
            for rxn2 in PCCdict.get(rxn1):
                if(PCCdict.get(rxn1).get(rxn2)[0] > 0.9999 and PCCdict.get(rxn1).get(rxn2)[1] < 0.0001):
                    usedSet.add(rxn2)
                    tempset.add(rxn2)
            clustered.add(frozenset(tempset))

    grouping = dict()
    i=0
    for s in clustered:
        if(len(s)>1):
            grouping[i]=s
            i+=1

    return rxnvalues, clustered, grouping

def mapSource(cobra_model,sol,current_met='SUCROSE_c_SE_d',i_max=2,weighting=50,trav_mets=[],ignore=1,arrow_labels=1,cell_trans=None,transcriptdata=None):
    now = datetime.now().strftime('%Y-%m-%d-%H-%M')
    f=Digraph('Sources:'+current_met+now,filename='Sources/'+current_met+now+'.gv',strict=True)
    f.attr(rankdir='LR',size='8,5')
    f.attr('node', shape='diamond',style='filled',fillcolor='darkgoldenrod1')
    f.node(current_met)
    f.attr('node', shape='circle',style='filled',fillcolor='white')
    if ignore == 'PROTON':
        ignore_mets=['Protein_processing_cost','Protein_polymerisation_cost','Protein_translocation_cost',
                     'ATPase_NADPHoxidase_constraint_','aPi_','aADP_','aATP_']+trav_mets
    elif ignore == 'ATP':
        ignore_mets=['Protein_processing_cost','Protein_polymerisation_cost','Protein_translocation_cost',
                     'ATPase_NADPHoxidase_constraint_','PROTON_','aPi_','aADP_','aATP_']+trav_mets
    elif ignore:
        ignore_mets=['ATP_','aATP_','ADP_','aADP_','WATER_','PROTON_','NADPH_','NADH_','NADP_','NAD_',
                     'ATPase_NADPHoxidase_constraint_','Pi_','aPi_','Protein_processing_cost',
                     'Protein_polymerisation_cost','Protein_translocation_cost','CARBON_CC','CARBON_MC']+trav_mets
    else:
        ignore_mets=['Protein_processing_cost','Protein_polymerisation_cost','Protein_translocation_cost',
                     'ATPase_NADPHoxidase_constraint_','PROTON_','Pi_','aPi_','aATP_','aADP_']+trav_mets
    ignore_rxns=[]
    trav_mets=mapSourceRec(f,current_met,cobra_model,sol,0,ignore_mets,ignore_rxns,i_max,weighting,arrow_labels,cell_trans,transcriptdata)
    print(trav_mets)
    f.view()
    return f

def mapSourceRec(f,current_met,cobra_model,sol,i=0,trav_mets=[],trav_rxns=[],i_max=2,weight_mult=50,arrow_labels=1,cell_trans=None,transcriptdata=None):
    if i>i_max-1:
        return ''
    i=i+1
    ig_rxns=[]
    colours=['azure3','royalblue','deeppink4','red','forestgreen','darkorchid','lightseagreen','orange']
    # weight_mult=100
    for rxn in [x for x in cobra_model.metabolites.get_by_id(current_met).reactions if sol[x.id]*x.get_coefficient(current_met)>0]:
        if abs(sol[rxn.id])>1e-5:
            if transcriptdata and str(transcriptdata[rxn.id])!='nan':
                datathere='*'
                rxnsuff='_'.join(rxn.id.split('_')[-2:])
                mcver=rxn.id.replace(rxnsuff,'MC_l')
                if 'MC_l' not in rxn.id and mcver in transcriptdata.keys():
                    if 20*transcriptdata[rxn.id]<transcriptdata[mcver]:
                        datathere+='*'
                    elif 20*transcriptdata[rxn.id]>transcriptdata[mcver]:
                        datathere+='^'
                datathere+=str(transcriptdata[rxn.id])
            else:
                datathere=''
            if rxn.get_coefficient(current_met)<0:
                for met in rxn.products:
                    if not any(x == met.id[:min(len(x),len(met.id))] for x in trav_mets):
                        if not any(x in rxn.id for x in trav_rxns):
                            this_width=min(weight_mult*abs(sol[rxn.id]*rxn.get_coefficient(current_met)),60)
                            this_edge=str(max(this_width,.5))
                            this_arrow=str(min(max(this_width/20,0.5),10))
                            if 'dielTransfer' in rxn.id:
                                this_col=1
                            elif any(x in rxn.id for x in ['_CC_SE_','_SE_pSE_','_ec_']):
                                this_col=2
                                if cell_trans:
                                    continue
                            elif any(['ATP' in met.id for met in rxn.products]):
                                this_col=3
                            elif any(['ATP' in met.id for met in rxn.reactants]):
                                this_col=4
                            elif any(['NADH' in met.id for met in rxn.products]) or any(['NADPH' in met.id for met in rxn.products]):
                                this_col=5
                            elif any(['NADH' in met.id for met in rxn.reactants]) or any(['NADPH' in met.id for met in rxn.reactants]):
                                this_col=6
                            elif any(['PPI' in met.id for met in rxn.products]):
                                this_col=7
                            else:
                                this_col=0
                            # f.attr(size=str((sol[rxn.id]*rxn.get_coefficient(current_met))/50))
                            if arrow_labels:
                                f.edge(met.id,current_met, label=rxn.id+datathere+': '+str(rxn.get_coefficient(current_met))+'*'+str(sol[rxn.id])+'\n'+rxn.reaction,arrowsize=this_arrow,penwidth=this_edge,color=colours[this_col])
                            else:
                                f.edge(met.id,current_met,label=str(sol[rxn.id]),arrowsize=this_arrow,penwidth=this_edge,color=colours[this_col])
                            # f.edge_attr.update(arrowhead='vee', arrowsize=this_width,penwidth=this_width)
#                                 print(this_edge)
                            ig_rxns=trav_rxns.copy()
                            ig_rxns.append(rxn.id)
                            temp=mapSourceRec(f,met.id,cobra_model,sol,i,trav_mets,ig_rxns,i_max,weight_mult,arrow_labels,cell_trans,transcriptdata)
            else:
                for met in rxn.reactants:
                    if not any(x == met.id[:min(len(x),len(met.id))] for x in trav_mets):
                        if not any(x in rxn.id for x in trav_rxns):
                            this_width=min(weight_mult*abs(sol[rxn.id]*rxn.get_coefficient(current_met)),60)
                            this_edge=str(max(this_width,.5))
                            this_arrow=str(min(max(this_width/20,0.5),10))
                            if 'dielTransfer' in rxn.id:
                                this_col=1
                            elif any(x in rxn.id for x in ['_CC_SE_','_SE_pSE_','_ec_']):
                                this_col=2
                                if cell_trans:
                                    continue
                            elif any(['ATP' in met.id for met in rxn.products]):
                                this_col=3
                            elif any(['ATP' in met.id for met in rxn.reactants]):
                                this_col=4
                            elif any(['NADH' in met.id for met in rxn.products]) or any(['NADPH' in met.id for met in rxn.products]):
                                this_col=5
                            elif any(['NADH' in met.id for met in rxn.reactants]) or any(['NADPH' in met.id for met in rxn.reactants]):
                                this_col=6
                            elif any(['PPI' in met.id for met in rxn.products]):
                                this_col=7
                            else:
                                this_col=0
                            if arrow_labels:
                                f.edge(met.id,current_met, label=rxn.id+datathere+': '+str(rxn.get_coefficient(current_met))+'*'+str(sol[rxn.id])+'\n'+rxn.reaction,arrowsize=this_arrow,penwidth=this_edge,color=colours[this_col])
                            else:
                                f.edge(met.id,current_met, label=str(sol[rxn.id]) ,arrowsize=this_arrow,penwidth=this_edge,color=colours[this_col])
                            # f.edge_attr.update(arrowhead='vee', arrowsize=this_width,penwidth=this_width)
#                                 print(this_edge)
                            ig_rxns=trav_rxns.copy()
                            ig_rxns.append(rxn.id)
                            temp=mapSourceRec(f,met.id,cobra_model,sol,i,trav_mets,ig_rxns,i_max,weight_mult,arrow_labels,cell_trans,transcriptdata)
    return ig_rxns

def mapSink(cobra_model,sol,current_met='SUCROSE_c_SE_d',i_max=2,weighting=50,trav_mets=[],cell_trans=None,ends=None,arrow_labels=1,data=None,ignore=None):
    now = datetime.now().strftime('%Y-%m-%d-%H-%M')
    f=Digraph('Sources:'+current_met+now,filename='Sources/'+current_met+now+'.gv',strict=True)
    f.attr(rankdir='LR',size='8,5')
    f.attr('node', shape='diamond',style='filled',fillcolor='darkgoldenrod1')
    f.node(current_met)
    f.attr('node', shape='circle',style='filled',fillcolor='white')
    if ignore=='PROTON':
        ignore_mets=['Protein_processing_cost','Protein_polymerisation_cost','Protein_translocation_cost',
                     'ATPase_NADPHoxidase_constraint_','aPi_','aADP_','aATP_']+trav_mets
    else:
        ignore_mets=['ATP_','aATP_','ADP_','aADP_','AMP_','WATER_','PROTON_','NADPH_','NADH_','NADP_','NAD_','ATPase_NADPHoxidase_constraint_','Pi_','aPi_','Protein_processing_cost','Protein_polymerisation_cost','Protein_translocation_cost','CARBON_CC','CARBON_MC','UBIQUIN','AMMONIUM_','PPI_','Phloem_contri','sucrose_to_be','CARBON_DIOXIDE_']+trav_mets
    ignore_rxns=[]
    trav_mets=mapSinkRec(f,current_met,cobra_model,sol,0,ignore_mets,ignore_rxns,i_max,weighting,cell_trans,ends,arrow_labels,data)
    print(trav_mets)
    f.view()
    return f

def mapSinkRec(f,current_met,cobra_model,sol,i=0,trav_mets=[],trav_rxns=[],i_max=2,weight_mult=50,cell_trans=None,ends=None,arrow_labels=1,data=None):
    if i>i_max-1:
        return ''
    i=i+1
    ig_rxns=[]
    colours=['azure3','royalblue','deeppink4','red','forestgreen','darkorchid','lightseagreen','orange']
    gltmets=['GLT_','2_KETOGLUTARATE_']
    for rxn in [x for x in cobra_model.metabolites.get_by_id(current_met).reactions if sol[x.id]*x.get_coefficient(current_met)<0]:
        if abs(sol[rxn.id])>1e-5:
#                 print(str(i)+': '+rxn.id)
            if rxn.get_coefficient(current_met)<0:
#                 print(rxn.id,rxn.reaction)
                mets=rxn.products
                if any([x in current_met for x in gltmets]) and any([x in met.id for x in gltmets for met in mets]):
                    mets=[met for met in mets if any([x in met.id for x in gltmets])]
                for met in mets:
                    cont = 1
                    if not (any([x in rxn.id for x in ['_pc_','_mc_','_xc_','_vc_']]) and met.id.split('_')[-3]==current_met.split('_')[-3]):
                        if any(x == met.id[:min(len(x),len(met.id))] for x in trav_mets):
                            if ends:
                                cont=0
                            else:
                                continue
                        if any(x in rxn.id for x in trav_rxns):
                            cont=0
#                             print(met.id,'end')
                        this_width=min(weight_mult*abs(sol[rxn.id]*rxn.get_coefficient(current_met)),60)
                        this_edge=str(max(this_width,.5))
                        this_arrow=str(min(max(this_width/20,0.5),10))
                        if 'dielTransfer' in rxn.id:
                            this_col=1
                        elif any(x in rxn.id for x in ['_CC_SE_','_SE_pSE_','_ec_']):
                            this_col=2
                            if cell_trans:
                                continue
                        elif any(['ATP' in met.id for met in rxn.products]):
                            this_col=3
                        elif any(['ATP' in met.id for met in rxn.reactants]):
                            this_col=4
                        elif any(['NADH' in met.id for met in rxn.products]) or any(['NADPH' in met.id for met in rxn.products]):
                            this_col=5
                        elif any(['NADH' in met.id for met in rxn.reactants]) or any(['NADPH' in met.id for met in rxn.reactants]):
                            this_col=6
                        elif any(['PPI' in met.id for met in rxn.products]):
                            this_col=7
                        else:
                            this_col=0
                        # f.attr(size=str((sol[rxn.id]*rxn.get_coefficient(current_met))/50))
                        tailtype=''
                        if data:
#                             print('1a',rxn.id)
                            if data[rxn.id]!='nan' and data[rxn.id]!='':
                                tailtype='dashed'
#                                 print('2a')
                        if arrow_labels:
                            f.edge(current_met,met.id,label=rxn.id+': '+str(rxn.get_coefficient(current_met))+'*'+str(sol[rxn.id])+'\n'+rxn.reaction,arrowsize=this_arrow,penwidth=this_edge,color=colours[this_col],style=tailtype)
                        else:
                            f.edge(current_met,met.id,label=str(sol[rxn.id]),arrowsize=this_arrow,penwidth=this_edge,color=colours[this_col],style=tailtype)
                        # f.edge_attr.update(arrowhead='vee', arrowsize=this_width,penwidth=this_width)
    #                                 print(this_edge)
                        ig_rxns=trav_rxns.copy()
                        ig_rxns.append(rxn.id)
                        if cont:
                            temp=mapSinkRec(f,met.id,cobra_model,sol,i,trav_mets,ig_rxns,i_max,weight_mult,cell_trans,ends,arrow_labels,data)
            else:
                mets=rxn.reactants
                if any([x in current_met for x in gltmets]) and any([x in met.id for x in gltmets for met in mets]):
                    mets=[met for met in mets if any([x in met.id for x in gltmets])]
                for met in mets:
                    cont = 1
                    if not (any([x in rxn.id for x in ['_pc_','_mc_','_xc_','_vc_']]) and met.id.split('_')[-3]==current_met.split('_')[-3]):
                        if any(x == met.id[:min(len(x),len(met.id))] for x in trav_mets):
                            if ends:
                                cont=0
                            else:
                                continue
                        if any(x in rxn.id for x in trav_rxns):
                            cont = 0
#                             print(met.id,'end')
                        this_width=min(weight_mult*abs(sol[rxn.id]*rxn.get_coefficient(current_met)),60)
                        this_edge=str(max(this_width,.5))
                        this_arrow=str(min(max(this_width/20,0.5),10))
                        if 'dielTransfer' in rxn.id:
                            this_col=1
                        elif any(x in rxn.id for x in ['_CC_SE_','_SE_pSE_','_ec_']):
                            if cell_trans:
                                continue
                            this_col=2
                        elif any(['ATP' in met.id for met in rxn.products]):
                            this_col=4
                        elif any(['ATP' in met.id for met in rxn.reactants]):
                            this_col=3
                        elif any(['NADH' in met.id for met in rxn.products]) or any(['NADPH' in met.id for met in rxn.products]):
                            this_col=6
                        elif any(['NADH' in met.id for met in rxn.reactants]) or any(['NADPH' in met.id for met in rxn.reactants]):
                            this_col=5
                        elif any(['PPI' in met.id for met in rxn.reactants]):
                            this_col=7
                        else:
                            this_col=0
                        tailtype=''
                        if data:
#                             print('1b')
                            if data[rxn.id]!='nan' and data[rxn.id]!='':
#                                 print('2b')
                                tailtype='dashed'
                        if arrow_labels:
                            f.edge(current_met,met.id,label=rxn.id+': '+str(rxn.get_coefficient(current_met))+'*'+str(sol[rxn.id])+'\n'+rxn.reaction,arrowsize=this_arrow,penwidth=this_edge,color=colours[this_col],style=tailtype)
                        else:
                            f.edge(current_met,met.id,label=str(sol[rxn.id]),arrowsize=this_arrow,penwidth=this_edge,color=colours[this_col],style=tailtype)
                        # f.edge_attr.update(arrowhead='vee', arrowsize=this_width,penwidth=this_width)
    #                                 print(this_edge)
                        ig_rxns=trav_rxns.copy()
                        ig_rxns.append(rxn.id)
                        if cont==1:
                            temp=mapSinkRec(f,met.id,cobra_model,sol,i,trav_mets,ig_rxns,i_max,weight_mult,cell_trans,ends,arrow_labels,data)
    return ig_rxns



def generateXbudget(model,solution,budg_met='PROTON',outfile="",show_plot=True,percentage=False,day_or_night_tag="_MC_l",compartments=["c","p","m","x"],save_plot_to="temp.png",colourDict={},prevlabels=None):
    if outfile!="":
        fout = open(outfile,"w")
    ATPdict = dict()
    total = 0
    for p in compartments:
        met=model.metabolites.get_by_id(budg_met+"_"+p+day_or_night_tag)
        for rxn in met.reactions:
            sto=rxn.metabolites.get(met)
            if outfile!="":
                fout.write(rxn.id+"\t"+rxn.reaction+"\t"+str(solution[rxn.id]*(sto))+"\t"+met.compartment+"\n")
            ATPdict[rxn.id]=solution[rxn.id]*(sto)
            if solution[rxn.id]*(sto) > 0:
                total = total + (solution[rxn.id]*(sto))
    if outfile!="":
        fout.close()

    tempDict = dict()
    for rxn in ATPdict.keys():
        tempDict[rxn]=abs(ATPdict[rxn])

    #sort ATPdict by values
    sorted_by_value = sorted(tempDict.items(), key= lambda x:x[1],reverse=True)

    ATPdict2 = dict()
    ATPdict2["Others-pos"]=0
    ATPdict2["Others-neg"]=0
    baseline = dict()
    pos_base=0
    neg_base=0
    i=0
#     print(sorted_by_value)
#     print(([y for y in sorted_by_value if y[1]<sorted_by_value[2][1]]))
#     print([sorted_by_value[x][1] for x in range(18,1,-1)])
    other_thresh=[sorted_by_value[x][1] for x in range(17,1,-1) if len([y for y in sorted_by_value if y[1]>sorted_by_value[x][1]])<=17]
#     print(other_thresh)
    other_thresh=max(other_thresh+[5e-5])
#     print(other_thresh)
#     other_thresh=total*0.0001
    for TEMP in sorted_by_value:
        rxn = TEMP[0]
        if ATPdict[rxn]>0:
            if ATPdict[rxn] < other_thresh:
                if percentage:
                    ATPdict2["Others-pos"]=ATPdict2["Others-pos"]+float(ATPdict[rxn]*100)/total
                else:
                    ATPdict2["Others-pos"]=ATPdict2["Others-pos"]+ATPdict[rxn]
                continue
            base = pos_base
            if percentage:
                ATPdict2[rxn]=float(ATPdict[rxn]*100)/total
                pos_base = pos_base + float(ATPdict[rxn]*100)/total
            else:
                pos_base = pos_base + ATPdict[rxn]
                ATPdict2[rxn]=ATPdict[rxn]
        else:
            if abs(ATPdict[rxn]) < other_thresh:
                if percentage:
                    ATPdict2["Others-neg"]=ATPdict2["Others-neg"]+float(ATPdict[rxn]*100)/total
                else:
                    ATPdict2["Others-neg"]=ATPdict2["Others-neg"]+ATPdict[rxn]
                continue
            base = neg_base
            if percentage:
                ATPdict2[rxn]=float(ATPdict[rxn]*100)/total
                neg_base = neg_base + float(ATPdict[rxn]*100)/total
            else:
                neg_base = neg_base + ATPdict[rxn]
                ATPdict2[rxn]=ATPdict[rxn]
        i=i+1
        baseline[rxn]=base
    baseline["Others-pos"]=pos_base
    baseline["Others-neg"]=neg_base

    if show_plot:
        plt.rcParams.update({'font.size': 10}) #sets a global fontsize
        plt.rcParams['xtick.major.size'] = 5 # adjusts tick line length and width
        plt.rcParams['xtick.major.width'] = 1
        plt.rcParams['ytick.major.size'] = 5
        plt.rcParams['ytick.major.width'] = 1
        plt.rcParams['axes.linewidth']=2 # makes axes line thicker
        thisplot = plt.figure(figsize=(3,4))
        if prevlabels:
            atpdictkeys=list(ATPdict2.keys())
            for it,lb in enumerate(prevlabels):
                if lb in atpdictkeys and lb!= atpdictkeys[it]:
                    temp=atpdictkeys[it]
                    lbpl = atpdictkeys.index(lb)
                    atpdictkeys[it]=lb
                    atpdictkeys[lbpl]=temp
        NUM_COLORS=len(ATPdict2.keys())
        if NUM_COLORS<=20:
            cm = plt.get_cmap('tab20')
        else:
            cm = plt.get_cmap('gist_ncar')
        for i,rxn in enumerate(ATPdict2.keys()):
            if colourDict.keys().__contains__(rxn):
                plt.bar(1,ATPdict2[rxn],width=0.1,bottom=baseline[rxn],label=rxn,color=colourDict[rxn])
            else:
                plt.bar(1,ATPdict2[rxn],width=0.1,bottom=baseline[rxn],label=rxn,color=cm(1.*i/NUM_COLORS))
        plt.xlim(0.8,1.2)
        if percentage:
            plt.ylabel(budg_met+" produced/consumed (%)")
        else:
            plt.ylabel(budg_met+" produced/consumed (in moles)")
        handles, labels = plt.gca().get_legend_handles_labels()
        labels2=list(set(labels)-set(["Others-neg","Others-pos"]))+list(["Others-neg","Others-pos"])
        handles2=[handles[labels.index(i)] for i in labels2]
        lgd=plt.legend(handles2,labels2,bbox_to_anchor=(1,1))
        plt.axhline(0,linestyle="--",color="black")
        plt.tight_layout
        plt.savefig(save_plot_to, bbox_extra_artists=(lgd,), bbox_inches='tight')
        return thisplot,labels
     
def plot_carbon_uptake(stem_model,psol,thresh=1,excl_cells='',scaling=0,overall=0):
    carbon_dict={}
    for rxn in stem_model.reactions:
        if any([('CARBON_DIOXIDE' in x.id or 'HCO3' in x.id) for x in rxn.metabolites])and abs(psol[rxn.id])>0 and not any(x in rxn.id for x in excl_cells):
            carbon_dict[rxn.id]=0
            for met in rxn.metabolites:
                if ('CARBON_DIOXIDE' in met.id or 'HCO3' in met.id):
                    carbon_dict[rxn.id]+=psol[rxn.id]*rxn.get_coefficient(met.id)
    if scaling:
        for rxn in carbon_dict.keys():
            if '_CC_' in rxn:
                carbon_dict[rxn]=carbon_dict[rxn]*20
            if '_SE_' in rxn:
                carbon_dict[rxn]=carbon_dict[rxn]*4
            if '_pSE_' in rxn:
                carbon_dict[rxn]=carbon_dict[rxn]*10/3
    if overall:
        new_carbon_dict={}
        new_carbon_dict['MC']=0
        new_carbon_dict['CC']=0
        new_carbon_dict['SE']=0
        new_carbon_dict['pSE']=0
        for rxn in carbon_dict.keys():
            if '_MC_' in rxn:
                new_carbon_dict['MC']+=carbon_dict[rxn]
            elif '_CC_' in rxn:
                new_carbon_dict['CC']+=carbon_dict[rxn]
            elif '_SE_' in rxn:
                new_carbon_dict['SE']+=carbon_dict[rxn]
            elif '_pSE_' in rxn:
                new_carbon_dict['pSE']+=carbon_dict[rxn]
        carbon_dict=new_carbon_dict
    labels=[x for x in carbon_dict.keys()if abs(carbon_dict[x])>thresh]
    data1=[carbon_dict[x] if '_tx' in x else -carbon_dict[x] for x in labels]
    width=0.6
    plt.figure(figsize=(len(labels)*6.4/20,4.8))
    plt.xticks(range(len(data1)), labels,rotation='vertical')
    plt.bar(np.arange(len(data1)), data1, width=width)
    # plt.bar(np.arange(len(data2))+ width, data2, width=width)
    plt.legend(['carbon -> model'])
    # plt.savefig('RegrEx_plots/'+'1_transfer_large_day.png', bbox_inches='tight')
    print('Carbon in: ',sum([x for x in carbon_dict.values()if x>0 ]))
    print('Carbon out: ',sum([x for x in carbon_dict.values()if x>0 ]))
    plt.show()

def plot_sym_transport(stem_model,psol,light_phase='_l',threshmin=5e-5,threshmax=1001,ybreak=None,fva_sol=None):
    ccse_rxns=[x for x in stem_model.reactions if '_CC_SE'+light_phase in x.id]
    sepse_rxns=[x for x in stem_model.reactions if '_SE_pSE'+light_phase in x.id]
    if fva_sol:
        range1=[x for x in ccse_rxns+sepse_rxns if abs(fva_sol['minimum'][x.id])+abs(fva_sol['maximum'][x.id])>threshmin]
    else:
        range1=[x for x in ccse_rxns+sepse_rxns if abs(psol[x.id])>threshmin and abs(psol[x.id])<threshmax]
    range_all=[]
    for ii in range1:
        if ii.id.replace('_SE_pSE'+light_phase,'_CC_SE'+light_phase) not in range_all:
            range_all.append(ii.id)
    data1=[psol[x] for x in range_all]
    data2=[psol[x] if x in stem_model.reactions else 0 for x in [y.replace('_CC_SE'+light_phase,'_SE_pSE'+light_phase) for y in range_all]]
    if fva_sol:
        data1_min=[fva_sol['minimum'][x] for x in range_all]
        data1_max=[fva_sol['maximum'][x] for x in range_all]
        data2_min=[fva_sol['minimum'][x] if x in list(psol.keys()) else 0 for x in [y.replace('_CC_SE'+light_phase,'_SE_pSE'+light_phase) for y in range_all]]
        data2_max=[fva_sol['maximum'][x] if x in list(psol.keys()) else 0 for x in [y.replace('_CC_SE'+light_phase,'_SE_pSE'+light_phase) for y in range_all]]
    labels=range_all.copy()
    for x in range(len(labels)):
        temp=labels[x].split('_')
        labels[x]='_'.join(temp[:-4]+[temp[-1]])
    width=0.3
    if ybreak:
        fig,(ax,ax2,ax3) = plt.subplots(3,1, sharex=True,figsize=(len(labels)*6.4/20,4.8))
    else:
        plt.figure(figsize=(len(labels)*6.4/20,4.8))
    plt.xticks(range(len(data1)), labels,rotation='vertical')
    if ybreak:
        breakthresh=ybreak #0.005
        ax.bar(np.arange(len(data1)), data1, width=width)
        ax.bar(np.arange(len(data2))+ width, data2, width=width)
        ax2.bar(np.arange(len(data1)), data1, width=width)
        ax2.bar(np.arange(len(data2))+ width, data2, width=width)
        ax3.bar(np.arange(len(data1)), data1, width=width)
        ax3.bar(np.arange(len(data2))+ width, data2, width=width)
        ax2.set_ylim(-breakthresh,breakthresh)
        if max(data1)>breakthresh:
            min2 = max(0.95*min([x for x in data1 if x>breakthresh]),breakthresh)
        else:
            min2=breakthresh
        min3 = 1.05*min(data1)
        max2 = 1.05*max(data1)
        if min(data1)<-breakthresh:
            max3 = max(0.95*max([x for x in data1 if x<-breakthresh]),-breakthresh)
        else:
            max3=-breakthresh
        ax.set_ylim(min2,max2)
        ax3.set_ylim(min3,max3)
        ax2.spines['top'].set_visible(False)
        ax2.spines['bottom'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
        ax3.spines['top'].set_visible(False)
        ax2.tick_params(bottom=False)
        ax.xaxis.tick_top()
        ax.tick_params(labeltop=False)
        ax3.xaxis.tick_bottom()
        ax.yaxis.set_major_formatter(FormatStrFormatter('%.3f'))
        ax2.yaxis.set_major_formatter(FormatStrFormatter('%.3f'))
        ax3.yaxis.set_major_formatter(FormatStrFormatter('%.3f'))

        d = .5  # proportion of vertical to horizontal extent of the slanted line
        kwargs = dict(marker=[(-1, -d), (1, d)], markersize=12,
                      linestyle="none", color='k', mec='k', mew=1, clip_on=False)
        ax.plot([0, 1], [0, 0], transform=ax.transAxes, **kwargs, label='_nolegend_')
        ax2.plot([0, 1], [1, 1], transform=ax2.transAxes, **kwargs, label='_nolegend_')
        ax2.plot([0, 1], [0, 0], transform=ax2.transAxes, **kwargs, label='_nolegend_')
        ax3.plot([0, 1], [1, 1], transform=ax2.transAxes, **kwargs, label='_nolegend_')
    else:
        plt.bar(np.arange(len(data1)), data1, width=width)
        plt.bar(np.arange(len(data2))+ width, data2, width=width)
    plt.legend(['CC -> SE','SE -> pSE'])
    # plt.savefig('RegrEx_plots/'+'1_transfer_large_day.png', bbox_inches='tight')
    plt.show()

def plot_apo_transport(stem_model,psol,cell_type='CC_l',thresh=5e-5,yscale=None,fva_sol_min=None,fva_sol_max=None,rangeplot=None,exclude_rxns=[],ybreak=None,savebool=None,ymax=0,plttitle=None,exp=None):
    if not rangeplot:
        labels=[rxn.id for rxn in stem_model.reactions if (('_ec_'+cell_type in rxn.id)and abs(psol[rxn.id])>thresh)and not any([x in rxn.id for x in exclude_rxns])]
    else:
        labels=[rxn for rxn in list(fva_sol_min.keys()) if (('_ec_'+cell_type in rxn)and (abs(fva_sol_min[rxn])>thresh or abs(fva_sol_max[rxn])>thresh))and not any([x in rxn for x in exclude_rxns])]
    labels = sorted(labels, key=lambda x: x,reverse=False)
    data1=[psol[rxn] for rxn in labels]
    if rangeplot:
        datamax = [round(fva_sol_max[rxn]-psol[rxn],4) for rxn in labels]
        datamin = [round(psol[rxn]-fva_sol_min[rxn],4) for rxn in labels]
#         print(labels)
#         print(data1)
#         print(datamax)
#         print(datamin)
    labels = ['_'.join(x.split('_')[:-3]) for x in labels]
    width=0.6
    
    if ybreak:
        fig=plt.figure(dpi=1000)
        axb = fig.add_subplot(111)
        axb.spines['top'].set_color('none')
        axb.spines['bottom'].set_color('none')
        axb.spines['left'].set_color('none')
        axb.spines['right'].set_color('none')
        axb.tick_params(labelcolor='w', top=False, bottom=False, left=False, right=False)
        if plttitle:
            plt.title(plttitle)
        # fig,(ax,ax2) = plt.subplots(2,1, sharex=True)
        (ax,ax2) = fig.subplots(2,1, sharex=True)
    if exp:
        for ii in range(len(labels)):
            labels[ii]=labels[ii].replace('GLC','Glu ')
            labels[ii]=labels[ii].replace('FRU','Fru ')
            labels[ii]=labels[ii].replace('SUCROSE','Suc ')
            
    plt.xticks(range(len(data1)), labels,rotation='vertical')
    if ybreak:
        ax.bar(np.arange(len(data1)), data1, width=width)
        ax2.bar(np.arange(len(data1)), data1, width=width)
        ax2.set_ylim(0,ybreak)
        min2 = min([x for x in data1 if x>ybreak])
        max2 = max(1.05*max(data1),ymax)
        ax.set_ylim(min2,max2)
        ax2.spines['top'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
        ax.xaxis.tick_top()
        ax.tick_params(labeltop=False)
        ax2.xaxis.tick_bottom()
        
        d = .5  # proportion of vertical to horizontal extent of the slanted line
        kwargs = dict(marker=[(-1, -d), (1, d)], markersize=12,
                      linestyle="none", color='k', mec='k', mew=1, clip_on=False)
        ax.plot([0, 1], [0, 0], transform=ax.transAxes, **kwargs)
        ax2.plot([0, 1], [1, 1], transform=ax2.transAxes, **kwargs)
        ax.yaxis.set_major_formatter(FormatStrFormatter('%.3f'))
        ax2.yaxis.set_major_formatter(FormatStrFormatter('%.3f'))
        if rangeplot:
            ax.errorbar(np.arange(len(data1)),
                data1,
                yerr=[datamin,datamax], 
                fmt="o", markersize=1,color="grey",label='_nolegend_',capsize=5)
            ax2.errorbar(np.arange(len(data1)),
                data1,
                yerr=[datamin,datamax], 
                fmt="o", markersize=1,color="grey",label='_nolegend_',capsize=5)
            # print(data1)
            # print(datamin)
            # print(datamax)
    else:
        plt.bar(np.arange(len(data1)), data1, width=width)
        plt.ylim(0,max(ymax,max(data1)))
        if rangeplot:
            plt.errorbar(np.arange(len(data1)),
                data1,
                yerr=[datamin,datamax], 
                fmt="o", markersize=1, color="grey",label='_nolegend_',capsize=5)
    # plt.bar(np.arange(len(data2))+ width, data2, width=width)
#     plt.legend(['apo ->'+cell_type])
    if ybreak:
        axb.set_ylabel(r'flux ($\mu$mol/m$^2$/s)',labelpad=3)
        axb.yaxis.set_label_coords(-.2, 0.5)
    else:
        plt.ylabel(r'flux ($\mu$mol/m$^2$/s)')
    plt.rc('font', size=18)
    if yscale:
        if not rangeplot and not any([x<-5e-5 for x in data1]):
            plt.ylim(0,min(100,yscale))
        else:
            plt.ylim(max(-80,-yscale),min(100,yscale))
    if savebool:
        now = datetime.now().strftime('%Y_%m_%d')
        plt.savefig('Figures/'+'apo_transport_'+savebool+now+'.jpeg', bbox_inches='tight')
    plt.show()

def plot_diel(cons_model, csol,thresh=5e-5,ymax=None,pc=0,fva_sol=None,rangeplot=0,ybreak=None,savebool=None):
    cell_list=['MC','CC','SE','pSE']
    aa_syn_dict={'MC':[],'CC':[],'SE':[],'pSE':[]}
    if rangeplot:
        aa_syn_dictmin={'MC':[],'CC':[],'SE':[],'pSE':[]}
        aa_syn_dictmax={'MC':[],'CC':[],'SE':[],'pSE':[]}
    if rangeplot:
        tot_synmin=[]
        tot_synmax=[]
    syn_rxns=[rxn.id for rxn in cons_model.reactions if ('dielTransfer' in rxn.id and abs(csol[rxn.id])>thresh)]
    cell_syn={}
    if rangeplot:
        cell_synmin={}
        cell_synmax={}
    for y in cell_list:
        syn_rxns=[x.replace('_'+y,'_'+'MC') for x in syn_rxns]
    syn_rxns=list(set(syn_rxns))
    syn_rxns = sorted(syn_rxns, key=lambda x: x,reverse=True)
    for ii in range(4):
        for it,jj in enumerate(syn_rxns):
            rxnid=jj.replace('_MC','_'+cell_list[ii])
            if rxnid in list(csol.keys()):
                aa_syn_dict[cell_list[ii]]+=[csol[rxnid]]
            else:
                aa_syn_dict[cell_list[ii]]+=[0]
            if rangeplot:
                if rxnid in list(csol.keys()):
                    aa_syn_dictmin[cell_list[ii]]+=[fva_sol['minimum'][rxnid]]
                    aa_syn_dictmax[cell_list[ii]]+=[fva_sol['maximum'][rxnid]]
                else:
                    aa_syn_dictmin[cell_list[ii]]+=[0]
                    aa_syn_dictmax[cell_list[ii]]+=[0]
    
    if pc=='rat':
        ratios=[1,20,5,10/3]
        mc_syn = aa_syn_dict[cell_list[0]]
        for ii in range(4):
            aa_syn_dict[cell_list[ii]]=[aa_syn_dict[cell_list[ii]][x]*ratios[ii] for x in range(len(aa_syn_dict[cell_list[ii]]))]
            if rangeplot:
                aa_syn_dictmin[cell_list[ii]]=[aa_syn_dictmin[cell_list[ii]][x]*ratios[ii] for x in range(len(aa_syn_dictmin[cell_list[ii]]))]
                aa_syn_dictmin[cell_list[ii]] = [abs(aa_syn_dictmin[cell_list[ii]][x]-aa_syn_dict[cell_list[ii]][x]) for x in range(len(aa_syn_dictmin[cell_list[ii]]))]
                aa_syn_dictmax[cell_list[ii]]=[aa_syn_dictmax[cell_list[ii]][x]*ratios[ii] for x in range(len(aa_syn_dictmax[cell_list[ii]]))]
                aa_syn_dictmax[cell_list[ii]] = [aa_syn_dictmax[cell_list[ii]][x]-aa_syn_dict[cell_list[ii]][x] for x in range(len(aa_syn_dictmax[cell_list[ii]]))]
                
    labels=['_'.join(x.split('_')[:-3]) for x in syn_rxns]

    width=0.2
    if ybreak:
        ax3b=0
        if any([x<0 for x in aa_syn_dict[cell_list[0]]+aa_syn_dict[cell_list[1]]]):
            fig,(ax,ax2,ax3) = plt.subplots(3,1, sharex=True,figsize=(len(labels)*6.4/12,4.8))
            ax3b=1
        else:
            fig,(ax,ax2) = plt.subplots(2,1, sharex=True,figsize=(len(labels)*6.4/12,4.8))
    else:
        plt.figure(figsize=(len(labels)*6.4/6,4.8))
    plt.xticks(range(len(aa_syn_dict[cell_list[ii]])), labels,rotation='vertical')
    for ii in range(len(cell_list)):
        if ybreak:
            ax.bar(np.arange(len(aa_syn_dict[cell_list[ii]]))+(ii-1.5)*width,
                    aa_syn_dict[cell_list[ii]], width=width)
            ax2.bar(np.arange(len(aa_syn_dict[cell_list[ii]]))+(ii-1.5)*width,
                    aa_syn_dict[cell_list[ii]], width=width)
            if ax3b:
                ax3.bar(np.arange(len(aa_syn_dict[cell_list[ii]]))+(ii-1.5)*width,
                    aa_syn_dict[cell_list[ii]], width=width)
            if ii==0:
                breakthresh=ybreak
                if ax3b:
                    ax2.set_ylim(-breakthresh,breakthresh)
                    min2 = max((1-ymax)*min([x for x in aa_syn_dict[cell_list[ii]] if x>breakthresh]),breakthresh)
                    min3 = (1+ymax)*min(aa_syn_dict[cell_list[ii]])
                    max2 = (1+ymax)*max(aa_syn_dict[cell_list[ii]])
                    max3 = max((1-ymax)*max([x for x in aa_syn_dict[cell_list[ii]] if x<-breakthresh]),-breakthresh)
                    ax.set_ylim(min2,max2)
                    ax3.set_ylim(min3,max3)
                    ax2.spines['top'].set_visible(False)
                    ax2.spines['bottom'].set_visible(False)
                    ax.spines['bottom'].set_visible(False)
                    ax3.spines['top'].set_visible(False)
                    ax2.tick_params(bottom=False)
                    ax.xaxis.tick_top()
                    ax.tick_params(labeltop=False)
                    ax3.xaxis.tick_bottom()
                    ax.yaxis.set_major_formatter(FormatStrFormatter('%.3f'))
                    ax2.yaxis.set_major_formatter(FormatStrFormatter('%.3f'))
                    ax3.yaxis.set_major_formatter(FormatStrFormatter('%.3f'))
                    ax2.set_ylabel(r'flux ($\mu$mol/m$^2$/s)')
                    if rangeplot:
                        ax.errorbar(np.arange(len(aa_syn_dict[cell_list[ii]]))+(ii-1.5)*width,
                            aa_syn_dict[cell_list[ii]],
                                    yerr=[aa_syn_dictmin[cell_list[ii]],aa_syn_dictmax[cell_list[ii]]], 
                                    fmt="o", color="grey",label='_nolegend_',markersize=1,capsize=5)
                        ax2.errorbar(np.arange(len(aa_syn_dict[cell_list[ii]]))+(ii-1.5)*width,
                            aa_syn_dict[cell_list[ii]],
                                    yerr=[aa_syn_dictmin[cell_list[ii]],aa_syn_dictmax[cell_list[ii]]], 
                                    fmt="o", color="grey",label='_nolegend_',markersize=1,capsize=5)
                        ax3.errorbar(np.arange(len(aa_syn_dict[cell_list[ii]]))+(ii-1.5)*width,
                            aa_syn_dict[cell_list[ii]],
                                    yerr=[aa_syn_dictmin[cell_list[ii]],aa_syn_dictmax[cell_list[ii]]], 
                                    fmt="o", color="grey",label='_nolegend_',markersize=1,capsize=5)
                else:
                    ax2.set_ylim(0,0.01)
                    min2 = min([x for x in aa_syn_dict[cell_list[ii]] if x>0.01])
                    max2 = (1+ymax)*max(aa_syn_dict[cell_list[ii]])
                    ax.set_ylim(min2,max2)
                    ax2.spines['top'].set_visible(False)
                    ax.spines['bottom'].set_visible(False)
                    ax.xaxis.tick_top()
                    ax.tick_params(labeltop=False)
                    ax2.xaxis.tick_bottom()
                    ax.yaxis.set_major_formatter(FormatStrFormatter('%.3f'))
                    ax2.yaxis.set_major_formatter(FormatStrFormatter('%.3f'))

                    d = .5  # proportion of vertical to horizontal extent of the slanted line
                    kwargs = dict(marker=[(-1, -d), (1, d)], markersize=12,
                                  linestyle="none", color='k', mec='k', mew=1, clip_on=False)
                    ax.plot([0, 1], [0, 0], transform=ax.transAxes, **kwargs, label='_nolegend_')
#                     plt.ylabel(r'flux ($\mu$mol/m$^2$/s)')
                    ax2.plot([0, 1], [1, 1], transform=ax2.transAxes, **kwargs, label='_nolegend_')
                    ax.set_ylabel(r'flux ($\mu$mol/m$^2$/s)')
                    if rangeplot:
                        ax.errorbar(np.arange(len(aa_syn_dict[cell_list[ii]]))+(ii-1.5)*width,
                            aa_syn_dict[cell_list[ii]],
                                    yerr=[aa_syn_dictmin[cell_list[ii]],aa_syn_dictmax[cell_list[ii]]], 
                                    fmt="o", color="grey",label='_nolegend_',markersize=1,capsize=5)
                        ax2.errorbar(np.arange(len(aa_syn_dict[cell_list[ii]]))+(ii-1.5)*width,
                            aa_syn_dict[cell_list[ii]],
                                    yerr=[aa_syn_dictmin[cell_list[ii]],aa_syn_dictmax[cell_list[ii]]], 
                                    fmt="o", color="grey",label='_nolegend_',markersize=1,capsize=5)
        else:
            plt.bar(np.arange(len(aa_syn_dict[cell_list[ii]]))+(ii-1.5)*width,
                    aa_syn_dict[cell_list[ii]], width=width)
            plt.ylabel(r'flux ($\mu$mol/m$^2$/s)')
            if rangeplot:
                plt.errorbar(np.arange(len(aa_syn_dict[cell_list[ii]]))+(ii-1.5)*width,
                    aa_syn_dict[cell_list[ii]],
                            yerr=[aa_syn_dictmin[cell_list[ii]],aa_syn_dictmax[cell_list[ii]]], 
                            fmt="o", color="grey",label='_nolegend_',markersize=1,capsize=5)
    plt.legend(cell_list)
    if pc=='rat':
        tadd = ' (scaled by cell ratios)'
    else:
        tadd=''
#     plt.title('dielTransfer reactions'+tadd)
    # plt.savefig('RegrEx_plots/'+'1_transfer_large_day.png', bbox_inches='tight')
    plt.show()
    if savebool:
        now = datetime.now().strftime('%Y_%m_%d')
        fig.savefig('Figures/diel_'+savebool+now+'.png',bbox_inches='tight')
    
def plot_subsys_btwn_cells(cons_model, csol,ss=['glycolysis'],thresh=0,plttitle='',phases=['_l','_d'],percell=0,comps=None):
    if not plttitle:
        plttitle=ss[0]
    ss_rxns=[x.id for x in cons_model.reactions if any([y in str(x.notes)for y in ss])]
    if comps:
        mc_rxns=[x for x in ss_rxns if any(['_MC'+y in x for y in phases]) and any([y in x for y in comps])]
    else:
        mc_rxns=[x for x in ss_rxns if any(['_MC'+y in x for y in phases])]
    cc_rxns=[x.replace('_MC_','_CC_') if x.replace('_MC_','_CC_') in ss_rxns else '' for x in mc_rxns]
    se_rxns=[x.replace('_MC_','_SE_') if x.replace('_MC_','_SE_') in ss_rxns else '' for x in mc_rxns]
    pse_rxns=[x.replace('_MC_','_pSE_') if x.replace('_MC_','_pSE_') in ss_rxns else '' for x in mc_rxns]
    datamc=[csol[x] if x else 0 for x in mc_rxns]
    datacc=[csol[x] if x else 0 for x in cc_rxns]
    datase=[csol[x] if x else 0 for x in se_rxns]
    datapse=[csol[x] if x else 0 for x in pse_rxns]
    no_flux=[]
    for ii in range(len(datamc)):
        if abs(datamc[ii])<=thresh and abs(datacc[ii])<=thresh and abs(datase[ii])<=thresh and abs(datapse[ii])<=thresh:
            no_flux.append(ii)
    mc_rxns=[mc_rxns[x] for x in range(len(mc_rxns)) if x not in no_flux]
    cc_rxns=[cc_rxns[x] for x in range(len(cc_rxns)) if x not in no_flux]
    se_rxns=[se_rxns[x] for x in range(len(se_rxns)) if x not in no_flux]
    pse_rxns=[pse_rxns[x] for x in range(len(pse_rxns)) if x not in no_flux]
    datamc=[datamc[x] for x in range(len(datamc)) if x not in no_flux]
    datacc=[datacc[x] for x in range(len(datacc)) if x not in no_flux]
    datase=[datase[x] for x in range(len(datase)) if x not in no_flux]
    datapse=[datapse[x] for x in range(len(datapse)) if x not in no_flux]
    if percell:
        datacc = [20*x for x in datacc]
        datase = [4*x for x in datase]
        datapse = [20*x/6 for x in datapse]
    labels=mc_rxns.copy()
    for x in range(len(labels)):
        temp=labels[x].split('_')
        labels[x]='_'.join(temp[:-2]+[temp[-1]])
    width=0.2
    plt.figure(figsize=(len(labels)*6.4/20,4.8))
    plt.xticks(range(len(datamc)), labels,rotation='vertical')
    plt.bar(np.arange(len(datamc)), datamc, width=width)
    plt.bar(np.arange(len(datacc))+ width, datacc, width=width)
    plt.bar(np.arange(len(datase))+ 2*width, datase, width=width)
    plt.bar(np.arange(len(datapse))+ 3*width, datapse, width=width)
    plt.legend(['MC','CC','SE','pSE'])
    plt.title(plttitle)
    # plt.savefig('RegrEx_plots/'+'1_transfer_large_day.png', bbox_inches='tight')
    plt.show()
                                       
def plot_aa_synth(cons_model, csolin,ss=['glycolysis'],thresh=1e-5,phases=['_l','_d'],ymax=0,pc=0,fva_solin=[],rangeplot=0,ybreak=None,room=0.05,deg=0,comps=0):
    if deg:
        csol = {}
        for rxn in cons_model.reactions:
            csol[rxn.id]=-csolin[rxn.id]
        if fva_solin:
            fva_sol={}
            fva_sol['maximum']={}
            fva_sol['minimum']={}
            for rxn in cons_model.reactions:
                fva_sol['maximum'][rxn.id]=-fva_solin['maximum'][rxn.id]
                fva_sol['minimum'][rxn.id]=-fva_solin['minimum'][rxn.id]
    else:
        csol=csolin
        fva_sol=fva_solin
    aa_list=['ARG_', 'ASN_','GLN_','GLT_','GLY_','HIS_','ILE_','L_ALPHA_ALANINE_',
             'L_ASPARTATE_','LEU_','LYS_','MET_','PHE_','SER_','THR_','TRP_','TYR_','VAL_','CYS_','PRO_']
    if comps:

        aa_list = [x+y for x in aa_list for y in ['c_','p_','m_','x_']]
        aa_list=sorted(aa_list,key=lambda x:x)
        aa_list=sorted(aa_list,key=lambda x:x[-1])
    cell_list=['MC','CC','SE','pSE']
    aa_syn_dict={'MC':[],'CC':[],'SE':[],'pSE':[]}
    aa_deg_dict={'MC':[],'CC':[],'SE':[],'pSE':[]}
    if rangeplot:
        aa_syn_dictmin={'MC':[],'CC':[],'SE':[],'pSE':[]}
        aa_syn_dictmax={'MC':[],'CC':[],'SE':[],'pSE':[]}
        aa_deg_dictmin={'MC':[],'CC':[],'SE':[],'pSE':[]}
        aa_deg_dictmax={'MC':[],'CC':[],'SE':[],'pSE':[]}
    tot_syn=[]
    tot_deg=[]
    aa_syn={}
    aa_deg={}
    if rangeplot:
        tot_synmin=[]
        tot_degmin=[]
        tot_synmax=[]
        tot_degmax=[]
    # List reactions that produce/destroy each AA
    for aa in aa_list:
        if rangeplot:
            aa_syn[aa]=[rxn for rxn in cons_model.reactions if any([((y.id.startswith(aa)) and (fva_sol['minimum'][rxn.id]*rxn.get_coefficient(y)>thresh or fva_sol['maximum'][rxn.id]*rxn.get_coefficient(y)>thresh) and 'RNA' not in y.id) for y in rxn.metabolites])]
            aa_deg[aa]=[rxn for rxn in cons_model.reactions if any([((y.id.startswith(aa)) and (fva_sol['minimum'][rxn.id]*rxn.get_coefficient(y)<-thresh or fva_sol['maximum'][rxn.id]*rxn.get_coefficient(y)<-thresh) and 'RNA' not in y.id) for y in rxn.metabolites])]
        else:
            aa_syn[aa]=[rxn for rxn in cons_model.reactions if any([((y.id.startswith(aa)) and csol[rxn.id]*rxn.get_coefficient(y)>thresh and 'RNA' not in y.id) for y in rxn.metabolites])]
            aa_deg[aa]=[rxn for rxn in cons_model.reactions if any([((y.id.startswith(aa)) and csol[rxn.id]*rxn.get_coefficient(y)<-thresh and 'RNA' not in y.id) for y in rxn.metabolites])]
        aa_syn[aa]=[x for x in aa_syn[aa] if any([y in x.id for y in phases]) and len(x.id.split('_')[-3])<2]
        aa_deg[aa]=[x for x in aa_deg[aa] if any([y in x.id for y in phases]) and len(x.id.split('_')[-3])<2]
        
        cell_syn={}
        cell_deg={}
        if rangeplot:
            cell_synmin={}
            cell_synmax={}
            cell_degmin={}
            cell_degmax={}
        # For each cell, list...
        for ii in range(4):
            cell_syn[cell_list[ii]]=[csol[rxn.id]*rxn.get_coefficient(y) for rxn in aa_syn[aa] for y in rxn.metabolites if any([(y.id.startswith(aa))])and '_'+cell_list[ii]+'_' in rxn.id]
            cell_deg[cell_list[ii]]=[csol[rxn.id]*rxn.get_coefficient(y) for rxn in aa_deg[aa] for y in rxn.metabolites if any([(y.id.startswith(aa))])and '_'+cell_list[ii]+'_' in rxn.id]
            aa_syn_dict[cell_list[ii]].append(sum(cell_syn[cell_list[ii]]))
            aa_deg_dict[cell_list[ii]].append(-sum(cell_deg[cell_list[ii]]))
            if rangeplot:
                cell_synmin[cell_list[ii]]=[fva_sol['minimum'][rxn.id]*rxn.get_coefficient(y) if rxn.get_coefficient(y)>=0 else fva_sol['maximum'][rxn.id]*rxn.get_coefficient(y) for rxn in aa_syn[aa] for y in rxn.metabolites if any([(y.id.startswith(aa))])and '_'+cell_list[ii]+'_' in rxn.id]
                cell_synmax[cell_list[ii]]=[fva_sol['maximum'][rxn.id]*rxn.get_coefficient(y) if rxn.get_coefficient(y)>=0 else fva_sol['minimum'][rxn.id]*rxn.get_coefficient(y) for rxn in aa_syn[aa] for y in rxn.metabolites if any([(y.id.startswith(aa))])and '_'+cell_list[ii]+'_' in rxn.id]
                cell_degmin[cell_list[ii]]=[fva_sol['minimum'][rxn.id]*rxn.get_coefficient(y) if rxn.get_coefficient(y)>=0 else fva_sol['maximum'][rxn.id]*rxn.get_coefficient(y) for rxn in aa_deg[aa] for y in rxn.metabolites if any([(y.id.startswith(aa))])and '_'+cell_list[ii]+'_' in rxn.id]
                cell_degmax[cell_list[ii]]=[fva_sol['maximum'][rxn.id]*rxn.get_coefficient(y) if rxn.get_coefficient(y)>=0 else fva_sol['minimum'][rxn.id]*rxn.get_coefficient(y) for rxn in aa_deg[aa] for y in rxn.metabolites if any([(y.id.startswith(aa))])and '_'+cell_list[ii]+'_' in rxn.id]
                aa_syn_dictmin[cell_list[ii]].append(sum(cell_synmin[cell_list[ii]]))
                aa_syn_dictmax[cell_list[ii]].append(sum(cell_synmax[cell_list[ii]]))
                aa_deg_dictmin[cell_list[ii]].append(-sum(cell_degmin[cell_list[ii]]))
                aa_deg_dictmax[cell_list[ii]].append(-sum(cell_degmax[cell_list[ii]]))
#         print(cell_syn)
        tot_syn.append(sum([aa_syn_dict[cell_list[x]][-1] for x in range(4)]))
        tot_deg.append(sum([aa_deg_dict[cell_list[x]][-1] for x in range(4)]))
        if rangeplot:
            tot_synmin.append(sum([aa_syn_dictmin[cell_list[x]][-1] for x in range(4)]))
            tot_degmin.append(sum([aa_deg_dictmin[cell_list[x]][-1] for x in range(4)]))
            tot_synmax.append(sum([aa_syn_dictmax[cell_list[x]][-1] for x in range(4)]))
            tot_degmax.append(sum([aa_deg_dictmax[cell_list[x]][-1] for x in range(4)]))
    if pc=='%':
        for ii in range(4):
#         print(aa_syn_dict,tot_syn)
            aa_syn_dict[cell_list[ii]]=[aa_syn_dict[cell_list[ii]][x]/tot_syn[x] if tot_syn[x]!=0 else 0 for x in range(len(aa_syn_dict[cell_list[ii]]))]
    #         print(aa_syn_dict)
            aa_deg_dict[cell_list[ii]]=[aa_deg_dict[cell_list[ii]][x]/tot_deg[x] if tot_deg[x]!=0 else 0 for x in range(len(aa_deg_dict[cell_list[ii]]))]
    elif pc=='rat':
        ratios=[1,20,5,6]
        mc_syn = aa_syn_dict[cell_list[0]]
        mc_deg = aa_deg_dict[cell_list[0]]
        for ii in range(4):
            aa_syn_dict[cell_list[ii]]=[aa_syn_dict[cell_list[ii]][x]*ratios[ii]/mc_syn[x] if mc_syn[x]!=0 else aa_syn_dict[cell_list[ii]][x]*ratios[ii]/max(tot_syn[x],0.01) for x in range(len(aa_syn_dict[cell_list[ii]]))]
            aa_deg_dict[cell_list[ii]]=[aa_deg_dict[cell_list[ii]][x]*ratios[ii]/mc_deg[x] if mc_deg[x]!=0 else aa_deg_dict[cell_list[ii]][x]*ratios[ii]/max(tot_deg[x],0.01) for x in range(len(aa_deg_dict[cell_list[ii]]))]
            if rangeplot:
                aa_syn_dictmin[cell_list[ii]]=[aa_syn_dictmin[cell_list[ii]][x]*ratios[ii]/mc_syn[x] if mc_syn[x]!=0 else aa_syn_dictmin[cell_list[ii]][x]*ratios[ii]/max(tot_syn[x],0.01) for x in range(len(aa_syn_dictmin[cell_list[ii]]))]
                aa_syn_dictmin[cell_list[ii]] = [abs(aa_syn_dictmin[cell_list[ii]][x]-aa_syn_dict[cell_list[ii]][x]) for x in range(len(aa_syn_dictmin[cell_list[ii]]))]
                aa_deg_dictmin[cell_list[ii]]=[aa_deg_dictmin[cell_list[ii]][x]*ratios[ii]/mc_deg[x] if mc_deg[x]!=0 else aa_deg_dictmin[cell_list[ii]][x]*ratios[ii]/max(tot_deg[x],0.01) for x in range(len(aa_deg_dictmin[cell_list[ii]]))]
                aa_deg_dictmin[cell_list[ii]] = [abs(aa_deg_dictmin[cell_list[ii]][x]-aa_deg_dict[cell_list[ii]][x]) for x in range(len(aa_deg_dictmin[cell_list[ii]]))]
                aa_syn_dictmax[cell_list[ii]]=[aa_syn_dictmax[cell_list[ii]][x]*ratios[ii]/mc_syn[x] if mc_syn[x]!=0 else aa_syn_dictmax[cell_list[ii]][x]*ratios[ii]/max(tot_syn[x],0.01) for x in range(len(aa_syn_dictmax[cell_list[ii]]))]
                aa_syn_dictmax[cell_list[ii]] = [aa_syn_dictmax[cell_list[ii]][x]-aa_syn_dict[cell_list[ii]][x] for x in range(len(aa_syn_dictmax[cell_list[ii]]))]
                aa_deg_dictmax[cell_list[ii]]=[aa_deg_dictmax[cell_list[ii]][x]*ratios[ii]/mc_deg[x] if mc_deg[x]!=0 else aa_deg_dictmax[cell_list[ii]][x]*ratios[ii]/max(tot_deg[x],0.01) for x in range(len(aa_deg_dictmax[cell_list[ii]]))]
                aa_deg_dictmax[cell_list[ii]] = [aa_deg_dictmax[cell_list[ii]][x]-aa_deg_dict[cell_list[ii]][x] for x in range(len(aa_deg_dictmax[cell_list[ii]]))]
    elif pc=='bycell':
        ratios=[1,20,4,10/3]
        for ii in range(4):
            aa_syn_dict[cell_list[ii]]=[aa_syn_dict[cell_list[ii]][x]*ratios[ii] for x in range(len(aa_syn_dict[cell_list[ii]]))]
            aa_deg_dict[cell_list[ii]]=[aa_deg_dict[cell_list[ii]][x]*ratios[ii] for x in range(len(aa_deg_dict[cell_list[ii]]))]
            if rangeplot:
                aa_syn_dictmin[cell_list[ii]]=[aa_syn_dictmin[cell_list[ii]][x]*ratios[ii] for x in range(len(aa_syn_dictmin[cell_list[ii]]))]
                aa_syn_dictmin[cell_list[ii]] = [abs(aa_syn_dictmin[cell_list[ii]][x]-aa_syn_dict[cell_list[ii]][x]) for x in range(len(aa_syn_dictmin[cell_list[ii]]))]
                aa_syn_dictmax[cell_list[ii]]=[aa_syn_dictmax[cell_list[ii]][x]*ratios[ii] for x in range(len(aa_syn_dictmax[cell_list[ii]]))]
                aa_syn_dictmax[cell_list[ii]] = [aa_syn_dictmax[cell_list[ii]][x]-aa_syn_dict[cell_list[ii]][x] for x in range(len(aa_syn_dictmax[cell_list[ii]]))]
    else:
        # print('blank pc')
        if rangeplot:
            for ii in range(4):
                # tempdict=aa_syn_dictmin[cell_list[ii]]
                aa_syn_dictmin[cell_list[ii]] = [abs(aa_syn_dictmin[cell_list[ii]][x]-aa_syn_dict[cell_list[ii]][x]) for x in range(len(aa_syn_dictmin[cell_list[ii]]))]
                aa_deg_dictmin[cell_list[ii]] = [abs(aa_deg_dictmin[cell_list[ii]][x]-aa_deg_dict[cell_list[ii]][x]) for x in range(len(aa_deg_dictmin[cell_list[ii]]))]
                aa_syn_dictmax[cell_list[ii]] = [aa_syn_dictmax[cell_list[ii]][x]-aa_syn_dict[cell_list[ii]][x] for x in range(len(aa_syn_dictmax[cell_list[ii]]))]
                aa_deg_dictmax[cell_list[ii]] = [aa_deg_dictmax[cell_list[ii]][x]-aa_deg_dict[cell_list[ii]][x] for x in range(len(aa_deg_dictmax[cell_list[ii]]))]
            # print(tempdict)
            # print(aa_syn_dictmin[cell_list[ii]])

    labels=[x[:-1] for x in aa_list]
    width=0.2
    if ybreak:
        fig,(ax,ax2) = plt.subplots(2,1, sharex=True,figsize=(len(labels)/2,4.8))
    else:
        plt.figure(figsize=(len(labels)/2,4.8))
    plt.xticks(range(len(aa_syn_dict[cell_list[ii]])), labels,rotation='vertical')
    for ii in range(len(cell_list)):
        if ybreak:
            breakthresh=ybreak
            ax.bar(np.arange(len(aa_syn_dict[cell_list[ii]]))+(ii-1.5)*width,
                aa_syn_dict[cell_list[ii]], width=width)
            ax2.bar(np.arange(len(aa_syn_dict[cell_list[ii]]))+(ii-1.5)*width,
                aa_syn_dict[cell_list[ii]], width=width)
            if ii==0:
                ax2.set_ylim(0,breakthresh)
                # print(aa_syn_dict[cell_list[ii]])
                min2 = (1-room)*min([x for x in aa_syn_dict[cell_list[ii]] if x>breakthresh])
                max2 = (1+room)*max(aa_syn_dict[cell_list[ii]])
                ax.set_ylim(min2,max(max2,ymax))
                ax2.spines['top'].set_visible(False)
                ax.spines['bottom'].set_visible(False)
                ax.xaxis.tick_top()
                ax.tick_params(labeltop=False)
                ax.yaxis.set_major_formatter(FormatStrFormatter('%.3f'))
                ax2.yaxis.set_major_formatter(FormatStrFormatter('%.3f'))

                d = .5  # proportion of vertical to horizontal extent of the slanted line
                kwargs = dict(marker=[(-1, -d), (1, d)], markersize=12,
                              linestyle="none", color='k', mec='k', mew=1, clip_on=False)
                ax.plot([0, 1], [0, 0], transform=ax.transAxes, **kwargs, label='_nolegend_')
                ax2.plot([0, 1], [1, 1], transform=ax2.transAxes, **kwargs, label='_nolegend_')
            if rangeplot:
                # if any([x<-1e-3 for x in aa_syn_dictmax[cell_list[ii]]]):
                #     print('problem with ',[labels[x] for x in range(len(aa_syn_dictmax[cell_list[ii]])) if aa_syn_dictmax[cell_list[ii]][x]<-1e-3],
                #     [aa_syn_dictmax[cell_list[ii]][x] for x in range(len(aa_syn_dictmax[cell_list[ii]])) if aa_syn_dictmax[cell_list[ii]][x]<-1e-3])
                # if any([x<0 for x in aa_syn_dictmin[cell_list[ii]]]):
                #     print('problem with ',[labels[x] for x in range(len(aa_syn_dictmin[cell_list[ii]])) if aa_syn_dictmin[cell_list[ii]][x]<-1e-3],
                #     [aa_syn_dictmin[cell_list[ii]][x] for x in range(len(aa_syn_dictmin[cell_list[ii]])) if aa_syn_dictmin[cell_list[ii]][x]<-1e-3])
                ax.errorbar(np.arange(len(aa_syn_dict[cell_list[ii]]))+(ii-1.5)*width,
                    aa_syn_dict[cell_list[ii]],
                    yerr=[aa_syn_dictmin[cell_list[ii]],aa_syn_dictmax[cell_list[ii]]],
                    fmt="o", color="grey",label='_nolegend_',markersize=1,capsize=5)
                ax2.errorbar(np.arange(len(aa_syn_dict[cell_list[ii]]))+(ii-1.5)*width,
                    aa_syn_dict[cell_list[ii]],
                    yerr=[aa_syn_dictmin[cell_list[ii]],aa_syn_dictmax[cell_list[ii]]], 
                    fmt="o", color="grey",label='_nolegend_',markersize=1,capsize=5)
        else:
            plt.bar(np.arange(len(aa_syn_dict[cell_list[ii]]))+(ii-1.5)*width,
                    aa_syn_dict[cell_list[ii]], width=width)
            if rangeplot:
                plt.errorbar(np.arange(len(aa_syn_dict[cell_list[ii]]))+(ii-1.5)*width,
                    aa_syn_dict[cell_list[ii]],
                             yerr=[aa_syn_dictmin[cell_list[ii]],aa_syn_dictmax[cell_list[ii]]], 
                             fmt="o", color="black",label='_nolegend_')
    if ybreak:
        ax.legend(cell_list)
    else:
        plt.legend(cell_list)
#     plt.title('AA synthesis')
    if ymax and not ybreak:
        plt.ylim([0,ymax])
    # plt.savefig('RegrEx_plots/'+'1_transfer_large_day.png', bbox_inches='tight')
    plt.ylabel(r'flux ($\mu$mol/m$^2$/s)')
    plt.show()
    # if deg:
    #     plt.figure(figsize=(len(labels)/2,4.8))
    #     plt.xticks(range(len(aa_deg_dict[cell_list[ii]])), labels,rotation='vertical')
    #     for ii in range(len(cell_list)):
    #         plt.bar(np.arange(len(aa_deg_dict[cell_list[ii]]))+(ii-1.5)*width,
    #                 aa_deg_dict[cell_list[ii]], width=width)
    #     plt.legend(cell_list)
    #     plt.ylabel(r'flux ($\mu$mol/m$^2$/s)')
    #     plt.title('AA degradation')
    #     if ymax:
    #         plt.ylim([0,ymax])
    #     # plt.savefig('RegrEx_plots/'+'1_transfer_large_day.png', bbox_inches='tight')
    #     plt.show()
    return aa_syn,aa_deg,aa_syn_dict,labels

def plot_aa_synth_presentation(cons_model, csol,ss=['glycolysis'],thresh=1e-5,phases=['_l','_d'],ymax=0,pc=0,fva_sol=[],rangeplot=0,ybreak=None,room=0.05,deg=0,comps=None):
    aa_list=['ARG_', 'ASN_','GLN_','GLT_','GLY_','HIS_','ILE_','L_ALPHA_ALANINE_',
             'L_ASPARTATE_','LEU_','LYS_','MET_','PHE_','SER_','THR_','TRP_','TYR_','VAL_','CYS_']
    if comps:
        aa_list = [x+y for x in aa_list for y in ['c_','p_','m_','x_']]
        aa_list=sorted(aa_list,key=lambda x:x)
        aa_list=sorted(aa_list,key=lambda x:x[-1])
    cell_list=['MC','CC','pSE']
    aa_syn_dict={'MC':[],'CC':[],'pSE':[]}
    if rangeplot:
        aa_syn_dictmin={'MC':[],'CC':[],'pSE':[]}
        aa_syn_dictmax={'MC':[],'CC':[],'pSE':[]}
    tot_syn=[]
    aa_syn={}
    if rangeplot:
        tot_synmin=[]
        tot_synmax=[]
    # List reactions that produce/destroy each AA
    for aa in aa_list:
        if rangeplot:
            aa_syn[aa]=[rxn for rxn in cons_model.reactions if any([((y.id.startswith(aa)) and (fva_sol['minimum'][rxn.id]*rxn.get_coefficient(y)>thresh or fva_sol['maximum'][rxn.id]*rxn.get_coefficient(y)>thresh) and 'RNA' not in y.id) for y in rxn.metabolites])]
        else:
            aa_syn[aa]=[rxn for rxn in cons_model.reactions if any([((y.id.startswith(aa)) and csol[rxn.id]*rxn.get_coefficient(y)>thresh and 'RNA' not in y.id) for y in rxn.metabolites])]
        aa_syn[aa]=[x for x in aa_syn[aa] if any([y in x.id for y in phases]) and len(x.id.split('_')[-3])<2]
        
        cell_syn={}
        if rangeplot:
            cell_synmin={}
            cell_synmax={}
        # For each cell, list...
        for ii in range(len(cell_list)):
            cell_syn[cell_list[ii]]=[csol[rxn.id]*rxn.get_coefficient(y) for rxn in aa_syn[aa] for y in rxn.metabolites if any([(y.id.startswith(aa))])and '_'+cell_list[ii]+'_' in rxn.id]
            aa_syn_dict[cell_list[ii]].append(sum(cell_syn[cell_list[ii]]))
            if rangeplot:
                cell_synmin[cell_list[ii]]=[fva_sol['minimum'][rxn.id]*rxn.get_coefficient(y) if rxn.get_coefficient(y)>=0 else fva_sol['maximum'][rxn.id]*rxn.get_coefficient(y) for rxn in aa_syn[aa] for y in rxn.metabolites if any([(y.id.startswith(aa))])and '_'+cell_list[ii]+'_' in rxn.id]
                cell_synmax[cell_list[ii]]=[fva_sol['maximum'][rxn.id]*rxn.get_coefficient(y) if rxn.get_coefficient(y)>=0 else fva_sol['minimum'][rxn.id]*rxn.get_coefficient(y) for rxn in aa_syn[aa] for y in rxn.metabolites if any([(y.id.startswith(aa))])and '_'+cell_list[ii]+'_' in rxn.id]
                aa_syn_dictmin[cell_list[ii]].append(sum(cell_synmin[cell_list[ii]]))
                aa_syn_dictmax[cell_list[ii]].append(sum(cell_synmax[cell_list[ii]]))
            
        tot_syn.append(sum([aa_syn_dict[cell_list[x]][-1] for x in range(len(cell_list))]))
        if rangeplot:
            tot_synmin.append(sum([aa_syn_dictmin[cell_list[x]][-1] for x in range(len(cell_list))]))
            tot_synmax.append(sum([aa_syn_dictmax[cell_list[x]][-1] for x in range(len(cell_list))]))
    if pc=='%':
        for ii in range(len(cell_list)):
#         print(aa_syn_dict,tot_syn)
            aa_syn_dict[cell_list[ii]]=[aa_syn_dict[cell_list[ii]][x]/tot_syn[x] if tot_syn[x]!=0 else 0 for x in range(len(aa_syn_dict[cell_list[ii]]))]
    #         print(aa_syn_dict)
    elif pc=='rat':
        ratios=[1,20,5,6]
        mc_syn = aa_syn_dict[cell_list[0]]
        for ii in range(len(cell_list)):
            aa_syn_dict[cell_list[ii]]=[aa_syn_dict[cell_list[ii]][x]*ratios[ii]/mc_syn[x] if mc_syn[x]!=0 else aa_syn_dict[cell_list[ii]][x]*ratios[ii]/max(tot_syn[x],0.01) for x in range(len(aa_syn_dict[cell_list[ii]]))]
            if rangeplot:
                aa_syn_dictmin[cell_list[ii]]=[aa_syn_dictmin[cell_list[ii]][x]*ratios[ii]/mc_syn[x] if mc_syn[x]!=0 else aa_syn_dictmin[cell_list[ii]][x]*ratios[ii]/max(tot_syn[x],0.01) for x in range(len(aa_syn_dictmin[cell_list[ii]]))]
                aa_syn_dictmin[cell_list[ii]] = [abs(aa_syn_dictmin[cell_list[ii]][x]-aa_syn_dict[cell_list[ii]][x]) for x in range(len(aa_syn_dictmin[cell_list[ii]]))]
                aa_syn_dictmax[cell_list[ii]]=[aa_syn_dictmax[cell_list[ii]][x]*ratios[ii]/mc_syn[x] if mc_syn[x]!=0 else aa_syn_dictmax[cell_list[ii]][x]*ratios[ii]/max(tot_syn[x],0.01) for x in range(len(aa_syn_dictmax[cell_list[ii]]))]
                aa_syn_dictmax[cell_list[ii]] = [aa_syn_dictmax[cell_list[ii]][x]-aa_syn_dict[cell_list[ii]][x] for x in range(len(aa_syn_dictmax[cell_list[ii]]))]
    else:
        # print('blank pc')
        if rangeplot:
            for ii in range(len(cell_list)):
                # tempdict=aa_syn_dictmin[cell_list[ii]]
                aa_syn_dictmin[cell_list[ii]] = [abs(aa_syn_dictmin[cell_list[ii]][x]-aa_syn_dict[cell_list[ii]][x]) for x in range(len(aa_syn_dictmin[cell_list[ii]]))]
                aa_syn_dictmax[cell_list[ii]] = [aa_syn_dictmax[cell_list[ii]][x]-aa_syn_dict[cell_list[ii]][x] for x in range(len(aa_syn_dictmax[cell_list[ii]]))]
                
    labels=[x[:-1] for x in aa_list]
    width=0.6
    if ybreak:
        fig,axes = plt.subplots(2,3, sharex=True,figsize=(len(labels)/2*3,4.8))
    else:
        fig,axes = plt.subplots(1,3, sharex=True,figsize=(len(labels)/2*3,4.8))

    # colours = ['blue','orange','green','red']
    colours = plt.rcParams['axes.prop_cycle'].by_key()['color']
    handles=[]
    for ii in [0,1,2]:
        if ybreak:
            breakthresh=ybreak
            if ii==0:
                min2 = (1-room)*min([x for x in aa_syn_dict[cell_list[0]] if x>breakthresh])
                max2 = (1+room)*max(aa_syn_dict[cell_list[0]])
            handles+=[axes[0,ii].bar(np.arange(len(aa_syn_dict[cell_list[ii]])),
                aa_syn_dict[cell_list[ii]], width=width,color=colours[ii])]
            axes[1,ii].bar(np.arange(len(aa_syn_dict[cell_list[ii]])),
                aa_syn_dict[cell_list[ii]], width=width,color=colours[ii])
            axes[1,ii].set_ylim(0,breakthresh)
            axes[0,ii].set_ylim(min2,max(max2,ymax))
            axes[1,ii].spines['top'].set_visible(False)
            axes[0,ii].spines['bottom'].set_visible(False)
            axes[0,ii].xaxis.tick_top()
            axes[0,ii].tick_params(labeltop=False)
            axes[0,ii].yaxis.set_major_formatter(FormatStrFormatter('%.3f'))
            axes[1,ii].yaxis.set_major_formatter(FormatStrFormatter('%.3f'))
            d = .5  # proportion of vertical to horizontal extent of the slanted line
            kwargs = dict(marker=[(-1, -d), (1, d)], markersize=12,
                            linestyle="none", color='k', mec='k', mew=1, clip_on=False)
            axes[0,ii].plot([0, 1], [0, 0], transform=axes[0,ii].transAxes, **kwargs, label='_nolegend_')
            axes[1,ii].plot([0, 1], [1, 1], transform=axes[1,ii].transAxes, **kwargs, label='_nolegend_')
            plt.subplot(2,3,ii+4)
            plt.xticks(range(len(aa_syn_dict[cell_list[ii]])), labels,rotation='vertical')
            plt.rc('font', size=18)
            axes[1,ii].set_ylabel(r'flux')
            axes[0,ii].set_ylabel(r'($\mu$mol/m$^2$/s)')
            if rangeplot:
                # if any([x<-1e-3 for x in aa_syn_dictmax[cell_list[ii]]]):
                #     print('problem with ',[labels[x] for x in range(len(aa_syn_dictmax[cell_list[ii]])) if aa_syn_dictmax[cell_list[ii]][x]<-1e-3],
                #     [aa_syn_dictmax[cell_list[ii]][x] for x in range(len(aa_syn_dictmax[cell_list[ii]])) if aa_syn_dictmax[cell_list[ii]][x]<-1e-3])
                # if any([x<0 for x in aa_syn_dictmin[cell_list[ii]]]):
                #     print('problem with ',[labels[x] for x in range(len(aa_syn_dictmin[cell_list[ii]])) if aa_syn_dictmin[cell_list[ii]][x]<-1e-3],
                #     [aa_syn_dictmin[cell_list[ii]][x] for x in range(len(aa_syn_dictmin[cell_list[ii]])) if aa_syn_dictmin[cell_list[ii]][x]<-1e-3])
                axes[0,ii].errorbar(np.arange(len(aa_syn_dict[cell_list[ii]]))+(ii-1.5)*width,
                    aa_syn_dict[cell_list[ii]],
                    yerr=[aa_syn_dictmin[cell_list[ii]],aa_syn_dictmax[cell_list[ii]]],
                    fmt="o", color="grey",label='_nolegend_',markersize=1,capsize=5)
                axes[0,ii].errorbar(np.arange(len(aa_syn_dict[cell_list[ii]]))+(ii-1.5)*width,
                    aa_syn_dict[cell_list[ii]],
                    yerr=[aa_syn_dictmin[cell_list[ii]],aa_syn_dictmax[cell_list[ii]]], 
                    fmt="o", color="grey",label='_nolegend_',markersize=1,capsize=5)
        else:
            axes[1,ii].bar(np.arange(len(aa_syn_dict[cell_list[ii]]))+(ii-1.5)*width,
                    aa_syn_dict[cell_list[ii]], width=width,color=colours[ii])
            plt.xticks(range(len(aa_syn_dict[cell_list[ii]])), labels,rotation='vertical')
            if rangeplot:
                axes[1,ii].errorbar(np.arange(len(aa_syn_dict[cell_list[ii]]))+(ii-1.5)*width,
                    aa_syn_dict[cell_list[ii]],
                             yerr=[aa_syn_dictmin[cell_list[ii]],aa_syn_dictmax[cell_list[ii]]], 
                             fmt="o", color="black",label='_nolegend_')
    legendlist = ['Mesophyll', 'Companion Cells','Sieve Elements']
    if ybreak:
        axes[0,ii].legend(handles,legendlist)
    else:
        plt.legend(legendlist)
    if ymax and not ybreak:
        plt.ylim([0,ymax])
    # plt.savefig('RegrEx_plots/'+'1_transfer_large_day.png', bbox_inches='tight')
    plt.show()
    return aa_syn


def plot_met_catabolism(cons_model, csol,mets=['SUCROSE'],thresh=1e-5,phases=['_l','_d'],ymax=0,pc=0,fva_sol=[],rangeplot=0,ybreak=None,room=0.05,deg=0):
    aa_list=['ARG_', 'ASN_','GLN_','GLT_','GLY_','HIS_','ILE_','L_ALPHA_ALANINE_',
             'L_ASPARTATE_','LEU_','LYS_','MET_','PHE_','SER_','THR_','TRP_','TYR_','VAL_','CYS_']
    aa_list=mets+aa_list
    cell_list=['MC','CC','SE','pSE']
    met_cat_dict={'MC':[],'CC':[],'SE':[],'pSE':[]}
    if rangeplot:
        met_cat_dictmin={'MC':[],'CC':[],'SE':[],'pSE':[]}
        met_cat_dictmax={'MC':[],'CC':[],'SE':[],'pSE':[]}
    tot_cat=[]
    met_cat={}
    if rangeplot:
        tot_catmin=[]
        tot_catmax=[]
    # List reactions that produce/destroy each AA
    for aa in aa_list:
        if rangeplot:
            met_cat[aa]=[rxn for rxn in cons_model.reactions if any([((y.id.startswith(aa)) and (fva_sol['minimum'][rxn.id]*rxn.get_coefficient(y.id)<-thresh or fva_sol['maximum'][rxn.id]*rxn.get_coefficient(y.id)<-thresh) and 'RNA' not in y.id) for y in rxn.metabolites])]
        else:
            met_cat[aa]=[rxn for rxn in cons_model.reactions if any([((y.id.startswith(aa)) and csol[rxn.id]*rxn.get_coefficient(y)<-thresh and 'RNA' not in y.id) for y in rxn.metabolites])]
        met_cat[aa]=[x for x in met_cat[aa] if any([y in x.id for y in phases]) and len(x.id.split('_')[-3])<2]
        
        cell_cat={}
        if rangeplot:
            cell_catmin={}
            cell_catmax={}
        # For each cell, list...
        for ii in range(len(cell_list)):
            cell_cat[cell_list[ii]]=[csol[rxn.id]*rxn.get_coefficient(y) for rxn in met_cat[aa] for y in rxn.metabolites if any([(y.id.startswith(aa))])and '_'+cell_list[ii]+'_' in rxn.id]
            met_cat_dict[cell_list[ii]].append(-sum(cell_cat[cell_list[ii]]))
            if rangeplot:
                cell_catmin[cell_list[ii]]=[fva_sol['minimum'][rxn.id]*rxn.get_coefficient(y) if rxn.get_coefficient(y)<0 else fva_sol['maximum'][rxn.id]*rxn.get_coefficient(y) for rxn in met_cat[aa] for y in rxn.metabolites if any([(y.id.startswith(aa))])and '_'+cell_list[ii]+'_' in rxn.id]
                cell_catmax[cell_list[ii]]=[fva_sol['maximum'][rxn.id]*rxn.get_coefficient(y) if rxn.get_coefficient(y)<=0 else fva_sol['minimum'][rxn.id]*rxn.get_coefficient(y) for rxn in met_cat[aa] for y in rxn.metabolites if any([(y.id.startswith(aa))])and '_'+cell_list[ii]+'_' in rxn.id]
                met_cat_dictmin[cell_list[ii]].append(-sum(cell_catmin[cell_list[ii]]))
                met_cat_dictmax[cell_list[ii]].append(-sum(cell_catmax[cell_list[ii]]))
        tot_cat.append(sum([met_cat_dict[cell_list[x]][-1] for x in range(len(cell_list))]))
        if rangeplot:
            tot_catmin.append(sum([met_cat_dictmin[cell_list[x]][-1] for x in range(len(cell_list))]))
            tot_catmax.append(sum([met_cat_dictmax[cell_list[x]][-1] for x in range(len(cell_list))]))
    if pc=='%':
        for ii in range(len(cell_list)):
            met_cat_dict[cell_list[ii]]=[met_cat_dict[cell_list[ii]][x]/tot_cat[x] if tot_cat[x]!=0 else 0 for x in range(len(met_cat_dict[cell_list[ii]]))]
    elif pc=='rat':
        ratios=[1,20,5,6]
        mc_cat = met_cat_dict[cell_list[0]]
        for ii in range(len(cell_list)):
            met_cat_dict[cell_list[ii]]=[met_cat_dict[cell_list[ii]][x]*ratios[ii]/mc_cat[x] if mc_cat[x]!=0 else met_cat_dict[cell_list[ii]][x]*ratios[ii]/max(tot_cat[x],0.01) for x in range(len(met_cat_dict[cell_list[ii]]))]
            if rangeplot:
                met_cat_dictmin[cell_list[ii]]=[met_cat_dictmin[cell_list[ii]][x]*ratios[ii]/mc_cat[x] if mc_cat[x]!=0 else met_cat_dictmin[cell_list[ii]][x]*ratios[ii]/max(tot_cat[x],0.01) for x in range(len(met_cat_dictmin[cell_list[ii]]))]
                met_cat_dictmin[cell_list[ii]] = [abs(met_cat_dictmin[cell_list[ii]][x]-met_cat_dict[cell_list[ii]][x]) for x in range(len(met_cat_dictmin[cell_list[ii]]))]
                met_cat_dictmax[cell_list[ii]]=[met_cat_dictmax[cell_list[ii]][x]*ratios[ii]/mc_cat[x] if mc_cat[x]!=0 else met_cat_dictmax[cell_list[ii]][x]*ratios[ii]/max(tot_cat[x],0.01) for x in range(len(met_cat_dictmax[cell_list[ii]]))]
                met_cat_dictmax[cell_list[ii]] = [met_cat_dictmax[cell_list[ii]][x]-met_cat_dict[cell_list[ii]][x] for x in range(len(met_cat_dictmax[cell_list[ii]]))]
    else:
        # print('blank pc')
        if rangeplot:
            for ii in range(len(cell_list)):
                # tempdict=met_cat_dictmin[cell_list[ii]]
                met_cat_dictmin[cell_list[ii]] = [abs(met_cat_dictmin[cell_list[ii]][x]-met_cat_dict[cell_list[ii]][x]) for x in range(len(met_cat_dictmin[cell_list[ii]]))]
                met_cat_dictmax[cell_list[ii]] = [met_cat_dictmax[cell_list[ii]][x]-met_cat_dict[cell_list[ii]][x] for x in range(len(met_cat_dictmax[cell_list[ii]]))]
            # print(tempdict)
            # print(met_cat_dictmin[cell_list[ii]])

    labels=[x[:-1] for x in aa_list]
    width=0.2
    if ybreak:
        fig,(ax,ax2) = plt.subplots(2,1, sharex=True,figsize=(len(labels)/2,4.8))
    else:
        plt.figure(figsize=(len(labels)/2,4.8))
    plt.xticks(range(len(met_cat_dict[cell_list[ii]])), labels,rotation='vertical')
    for ii in range(len(cell_list)):
        if ybreak:
            breakthresh=ybreak
            ax.bar(np.arange(len(met_cat_dict[cell_list[ii]]))+(ii-1.5)*width,
                met_cat_dict[cell_list[ii]], width=width)
            ax2.bar(np.arange(len(met_cat_dict[cell_list[ii]]))+(ii-1.5)*width,
                met_cat_dict[cell_list[ii]], width=width)
            if ii==0:
                ax2.set_ylim(0,breakthresh)
                min2 = (1-room)*min([x for x in met_cat_dict[cell_list[ii]] if x>breakthresh]+[breakthresh])
                max2 = (1+room)*max(met_cat_dict[cell_list[ii]])
                ax.set_ylim(min2,max(max2,ymax))
                ax2.spines['top'].set_visible(False)
                ax.spines['bottom'].set_visible(False)
                ax.xaxis.tick_top()
                ax.tick_params(labeltop=False)
                ax.yaxis.set_major_formatter(FormatStrFormatter('%.3f'))
                ax2.yaxis.set_major_formatter(FormatStrFormatter('%.3f'))

                d = .5  # proportion of vertical to horizontal extent of the slanted line
                kwargs = dict(marker=[(-1, -d), (1, d)], markersize=12,
                              linestyle="none", color='k', mec='k', mew=1, clip_on=False)
                ax.plot([0, 1], [0, 0], transform=ax.transAxes, **kwargs, label='_nolegend_')
                ax2.plot([0, 1], [1, 1], transform=ax2.transAxes, **kwargs, label='_nolegend_')
            if rangeplot:
                # if any([x<-1e-3 for x in met_cat_dictmax[cell_list[ii]]]):
                #     print('problem with ',[labels[x] for x in range(len(met_cat_dictmax[cell_list[ii]])) if met_cat_dictmax[cell_list[ii]][x]<-1e-3],
                #     [met_cat_dictmax[cell_list[ii]][x] for x in range(len(met_cat_dictmax[cell_list[ii]])) if met_cat_dictmax[cell_list[ii]][x]<-1e-3])
                # if any([x<0 for x in met_cat_dictmin[cell_list[ii]]]):
                #     print('problem with ',[labels[x] for x in range(len(met_cat_dictmin[cell_list[ii]])) if met_cat_dictmin[cell_list[ii]][x]<-1e-3],
                #     [met_cat_dictmin[cell_list[ii]][x] for x in range(len(met_cat_dictmin[cell_list[ii]])) if met_cat_dictmin[cell_list[ii]][x]<-1e-3])
                ax.errorbar(np.arange(len(met_cat_dict[cell_list[ii]]))+(ii-1.5)*width,
                    met_cat_dict[cell_list[ii]],
                    yerr=[met_cat_dictmin[cell_list[ii]],met_cat_dictmax[cell_list[ii]]],
                    fmt="o", color="grey",label='_nolegend_',markersize=1,capsize=5)
                ax2.errorbar(np.arange(len(met_cat_dict[cell_list[ii]]))+(ii-1.5)*width,
                    met_cat_dict[cell_list[ii]],
                    yerr=[met_cat_dictmin[cell_list[ii]],met_cat_dictmax[cell_list[ii]]], 
                    fmt="o", color="grey",label='_nolegend_',markersize=1,capsize=5)
        else:
            plt.bar(np.arange(len(met_cat_dict[cell_list[ii]]))+(ii-1.5)*width,
                    met_cat_dict[cell_list[ii]], width=width)
            if rangeplot:
                plt.errorbar(np.arange(len(met_cat_dict[cell_list[ii]]))+(ii-1.5)*width,
                    met_cat_dict[cell_list[ii]],
                             yerr=[met_cat_dictmin[cell_list[ii]],met_cat_dictmax[cell_list[ii]]], 
                             fmt="o", color="black",label='_nolegend_')
    if ybreak:
        ax.legend(cell_list)
    else:
        plt.legend(cell_list)
    if ymax and not ybreak:
        plt.ylim([0,ymax])
    # plt.savefig('RegrEx_plots/'+'1_transfer_large_day.png', bbox_inches='tight')
    plt.ylabel(r'flux ($\mu$mol/m$^2$/s)')
    plt.show()
    return met_cat


def plot_pwy_btwn_cells(cons_model, csol,pwy_dict=[],ss=['glycolysis'],exact=0,thresh=0,plttitle='',phases=['_l','_d'],percell=0):
#     aa_syn_dict={'ARGI':'ARGSYN-PWY','ARGII':'ARGSYNBSUB-PWY','ARGIII':'PWY-5154','ARGIV':'PWY-7400',
#              'CYST':'CYSTSYN-PWY','GLN':'GLNSYN-PWY','GLY':'GLYSYN-THR-PWY','HIS':'HISTSYN-PWY',
#              'ILEI':'ILEUSYN-PWY','ILEII':'PWY-5101','ILEIII':'PWY-5103','ILEIV':'PWY-5104','ILEV':'PWY-5108'
#              'METI':'HOMOSER-METSYN-PWY','METII':'PWY-702','METIII':'HSERMETANA-PWY','METIV':'PWY-7977'
#              'THR':'HOMOSER-THRESYN-PWY','TRP':'TRPSYN-PWY','CIT':'CITRULBIO-PWY',
#              'LEU':'LEUSYN-PWY','LYSVI':'PWY-5097',
#              'PHEI':'PHESYN','PHEII':'PWY-3462','PHEIII':'PWY-7432'
#              'PROIII':'PWY-3341','PROI':'PROSYN-PWY','SERI':'SERSYN-PWY'
#              'VAL':'VALSYN-PWY'}
#     aa_deg_dict={'ALAI':'ALADEG-PWY','ALAII':'ALACAT2-PWY','ALAIII':'ALANINE-DEG3-PWY',
#                  'ALAIV':'PWY1-2','ALAV':'PWY-8189','ALAVI':'PWY-8188',
#                  'ASPI':'ASPARTATE-DEG1-PWY'
#                  'ARGVI':'ARG-PRO-PWY','ARGI':'ARGASEDEG-PWY','ARGXIII':'PWY-8187',
#                  'GLNI':'GLUTAMINDEG-PWY'
#                  'ILEI':'ILEUDEG-PWY','LEUI':'LEU-DEG2-PWY'
#                  'METI':'METHIONINE-DEG1-PWY','METII':'PWY-701','METIII':'PWY-5082'
#                  'THRIV':'PWY-5436','VALI':'VALDEG-PWY'}
#     # glycolysisI from glucose-6P, II from fructose-6P, II from glucose
#     other_pwys={'ethanolfermentationI':'PWY-5480','ethanolfermentationII':'PWY-5486','ethanolfermentationIII':'PWY-6587',
#                'lactatefermentation(S)':'PWY-5481','lactatefermentation(R)':'PWY-8274','heterolacticfermentation':'P122-PWY',
#                'glycolysisI':'GLYCOLYSIS','glycolysisII':'PWY-5484','glycolysisIII':'ANAGLYCOLYSIS',#'glycolysisIV':'PWY-1042','glycolysisV':'P341-PWY'
#                'TCA':'PWY-5690','Calvin-Benson-BasshamCycle':'CALVIN-PWY','RubiscoShunt':'PWY-5723',
#               'NonoxidativePPPI':'NONOXIPENT-PWY','NonoxidativePPPI':'PWY-8178','OPPP':'OXIDATIVEPENT-PWY',
#               'photosynthesislightreactions':'PWY-101','photorespiration':'PWY-181'}
    if not pwy_dict:
        pwy_dict = {'ethanolfermentationI':'PWY-5480','ethanolfermentationII':'PWY-5486','ethanolfermentationIII':'PWY-6587',
               'lactatefermentation(S)':'PWY-5481','lactatefermentation(R)':'PWY-8274','heterolacticfermentation':'P122-PWY',
               'glycolysisI':'GLYCOLYSIS','glycolysisII':'PWY-5484','glycolysisIII':'ANAGLYCOLYSIS',#'glycolysisIV':'PWY-1042','glycolysisV':'P341-PWY'
               'TCA':'PWY-5690','Calvin-Benson-BasshamCycle':'CALVIN-PWY','RubiscoShunt':'PWY-5723',
              'NonoxidativePPPI':'NONOXIPENT-PWY','NonoxidativePPPI':'PWY-8178','OPPP':'OXIDATIVEPENT-PWY',
              'photosynthesislightreactions':'PWY-101','photorespiration':'PWY-181'}
    if not plttitle:
        plttitle=ss[0]
    if exact:
        ss_keys = [x for x in pwy_dict.keys() if any([y == x for y in ss])]
    else:
        ss_keys = [x for x in pwy_dict.keys() if any([y in x for y in ss])]
    ss = []
    for ii in ss_keys:
        ss+=[pwy_dict[ii]]
    ss_rxns=[x.id for x in cons_model.reactions if (x.notes and 'PWY' in x.notes.keys() and any([y in str(x.notes['PWY'])for y in ss]))]
    mc_rxns=[x for x in ss_rxns if any(['_MC'+y in x for y in phases])]
    cc_rxns=[x.replace('_MC_','_CC_') if x.replace('_MC_','_CC_') in ss_rxns else '' for x in mc_rxns]
    se_rxns=[x.replace('_MC_','_SE_') if x.replace('_MC_','_SE_') in ss_rxns else '' for x in mc_rxns]
    pse_rxns=[x.replace('_MC_','_pSE_') if x.replace('_MC_','_pSE_') in ss_rxns else False for x in mc_rxns]
    datamc=[csol[x] if x else 0 for x in mc_rxns]
    datacc=[csol[x] if x else 0 for x in cc_rxns]
    datase=[csol[x] if x else 0 for x in se_rxns]
    datapse=[csol[x] if x else 0 for x in pse_rxns]
    if percell:
        datacc = [20*x for x in datacc]
        datase = [4*x for x in datase]
        datapse = [20*x/6 for x in datapse]
    no_flux=[]
    for ii in range(len(datamc)):
        if abs(datamc[ii])<=thresh and abs(datacc[ii])<=thresh and abs(datase[ii])<=thresh and abs(datapse[ii])<=thresh:
            no_flux.append(ii)
    mc_rxns=[mc_rxns[x] for x in range(len(mc_rxns)) if x not in no_flux]
    cc_rxns=[cc_rxns[x] for x in range(len(cc_rxns)) if x not in no_flux]
    se_rxns=[se_rxns[x] for x in range(len(se_rxns)) if x not in no_flux]
    pse_rxns=[pse_rxns[x] for x in range(len(pse_rxns)) if x not in no_flux]
    datamc=[datamc[x] for x in range(len(datamc)) if x not in no_flux]
    datacc=[datacc[x] for x in range(len(datacc)) if x not in no_flux]
    datase=[datase[x] for x in range(len(datase)) if x not in no_flux]
    datapse=[datapse[x] for x in range(len(datapse)) if x not in no_flux]
    labels=mc_rxns.copy()
    for x in range(len(labels)):
        temp=labels[x].split('_')
        labels[x]='_'.join(temp[:-2]+[temp[-1]])
    width=0.2
    plt.figure(figsize=(len(labels)*6.4/20,4.8))
    plt.xticks(range(len(datamc)), labels,rotation='vertical')
    plt.bar(np.arange(len(datamc)), datamc, width=width)
    plt.bar(np.arange(len(datacc))+ width, datacc, width=width)
    plt.bar(np.arange(len(datase))+ 2*width, datase, width=width)
    plt.bar(np.arange(len(datapse))+ 3*width, datapse, width=width)
    plt.legend(['MC','CC','SE','pSE'])
    plt.title(plttitle)
    # plt.savefig('RegrEx_plots/'+'1_transfer_large_day.png', bbox_inches='tight')
    plt.show()

def plot_met_node_btwn_cells(cons_model, csol,met='PYRUVATE',thresh=0,threshmax=1000,plttitle=''):
    if not plttitle:
        plttitle=met
    met_rxns=[rxn.id for rxn in cons_model.reactions if (any([met in x.id for x in rxn.reactants]))^(any([met in x.id for x in rxn.products]))]
    mc_rxns=[x for x in met_rxns if '_MC_l' in x]
    cc_rxns=[x.replace('_MC_l','_CC_l') if x.replace('_MC_l','_CC_l') in met_rxns else '' for x in mc_rxns]
    se_rxns=[x.replace('_MC_l','_SE_l') if x.replace('_MC_l','_SE_l') in met_rxns else '' for x in mc_rxns]
    pse_rxns=[x.replace('_MC_l','_pSE_l') if x.replace('_MC_l','_pSE_l') in met_rxns else '' for x in mc_rxns]
    datamc=[csol[x] if x else 0 for x in mc_rxns]
    datacc=[csol[x] if x else 0 for x in cc_rxns]
    datase=[csol[x] if x else 0 for x in se_rxns]
    datapse=[csol[x] if x else 0 for x in pse_rxns]
    no_flux=[]
    for ii in range(len(datamc)):
        if abs(datamc[ii])<=thresh and abs(datacc[ii])<=thresh and abs(datase[ii])<=thresh and abs(datapse[ii])<=thresh:
            no_flux.append(ii)
    mc_rxns=[mc_rxns[x] for x in range(len(mc_rxns)) if x not in no_flux]
    cc_rxns=[cc_rxns[x] for x in range(len(cc_rxns)) if x not in no_flux]
    se_rxns=[se_rxns[x] for x in range(len(se_rxns)) if x not in no_flux]
    pse_rxns=[pse_rxns[x] for x in range(len(pse_rxns)) if x not in no_flux]
    datamc=[datamc[x] for x in range(len(datamc)) if x not in no_flux]
    datacc=[datacc[x] for x in range(len(datacc)) if x not in no_flux]
    datase=[datase[x] for x in range(len(datase)) if x not in no_flux]
    datapse=[datapse[x] for x in range(len(datapse)) if x not in no_flux]
    labels=mc_rxns.copy()
    for x in range(len(labels)):
        temp=labels[x].split('_')
        labels[x]='_'.join(temp[:-2]+[temp[-1]])
    width=0.2
    plt.xticks(range(len(datamc)), labels,rotation='vertical')
    plt.bar(np.arange(len(datamc)), datamc, width=width)
    plt.bar(np.arange(len(datacc))+ width, datacc, width=width)
    plt.bar(np.arange(len(datase))+ 2*width, datase, width=width)
    plt.bar(np.arange(len(datapse))+ 3*width, datapse, width=width)
    plt.legend(['MC','CC','SE','pSE'])
    plt.title(plttitle)
    # plt.savefig('RegrEx_plots/'+'1_transfer_large_day.png', bbox_inches='tight')
    plt.show()
        

def gaussian(x, a, b, c, d=0):
    return a * math.exp(-(x - b)**2 / (2 * c**2)) + d


    
def mapSinkC(cobra_model,sol,current_met='SUCROSE_c_SE_d',i_max=2,weighting=50,trav_mets=[],gradations=21):
    now = datetime.now().strftime('%Y-%m-%d-%H-%M')
    f=Digraph('Sources:'+current_met+now,filename='Sources/'+current_met+now+'.gv',strict=True)
    f.attr(rankdir='LR',size='8,5')
    f.attr('node', shape='diamond',style='filled',fillcolor='darkgoldenrod1')
    f.node(current_met)
    f.attr('node', shape='circle',style='filled',fillcolor='white')
    ignore_mets=['ATP_','aATP_','ADP_','aADP_','WATER_','PROTON_','NADPH_','NADH_','NADP_','NAD_','ATPase_NADPHoxidase_constraint_','Pi_','aPi_','Protein_processing_cost','Protein_polymerisation_cost','Protein_translocation_cost','CARBON_CC','CARBON_MC','UBIQUIN']
    ignore_rxns=[]
    abs_max=max(abs(min(sol)),max(sol))
    abs_min=min([abs(sol[x]) for x in sol.keys()])
    gradation = np.logspace(abs_min,abs_max,gradations)
    # create rgb from green to white for pos and red to white for neg
    for x in range(gradations):
        r = int(gaussian(x, 158.8242, 201, 87.0739) + gaussian(x, 158.8242, 402, 87.0739))
        g = int(gaussian(x, 129.9851, 157.7571, 108.0298) + gaussian(x, 200.6831, 399.4535, 143.6828))
        b = int(gaussian(x, 231.3135, 206.4774, 201.5447) + gaussian(x, 17.1017, 395.8819, 39.3148))
    
    trav_mets=mapSinkRecC(f,current_met,cobra_model,sol,0,ignore_mets,ignore_rxns,i_max,weighting)
    print(trav_mets)
    f.view()
    return f

def mapSinkRecC(f,current_met,cobra_model,sol,i=0,trav_mets=[],trav_rxns=[],i_max=2,weight_mult=50):
    if i>i_max-1:
        return ''
    i=i+1
    ig_rxns=[]
    colours=['azure3','royalblue','deeppink4']
    # weight_mult=100
    for rxn in cobra_model.reactions:
        if cobra_model.metabolites.get_by_id(current_met) in rxn.metabolites:
            if sol[rxn.id]*rxn.get_coefficient(current_met)<-1e-5:
#                 print(str(i)+': '+rxn.id)
                if rxn.get_coefficient(current_met)<0:
                    for met in rxn.products:
                        if any([x in rxn.id for x in ['_pc_','_mc_','_xc_','_vc_']]) and met.id.split('_')[-3]==current_met.split('_')[-3]:
                            continue
                        if not any(x in met.id for x in trav_mets):
                            if not any(x in rxn.id for x in trav_rxns):
                                this_width=min(weight_mult*abs(sol[rxn.id]*rxn.get_coefficient(current_met)),60)
                                this_edge=str(max(this_width,.5))
                                this_arrow=str(min(max(this_width/20,0.5),10))
                                if 'dielTransfer' in rxn.id:
                                    this_col=1
                                elif any(x in rxn.id for x in ['_CC_SE_','_SE_pSE_','_ec_']):
                                    this_col=2
                                else:
                                    this_col=0
                                # f.attr(size=str((sol[rxn.id]*rxn.get_coefficient(current_met))/50))
                                f.edge(current_met,met.id,label=rxn.id+': '+str(rxn.get_coefficient(current_met))+'*'+str(sol[rxn.id])+'\n'+rxn.reaction,arrowsize=this_arrow,penwidth=this_edge,color=colours[this_col])
                                # f.edge_attr.update(arrowhead='vee', arrowsize=this_width,penwidth=this_width)
#                                 print(this_edge)
                                ig_rxns=trav_rxns.copy()
                                ig_rxns.append(rxn.id)
                                temp=mapSinkRecC(f,met.id,cobra_model,sol,i,trav_mets,ig_rxns,i_max,weight_mult)
                else:
                    for met in rxn.reactants:
                        if any([x in rxn.id for x in ['_pc_','_mc_','_xc_','_vc_']]) and met.id.split('_')[-3]==current_met.split('_')[-3]:
                            continue
                        if not any(x in met.id for x in trav_mets):
                            if not any(x in rxn.id for x in trav_rxns):
                                this_width=min(weight_mult*abs(sol[rxn.id]*rxn.get_coefficient(current_met)),60)
                                this_edge=str(max(this_width,.5))
                                this_arrow=str(min(max(this_width/20,0.5),10))
                                if 'dielTransfer' in rxn.id:
                                    this_col=1
                                elif any(x in rxn.id for x in ['_CC_SE_','_SE_pSE_','_ec_']):
                                    this_col=2
                                else:
                                    this_col=0
                                f.edge(current_met,met.id,label=rxn.id+': '+str(rxn.get_coefficient(current_met))+'*'+str(sol[rxn.id])+'\n'+rxn.reaction,arrowsize=this_arrow,penwidth=this_edge,color=colours[this_col])
                                # f.edge_attr.update(arrowhead='vee', arrowsize=this_width,penwidth=this_width)
#                                 print(this_edge)
                                ig_rxns=trav_rxns.copy()
                                ig_rxns.append(rxn.id)
                                temp=mapSinkRecC(f,met.id,cobra_model,sol,i,trav_mets,ig_rxns,i_max,weight_mult)
    return ig_rxns

def plot_proton_exchange(stem_model,psol,cell_type='CC_l',thresh=0,met='PROTON_e_'):
    if met+cell_type in stem_model.metabolites:
        He=met+cell_type
    elif met+cell_type[:-2] in stem_model.metabolites:
        He=met+cell_type[:-2]
    elif met+cell_type[-1] in stem_model.metabolites:
        He=met+cell_type[-1]
    else:
        print('He = '+met+cell_type[-1]+'?')
    labels=[rxn.id for rxn in stem_model.reactions if (any(He in x.id for x in rxn.metabolites)and any('_'+cell_type.split('_')[0]+'_' in x.id for x in rxn.metabolites) and abs(psol[rxn.id])>0)]
    labels = sorted(labels, key=lambda x: x[-1:-5:-1],reverse=True)
    print(cell_type.split('_')[0],He,labels)
    data1=[-psol[rxn]*stem_model.reactions.get_by_id(rxn).get_coefficient(He) for rxn in labels]
    print('Net H+ transport: '+str(sum(data1)))
    labels = ['_'.join(x.split('_')[:-3]+[x.split('_')[-1]]) for x in labels]
    width=0.6
    plt.xticks(range(len(data1)), labels,rotation='vertical')
    plt.bar(np.arange(len(data1)), data1, width=width)
    # plt.bar(np.arange(len(data2))+ width, data2, width=width)
    plt.legend(['apo -> '+cell_type])
    plt.ylabel('flux')
#     plt.ylim(-80,100)
    # plt.savefig('RegrEx_plots/'+'1_transfer_large_day.png', bbox_inches='tight')
    plt.show()
    
def compareModels(model1,model2):
    m1 = io.read_sbml_model(model1)
    m2 = io.read_sbml_model(model2)
    rxn_components=['bounds','compartments','id','reaction']#,'notes']
    for rxn1 in m1.reactions:
        if rxn1.id in m2.reactions:
            rxn2=m2.reactions.get_by_id(rxn1.id)
            for bit in rxn_components:
                if(getattr(rxn1,bit)!=getattr(rxn2,bit)):
                    print(rxn1.id,bit,getattr(rxn1,bit),getattr(rxn2,bit),'\n')
    print('Reactions not in m2:',[rxn.id for rxn in m1.reactions if rxn.id not in m2.reactions])
    print('Reactions not in m1:',[rxn.id for rxn in m2.reactions if rxn.id not in m1.reactions])
    
def constraint_scan(model,cst_reaction='diel_biomass',save_file_name='constraint_scan',start_val=None,step_size=0.05,number_of_steps=5,step_dir=-1): 
    ub=model.reactions.get_by_id(cst_reaction).upper_bound
    lb=model.reactions.get_by_id(cst_reaction).lower_bound
    fva_sol=flux_variability_analysis(model)
    psol=pfba(model)
    data={'reaction_id':[x.id for x in model.reactions],
          'cell':['_'.join(x.id.split('_')[-2:]) if ('_l' in x.id or '_d' in x.id) else '' for x in model.reactions],
          'subsystem':[x.notes['SUBSYSTEM']  if 'SUBSYSTEM' in str(x.notes) else '' for x in model.reactions],
          'reaction':[x.reaction for x in model.reactions],
         'fva_min':[fva_sol.minimum[x.id] for x in model.reactions],
         'fva_max':[fva_sol.maximum[x.id] for x in model.reactions],
          '1':[psol[x.id] for x in model.reactions]}
    if start_val==None:
        start_val=psol[cst_reaction]
    for ii in range(1,number_of_steps):
        rxn_val = (1+step_dir*step_size*ii)
        model.reactions.get_by_id(cst_reaction).upper_bound=rxn_val*start_val
        model.reactions.get_by_id(cst_reaction).lower_bound=rxn_val*start_val
        this_sol=pfba(model)
        data[rxn_val]=[this_sol[x.id] for x in model.reactions]
    model.reactions.get_by_id(cst_reaction).upper_bound=ub
    model.reactions.get_by_id(cst_reaction).lower_bound=lb
    df=pd.DataFrame(data)

    corr=list()
    m=[]
    start_i=[x for x in range(len(df.columns)) if df.columns[x]=='1'][0]
    x = np.array(df.columns[start_i:(start_i+number_of_steps-1)]).reshape((-1,1))
    for ii in range(len(df)):
        y = np.array(df.loc[ii][start_i:(start_i+number_of_steps-1)])
        if len(set(y))<=1:
            m.append(0)
            if y[2]!=0:
                corr.append(2)
            else:
                corr.append(0)
        else:
    #         print(x,y)
            m.append(LinearRegression().fit(x,y).coef_)
            corr.append(LinearRegression().fit(x,y).score(x,y))
#             corr.append(max(LinearRegression().fit(x,y).score(x,y),2))
    print(corr[:10])
    df['correlation']=corr
    df['slope']=m

    now = datetime.now().strftime('%Y_%m_%d')
    df.to_csv(save_file_name+now+'.csv')

def showReactions(stem_model,sols,met,flux=None,transcripts=None): 
    printdict={}
    if str(type(sols))=="<class 'dict'>":
        if (str(type(list(sols.values())[0]))=="<class 'dict'>" or type(list(sols.values())[0])==cobra.core.solution.Solution):
            for rxn in stem_model.metabolites.get_by_id(met).reactions:
                if any([sols[key][rxn.id]!=0 for key in list(sols.keys())]) or flux:
                    print([str(key)+': '+str(sols[key][rxn.id]*rxn.get_coefficient(met)) for key in list(sols.keys())],rxn.id)
                    print('   ',rxn.reaction)
        else:
            for rxn in stem_model.metabolites.get_by_id(met).reactions:
                if sols[rxn.id]!=0 or flux:
                    printdict[rxn.id]=sols[rxn.id]*rxn.get_coefficient(met)
                    # print(sols[rxn.id]*rxn.get_coefficient(met),rxn.id)
                    # print('   ',rxn.reaction)
            newdict=sorted(printdict.items(),key=lambda x:abs(x[1]), reverse=True)
            for x in newdict:
                print(x[1],x[0])
                if transcripts and str(transcripts[x[0]])!='nan':
                    if '_MC_' in x[0]:
                        print(transcripts[x[0]],20*transcripts[x[0].replace('_MC_','_CC_')])
                    elif x[0].replace(x[0].split('_')[-2],'MC') in transcripts.keys():
                        print(20*transcripts[x[0]],transcripts[x[0].replace(x[0].split('_')[-2],'MC')])
                    else:
                        print(transcripts[x[0]])
                print('   ',stem_model.reactions.get_by_id(x[0]).reaction)
    elif type(sols)==cobra.core.solution.Solution:
        print('opt21')
        for rxn in stem_model.metabolites.get_by_id(met).reactions:
            if sols[rxn.id]!=0 or flux:
                print(sols[rxn.id]*rxn.get_coefficient(met),rxn.id)
                print('   ',rxn.reaction)


def showReactionsWith(stem_model,sols,rxn_frag,flux=None,transcripts=None): 
    printdict={}
    if str(type(sols))=="<class 'dict'>":
        if (str(type(list(sols.values())[0]))=="<class 'dict'>" or type(list(sols.values())[0])==cobra.core.solution.Solution):
            for rxn in [x for x in stem_model.reactions if rxn_frag in x.id]:
                if any([sols[key][rxn.id]!=0 for key in list(sols.keys())]) or flux:
                    print([str(key)+': '+str(sols[key][rxn.id]) for key in list(sols.keys())],rxn.id)
                    print('   ',rxn.reaction)
        else:
            for rxn in [x for x in list(sols.keys()) if rxn_frag in x]:
                if sols[rxn]!=0 or flux:
                    printdict[rxn]=sols[rxn]
                    # print(sols[rxn.id],rxn.id)
                    # print('   ',rxn.reaction)
            newdict=sorted(printdict.items(),key=lambda x:abs(x[1]), reverse=True)
            for x in newdict:
                print(x[1],x[0])
                print('   ',stem_model.reactions.get_by_id(x[0]).reaction)
    elif type(sols)==cobra.core.solution.Solution:
        for rxn in [x for x in list(sols.fluxes.keys()) if rxn_frag in x]:
            if sols[rxn]!=0 or flux:
                print(sols[rxn],rxn)
                if transcripts and str(transcripts[rxn])!='nan':
                    if '_MC_' in rxn:
                        print(transcripts[rxn],20*transcripts[rxn.replace('_MC_','_CC_')])
                    elif rxn.replace(rxn.split('_')[-2],'MC') in transcripts.keys():
                        print(20*transcripts[rxn],transcripts[rxn.replace(rxn.split('_')[-2],'MC')])
                    else:
                        print(transcripts[rxn])
                print('   ',stem_model.reactions.get_by_id(rxn).reaction)
    elif type(sols)==list:
        for rxn in [x for x in list(sols.keys()) if rxn_frag in x]:
            if any([x[rxn.id]!=0 for x in sols]) or flux:
                print([x[rxn.id] for x in sols],rxn.id)
                print('   ',rxn.reaction)