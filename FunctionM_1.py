# -*- coding: utf-8 -*-
"""
The distribution and transformation of PFAS in wastewater treatment plant (WWTP), DT_PFAS_investigator.version 0

Function definition module file 1.
Define functions used for PFAS mass flow calculation and the investigation of mechanismes driving PFAS concentration changes in WWTP.

ENV: BeeWare

Created on Wed Sep 17 09:35:06 2025

@author: Xiangui Huang, Ben-Gurion University, Israel.
"""


import pandas as pd
from importlib import resources #used in native app.



#-------------------------------------------------------------------------------------------------------------------
#Here, some global variables are defined as assigned as default value,
# they can be used in any fucntion in this module also can be imported into other script and changed their values.
# In this app development, they will be imported into my app script, and will be changed through GUI by user.
#---------------------------------------------------------------------------------------------------------------
consider_SD_PFASmass_CalC=True

















#---------------------
#define functions for reading basic information from file and other database.
#---------------------
def system_infor_reader(inputdata_file):
    """
    inputdata_file-the dictionary of the input data file.
    ----------
    inputdata_file : csv file
    It contains two type of information: the PFAS concentration in each compartment; the basic information of each compartment
    and the connection among compartment.

    Returns
    -------
    A dictionary object, containing the basic information of the system
    {'label': , 'temp': , "pascal": , 'condition':, 'flux_unit':, 'input_source':, 'input_flux':,
    'output_source':, 'output_flux':}.

    """
    input_data_df = pd.read_csv(inputdata_file, header=0, index_col=0)
    properties_df=input_data_df[input_data_df['class']=='properties']
    properties_df.set_index('parameter', inplace=True)
    
    column_list=properties_df.columns.tolist()
    compt_list=[compt_i for compt_i in column_list if 'compt_' in compt_i]
    
    system_infor={}
    for compt_i in compt_list:
        label=properties_df.loc['label', compt_i]
        temp=properties_df.loc['temp', compt_i]
        pascal=properties_df.loc['pascal', compt_i]
        condition=properties_df.loc['condition', compt_i]
        flux_unit=properties_df.loc['flux_unit', compt_i]
        
        index_list=properties_df.index.tolist()
        input_list=[input_i for input_i in index_list if 'receive_from' in input_i]
        input_flux=[input_i for input_i in index_list if 'receive_flux' in input_i]
        input_sources=[properties_df.loc[input_i, compt_i] for input_i in input_list]
        source_flux=[properties_df.loc[flux_i, compt_i] for flux_i in input_flux]
        
        if 'NO' in input_sources: input_sources.remove('NO')
        if 'NO' in source_flux: source_flux.remove('NO')
        
        output_list=[output_i for output_i in index_list if 'discharge_to' in output_i]
        output_flux=[output_i for output_i in index_list if 'discharge_flux' in output_i]
        output_sources=[properties_df.loc[output_i, compt_i] for output_i in output_list]
        output_flux=[properties_df.loc[flux_i, compt_i] for flux_i in output_flux]
        
        if 'NO' in output_sources: output_sources.remove('NO')
        if 'NO' in output_flux: output_flux.remove('NO')
        
        system_infor[compt_i]={'label': label, 'temp': temp, "pascal": pascal, 'condition':condition, 'flux_unit':flux_unit, 'input_source':input_sources, 'input_flux':source_flux,
                               'output_source':output_sources, 'output_flux':output_flux}
        
    
    return system_infor
        







#Get the PFAS_trans_lib file path, which will be used in Step 3.

PFAS_Trans_lib_file = resources.files("pfas_dt_investigator.resources").joinpath("Pre_Trans_Lib.csv") #used in native app

#PFAS_Trans_lib_file = "C:\Aphd_PFAS in WWTP\DataAnalysis\DT_PFAS_Investigator\Pre_Trans_Lib.csv"  #used in script

def reactant_product_lookup(PFAS, condition, PFAS_Trans_lib_file=PFAS_Trans_lib_file):
    """
    This function read the lib file path, the PFAS species, and transformation condition,
    and to extract the possibile precursors and products.
    PFAS: a PFAS species
    
    condition: the condition interms of oxygen or other reaction condition, check the lib file for specific name.
    
    PFAS_Trans_lib_file: the file path of PFAS precursors transformation lib.
    Notes: input file is the transformation lib, which is set as the defult lib file. User can modify or define their own lib file,
    and employed it to this function.
    
    return: a list contain dict elements, each dict contain the precursor and product of the received PFAS species.
    """
    
    with PFAS_Trans_lib_file.open("r", encoding="utf-8") as f: Pre_Trans_lib_v0_df = pd.read_csv(f, header=0, index_col=0) #used in native app
    #Pre_Trans_lib_v0_df = pd.read_csv(PFAS_Trans_lib_file, header=0, index_col=0) #used in script.
    
    
    Pre_Trans_lib_df=Pre_Trans_lib_v0_df[Pre_Trans_lib_v0_df['condition']==condition]
    tree_No_list=Pre_Trans_lib_df['tree_No'].unique()
    
    extract=[]
    for tree_i in tree_No_list:
        tree_i_df=Pre_Trans_lib_df[Pre_Trans_lib_df['tree_No']==tree_i]
        index_list=tree_i_df.index.tolist()
        
        precursors=[]
        products=[]
        
        isin_reactant = tree_i_df['reactant'].isin([PFAS]).any()
        isin_product=tree_i_df['product'].isin([PFAS]).any()
        
        #case 1: this PFAS species is a reactant.
        if isin_reactant and not isin_product:
            position='reactant'
            product_i=PFAS
            for index_i in index_list:
                if tree_i_df.loc[index_i, 'reactant'] == product_i or tree_i_df.loc[index_i, 'reactant'] in products:
                    product_i=tree_i_df.loc[index_i, 'product']
                    products=products+[product_i]
            #save the results into a dir.
            PFAS_dir={'tree_i':tree_i, 'position':position, 'precursor':precursors, 'product':products}
            extract=extract+[PFAS_dir]
        
            
        #case 2: this PFAS species is a intermediate.
        elif isin_reactant and isin_product:
            position='intermediate'
            product_i=PFAS
            precursor_i=PFAS
            
            #find the products
            for index_i in index_list:
                if tree_i_df.loc[index_i, 'reactant'] == product_i or tree_i_df.loc[index_i, 'reactant'] in products:
                    product_i=tree_i_df.loc[index_i, 'product']
                    products=products+[product_i]
            
            #find the precursors
            for index_i in  list(reversed(index_list)): #reverse the order, look up the precursor from bottom to top in the lib.
                #print(index_i)
                if tree_i_df.loc[index_i, 'product'] in precursors or tree_i_df.loc[index_i, 'product'] == PFAS:
                    precursor_i=tree_i_df.loc[index_i, 'reactant']
                    precursors=precursors+[precursor_i]
            
            PFAS_dir={'tree_i':tree_i, 'position':position, 'precursor':precursors, 'product':products}
            extract=extract+[PFAS_dir]
                    
        
        #case 3: this PFAS is a terminal product.
        elif not isin_reactant and isin_product:
            position='terminal'
            precursor_i=PFAS
            for index_i in  list(reversed(index_list)): #reverse the order, look up the precursor from bottom to top in the lib.
                if tree_i_df.loc[index_i, 'product'] in precursors or tree_i_df.loc[index_i, 'product'] == PFAS:
                    precursor_i=tree_i_df.loc[index_i, 'reactant']
                    precursors=precursors+[precursor_i]
            
            PFAS_dir={'tree_i':tree_i, 'position':position, 'precursor':precursors, 'product':products}
            extract=extract+[PFAS_dir]
        
    #in case the required PFAS species is not included in the PFAS trans lib.
    if extract ==[]: extract=['not_found']
                   
    return extract







def compt_label_conveter(compt_label, inputdata_file):
    """
    This is function is a conveter, providing compt_number it will give you its label, or vice versa.

    Parameters
    ----------
    compt_label : text
        a compt number, compt_i or a label of a compt.
    inputdata_file : a csv file
        It contains two type of information: the PFAS concentration in each compartment; the basic information of each compartment
        and the connection among compartment.

    Returns
    -------
    compt numner or a label of the corresponding compt_i.

    """
    system_infor_dir=system_infor_reader(inputdata_file)
    compt_list=list(system_infor_dir.keys())
    label_list=[system_infor_dir[compt_i]['label'] for compt_i in compt_list]
    
    compt_label_dir={}
    for i in range(len(compt_list)):
        compt_label_dir[compt_list[i]]=label_list[i]
        compt_label_dir[label_list[i]]=compt_list[i]
        
    if 'Influent_' in compt_label or 'Effluent_' in compt_label:
        value=compt_label
    else:
        value=compt_label_dir[compt_label]
    
    return value
    



def SD_column(compt_i):
    """
    This function is dedicated to give the SD column name for the provided compt_i column 

    Parameters
    ----------
    compt_i: str object
        a name of compartment or Influent, such as compt_1, Influent_1.

    Returns
    -------
    SD column name, such as SD_Influent(1)

    """
    compt_i_split=compt_i.split("_")
    SD_compt_i=f"SD_{compt_i_split[0]}({compt_i_split[1]})"
    return SD_compt_i





def PFAS_species_reader(inputdata_file):
    """
    Parameters
    ----------
    inputdata_file :a csv file path.
        a file contains the PFAS data for each compartment and the basic information of each compartment
        and the connection among compartment.
        
    Returns
    -------
    A dir contain terminal and precursors's PFAS species.

    """
    input_data_df = pd.read_csv(inputdata_file, header=0, index_col=0)
    terminal_short_species=input_data_df[input_data_df['class']=='terminal_SPFAS']
    terminal_long_species=input_data_df[input_data_df['class']=='terminal_LPFAS']
    precursor_species=input_data_df[input_data_df['class']=='precursor_PFAS']
    
    terminal_short_list= terminal_short_species['parameter'].tolist()
    terminal_long_list= terminal_long_species['parameter'].tolist()
    precursor_list= precursor_species['parameter'].tolist()
    
    PFAS_species_dir={}
    PFAS_species_dir['terminal short']=terminal_short_list
    PFAS_species_dir['terminal long']=terminal_long_list
    PFAS_species_dir['precursor']=precursor_list
    
    return PFAS_species_dir
    
    
    
       
    
def PFAS_data_reader(inputdata_file, PFAS, compt_i):
    """
    This function is dedicated to calculate the input and output mass of a PFAS species for a given compartment.
    Parameters
    ----------
    inputdata_file: a csv file path.
        a file contains the PFAS data for each compartment and the basic information of each compartment
        and the connection among compartment.
    PFAS : text
        a name of PFAS species.
    compt_i : text
        a number of a compartment in the system.

    Returns
    -------
    a dir object contain the total input and output PFAS mass {'input': sum_input_mass, 'input_dir':input_mass_dir, 'output': sum_output_mass,'output_dir':output_mass_dir, 'unit': ug/d or mg/d}.

    """
    input_data_df = pd.read_csv(inputdata_file, header=0, index_col=0)
    PFAS_df=input_data_df[input_data_df['class'].str.contains('PFAS', na=False)]
    PFAS_df.set_index('parameter', inplace=True)
    
    compt_propert_dir=system_infor_reader(inputdata_file)
    compt_i_dir=compt_propert_dir[compt_i]
    
    sum_PFAS_input=[]
    sum_PFAS_output=[]
    
    PFAS_unit=PFAS_df.loc['PFAS_unit', compt_i]
    flux_unit=compt_i_dir['flux_unit']
    
    if (PFAS_unit == 'ng/L' and flux_unit=='m3/d') or (PFAS_unit == 'ng/kg' and flux_unit=='ton/d') or (PFAS_unit == 'ug/kg' and flux_unit=='ton/d') or (PFAS_unit == 'ug/L' and flux_unit=='m3/d'):
        for i, source_i in enumerate( compt_i_dir['input_source'] ):
            PFAS_conc_input=PFAS_df.loc[PFAS, source_i]
            if PFAS_conc_input == 'ND' or PFAS_conc_input == 'D':
                sum_PFAS_input=sum_PFAS_input+[PFAS_conc_input]
            else:
                input_mass=float(PFAS_conc_input)*float(compt_i_dir['input_flux'][i]) #calculate as ug/d or mg/d
                sum_PFAS_input=sum_PFAS_input+[input_mass]
        
        for i, output_i in enumerate( compt_i_dir['output_source'] ):
            PFAS_conc_output=PFAS_df.loc[PFAS, compt_i]
            
            if PFAS_conc_output == 'ND' or PFAS_conc_output == 'D':
                sum_PFAS_output=sum_PFAS_output+[PFAS_conc_output]
            else:
                output_mass=float(PFAS_conc_output)*float(compt_i_dir['output_flux'][i]) #calculate as ug/d or mg/d
                sum_PFAS_output=sum_PFAS_output+[output_mass]
    else:
        print('Error, the unit between PFAS conc and flux are not consistent.')
    
    
    
    sum_PFAS_mass_flow={}
    input_dir={}
    output_dir={}
    
    for i in range(len(compt_i_dir['input_source'])):
        input_source=compt_i_dir['input_source'][i]
        if 'compt_' in input_source:
            input_source=compt_propert_dir[input_source]['label']
        input_dir[input_source]=sum_PFAS_input[i]
        
    for j in range(len(compt_i_dir['output_source'])):
        output_source=compt_i_dir['output_source'][j]
        if 'compt_' in output_source:
            output_source=compt_propert_dir[output_source]['label']
        output_dir[output_source]=sum_PFAS_output[j]
        
        
    sum_PFAS_mass_flow['input_dir']=input_dir
    sum_PFAS_mass_flow['output_dir']=output_dir
    
    
    
    sum_PFAS_inputmass='ND'
    sum_PFAS_outputmass='ND'
    
    PFAS_mass_input=0
    for i in sum_PFAS_input:
        if i == 'D':
            sum_PFAS_inputmass='D'
        elif isinstance(i, (int, float)):
            PFAS_mass_input +=i
    
    if PFAS_mass_input ==0:
        sum_PFAS_mass_flow['input']= sum_PFAS_inputmass
    else:
        sum_PFAS_mass_flow['input']= PFAS_mass_input
    
    
    PFAS_mass_output=0
    for i in sum_PFAS_output:
        if i == 'D':
            sum_PFAS_outputmass='D'
        elif isinstance(i, (int, float)):
            PFAS_mass_output+=i
    
    if PFAS_mass_output ==0:
        sum_PFAS_mass_flow['output']=sum_PFAS_outputmass
    else:
        sum_PFAS_mass_flow['output']= PFAS_mass_output
    
    
    #define the unit.
    if (PFAS_unit == 'ng/L' and flux_unit=='m3/d') or (PFAS_unit == 'ng/kg' and flux_unit=='ton/d'): unit='ug/d'
    if (PFAS_unit == 'ug/kg' and flux_unit=='ton/d') or (PFAS_unit == 'ug/L' and flux_unit=='m3/d'): unit='mg/d'
    
    sum_PFAS_mass_flow['unit']= unit
    
    
    return sum_PFAS_mass_flow





def PFAS_data_SD_reader(inputdata_file, PFAS, compt_i):
    """
    This is a twin function with PFAS_data_reader(), which dedicated to read and calculate the standard deviation (SD).

    ----------
    inputdata_file: a csv file path.
        a file contains the PFAS data for each compartment and the basic information of each compartment
        and the connection among compartment.
    PFAS : text
        a name of PFAS species.
    compt_i : text
        a number of a compartment in the system.

    Returns
    -------
    a dir object contain the total input and output PFAS mass {'input': sum_input_mass, 'input_dir':input_mass_dir, 'output': sum_output_mass,'output_dir':output_mass_dir, 'unit': ug/d or mg/d}.

    """
    input_data_df = pd.read_csv(inputdata_file, header=0, index_col=0)
    PFAS_df=input_data_df[input_data_df['class'].str.contains('PFAS', na=False)]
    PFAS_df.set_index('parameter', inplace=True)
    
    compt_propert_dir=system_infor_reader(inputdata_file)
    compt_i_dir=compt_propert_dir[compt_i]
    
    sum_PFAS_input=[]
    sum_PFAS_output=[]
    
    PFAS_unit=PFAS_df.loc['PFAS_unit', compt_i]
    flux_unit=compt_i_dir['flux_unit']
    
    if (PFAS_unit == 'ng/L' and flux_unit=='m3/d') or (PFAS_unit == 'ng/kg' and flux_unit=='ton/d') or (PFAS_unit == 'ug/kg' and flux_unit=='ton/d') or (PFAS_unit == 'ug/L' and flux_unit=='m3/d'):
        for i, source_i in enumerate( compt_i_dir['input_source'] ):
            source_i_SD=SD_column(source_i)
            PFAS_conc_input=PFAS_df.loc[PFAS, source_i]
            PFAS_conc_input_SD=PFAS_df.loc[PFAS, source_i_SD]
            
            if PFAS_conc_input_SD == 'ND' or PFAS_conc_input_SD == 'D':
                sum_PFAS_input=sum_PFAS_input+[PFAS_conc_input_SD]
            else:
                input_mass=float(PFAS_conc_input)*float(compt_i_dir['input_flux'][i]) #calculate as ug/d or mg/d
                input_mass_SD=float(PFAS_conc_input_SD)*input_mass/float(PFAS_conc_input)
                sum_PFAS_input=sum_PFAS_input+[input_mass_SD]
        
        for i, output_i in enumerate( compt_i_dir['output_source'] ):
            compt_i_SD=SD_column(compt_i)
            PFAS_conc_output=PFAS_df.loc[PFAS, compt_i]
            PFAS_conc_output_SD=PFAS_df.loc[PFAS, compt_i_SD]
            
            if PFAS_conc_output == 'ND' or PFAS_conc_output == 'D':
                sum_PFAS_output=sum_PFAS_output+[PFAS_conc_output_SD]
            else:
                output_mass=float(PFAS_conc_output)*float(compt_i_dir['output_flux'][i]) #calculate as ug/d or mg/d
                output_mass_SD=float(PFAS_conc_output_SD)*output_mass/float(PFAS_conc_output)
                sum_PFAS_output=sum_PFAS_output+[output_mass_SD]
    else:
        print('Error, the unit between PFAS conc and flux are not consistent.')
    
    
    
    sum_PFAS_mass_flow={}
    input_dir={}
    output_dir={}
    
    for i in range(len(compt_i_dir['input_source'])):
        input_source=compt_i_dir['input_source'][i]
        if 'compt_' in input_source:
            input_source=compt_propert_dir[input_source]['label']
        input_dir[input_source]=sum_PFAS_input[i]
        
    for j in range(len(compt_i_dir['output_source'])):
        output_source=compt_i_dir['output_source'][j]
        if 'compt_' in output_source:
            output_source=compt_propert_dir[output_source]['label']
        output_dir[output_source]=sum_PFAS_output[j]
        
        
    sum_PFAS_mass_flow['input_dir']=input_dir
    sum_PFAS_mass_flow['output_dir']=output_dir
    
    
    
    sum_PFAS_inputmass='ND'
    sum_PFAS_outputmass='ND'
    
    PFAS_mass_input=0
    for i in sum_PFAS_input:
        if i == 'D':
            sum_PFAS_inputmass='D'
        elif isinstance(i, (int, float)):
            PFAS_mass_input +=i**2
    
    if PFAS_mass_input ==0:
        sum_PFAS_mass_flow['input']= sum_PFAS_inputmass
    else:
        sum_PFAS_mass_flow['input']= PFAS_mass_input**0.5
    
    
    PFAS_mass_output=0
    for i in sum_PFAS_output:
        if i == 'D':
            sum_PFAS_outputmass='D'
        elif isinstance(i, (int, float)):
            PFAS_mass_output+=i**2
    
    if PFAS_mass_output ==0:
        sum_PFAS_mass_flow['output']=sum_PFAS_outputmass
    else:
        sum_PFAS_mass_flow['output']= PFAS_mass_output**0.5
    
    
    #define the unit.
    if (PFAS_unit == 'ng/L' and flux_unit=='m3/d') or (PFAS_unit == 'ng/kg' and flux_unit=='ton/d'): unit='ug/d'
    if (PFAS_unit == 'ug/kg' and flux_unit=='ton/d') or (PFAS_unit == 'ug/L' and flux_unit=='m3/d'): unit='mg/d'
    
    sum_PFAS_mass_flow['unit']= unit
    
    
    return sum_PFAS_mass_flow






def thelta_PFAS(sum_input, sum_output, sum_input_SD, sum_output_SD, consider_SD=True):
    """This function is define with the goal to judge the PFAS mass change across a step.
     The input value could be a number, a text indicating the level of PFAS: D and ND.
     The related delta exceed 0.05 is reconginzd as there is PFAS mass change, it will be equilibrium otherwise.
     return: increased, decreased, or equilibrium.
     """
   
    isnumber_input=isinstance(sum_input, (int, float))
    isnumber_output=isinstance(sum_output, (int, float))
    isnumber_input_SD=isinstance(sum_input_SD, (int, float))
    isnumber_output_SD=isinstance(sum_output_SD, (int, float))
    
    thelta='invalid_value'
    consider_SD=consider_SD_PFASmass_CalC
    
    if (isnumber_input and isnumber_output) and (isnumber_input_SD and  isnumber_output_SD):
        if consider_SD:
            sum_input_range=[sum_input-sum_input_SD, sum_input+sum_input_SD]
            sum_output_range=[sum_output-sum_output_SD, sum_output+sum_output_SD]
            
            if (sum_input_range[0] <= sum_output_range[0] and sum_input_range[1] >= sum_output_range[1]) or (sum_input_range[0] >= sum_output_range[0] and sum_input_range[1] <= sum_output_range[1]):
                thelta='equilibrium'
            elif (sum_input_range[0] > sum_output_range[0]) and (sum_input_range[1] > sum_output_range[1]):
                thelta='decreased'
            elif sum_input_range[0] < sum_output_range[0] and sum_input_range[1] < sum_output_range[1]:
                thelta='increased'
        if not consider_SD:
            delta=sum_input - sum_output
            if delta > 0.0: thelta='decreased'
                        
            elif delta < 0.0: thelta = 'increased'
               
            else: thelta='equilibrium'
    

            
    elif isnumber_input and (sum_output=='D' or sum_output=='ND'):
        thelta='decreased'
    
    elif (sum_input=='ND' or sum_input=='D') and isnumber_output:
        thelta='increased'
    
    elif not isnumber_input and not isnumber_output:
        if sum_input=='D' and sum_output=='ND':
            thelta='decreased'
        elif sum_input=='ND' and sum_output=='D':
            thelta='increased'
        elif sum_input == sum_output:
            thelta='equilibrium'
    
    return thelta












#------------------------------------
# defing functions for data analysis.
#--------------------------------------

def DT_increase_analyzer(PFAS_species, PFAS_species_dir, PFASdata_file, compt_i, compt_condition):
    """
    This function performe the data analysis, which return the possibile reason why the PFAS_species increase across a compartment.

    Parameters
    ----------
    PFAS_species : text
        The name of a PFAS species.
    PFAS_species_dir : dirctionary object
        A dirctionary contains all the terminal species and precursors species.
    PFASdata_file: a csv file path
        The file path to a file contains the PFAS data for each compartment and the basic information of each compartment.
    compt_i: the number of the current compartment.
    compt_condition: the condition of the current compartment.

    Returns
    -------
    Return a list, which contain the reason why the PFAS species change in the compartment,
    in which key words: ND means it was found in Pre_Trans_lib but not detected; NF_lib denotes it was not found in Pre_Trans_lib.

    """
    reason_list=['Plus > Minus:', 'Minus: adsorbed onto sludge']
    #case 1: if the species is a terminal product.
    if PFAS_species in PFAS_species_dir['terminal short']+PFAS_species_dir['terminal long']:
        #consider the ultra short-chain species, no adsorption/desorption.
        if PFAS_species not in ['TFA', 'PFPrA', 'PFEtS', 'PFPrS']:
            if PFAS_species in PFAS_species_dir['terminal short']: #the short-chain species can be desorbed from sludge in all the wastewater treatment steps.
                reason_list=reason_list+['Plus: desorbed from sludge']
            
            #consider the long-chain, which can be desorbed significantly in the Reactor.
            if PFAS_species in PFAS_species_dir['terminal long']:
                #compt_value=compt_label_conveter(compt_i, PFASdata_file)
                #if compt_value == 'Reactor':
                reason_list=reason_list+['Plus: desorbed from sludge']
        else:
            reason_list +=['Minus: volatilization']
        
        
        precursor_terminal=reactant_product_lookup(PFAS_species, compt_condition)
        if precursor_terminal != ['not_found']:
            for dir_i in precursor_terminal:
                precursor_list= dir_i['precursor']
                if precursor_list != []:
                    tree_i = dir_i['tree_i']
                    for precursor_i in precursor_list:
                        if precursor_i in PFAS_species_dir['precursor']:
                            precursor_i_conc=PFAS_data_reader(PFASdata_file, precursor_i, compt_i)
                            precursor_i_SD=PFAS_data_SD_reader(PFASdata_file, precursor_i, compt_i)
                            precursor_changes = thelta_PFAS(precursor_i_conc['input'], precursor_i_conc['output'], precursor_i_SD['input'], precursor_i_SD['output'])
                            
                            if precursor_changes == 'decreased':
                                reason_list=reason_list+['Plus: Tree_{0}_{1}_trans'.format(tree_i, precursor_i)]
           
            if not any('_trans' in s for s in reason_list):
               reason_list=reason_list+['precursor_ND']
                            
        
        elif precursor_terminal == ['not_found']:
            reason_list=reason_list+['precursor_NF_lib']
        
    
        
    
    
    #case 2: if the species is a parent precursor or a mediate precursor.
    elif PFAS_species in PFAS_species_dir['precursor']:
        reason_list=reason_list+['Plus: desorbed from sludge']
        
        precursor_terminal=reactant_product_lookup(PFAS_species, compt_condition)
        if precursor_terminal != ['not_found']:
            found_precursor=False
            found_product=False
            
            for dir_i in precursor_terminal:
                tree_i = dir_i['tree_i']
                precursor_list= dir_i['precursor']
                product_list=dir_i['product']
                
                if precursor_list !=[]:
                    found_precursor=True
                    #To check if there is initial precursor transformation.
                    for precursor_i in precursor_list:
                        if precursor_i in PFAS_species_dir['precursor']:
                            precursor_i_conc=PFAS_data_reader(PFASdata_file, precursor_i, compt_i)
                            precursor_i_SD=PFAS_data_SD_reader(PFASdata_file, precursor_i, compt_i)
                            precursor_changes = thelta_PFAS(precursor_i_conc['input'], precursor_i_conc['output'], precursor_i_SD['input'], precursor_i_SD['output'])
                            if precursor_changes == 'decreased':
                                reason_list=reason_list+['Plues: Tree_{0}_{1}_trans'.format(tree_i, precursor_i)]
                        
                
                #To check if the current precursor transformed to terminal products.
                if product_list !=[]:
                    found_product=True
                    for product_i in product_list:
                        if product_i in PFAS_species_dir['terminal short']+PFAS_species_dir['terminal long']:
                            product_i_conc=PFAS_data_reader(PFASdata_file, product_i, compt_i)
                            product_i_conc_SD=PFAS_data_SD_reader(PFASdata_file, product_i, compt_i)
                            input_mass=product_i_conc['input']
                            output_mass=product_i_conc['output']
                            input_mass_SD=product_i_conc_SD['input']
                            output_mass_SD=product_i_conc_SD['output']
                            thelta=thelta_PFAS(input_mass, output_mass, input_mass_SD, output_mass_SD)
                            
                            #Notes: here we only assume that the precusor transformation was detected only its products are dected to be increased.
                            #there are many other possibile reason, but we only consider this one sitution as it is more confident.
                            #thus there leave a space for this softwaare to further improve.
                            if thelta == 'increased':
                                reason_list=reason_list+['Minus: Tree_{0}_trans_to_{1}'.format(tree_i, product_i)]
                            
            if found_precursor and not any('_trans' in s for s in reason_list):
                reason_list=reason_list+['precursor_ND']
            if not found_precursor and not any('_trans' in s for s in reason_list):
                reason_list=reason_list+['precursor_NF_lib']
                
                
            if found_product and not any('trans_to' in s for s in reason_list):
                reason_list=reason_list+['trans_ND']
            if not found_product and not any('trans_to' in s for s in reason_list):
                reason_list=reason_list+['product_NF_lib']
            
                       
        
        elif precursor_terminal == ['not_found']:
            reason_list=reason_list+['pre_pro_NF_lib']
    
        
    
    return reason_list
        
    


def DT_decrease_analyzer(PFAS_species, PFAS_species_dir, PFASdata_file, compt_i, compt_condition):
    """
    This function performe the data analysis, which return the possibile reason why the PFAS_species decrease across a compartment.

    Parameters
    ----------
    PFAS_species: text
        The name of a PFAS species.
    PFAS_species_dir : dirctionary object
        A dirctionary contains all the terminal species and precursors species.
    PFASdata_file: a csv file path
        The file path to a file contains the PFAS data for each compartment and the basic information of each compartment.
    compt_i: the number of the current compartment.
    compt_condition: the condition of the current compartment.

    Returns
    -------
    Return a list, which contain the reason why the PFAS species change in the compartment,
    in which key words: ND means it was found in Pre_Trans_lib but not detected; NF_lib denotes it was not found in Pre_Trans_lib.

    """
    reason_list=['Plus < Minus:', 'Plus: desorbed from sludge']
    #case 1: if the species is a terminal product.
    if PFAS_species in PFAS_species_dir['terminal short']+PFAS_species_dir['terminal long']:
        
        #consider the ultra short-chain species, no adsorption/desorption.
        if PFAS_species in ['TFA', 'PFPrA', 'PFEtS', 'PFPrS']:
            reason_list=reason_list+['Minus: volatilization and adsorption']
        else:
            reason_list=reason_list+['Minus: adsorbed onto sludge']
    
    
    #case 2: if the species is a precursor.
    elif PFAS_species in PFAS_species_dir['precursor']:
        reason_list=reason_list+['Minus: adsorbed onto sludge']
        
        precursor_terminal=reactant_product_lookup(PFAS_species, compt_condition)
        if precursor_terminal != ['not_found']:
            found_precursor=False
            found_product=False
            
            for dir_i in precursor_terminal:
                tree_i = dir_i['tree_i']
                precursor_list= dir_i['precursor']
                product_list=dir_i['product']
                
                if precursor_list !=[]:
                    found_precursor=True
                    #To check if there is initial precursor transformation.
                    for precursor_i in precursor_list:
                        if precursor_i in PFAS_species_dir['precursor']:
                            precursor_i_conc=PFAS_data_reader(PFASdata_file, precursor_i, compt_i)
                            precursor_i_SD=PFAS_data_SD_reader(PFASdata_file, precursor_i, compt_i)
                            precursor_changes = thelta_PFAS(precursor_i_conc['input'], precursor_i_conc['output'], precursor_i_SD['input'], precursor_i_SD['output'])
                            if precursor_changes == 'decreased':
                                reason_list=reason_list+['Plus: Tree_{0}_{1}_trans'.format(tree_i, precursor_i)]
                        
                
                #To check if the current precursor transformed to terminal products.
                if product_list !=[]:
                    found_product=True
                    for product_i in product_list:
                        if product_i in PFAS_species_dir['terminal short']+PFAS_species_dir['terminal long']:
                            product_i_conc=PFAS_data_reader(PFASdata_file, product_i, compt_i)
                            product_i_conc_SD=PFAS_data_SD_reader(PFASdata_file, product_i, compt_i)
                            input_mass=product_i_conc['input']
                            output_mass=product_i_conc['output']
                            input_mass_SD=product_i_conc_SD['input']
                            output_mass_SD=product_i_conc_SD['output']
                            thelta=thelta_PFAS(input_mass, output_mass, input_mass_SD, output_mass_SD)
                            
                            #Notes: here we only assume that the precusor transformation was detected only its products are dected to be increased.
                            #there are many other possibile reason, but we only consider this one sitution as it is more confident.
                            #thus there leave a space for this softwaare to further improve.
                            if thelta == 'increased':
                                reason_list=reason_list+['Minus: Tree_{0}_trans_to_{1}'.format(tree_i, product_i)]
                            
            if found_precursor and not any('_trans' in s for s in reason_list):
                reason_list=reason_list+['precursor_ND']
            if not found_precursor and not any('_trans' in s for s in reason_list):
                reason_list=reason_list+['precursor_NF_lib']
                
                
            if found_product and not any('trans_to' in s for s in reason_list):
                reason_list=reason_list+['trans_ND']
            if not found_product and not any('trans_to' in s for s in reason_list):
                reason_list=reason_list+['product_NF_lib']
            
        
        elif precursor_terminal == ['not_found']:
            reason_list=reason_list+['pre_pro_NF_lib']
    return reason_list
 
                    
        
#--------------------------------------------------
#old version function: DT_equilibrium_analyzer.
#---------------------------------------------------
   
# def DT_equilibrium_analyzer(PFAS_species, PFAS_species_dir, compt_condition):
#     """
#     This function provides a general reason why the input and output PFAS is in equilibrium.

#     Parameters
#     ----------
#     PFAS_species: text
#         The name of a PFAS species.
#     PFAS_species_dir : dirctionary object
#         A dirctionary contains all the terminal species and precursors species.

#     Returns
#     -------
#     A list containing the general reason.

#     """
    
#     reason_list=['Plus = Minus:']
#     #case 1: if the species is a terminal product.
#     if PFAS_species in PFAS_species_dir['terminal short']+PFAS_species_dir['terminal long']:
#         if PFAS_species in ['TFA', 'PFPrA', 'PFEtS', 'PFPrS']:
#             reason_list +=['Plus: desorption and precursor trans', 'Minus: adsorption and volatilization']
#         else:
#             reason_list +=['Plus: desorption and precursor trans', 'Minus: adsorption']
    
    
#     #case 2: if the species is a precursor.
#     elif PFAS_species in PFAS_species_dir['precursor']:
#         precursor_terminal=reactant_product_lookup(PFAS_species, compt_condition)
#         if precursor_terminal != ['not_found']:
#             found_precursor=False
#             found_product=False
            
#             for dir_i in precursor_terminal:
#                 #tree_i = dir_i['tree_i']
#                 precursor_list= dir_i['precursor']
#                 product_list=dir_i['product']
                
#                 if precursor_list !=[]:
#                     found_precursor=True
#                 if product_list != []:
#                     found_product = True
            
#             if found_precursor and found_product:
#                 reason_list +=['Plus: desorption and precursor trans', 'Minus: adsorption and trans']
#             if not found_precursor and found_product:
#                 reason_list +=['Plus: desorption', 'Minus: adsorption and trans']
            
        
            
#     return reason_list




#-----------------------------------------------
#new version function: DT_equilibrium_analyzer.
#-----------------------------------------------

def DT_equilibrium_analyzer(PFAS_species, PFAS_species_dir, PFASdata_file, compt_i, compt_condition):
    """
    This function performe the data analysis, which return the possibile reason why the PFAS_species retain equilibrium across a compartment.

    Parameters
    ----------
    PFAS_species: text
        The name of a PFAS species.
    PFAS_species_dir : dirctionary object
        A dirctionary contains all the terminal species and precursors species.
    PFASdata_file: a csv file path
        The file path to a file contains the PFAS data for each compartment and the basic information of each compartment.
    compt_i: the number of the current compartment.
    compt_condition: the condition of the current compartment.

    Returns
    -------
    Return a list, which contain the reason why the PFAS species change in the compartment,
    in which key words: ND means it was found in Pre_Trans_lib but not detected; NF_lib denotes it was not found in Pre_Trans_lib.

    """
    reason_list=['Plus = Minus:', 'Plus: desorbed from sludge']
    #case 1: if the species is a terminal product.
    if PFAS_species in PFAS_species_dir['terminal short']+PFAS_species_dir['terminal long']:
        
        #consider the ultra short-chain species, no adsorption/desorption.
        if PFAS_species in ['TFA', 'PFPrA', 'PFEtS', 'PFPrS']:
            reason_list=reason_list+['Minus: volatilization and adsorption']
        else:
            reason_list=reason_list+['Minus: adsorbed onto sludge']
    
    
    #case 2: if the species is a precursor.
    elif PFAS_species in PFAS_species_dir['precursor']:
        reason_list=reason_list+['Minus: adsorbed onto sludge']
        
        precursor_terminal=reactant_product_lookup(PFAS_species, compt_condition)
        if precursor_terminal != ['not_found']:
            found_precursor=False
            found_product=False
            
            for dir_i in precursor_terminal:
                tree_i = dir_i['tree_i']
                precursor_list= dir_i['precursor']
                product_list=dir_i['product']
                
                if precursor_list !=[]:
                    found_precursor=True
                    #To check if there is initial precursor transformation.
                    for precursor_i in precursor_list:
                        if precursor_i in PFAS_species_dir['precursor']:
                            precursor_i_conc=PFAS_data_reader(PFASdata_file, precursor_i, compt_i)
                            precursor_i_SD=PFAS_data_SD_reader(PFASdata_file, precursor_i, compt_i)
                            precursor_changes = thelta_PFAS(precursor_i_conc['input'], precursor_i_conc['output'], precursor_i_SD['input'], precursor_i_SD['output'])
                            if precursor_changes == 'decreased':
                                reason_list=reason_list+['Plus: Tree_{0}_{1}_trans'.format(tree_i, precursor_i)]
                        
                
                #To check if the current precursor transformed to terminal products.
                if product_list !=[]:
                    found_product=True
                    for product_i in product_list:
                        if product_i in PFAS_species_dir['terminal short']+PFAS_species_dir['terminal long']:
                            product_i_conc=PFAS_data_reader(PFASdata_file, product_i, compt_i)
                            product_i_conc_SD=PFAS_data_SD_reader(PFASdata_file, product_i, compt_i)
                            input_mass=product_i_conc['input']
                            output_mass=product_i_conc['output']
                            input_mass_SD=product_i_conc_SD['input']
                            output_mass_SD=product_i_conc_SD['output']
                            thelta=thelta_PFAS(input_mass, output_mass, input_mass_SD, output_mass_SD)
                            
                            #Notes: here we only assume that the precusor transformation was detected only its products are dected to be increased.
                            #there are many other possibile reason, but we only consider this one sitution as it is more confident.
                            #thus there leave a space for this softwaare to further improve.
                            if thelta == 'increased':
                                reason_list=reason_list+['Minus: Tree_{0}_trans_to_{1}'.format(tree_i, product_i)]
                            
            if found_precursor and not any('_trans' in s for s in reason_list):
                reason_list=reason_list+['precursor_ND']
            if not found_precursor and not any('_trans' in s for s in reason_list):
                reason_list=reason_list+['precursor_NF_lib']
                
                
            if found_product and not any('trans_to' in s for s in reason_list):
                reason_list=reason_list+['trans_ND']
            if not found_product and not any('trans_to' in s for s in reason_list):
                reason_list=reason_list+['product_NF_lib']
            
        
        elif precursor_terminal == ['not_found']:
            reason_list=reason_list+['pre_pro_NF_lib']
    return reason_list






#-------------------------------------------------------------------------------------------------
#defing the functions for data analysis in solid (e.g., sludge and biosolid) treatment step.
#-------------------------------------------------------------------------------------------------
def PFAS_Temp_analyzer(PFAS, temp):
    """
    This function provids a general judgement for a given PFAS species, if it will be evaporated or been broken.

    Parameters
    ----------
    PFAS : test
        the name of a PFAS species.
    temp : number
        the temperature at celsius.

    Returns
    -------
    a text: no evaporation, evaporation, or thermal decomposition.

    """
    
    low_evap_group=["PFBA","PFBS","PFPeA","L-PFPeS","PFHxA","PFHxS","PFHpA","L-PFHpS", "4:2FTS","6:2FTS"]
    high_evap_group=["PFOA","PFOS","PFNA","L-PFNS","PFDA", "PFDS","PFUdA","PFDoA","PFTrDA","PFTeDA","8:2FTS","N-MeFOSAA","N-EtFOSAA","FOSA"]
    temp=float(temp)
    
    result='no judgement'
    if temp < 100.0:
        result='no evaporation'
        
    elif temp >= 100 and temp <=200.0:
        if PFAS in low_evap_group:
            result='evaporation'
        else:
            result='no evaporation'
            
    elif temp > 200.0 and temp < 600.0:
        if PFAS in low_evap_group or PFAS in high_evap_group:
            result='evaporation'
    elif temp >= 600.0:
        result='thermal decomposition'
    
    return result




def DT_increase_solid_analyzer(PFAS_species, PFAS_species_dir, PFASdata_file, compt_i, compt_condition):
    """
    This function performe the data analysis, which return the possibile reason why the PFAS_species increase in a compartment iwth pure solid matrix.

    Parameters
    ----------
    PFAS_species : text
        The name of a PFAS species.
    PFAS_species_dir : dirctionary object
        A dirctionary contains all the terminal species and precursors species.
    PFASdata_file: a csv file path
        The file path to a file contains the PFAS data for each compartment and the basic information of each compartment.
    compt_i: the number of the current compartment.
    compt_condition: the condition of the current compartment.

    Returns
    -------
    Return a list, which contain the reason why the PFAS species change in the compartment,
    in which key words: ND means it was found in Pre_Trans_lib but not detected; NF_lib denotes it was not found in Pre_Trans_lib.

    """
    
    reason_list=DT_increase_analyzer(PFAS_species, PFAS_species_dir, PFASdata_file, compt_i, compt_condition)
    reason_list.remove('desorbed from sludge')
    
    return reason_list

        
   
    
def DT_decrease_solid_analyzer(PFAS_species, PFAS_species_dir, PFASdata_file, compt_i, compt_condition, compt_temp):
    """
    This function performe the data analysis, which return the possibile reason why the PFAS_species decrease in a compartment.

    Parameters
    ----------
    PFAS_species: text
        The name of a PFAS species.
    PFAS_species_dir : dirctionary object
        A dirctionary contains all the terminal species and precursors species.
    PFASdata_file: a csv file path
        The file path to a file contains the PFAS data for each compartment and the basic information of each compartment.
    compt_i: the number of the current compartment.
    compt_condition: the condition of the compartment in terms of oxygen.
    compt_temp: the temp of the compartment.

    Returns
    -------
    Return a list, which contain the reason why the PFAS species change in the compartment,
    in which key words: ND means it was found in Pre_Trans_lib but not detected; NF_lib denotes it was not found in Pre_Trans_lib.

    """
    
    reason_list=[]
    
    temp_sensitive=PFAS_Temp_analyzer(PFAS_species, compt_temp)
    reason_list += [temp_sensitive]
    
    reason_list_DT=DT_decrease_analyzer(PFAS_species, PFAS_species_dir, PFASdata_file, compt_i, compt_condition)
    reason_list +=reason_list_DT
    
    reason_list.remove('adsorbed onto sludge')
    
    return reason_list
            
        

def DT_equilibrium_solid_analyzer(PFAS_species, PFAS_species_dir, compt_temp):
    """
    This function provides a general reason why the input and output PFAS is in equilibrium.

    Parameters
    ----------
    PFAS_species: text
        The name of a PFAS species.
    PFAS_species_dir : dirctionary object
        A dirctionary contains all the terminal species and precursors species.
    compt_temp: the temp of the compartment.

    Returns
    -------
    A list containing the general reason.

    """
    temp_sensitive=PFAS_Temp_analyzer(PFAS_species, compt_temp)
    
    reason_list=[]
    #case 1: if the species is a terminal product.
    if PFAS_species in PFAS_species_dir['terminal short'] + PFAS_species_dir['terminal long']:
        if 'no evaporation' == temp_sensitive:
            reason_list=['generation (pre_Trans) equal consumption (unknown sink)']
        elif 'evaporation' == temp_sensitive:
            reason_list=['generation (pre_Trans) equal consumption (evaporation)']
        elif 'thermal decomposition' == temp_sensitive:
            reason_list=['generation (pre_Trans) equal consumption (evaporation&thermal decomposition)']
    
    
    #case 2: if the species is a precursor.
    elif PFAS_species in PFAS_species_dir['precursor']:
        if 'no evaporation' == temp_sensitive:
            reason_list=['generation (pre_Trans) equal consumption (Trans_to)']
        elif 'evaporation' == temp_sensitive:
            reason_list=['generation (pre_Trans) equal consumption (Trans_to&evaporation)']
        elif 'thermal decomposition' == temp_sensitive:
            reason_list=['generation (pre_Trans) equal consumption (evaporation&thermal decomposition)']
                   
    return reason_list





#--------------------------------------------------------------------------------
#define function to performe the PFAS distribution and transformation analysis.
#---------------------------------------------------------------------------------
def DT_analyzer(system_infor, PFAS_species_list, compt_list, PFAS_species_dir, input_file, mass_unit='ug/d'):
    """
    This function is dedicated to performe the PFAS distribution and transformation based on the reuslt of mass balance judgement.

    Parameters
    ----------
    system_infor: a dir object
        a dir contain the basic information of the system.
    PFAS_species_list : list
        here we define which PFAS species are included in this analysis.
    compt_list : list
        a list comtain all the compartments of the system.
    mass_unit : text
        the unit of the mass flux, the default is ug/day.
    PFAS_species_dir : dir
        A dir object contain the terminal and precursors' species.
    input_file : a csv file.
        a path of the inputdata csv file.
    
    Returns
    -------
    Return a dict, it contains six columns: PFAS, compt, input, output, change, and reason.

    """
    #creating empty list to store the generated infor. convert this block as a function.
    PFAS_container=[]
    compt_container=[]
    inputmass_container=[]
    outputmass_container=[]
    mass_change_container=[]
    
    reason_container=[]
    
    
    for PFAS_i in PFAS_species_list:
        for compt_i in compt_list:
            PFAS_container += [PFAS_i]
            compt_container += [system_infor[compt_i]['label']]
            
            PFAS_mass_flow=PFAS_data_reader(input_file, PFAS_i, compt_i)
            PFAS_mass_flow_SD=PFAS_data_SD_reader(input_file, PFAS_i, compt_i)
            PFAS_mass_input=PFAS_mass_flow['input']
            PFAS_mass_output=PFAS_mass_flow['output']
            PFAS_mass_input_SD=PFAS_mass_flow_SD['input']
            PFAS_mass_output_SD=PFAS_mass_flow_SD['output']
            
            mass_unit=PFAS_mass_flow['unit']
            
            inputmass_container += [PFAS_mass_input]
            outputmass_container += [PFAS_mass_output]
            
            PFAS_change=thelta_PFAS(PFAS_mass_input, PFAS_mass_output, PFAS_mass_input_SD, PFAS_mass_output_SD)
            mass_change_container += [PFAS_change]
            
            reason_i=['no_reasoning']
            #The matrix of this compartment is in liquid phase.
            if system_infor[compt_i]['flux_unit'] == 'm3/d':
                if PFAS_change == 'increased':
                    reason_i=DT_increase_analyzer(PFAS_i, PFAS_species_dir, input_file, compt_i, system_infor[compt_i]['condition'])
                elif PFAS_change == 'decreased':
                    reason_i=DT_decrease_analyzer(PFAS_i, PFAS_species_dir, input_file, compt_i, system_infor[compt_i]['condition'])
                elif PFAS_change == 'equilibrium':
                    reason_i=DT_equilibrium_analyzer(PFAS_i, PFAS_species_dir, input_file, compt_i, system_infor[compt_i]['condition'])
             
            #The matrix of this compartment is in solid phase.
            compt_temp=system_infor[compt_i]['temp']
            if system_infor[compt_i]['flux_unit'] == 'ton/d':
                if PFAS_change == 'increased':
                    reason_i=DT_increase_solid_analyzer(PFAS_i, PFAS_species_dir, input_file, compt_i, system_infor[compt_i]['condition'])
                elif PFAS_change == 'decreased':
                    reason_i=DT_decrease_solid_analyzer(PFAS_i, PFAS_species_dir, input_file, compt_i, system_infor[compt_i]['condition'], compt_temp)
                elif PFAS_change == 'equilibrium':
                    reason_i=DT_equilibrium_solid_analyzer(PFAS_i, PFAS_species_dir, compt_temp)
            
            reason_container += [reason_i]
            
         
    
    
    #combining these container into a dataframe object.
    combinations_dir={'PFAS': PFAS_container, 'compt':compt_container, 'input ({})'.format(mass_unit):inputmass_container, 'output ({})'.format(mass_unit): outputmass_container,
                      'change':mass_change_container, 'reason': reason_container}
    
    #combinations_df=pd.DataFrame(combinations_dir, index=list(range(len(PFAS_container))))
    
    return combinations_dir







#------------------------------------------------------------
#define functions to performe the mass balance evaluation.
#------------------------------------------------------------
def compt_total_massflow_summary(PFAS_species, mass_unit, DT_analysis_results_file='./Output/PFAS_DT_analysis_results.csv'):
    """
    This function summary PFAS mass flow of input and output of each compartment. The number of PFAS species can be specifed in the input parameter - PFAS_species
    
    Notes:
    1. The 'total' in the function name indactes the total input and output sources for a compartment.
    2. This function can be employed only the PFAS_DT_analysis csv file are created.
    3. Only the number will be summaried, N and ND will be ignored.

    Parameters
    ----------
    PFAS_species: a list contain PFAS species.
        a single species means only calculate one species mass flow for each compartment, if all the species are included, it will calculate the sum_PFAS mass flow.
    
    mass_unit: ug/d or mg/d.
    
    DT_analysis_results_file :  a csv file containing the PFAS DT analysis results, optional
        DESCRIPTION. The default is './Output/PFAS_DT_analysis_results.csv'.

    Returns
    -------
    A dir object: {'compt_i': {'sum_input (unit)': xxx, 'sum_output (unit)': xxxx}}

    """
        
    reasoning_results_df=pd.read_csv(DT_analysis_results_file, header=0, index_col=0)
    
    #Extracting the columns required for the follwing calculation.
    extracted_columns_df=reasoning_results_df[['PFAS', 'compt', 'input ({})'.format(mass_unit), 'output ({})'.format(mass_unit)]]
    
    compt_massflow_df=extracted_columns_df[extracted_columns_df['PFAS'].isin(PFAS_species)]
    
    compt_massflow_dir={}
    compt_list=compt_massflow_df['compt'].unique()
    
    for compt_i in compt_list:
        compt_i_df= compt_massflow_df[compt_massflow_df['compt']==compt_i]
        
        compt_i_df_1=compt_i_df[compt_i_df['input ({})'.format(mass_unit)] != 'ND']
        compt_i_input_df=compt_i_df_1[compt_i_df_1['input ({})'.format(mass_unit)] != 'D']
        
        compt_i_df_2=compt_i_df[compt_i_df['output ({})'.format(mass_unit)] != 'ND']
        compt_i_output_df=compt_i_df_2[compt_i_df_2['output ({})'.format(mass_unit)] != 'D']
        
        input_mass_list=compt_i_input_df['input ({})'.format(mass_unit)].tolist()
        output_mass_list=compt_i_output_df['output ({})'.format(mass_unit)].tolist()
        
        input_mass_list=[float(i) for i in input_mass_list]
        output_mass_list=[float(i) for i in output_mass_list]
        
        sum_input_mass=round(sum(input_mass_list), 3)
        sum_output_mass=round(sum(output_mass_list), 3)
        
        if sum_input_mass ==0: sum_input_mass='<LOQ'
        if sum_output_mass == 0: sum_output_mass='<LOQ'
        
        compt_massflow_dir[compt_i]={'sum_input ({})'.format(mass_unit): sum_input_mass, 'sum_output ({})'.format(mass_unit): sum_output_mass}
    
    return compt_massflow_dir





def PFAS_massflow_reader(PFAS_species, from_compt, to_compt, inputdata_file):
    """
    This function is designed to read the PFAS mass for an edges among the compartment network, which will be employed in function PFAS_massflow_network.

    Parameters
    ----------
    PFAS_species : a list contain PFAS species.
        a single species means only calculate one species mass flow for each compartment, if all the species are included, it will calculate the sum_PFAS mass flow. 
        
    from_compt : text
        the name of a parent node in the compartment network.
    to_compt : text
        The name of a child node in the compartment network.
    inputdata_file : a csv file
       It contains two type of information: the PFAS concentration in each compartment; the basic information of each compartment
       and the connection among compartment.

    Returns
    -------
    [PFAS mass, PFAS mass SD], The PFAS mass, a number with unit ug/d or mg/d, or "<LOQ" if there is no quantitative results.

    """
    #step1: PFAS massflow
    
    PFAS_massflow_list=[]
    
    
    if 'compt_' in to_compt:
        for PFAS_i in PFAS_species:
            PFAS_massflow_compt_dir=PFAS_data_reader(inputdata_file, PFAS_i, to_compt)
            input_PFAS_mass_dir=PFAS_massflow_compt_dir['input_dir']
            
            from_compt_label=compt_label_conveter(from_compt, inputdata_file)
            input_PFAS_mass=input_PFAS_mass_dir[from_compt_label]
            
            PFAS_massflow_list +=[input_PFAS_mass]
    
    if 'Effluent' in to_compt:
        for PFAS_i in PFAS_species:
            PFAS_massflow_compt_dir=PFAS_data_reader(inputdata_file, PFAS_i, from_compt)
            input_PFAS_mass_dir=PFAS_massflow_compt_dir['output_dir']
            
            to_compt_label=compt_label_conveter(to_compt, inputdata_file)
            input_PFAS_mass=input_PFAS_mass_dir[to_compt_label]
            
            PFAS_massflow_list +=[input_PFAS_mass]
    
    
    PFAS_massflow_list_0=[i for i in PFAS_massflow_list if i !='D']
    PFAS_massflow_list_1 =[i for i in PFAS_massflow_list_0 if i !='ND']
    
    if len(PFAS_massflow_list_1) != 0:
        PFAS_massflow_list_2=[float(i) for i in PFAS_massflow_list_1]
        PFAS_mass=round(sum(PFAS_massflow_list_2), 3)
    
    elif len(PFAS_massflow_list_1) == 0: PFAS_mass='<LOQ'


    
    #step2: PFAS mass flow SD
    PFAS_massflow_SD_list=[]
    
    
    if 'compt_' in to_compt:
        for PFAS_i in PFAS_species:
            PFAS_massflow_SD_compt_dir=PFAS_data_SD_reader(inputdata_file, PFAS_i, to_compt)
            input_PFAS_mass_SD_dir=PFAS_massflow_SD_compt_dir['input_dir']
            
            from_compt_label=compt_label_conveter(from_compt, inputdata_file)
            input_PFAS_mass_SD=input_PFAS_mass_SD_dir[from_compt_label]
            
            PFAS_massflow_SD_list +=[input_PFAS_mass_SD]
    
    if 'Effluent' in to_compt:
        for PFAS_i in PFAS_species:
            PFAS_massflow_SD_compt_dir=PFAS_data_SD_reader(inputdata_file, PFAS_i, from_compt)
            input_PFAS_mass_SD_dir=PFAS_massflow_SD_compt_dir['output_dir']
            
            to_compt_label=compt_label_conveter(to_compt, inputdata_file)
            input_PFAS_mass_SD=input_PFAS_mass_SD_dir[to_compt_label]
            
            PFAS_massflow_SD_list +=[input_PFAS_mass_SD]
    
    
    PFAS_massflow_SD_list_0=[i for i in PFAS_massflow_SD_list if i !='D']
    PFAS_massflow_SD_list_1 =[i for i in PFAS_massflow_SD_list_0 if i !='ND']
    
    if len(PFAS_massflow_SD_list_1) !=0:
        PFAS_massflow_SD_list_2=[float(i) for i in PFAS_massflow_SD_list_1]
        PFAS_massflow_SD_list_3=[i**2 for i in PFAS_massflow_SD_list_2]
        PFAS_mass_SD=round( (sum(PFAS_massflow_SD_list_3))**0.5, 3)
    
    elif len(PFAS_massflow_SD_list_1) ==0: PFAS_mass_SD='<LOQ'
    
    
    return [PFAS_mass, PFAS_mass_SD]
        
    
    

def workflow_list_arrangment(from_list, to_list, inputdata_file):
    """
    This function used to divide the flow diagram into flow-chains, which is employed to show the workflow in GUI app PFAS_DT_Investigator.

    Parameters
    ----------
    from_list : list
        list contains parent nodes, in which each compartment represented as its label.
    to_list : list
        list contains child nodes, in which each compartment represented as its label.

    Returns
    -------
    A dir: {'chain_1': [node1, node2, ...]}, in which workflow from node1 >> node2 >>...

    """
    
    edges_number=[]
    for i in range(len(from_list)):
        parent_compt_i=compt_label_conveter(from_list[i], inputdata_file)
        if 'Influent_' in parent_compt_i:
            parent_label_i_to=to_list[i]
            parent_compt_i_to=compt_label_conveter(parent_label_i_to, inputdata_file)
            compt_number = ''.join(filter(lambda x: x.isdigit(), parent_compt_i_to))
            compt_number=float(compt_number)-0.5
            edges_number += [compt_number]
        
        else:
            compt_number = ''.join(filter(lambda x: x.isdigit(), parent_compt_i))
            compt_number=float(compt_number)
            edges_number += [compt_number]
    
    workflow_df=pd.DataFrame({'from': from_list, 'to': to_list, 'No.': edges_number}, index=list(range(len(from_list))))
    
    workflow_sorted_df=workflow_df.sort_values(by='No.', ascending=True)
    
    result_dict=workflow_sorted_df.to_dict(orient='list')
    
    return result_dict
    


def chain_seperation(from_list, to_list, PFAS_massflow, PFAS_massflow_SD):
    """
    Seperate the workflow list into different tables, in which each table represent a tree.

    Parameters
    ----------
    from_list : list
        list contains parent nodes, in which each compartment represented as its label.
    to_list : list
        list contains child nodes, in which each compartment represented as its label.

    Returns
    -------
    chain_dict : a dict object
        {'chain_1': {'from': [], 'to': [], 'tree': [], 'PFAS_massflow': [],  'PFAS_massflow_SD': []}}.

    """
    chain_dict={}
    tree_No=[]
    tree_number=0
    for i in list(range(len(from_list))):
        if 'Influent_' in from_list[i]: tree_number +=1
        tree_No += [tree_number]
    
    workflow_df=pd.DataFrame({'from': from_list, 'to': to_list, 'Tree_No': tree_No, 'PFAS_massflow':PFAS_massflow, 'PFAS_massflow_SD':PFAS_massflow_SD}, index=list(range(len(from_list))))
    unique_tree=list(set(tree_No))
    for tree_i in unique_tree:
        tree_i_df=workflow_df[workflow_df['Tree_No']==tree_i]
        tree_i_dict=tree_i_df.to_dict(orient='list')
        chain_dict['Tree_{}'.format(tree_i)]=tree_i_dict
    return chain_dict

        




def PFAS_massflow_network(system_infor_dir, PFAS_species, inputdata_file, mass_unit='ug/d'):
    """
    This function creat the network list, in which the PFAS mass flow (a single species or total) direction among the investigated compartments are presented.
    Based on this networklist, the flow chart of the investigated system can be plotted, which will be employed in GUI version development.

    Parameters
    ----------
    system_infor_dir : dir object
        a dir object contain the basic information of the investigated system, including the structure information.
    PFAS_species: a list of PFAS species.
        
    mass_unit: the ddefault is ug/day
    
    inputdata_file: the input data csv file

    Returns
    -------
    a network list, representing the arc from one compt to another compt, including the total PFAS mass and SD. 

    """
    
    compt_list=list(system_infor_dir.keys())
    
    from_col=[]
    to_col=[]
    massflow_col=[]
    massflow_SD_col=[]
    
    for compt_i in compt_list:
        input_source_list=system_infor_dir[compt_i]['input_source']
        output_source_list=system_infor_dir[compt_i]['output_source']
        
        compt_i_label=system_infor_dir[compt_i]['label']
        
        for input_source_i in input_source_list:
            input_source_i_label=compt_label_conveter(input_source_i, inputdata_file)
            from_col +=[input_source_i_label]
            to_col +=[compt_i_label]
        
        for output_source_i in output_source_list:
            output_source_i_label=compt_label_conveter(output_source_i, inputdata_file)
            from_col +=[compt_i_label]
            to_col +=[output_source_i_label]
    
        
    #remove replicate edges in the constructed network list.
    combined_intuple_list=[]
    for i in range(len(from_col)):
        from_i=from_col[i]
        to_i=to_col[i]
        
        combined_intuple_list +=[(from_i, to_i)]
              
    combined_intuple_list_1=list(set(combined_intuple_list)) #remove all replicate edges by converting it into a set object.
    
    from_col_1=[]
    to_col_1=[]
    for i in combined_intuple_list_1:
        from_col_1 +=[i[0]]
        to_col_1 +=[i[1]]
    
   #reorganize the list based on the sequence of the system workflow.
    organized_dict=workflow_list_arrangment(from_col_1, to_col_1, inputdata_file)
    from_col_1=organized_dict['from']
    to_col_1=organized_dict['to']
    
   
   #read the PFAS mass for each edges, from -> to.
    for j in range(len(from_col_1)):
       from_i = from_col_1[j]
       to_i = to_col_1[j]
       
       from_i_compt=compt_label_conveter(from_i, inputdata_file)
       to_i_compt=compt_label_conveter(to_i, inputdata_file)
       
       PFAS_mass=PFAS_massflow_reader(PFAS_species, from_i_compt, to_i_compt, inputdata_file)
       massflow_col +=[PFAS_mass[0]]
       massflow_SD_col +=[PFAS_mass[1]]
       
    #combining the results to build the network list.
    massflow_network_dir={}
    massflow_network_dir['from']=from_col_1
    massflow_network_dir['to']=to_col_1
    massflow_network_dir['PFAS_mass ({})'.format(mass_unit)]=massflow_col
    massflow_network_dir['SD']=massflow_SD_col
    
    
    return massflow_network_dir
        


        



    
            
            
        








