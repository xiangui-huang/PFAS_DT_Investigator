# -*- coding: utf-8 -*-
"""
The distribution and transformation of PFAS in wastewater treatment plant (WWTP), DT_PFAS_investigator.version 0

Function definition module file 2.
Define functions used to simulate the PFAS distribution in WWTP.

ENV: BeeWare

Created on Fri Jan 16 09:59:17 2026
Xiangui Huang, Ben-Gurion University, Israel.
@author: xiang
"""

import numpy as np
import pandas as pd

from scipy.optimize import least_squares

#-------------------------------------------------------------------------------------------------------------------
#Here, some global variables are defined as assigned as default value,
# they can be used in any fucntion in this module also can be imported into other script and changed their values.
# In this app development, they will be imported into my app script, and will be changed through GUI by user.
#------------------------------------------------------------------------------------------------------------------
RSL_Clarifier1=0.0001
RSL_Reactor_Clarifier2=0.001
optimization_method='dogbox'
optimization_loss='soft_l1'








#The function simulating the PFAS distribution and transforamtion in the Parimary Clarifier (compartment No. I).
def PFAS_DT_Clarifier1(V0=5000, Ci=100, Cj=100, Ck=10, #the boundary condition, volume of influent and PFAS concentration
                       RSL=RSL_Clarifier1, fj=0.05, fk=0.0025, #the parameters of matrixes and precursor transforamtion. 
                       PC_i=1000, PC_j=50, PC_k=1000, #the partition coefficient.
                       ):
    """
    This is a equilibrium model that simulating the PFAS distribution and transforamtion in the Parimary Clarifier.

    Parameters
    ----------
    #the boundary condition:
    Vo : int/float, optional
        The flux of influent (m3/day). The default is 5000.
    Ci : int/float, optional
        The precursor concentration in influent (ng/L). The default is 100.
    Cj : int/float, optional
        The short-chain PFAS in influent (ng/L). The default is 100.
    Ck : int/float, optional
        The long-chain PFAS in influent (ng/L). The default is 10.
     
    #the parameters of matrixes and precursor transforamtion
    RSL : float, optional
        The ratio of sludge (ton dry weight) to liquid (m3) in the parimary clarifier. The default is 0.0001.
    fj : float, optional
        The fraction of precursor transfer to short-chain PFAS. The default is 0.05.
    fk : float, optional
        The fraction of precursor transfer to long-chain PFAS. The default is 0.0025.
    
    #the partition coefficient.
    PC_i : int, optional
        The partition coefficient of precurosr. The default is 1000.
    PC_j : int, optional
        The partition coefficient of short-chain PFAS. The default is 50.
    PC_k : int, optional
        The partition coefficient of long-chain PFAS. The default is 1000.
   

    Returns
    -------
    dict object: The flux (V_I), precursor concentration (Cli_I), short-chain (Clj_I), long-chain (Clk_I) in the effluent of parimary clarifier.

    """
    V_I=V0
    ms=V0*RSL
    
    #total mass of precursors
    sum_mi_I=V0*Ci*(1-fj-fk)
    #total mass of short-chain
    sum_mj_I=V0*Cj+V0*Ci*fj
    #total mass of long-chain
    sum_mk_I=V0*Ck + V0*Ci*fk
    
    # Define the coefficient matrix A
    A_i = np.array([[ms, V0], [1, -PC_i]])
    A_j = np.array([[ms, V0], [1, -PC_j]])
    A_k = np.array([[ms, V0], [1, -PC_k]])

    # Define the constant vector b
    b_i = np.array([sum_mi_I, 0])
    b_j = np.array([sum_mj_I, 0])
    b_k = np.array([sum_mk_I, 0])
    
    
    # Solve the linear system
    try:
        x_i = np.linalg.solve(A_i, b_i)
        x_j = np.linalg.solve(A_j, b_j)
        x_k = np.linalg.solve(A_k, b_k)
    
    # Output: Solution: [2. 1.]
    except np.linalg.LinAlgError as e:
        print(f"Error solving linear system: {e}")
        print("This may be because the matrix A is singular (non-invertible).")
    
    
    x_i =x_i.tolist()
    x_j =x_j.tolist()
    x_k =x_k.tolist()
    
    results_dict={}
    
    results_dict['V_I']=V_I
    results_dict['Cli_I']=round(x_i[1], 3)
    results_dict['Clj_I']=round(x_j[1], 3)
    results_dict['Clk_I']=round(x_k[1], 3)
    
    return results_dict
    





    
#The function simulating the PFAS distribution and transforamtion in the Bio-reactor and Parimary Clarifier (compartment No. III).
def PFAS_DT_Ractor_Clarifier2(V_I, Cl_i, Cl_j, Cl_k, #the influent form the parimary clarifier.
                              Ci_RAS=10000, Cj_RAS=1000, Ck_RAS=10000, #the PFAS concentration in the recycled activited sludge (RAS).
                              RSL=RSL_Reactor_Clarifier2, fj=0.5, fk=0.025, f_RAS=0.8, #the parameters of matrixes and precursor transforamtion.
                              PC_i=1000, PC_j=50, PC_k=1000, #the partition coefficient.
        ):
    """
    This is a equilibrium model that simulating the PFAS distribution and transforamtion in the Bio-reactor and Parimary Clarifier.

    Parameters
    ----------
    #the influent form the parimary clarifier.
    V_I : int/float
        The wastewater flux from the parimary clarifier (m3/day).
    Cl_i : int/float
        The precursor's concentration in the V_I.
    Cl_j : int/float
        The short-chain PFAS's concentration in the V_I.
    Cl_k : int/float
        The long-chain PFAS's concentration in the V_I.
        
    #the PFAS concentration in the recycled activited sludge (RAS).
    Ci_RAS : int/float, optional
        The precursor's concentration in the RAS. The default is 10000.
    Cj_RAS : int/float, optional
        The short-chain PFAS's concentration in the RAS. The default is 1000.
    Ck_RAS : int/float, optional
        The long-chain PFAS's concentration in the RAS. The default is 10000.
    
    #the parameters of matrixes and precursor transforamtion.
    RSL : float, optional
        DESCRIPTION. The default is 0.001.
    fj : float, optional
        The fraction of precursor transfer to short-chain PFAS. The default is 0.65.
    fk : float, optional
        The fraction of precursor transfer to long-chain PFAS. The default is 0.25.
    f_RAS : float, optional
        The fraction of the RAS from the secondary clarifier. The default is 0.8.
    
    #the partition coefficient.
    PC_i : float, optional
        The partition coefficient of precurosr. The default is 1000.
    PC_j : float, optional
        The partition coefficient of short-chain PFAS. The default is 50.
    PC_k : float, optional
        The partition coefficient of long-chain PFAS. The default is 1000.
    
    
    Returns
    -------
    dict obj: the flux (V_III), PFAS concentration in the effluent and RAS:
        precursor (Cli_III and Csi_III), short-chain (Clj_III and Csj_III), long-chain (Clk_III and Csk_III).

    """
    V_III=V_I
    
    ms=V_III*RSL
    m_old_RAS=f_RAS*ms
    m_new_RAS=(1-f_RAS)*ms
    
    #total mass of precursors
    sum_mi_III=(V_III*Cl_i + m_old_RAS*Ci_RAS)*(1-fj-fk)
    #total mass of short-chain
    sum_mj_III=V_III*Cl_j+m_old_RAS*Cj_RAS + (V_III*Cl_i + m_old_RAS*Ci_RAS)*fj
    #total mass of long-chain
    sum_mk_III=V_III*Cl_k+ m_old_RAS*Ck_RAS + (V_III*Cl_i + m_old_RAS*Ci_RAS)*fk
    
    # Define the coefficient matrix A
    A_i = np.array([[ms, V_III], [1, -PC_i]])
    A_j = np.array([[ms, V_III], [1, -PC_j]])
    A_k = np.array([[ms, V_III], [1, -PC_k]])
    
    # Define the constant vector b
    b_i = np.array([sum_mi_III, 0])
    b_j = np.array([sum_mj_III, 0])
    b_k = np.array([sum_mk_III, 0])
    
    # Solve the linear system
    try:
        x_i = np.linalg.solve(A_i, b_i)
        x_j = np.linalg.solve(A_j, b_j)
        x_k = np.linalg.solve(A_k, b_k)
    
    # Output: Solution: [2. 1.]
    except np.linalg.LinAlgError as e:
        print(f"Error solving linear system: {e}")
        print("This may be because the matrix A is singular (non-invertible).")
    
    
    x_i =x_i.tolist()
    x_j =x_j.tolist()
    x_k =x_k.tolist()
    
    results_dict={}
    
    results_dict['V_III']=V_III
    results_dict['Cli_III']=round(x_i[1], 3)
    results_dict['Clj_III']=round(x_j[1], 3)
    results_dict['Clk_III']=round(x_k[1], 3)
    
    results_dict['Csi_III']=round(x_i[0], 3)
    results_dict['Csj_III']=round(x_j[0], 3)
    results_dict['Csk_III']=round(x_k[0], 3)
    
    
    
    return results_dict
    
    






#---------------------------------------------------------------------------------------
#Define the function for parameters calibration, given the input and output are observed.
#For this model, we calibrate six parameters: fj, fk, f_RAS, PC_i, PC_j, PC_k.
#calibration method: least square.
#---------------------------------------------------------------------------------------- 


#define the residuales function
def residules_Clarifier1(theta, X, Y_obs, sigma=None, use_log=True):
    """
    This function is dedicated to calculate the residules of the model simulating PFAS_DT_Ractor_Clarifier2.

    Parameters
    ----------
    theta : list
        The parameter intended to be calibrated [fj, fk, PC_i, PC_j, PC_k].
    X : list [N*4]
        The measured input data array with N observations [ [V0, C_i, C_j, C_k], ... []].
    Y_obs : list [N*4]
        The measured output data array [[V_I, Cli_I, Clj_I, Clk_I], ..., []].
    sigma : list, optional
        The weight for each calibrated parameters [w_fj, w_fk, w_PC_i, W_PC_j, W_PC_k]. The default is None.
    use_log : bool, optional
        re-scale the residules at log-scale. The default is True.
   

    Returns
    -------
    A 1D residual vector [res_Cli_I, res_Clj_I, res_Clk_I, ..., N].

    """
    res_list=[]
    
    #extract the parameter value.
    fj_0, fk_0, PC_i_0, PC_j_0, PC_k_0=theta
   
    
    #extract the input data and output data
    for i in range(len(X)):
        X_i=X[i]
        Y_obs_i=Y_obs[i]
        
        V_0, C_i, C_j, C_k =X_i       
        
        Y_hat_dict=PFAS_DT_Clarifier1(V0 = V_0, Ci = C_i, Cj = C_j, Ck = C_k, fj = fj_0, fk = fk_0, PC_i = PC_i_0, PC_j = PC_j_0, PC_k = PC_k_0)
    
        Y_hat_array=np.array([Y_hat_dict[i] for i in list(Y_hat_dict.keys()) ])
        Y_obs_array=np.array(Y_obs_i)
    
        if use_log:
            r_array=np.log( (abs(Y_hat_array-Y_obs_array)+1)/(Y_obs_array+1))
        else:
            r_array=abs(Y_hat_array-Y_obs_array)/(Y_obs_array+1)
        
        # Optional weighting by measurement uncertainty (preferred if you have it)
        if sigma is not None:
            r_array = r_array / sigma  # sigma shape (N,3)
            
        r_array_0=np.around(r_array, 6)
        r_array_1=np.delete(r_array_0, [0]) #here remove the flux, as it is known keep constant.
        
        res_list +=[list(r_array_1)]



    res_array=np.array(res_list)
    return res_array.ravel() 




def residules_Reactor_Clarifier2(theta, X, Y_obs, included_RAS= True, sigma=None, use_log=True):
    """
    This function is dedicated to calculate the residules of the model simulating PFAS_DT_Clarifier1.

    Parameters
    ----------
    theta : list
        The parameter intended to be calibrated [fj, fk, f_RAS, PC_i, PC_j, PC_k].
    X : list
        The measured input data array, the PFAS in the effluent from clarifier1 and in RAS, [V_I, Cli_I, Clj_I, Clk_I, Ci_RAS, Cj_RAS, Ck_RAS].
    Y_obs : list
        The measured output data array, the effluent from clarifier2:
            if included_RAS == True: [V_III, Cli_III, Clj_III, Clk_III, Csi_III, Csj_III, Csk_III]
            else: [V_III, Cli_III, Clj_III, Clk_III].
    included_RAS: bool, optional
        Here, the user can choice whether consider the PFAS conc in the RAS, the default is False.
    sigma : list, optional
        The weight for each calibrated parameters [w_fj, w_fk, w_f_RAS, w_PC_i, W_PC_j, W_PC_k]. The default is None.
    use_log : bool, optional
        re-scale the residules at log-scale. The default is True.
   

    Returns
    -------
    A 1D residual vector.

    """
    res_list=[]
    
    #extract the parameter value.
    fj_0, fk_0, f_RAS_0, PC_i_0, PC_j_0, PC_k_0=theta
   
    #extract the input data and output data
    for i in range(len(X)):
        X_i=X[i]
        Y_obs_i=Y_obs[i]
    
    
        V_I, C_i, C_j, C_k, Ci_RAS, Cj_RAS, Ck_RAS = X_i
        
        if included_RAS:
            V_III, Cli_III, Clj_III, Clk_III, Csi_III, Csj_III, Csk_III =Y_obs_i
        else:
            V_III, Cli_III, Clj_III, Clk_III = Y_obs_i
        
        
        Y_hat_dict=PFAS_DT_Ractor_Clarifier2(V_I = V_I, Cl_i= C_i, Cl_j=C_j, Cl_k=C_k, Ci_RAS=Ci_RAS, Cj_RAS=Cj_RAS, Ck_RAS=Ck_RAS, fj = fj_0, fk = fk_0, f_RAS= f_RAS_0, PC_i = PC_i_0, PC_j = PC_j_0, PC_k = PC_k_0)
        
        if included_RAS:
            Y_hat_array=np.array([Y_hat_dict[i] for i in list(Y_hat_dict.keys()) ])
        else:
            Y_hat_array=np.array([Y_hat_dict[i] for i in ['V_III', 'Cli_III', 'Clj_III', 'Clk_III'] ])
        
        Y_obs_array=np.array(Y_obs_i)
        
        if use_log:
            r_array=np.log( (abs(Y_hat_array-Y_obs_array)+1)/(Y_obs_array+1))
        else:
            r_array=abs(Y_hat_array-Y_obs_array)/(Y_obs_array+1)
        
        # Optional weighting by measurement uncertainty (preferred if you have it)
        if sigma is not None:
            r_array = r_array / sigma  # sigma shape (N,3)
            
        r_array_0=np.around(r_array, 6)
        r_array_1=np.delete(r_array_0, [0])
        
        res_list +=[list(r_array_1)]

    res_array=np.array(res_list)
    return res_array.ravel() 
     
    









#define the function to perform the parameter calibration by least_square method.
def residules_Clarifier1_optimization(theta_0, lb, ub, X, Y_obs, sigma=None, use_log=False):
    """
    This function is dedicated to perform the parameter optimization based on the given initaL parameter value and observations.
    Notice: five parameters are need to be calibrated, thus at least five observations are required.

    Parameters
    ----------
    theta_0 : list
        The initial parameters value [fj, fk, PC_i, PC_j, PC_k].
    lb : list
        The lower bound of the parameter [fj_lb, fk_lb, PC_i_lb, PC_j_lb, PC_k_lb].
    ub : list
        The lower bound of the parameter [fj_ub, fk_ub, PC_i_ub, PC_j_ub, PC_k_ub].
    X : list
        The measured input data array [V0, C_i, C_j, C_k].
    Y_obs : list
        The measured output data array [V_I, Vli_I, Vlj_I, Vlk_I].
    sigma : list, optional
        The weight for each calibrated parameters. The default is None.
    use_log : bool, optional
        re-scale the residules at log-scale. The default is False.

    Returns
    -------
    A object, contain all the infor of optimized results, e.g., optimized parameter and Minimalized residules.

    """
    
    res = least_squares(
    residules_Clarifier1,
    x0=theta_0,
    bounds=(lb, ub),
    args=(X, Y_obs, sigma, use_log),
    method=optimization_method,
    loss=optimization_loss,      
    f_scale=1.0,
    jac="3-point"      # or provide your own jacobian if you can
    )

    theta_hat = res.x
    optimized_residules = res.fun
    nfev= res.nfev
    cost = res.cost
        
    
    return {'theta_hat': np.around(theta_hat, 3), 'min_residules': np.around(optimized_residules, 3), 'nfev': np.around(nfev, 3), 'cost': round(cost, 3),'message': res.message}





def residules_Reactor_Clarifier2_optimization(theta_0, lb, ub, X, Y_obs, included_RAS=True, sigma=None, use_log=False):
    """
    This function is dedicated to perform the parameter optimization based on the given initaL parameter value and observations.
    Notice: six parameters are need to be calibrated, thus at least six observations are required.

    Parameters
    ----------
    theta_0 : list
        The initial parameters value [fj, fk, f_RAS, PC_i, PC_j, PC_k].
    lb : list
        The lower bound of the parameter [fj_lb, fk_lb, f_RAS_lb, PC_i_lb, PC_j_lb, PC_k_lb].
    ub : list
        The lower bound of the parameter [fj_ub, fk_ub, f_RAS_ub, PC_i_ub, PC_j_ub, PC_k_ub].
     X : list
         The measured input data array, the PFAS in the effluent from clarifier1 and in RAS, [V_I, Cli_I, Clj_I, Clk_I, Ci_RAS, Cj_RAS, Ck_RAS].
     Y_obs : list
         The measured output data array, the effluent from clarifier2:
             if included_RAS == True: [V_III, Cli_III, Clj_III, Clk_III, Csi_III, Csj_III, Csk_III]
             else: [V_III, Cli_III, Clj_III, Clk_III].
     included_RAS: bool, optional
         Here, the user can choice whether consider the PFAS conc in the RAS, the default is True.
         It should be noticed that the PFAS conc in RAS in the X should be equal to the one in the Y_obs as we assume the system is under equilibrium.
         
     sigma : list, optional
         The weight for each calibrated parameters [w_fj, w_fk, w_f_RAS, w_PC_i, W_PC_j, W_PC_k]. The default is None.
     use_log : bool, optional
         re-scale the residules at log-scale. The default is True.

    Returns
    -------
    A object, contain all the infor of optimized results, e.g., optimized parameter (theta) and Minimalized residules.

    """
    
    res = least_squares(
    residules_Reactor_Clarifier2,
    x0=theta_0,
    bounds=(lb, ub),
    args=(X, Y_obs, included_RAS, sigma, use_log),
    method=optimization_method,
    loss=optimization_loss,
    f_scale=1.0,
    jac="3-point"      # or provide your own jacobian if you can
    )

    theta_hat = res.x
    theta_hat = np.around(theta_hat, 3)
    theta_hat=theta_hat.tolist()
    optimized_residules = res.fun
    nfev= res.nfev
    cost = res.cost
    
    
    return {'theta_hat': theta_hat, 'min_residules': np.around(optimized_residules, 3), 'nfev': np.around(nfev, 3), 'cost': round(cost, 3),'message': res.message}








#---------------------------------------------------------------------------
#define the function read the observation and parameter set from csv file.
#---------------------------------------------------------------------------
def observation_reader(obervation_file, compt_name):
    """
    This function is dedicated to read the observation and parameter set from a csv file.

    Parameters
    ----------
    obervation_file : a csv file path
        a csv file contain the parameters setting and observations.
    compt_name : the compartment name: Clarifier1 and Reactor_Clarifier2.
        the compartment name allow the user to select which compoartment model will be used.

    Returns
    -------
    a dict object: {theta_0:, lb: , ub: , X: , Y_obs: }

    """
    if compt_name not in ['Clarifier1', 'Reactor_Clarifier2']:
        res='Error, the entered compt name not recongnized.'
    
    
    else:
    
        observation_df = pd.read_csv(obervation_file, header=0, index_col=0)
        observation_df=observation_df.fillna(0)
        
        observation_df_1=observation_df[observation_df.index.str.contains(compt_name)]
        observation_df_1.set_index('parameter', inplace=True)
        #observation_df_1.drop(columns=['compt'], inplace=True)
        
        theta_0=observation_df_1.loc['initial_value'].tolist()
        lb=observation_df_1.loc['lb'].tolist()
        ub=observation_df_1.loc['ub'].tolist()
        
        X_df = observation_df_1[observation_df_1.index.str.contains('X_obs')]
        Y_obs_df = observation_df_1[observation_df_1.index.str.contains('Y_obs')]
        
        
        if compt_name == 'Clarifier1':
            theta_0=theta_0[:5]
            lb=lb[:5]
            ub=ub[:5]
            
            X_df_1=X_df.iloc[:, :4]
            Y_obs_df_1=Y_obs_df.iloc[:, :4]
            
            
            
        elif compt_name == 'Reactor_Clarifier2':
            theta_0=theta_0[:6]
            lb=lb[:6]
            ub=ub[:6]
            
            X_df_1=X_df.iloc[:, :7]
            Y_obs_df_1=Y_obs_df.iloc[:, :7]
        
        theta_0=[float(i) for i in theta_0]
        lb=[float(i) for i in lb]
        ub=[float(i) for i in ub]
        
        
        X_array=X_df_1.values
        Y_obs_array=Y_obs_df_1.values
        
        X_array=X_array.astype(float)
        Y_obs_array=Y_obs_array.astype(float)
        
        res={'theta_0': theta_0, 'lb': lb, 'ub': ub, 'X': X_array, 'Y_obs': Y_obs_array}
    
           
    return res
    
    
        
    
    
    