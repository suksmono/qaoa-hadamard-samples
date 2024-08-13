# -*- coding: utf-8 -*-
"""
Created on Sun Apr 16 15:00:17 2023

@author: quanta
test coeffs setting on 4-order H matrix
H4=[s1,  1,  1,  1;
     1, s2,  1, -1;
     1,  1, s3, -1;
     1, -1, -1, s4]
correct solution: [1, -1, -1, 1] or "0110"
"""

##
from sympy import *
#
from qiskit import QuantumCircuit, ClassicalRegister, QuantumRegister
from qiskit import Aer, execute
from qiskit.circuit import Parameter
from scipy.optimize import minimize
from qiskit.visualization import plot_histogram
# ///////////////////
N_LAYER= 1 #6
##
#-------------------------------------
# convert binary {0,1} -> spin {-1, +1}
# -- variables
def b2s(x):
    if (x<1):
        s= 1.0
    if (x>0):
        s= -1.0
    return s
#-------------------------------------
############################################################
# Evaluasi fungsi secara simbolik
# Waktu eksekusi lama, kurang cocok untuk orde tinggi
############################################################
def obj_si(tv, Es):
    # evaluate objective function for a string tv
    # Es=sympify(E12)
    # //// possibly need bit-> spin variable conversion
    NQ=len(tv)
    ss=symbols('s0:%d'%NQ)
    #
    for n in range(0,NQ):
        Es=Es.subs(ss[n],tv[n])
    # --
    obj=Es.evalf()
    return obj
#--

# //// fungsi objective untuk kasus matriks Hadamard
# ///  di sini hanya untuk H-68
def hadamard_obj(vs, Es):
    # value of objective function
    # TO BE MODIFIED WITH HADAMARD MATRIX ??
    ##############################################
    """
    Given a bitstring as a solution, this function returns
    the number of edges shared between the two partitions
    of the graph.   
    Args:
        x: str
           solution bitstring
    Returns:
        obj: float
             Objective
    """
    ##///////////////////////
    obj=obj_H4(vs)
    #gIter=gIter+1
    F_Found=0
    if abs(obj)<0.00001:
        F_Found=1
        GEN_SOL.append(vs)
        print('Solution Found ...! Iter->', len(vObj), 
              'E->', obj, 'Sol->', vs)
    return obj, F_Found #, SOL

# /////////////////////////////////////////////////////
# hitung rata-rata energi dari semua solusi
# diam-diam diasumsikan, dengan menurunkan rata-rata
# energi semua solusi, maka solusi yang diinginkan
# atau E=0 dapat dicapai -->> belum tentu benar (?)
# /////////////////////////////////////////////////////
def compute_expectation(counts):  
    """
    Computes expectation value based on measurement results   
    Args:
        counts: dict
                key as bitstring, val as count                  
    Returns:
        avg: float
             expectation value
    """
    
    avg = 0
    sum_count = 0
    cnt_corr_sol=0
    for bitstring, count in counts.items():     
        #obj, SOL = hadamard_obj(bitstring)
        ### print('bitstring->', bitstring, 'Es>>', Es)
        obj, F_Found = hadamard_obj(bitstring, Es)
        if F_Found>0:
            cnt_corr_sol=cnt_corr_sol+count
            #print('Correct Sol->', count)
        avg += obj * count
        sum_count += count
    #print('Total Sample->',sum_count, ', correct sol->', cnt_corr_sol)
    ##
    vShots.append(sum_count)
    vCorrect.append(cnt_corr_sol)
    ##
    return avg/sum_count


# Finally we write a function that executes the circuit on the chosen backend
def get_expectation(p, NPARAMS, shots=1*512):  
    """
    Runs parametrized circuit
    Args:
        p: int,
           Number of repetitions of unitaries
    """
    
    backend = Aer.get_backend('qasm_simulator')
    backend.shots = shots
    #gIter=0
    #print('NPARAMS>>>',NPARAMS)
    def execute_circ(theta):
        
        # benchmark-> 
        #qc=create_qaoa_circ_44(NQ, theta)
        #/// sudah oke untuk order-44
        #qc=create_qaoa_circ_from_terms(qx,nqubits,theta)
        #qc=create_qaoa_circ_from_terms_multiparams(qx,nqubits,theta)
        qc=create_qaoa_circ_from_terms_multiparams_multilayers(qx, nqubits, NPARAMS, theta)
        #--
        #print(qc.draw())
        # --
        counts = backend.run(qc, seed_simulator=10, 
                             nshots=1*512).result().get_counts()
        mean_obj=compute_expectation(counts)
        vObj.append(mean_obj)
        #return compute_expectation(counts)
        return mean_obj
    #    
    return execute_circ
#---
#////////////////////////////////////////////////////
# order: H4
#////////////////////////////////////////////////////
def obj_H4(tv):
    # calculate the value of objective function
    # i.e. Hamiltonian's energy of a solution/bitstring
    # order 44--> 5 qubits
    # /////////// harusnya konversi [0,1] -> [-1, +1] ...??????
    ## >> revised
    ##-H44->> NQ=5 (non-hybrid)
    s0=b2s(int(tv[0]))
    s1=b2s(int(tv[1]))
    s2=b2s(int(tv[2]))
    #s3=b2s(int(tv[3]))
    #s4=b2s(int(tv[4]))
    ##################################### 
    # ------------------------------------------------------------------------
    #obj= (6*s0*s1 + 6*s0*s2 - 6*s0*s3 + 12*s0 + 2*s1*s2 - 2*s1*s3 + 4*s1 - 2*s2*s3 + 4*s2 - 4*s3 + 16)
    #obj=(24*s0*s1 + 24*s0*s2 - 24*s0*s3 + 8*s1*s2 - 8*s1*s3 - 8*s2*s3 + 48)
    #obj=(8*s0*s1 + 8*s0*s2*s3 + 24*s0 + 8*s1*s2*s3 + 24*s1 + 24*s2*s3 + 48)    
    obj=(8*s0*s1 - 8*s0*s2 + 24*s0 - 8*s1*s2 + 24*s1 - 24*s2 + 48)
    # ------------------------------------------------------------------------
    obj=abs(obj) ## since we seek for a minimum
    return obj
#////////////////////////////////////////////////////
#--
#////////////////////////////////////////////////////////////
#"""
import matplotlib.pyplot as plt
from collections import Counter
import random

#////////////////////////////////////////////////// 
# //import from propietary library from PRB folder
#////////////////////////////////////////////////// 
import sys
sys.path.append('../_turyn/PRB')
#--
from prb_qaoa import *
from prb_obj_turyn import *
#from prb_qaoa_auto_qcirc_cascaded import *
from prb_qaoa_auto_qcirc import *

#--
## define and initialize global variables
N_LAYER= 1 #6
NQ=21
nqubits=NQ
#H4 ="6*s0*s1 + 6*s0*s2 - 6*s0*s3 + 12*s0 + 2*s1*s2 - 2*s1*s3 + 4*s1 - 2*s2*s3 + 4*s2 - 4*s3 + 16"
#H4='24*s0*s1 + 24*s0*s2 - 24*s0*s3 + 8*s1*s2 - 8*s1*s3 - 8*s2*s3 + 48'
#H4='8*s0*s1 + 8*s0*s2*s3 + 24*s0 + 8*s1*s2*s3 + 24*s1 + 24*s2*s3 + 48'
H4='8*s0*s1 - 8*s0*s2 + 24*s0 - 8*s1*s2 + 24*s1 - 24*s2 + 48'
##################
vObj=[]
vShots=[]
vCorrect=[]
GEN_SOL=[]
if __name__ == '__main__':
    #fpath='../hamiltonians-williamson/'
    gIter=0
    # //// hanya order 68 -> N=6
    NQ=3  # number of qubits
    #print('order->', M, 'nqubits->', NQ)
    nqubits=NQ
    #-- text expression of H68
    Etxt=H4
    Es=sympify(Etxt) 
    ## 
    qx=efunct_2_qidx(Es)
    # determine n_gamma
    #//////////////
    n_gamma=0
    for n in range(0, len(qx)):
        kbody=len(qx[n][1])
        #print('k-body->', kbody)
        if kbody>0:
            n_gamma=n_gamma+1
        #
    #
    #print(n_gamma)
    n_beta=NQ
    #//////////////////////////////
    ##--create quantum circuits from terms
    #//// UBAH BAGIAN get_expectation
    ###############################
    # ///// Number of Layers /////
    NLAYER=NQ#NQ# NQ#NQ # int(1*NQ*NQ) #NQ*4 #NQ*5 #1 #10 #upto 10 layer is ok for orde-68
    ###############################
    # //// remark/remove seed for repeated experiments
    # //////////////////////////////////////
    RSEED=29 #2,3,5,7,11,13,17,19,23,29>>NQ
    # //////////////////////////////////////
    random.seed(RSEED)
    #theta0=[]
    #for i in range (0,2*NLAYER):
    #  theta0.append(random.uniform(-1.0*3.14, 1.0*3.14))
    theta0=[]
    NPARAMS=(n_gamma+n_beta) #*NLAYER
    for i in range (0,NPARAMS*NLAYER):
      theta0.append(random.uniform(-1.0, 1.0))
  
    # salah satu pilihan nilai terbaik
    # theta0=[0.1, -0.1]
    #
    #expectation = get_expectation(p=1)
    expectation = get_expectation(NLAYER, NPARAMS)
    ## minimize the expectation
    #res = minimize(expectation, 
    #                  theta0, 
    #                  method='COBYLA')
    
    res = minimize(expectation, 
                      theta0, 
                      method='COBYLA',
                      options={'xatol': 1e-18, 'disp': True}
                      )   

    print(res)
    print('/////////////////////////////////////////////////')
    print('N-LAYER = ', NLAYER)
    
    #///// RESULTS DURING ITERATION /////
    print('\n** results during iteration**')
    ttl_shots=0
    ttl_csol=0
    for m in range(0, len(vShots)):
        ttl_shots= ttl_shots + vShots[m]
        ttl_csol = ttl_csol + vCorrect[m]
    #///
    print('Total Shots->', ttl_shots, 'Total Correct->', ttl_csol)
    print('Prob of random algorithm (before acc all sol) >> ', 1/2**NQ)    
    print('Prob of QAOA>>', ttl_csol/ttl_shots)    
            
    print('\n **results after iteration**')
    #-- DISPLAY RESULTS    
    #backend = Aer.get_backend('aer_simulator')
    backend = Aer.get_backend('qasm_simulator')
    backend.shots = 4*512
    # --
    #qc_res=create_qaoa_circ_from_terms(qx,nqubits,res.x)
    qc_res=create_qaoa_circ_from_terms_multiparams_multilayers(qx, nqubits, NPARAMS, res.x)
    #--    
    counts = backend.run(qc_res, seed_simulator=RSEED).result().get_counts()
    #///
    print('Total distinct solutions->', len(counts))
    #--- cek hadamard
    # /// from a set of candidate
    cntr=Counter(counts)
       
    # //////////////////////////////////////////
    # ---analisis---
    # jumlah solusi dari NSOL- frekuensi tertinggi
    #//// selected solutions: most frequent ///
    NSOL= len(counts) # int(len(counts)/2.0) #20*25 # 10 #int(2**NQ/10)
    vbs=cntr.most_common(NSOL) # get solution candidates
    # /// scan all solution candidates
    cnt_sol=0  # solution counter
    cnt_freq=0 # frequency of solution
    TOT_SOL=0
    min_obj=1000*1000
    for m in range(0,NSOL): 
        tbs=vbs[m][0] # bitstring of solutions
        #---
        tfr=vbs[m][1] # frequency
        TOT_SOL=TOT_SOL+tfr
        vobj=obj_H4(tbs) #     obj_si_bs
        #---
        #print(tbs,':',tfr, '-> obj:', vobj)
        if vobj<min_obj:
          min_obj=vobj
          min_sol=tbs
        # ////////////
        #print(tbs,':',tfr, '-> obj:', vobj, 'min_obj->', min_obj) 
        # ////////////
        #
        if vobj<1:
          cnt_sol=+1
          cnt_freq=cnt_freq+tfr
          sol=tbs
    #--
    if cnt_sol<1:
      print('NO EXACT solution among the most frequent')
      #sol=tbs
      sol=min_sol
    else:
      print('Number of correct solution found->', cnt_freq, 'among', TOT_SOL)
    #--
    ##/// solution from GEN_SOL
    if len(GEN_SOL)>0:
        bs=GEN_SOL[0]
    else:
        bs=sol
    #--
    #bs='101110001101001111010'
    print('tbs->', bs)

    #////////////////
    # plot learning curve
    plt.plot(vObj)
    plt.xlabel("iteration")
    plt.ylabel("mean error")
    #---
    #/// plot histogram
    #plot_histogram(counts)