# -*- coding: utf-8 -*-
"""
Created on Mon May  1 09:39:29 2023

@author: quanta
"""
import numpy as np
#//////////////////////////////////////
# function, check whether already in PRB
#//////////////////////////////////////
#
def int2bin(x,ndig):
    # integer to n-digit binary conversion
    getbinary = lambda x, n: format(x, 'b').zfill(n)
    sb=getbinary(x, ndig)
    return sb
#--
#////////////////////////////////////////////////// 
# //import from propietary library from PRB folder
#////////////////////////////////////////////////// 
import sys
sys.path.append('../PRB')
#--
#from prb_qaoa import *
#from prb_obj_turyn import *
# -- variables
def b2s(x):
    if (x<1):
        s= 1.0
    if (x>0):
        s= -1.0
    return s
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
    s3=b2s(int(tv[3]))
    #s4=b2s(int(tv[4]))
    #s5=b2s(int(tv[5]))
    ##################################### 
    # ------------------------------------------------------------------------
    #obj=(s0+2*s1+3*s2 + 2*s0*s1 + 3*s1*s2 + 5*s1*s2*s3)
    #obj=(s0 + s1 + s2 + s0*s1 + s1*s2 + s1*s2*s3 + s0*s1*s2*s3 - 1)
    #obj=(s0 + 2*s1 + 3*s2 + 5*s0*s1 + 7*s1*s2 + 11*s1*s2*s3 + 13*s0*s1*s2*s3)
    #obj=(s0*s1 + 2*s0*s2 + 3*s0*s3 + 5*s1*s2 + 7*s1*s3+ 11*s2*s3-3)
    obj=(s0*s1 + s0*s2 + s0*s3 + s1*s2 + s1*s3+ s2*s3)
    # ------------------------------------------------------------------------
    obj=abs(obj) ## since we seek for a minimum
    return obj
#////////////////////////////////////////////////////
#/////////////////////////
if __name__ == '__main__':
    NQ=4 # number of qubits

    #print('Turyn>> k=', k, '; H-order->', M, '#qubits->', NQ)
    #--- check  --
    vMin=1e9
    vMax=0
    cntTOT=0
    cntSOL=0
    vObj=[]
    eTot=0
    for m in range(0,2**NQ):
        sb = int2bin(m,NQ)
        #//////////////////////////
        # objective function
        # order: 44, 68, 92, 116
        #//////////////////////////
        E=obj_H4(sb)
        vObj.append(E)
        eTot=eTot+E
        #//////////////////////////
        #=-----
        if E<vMin:
          vMin=E
        elif E>vMax:
          vMax=E
        #--
        if abs(E)<0.00001:
            cntSOL=cntSOL+1
            SOL=sb
            print('Solution Found ...!', E, SOL)
        #return obj #, SOL
        #--- endif
        cntTOT = cntTOT+1
    #---- end for
    #print('Turyn>> k=', k, '; H-order->', M, '#qubits->', NQ)
    print('Total bitstring= ', cntTOT, '>>', 'Min = ', vMin, 'Max = ', vMax)
    print('Total solution ', cntSOL, '>>', 100*cntSOL/cntTOT, ' %')
    vObj.sort()
    print('min Obj=', vObj[0], ';max Obj=', vObj[2**NQ-1], ';mean Obj=',eTot/cntTOT)
    print('first level >>', vObj[cntSOL])
    