# -*- coding: utf-8 -*-
"""
Created on Fri Dec  2 22:29:54 2016

@author: tomi
"""

#beads on string
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

aas = np.array(range(700))


def rand_skew_norm(fAlpha, fLocation, fScale):
    sigma = fAlpha / np.sqrt(1.0 + fAlpha**2) 

    afRN = np.random.randn(2)
    u0 = afRN[0]
    v = afRN[1]
    u1 = sigma*u0 + np.sqrt(1.0 -sigma**2) * v 

    if u0 >= 0:
        return u1*fScale + fLocation 
    return (-u1)*fScale + fLocation 

def randn_skew(N, skew=0.0):
    return [rand_skew_norm(skew, 0, 1) for x in range(N)]


def first_valve(aa):
    if aa in { 'W', 'Y', 'F', 'K', 'R' } : 
        alpha_skew = 3
    else:
        alpha_skew = 0
    
    step = randn_skew(1, alpha_skew)
#    print (step)
    return step[0]

"""
NUM_SAMPLES = 100000
SKEW_PARAMS = [-5, 0, 5]

plt.subplots(figsize=(7,4))
for alpha_skew in SKEW_PARAMS:
    p = randn_skew(NUM_SAMPLES, alpha_skew)
    sns.distplot(p)
"""    
    
position = 86 #go to the 86th aa ... this is first residue outside of secyeg cplx ... initiation will be solved later
positions1v = []
low, high = 0, len(aas)
counter = 0

while 1:
    counter += 1
    position = position + first_valve('W')
    positions1v.append(position)
#    print(position)
    
    if position > high:
        print ('Translocation has finished in ' + str(counter) + 'steps')
        break
    if position < low:
        print ('Peptide escaped in ' + str(counter) + 'steps')
        break


plt.plot(positions1v)
plt.show()
sns.distplot(np.diff(positions1v))
plt.show()







