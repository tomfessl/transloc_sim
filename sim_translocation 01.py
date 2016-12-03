# -*- coding: utf-8 -*-
"""
Created on Fri Dec  2 22:29:54 2016
@author: tomi
"""
from __future__ import division

def fig_seq_positions(pOA):
    dict={'proOmpA': pOA}
    aas = {'F','W','Y','R'}
    for motif in aas:
        counter=0
        for u, seq in dict.iteritems():
            i=-1 #start search at the beginning of the sequence=
            while True:
                try:
                    i= seq.index(motif, i+1) #get the index of the next occurrence
                    counter +=1
                    if motif == 'F': color = 'r'; ms = 10
                    if motif == 'W': color = 'g'; ms = 20
                    if motif == 'Y': color = 'b'; ms = 10
                    if motif == 'R': color = 'y'; ms = 10    
                    
                    plt.plot(-0.05,i,c=color,markersize=ms,marker='o',alpha = 0.5)                
#                    print "%s has been found in %d position of %s, %s"%(motif, i+1, u, counter)  
                except ValueError:
                    break #no more motifs found

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
#    if aa in { 'W', 'Y', 'F', 'K', 'R' } : 
#        alpha_skew = 3
    """
        locFlex.append(aa_props.local_flexibility(letter))
        volume.append(aa_props.volume(letter))
        area.append(aa_props.accessible_surface_area(letter))
        polarity.append(aa_props.polarity(letter))
        hydrophilicity.append(aa_props.hydrophilicity(letter))
        hydropathy.append(aa_props.hydropathy(letter))
    """
    if aa_props.volume(aa) >200: 
        alpha_skew = skew
    else:
        alpha_skew = 0
    
    step = randn_skew(1, alpha_skew)
#    print (step)
    return 0.2*step[0]

def sequence():
    pOA='MKKTAIAIAVALAGFATVAQAAPKDNTWYTGAKLGWSQYHDTGFINNNGPTHENQLGAGAFGGYQVNPYVGFEMGYDWLGRMPYKGSVENGAYKAQGVQLTAKLGYPITDDLDIYTRLGGMVWRADTKSNVYGKNHDTGVSPVFAGGVEYAITPEIATRLEYQWTNNIGDAHTIGTRPDNGMLSLGVSYRFGQGEAAPVVAPAPAPAPEVQTKHFTLKSDVLFNFNKATLKPEGQAALDQLYSQLSNLDPKDGSVVVLGYTDRIGSDAYNQGLSERRAQSVVDYLISKGIPADKISARGMGESNPVTGNTCDNVKQRAALIDCLAPDRRVEIEVKGIKDVVTQPQA'
    return pOA
    
def one_valve_process(pOA):    
    position = 10 #go to the 86th aa ... this is first residue outside of secyeg cplx ... initiation will be solved later
    positions1v = []
    low, high = 0, len(pOA)
    counter = 0
    while 1:
        counter += 1
        
        if position >= high:
            print ('Translocation has finished in ' + str(counter) + 'steps')
            break
        if position <= low:
            print ('Peptide escaped in ' + str(counter) + 'steps')
            break
        step = first_valve(pOA[int(position)])
        position = position + step
        positions1v.append(position)
    plot_rates(positions1v,pOA)
        
def plot_rates(positions1v,pOA):  
    length = len(pOA)
    steps = np.diff(positions1v)
    all_steps = []
    for i in range(length):
        try:
            nums = steps[np.where(np.round(positions1v) == i)[0]]
        except:
            pass
        all_steps.append(np.nanmean(nums))
    every_step = np.array(all_steps)
    
    N = 100
    x = range(len(positions1v))
    m_x = np.convolve(x, np.ones((N,))/N, mode='valid')  
    mov_av = np.convolve(positions1v, np.ones((N,))/N, mode='valid')  
    plt.figure(figsize=(12,7))
    plt.subplot(141)
    plt.ylim([0, len(pOA)])
    fig_seq_positions(pOA)
    plt.plot(m_x,mov_av,'r')
    plt.plot(x,positions1v)
    plt.xlim([-100, len(positions1v)])
    plt.xlabel('# of steps')
    plt.ylabel('position in sequence')

    plt.subplot(142)
    every_step[np.isnan(every_step)]=0
    plt.plot(every_step,range(len(every_step)), alpha = 0.5)
    M=5
    s_x = np.convolve(range(len(every_step)), np.ones((M,))/M, mode='valid')  
    s_y = np.convolve(every_step, np.ones((M,))/M, mode='valid')  
    plt.plot(s_y,s_x,'r')
    plt.ylim([0, len(pOA)])
    hist, bin_edge = np.histogram(positions1v,range(len(pOA)))
    rate = 1./hist
    rate[np.isinf(rate)]=0
#    
#    histA, bin_edgeA = np.histogram(mov_av,range(len(pOA)))
#    rateA = 1./histA
#    rateA[np.isinf(rateA)]=0
#    
    fig_seq_positions(pOA)
    plt.subplot(143)
    fig_seq_positions(pOA)
    plt.plot(rate, bin_edge[:-1],'b',alpha = 0.25)
#    plt.plot(rateA, bin_edgeA[:-1],'r', alpha = 0.5)
    N = 3    
    m_xr = np.convolve(bin_edge[:-1], np.ones((N,))/N, mode='valid')  
    m_yr = np.convolve(rate, np.ones((N,))/N, mode='valid')     
    plt.plot(m_yr, m_xr,'g', alpha = 0.5)
    plt.xlabel('|dx|/dt')
    plt.ylim([0, len(pOA)])
    plt.xlim([-0.1,1.1*np.max(rate[np.isfinite(rate)])])

    plt.subplot(248)
    for alpha_skew in [0,skew]:
        p = randn_skew(int(1e5), alpha_skew)
        sns.distplot(p)
    plt.xlabel('step size')
    plt.legend(['unbiased','biased'])    

if __name__ == "__main__":    
    import numpy as np
    import matplotlib.pyplot as plt
    import seaborn as sns
    import amino_acid as aa_props
    global skew
    skew = 30
    pOA = sequence()
    one_valve_process(pOA)
    plt.show()
