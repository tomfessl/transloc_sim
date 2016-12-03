# -*- coding: utf-8 -*-
"""
Created on Fri Jun 03 13:32:55 2016

@author: fbstf
"""
import numpy as np
import matplotlib.pyplot as plt
import amino_acid as aa_props


"""
Amino Acids with Hydrophobic Side Chain - Aromatic
Phenylalanine, Phe, F ... a van der Waals volume of 135 A3, and an accessible surface area of 218 A2
Tryptophan, Trp, W  ... a van der Waals volume of 163 A3, and an accessible surface area of 255 A2
Tyrosine, Tyr, Y ... a van der Waals volume of 105 A3, and an accessible surface area of 230 A2
SecYEG channel + SecA ~ 40AA length
"""

dict={'proOmpA': 'MKKTAIAIAVALAGFATVAQAAPKDNTWYTGAKLGWSQYHDTGFINNNGPTHENQLGAGAFGGYQVNPYVGFEMGYDWLGRMPYKGSVENGAYKAQGVQLTAKLGYPITDDLDIYTRLGGMVWRADTKSNVYGKNHDTGVSPVFAGGVEYAITPEIATRLEYQWTNNIGDAHTIGTRPDNGMLSLGVSYRFGQGEAAPVVAPAPAPAPEVQTKHFTLKSDVLFNFNKATLKPEGQAALDQLYSQLSNLDPKDGSVVVLGYTDRIGSDAYNQGLSERRAQSVVDYLISKGIPADKISARGMGESNPVTGNTCDNVKQRAALIDCLAPDRRVEIEVKGIKDVVTQPQA'}

#dict={'proOmpF': 'MMKRNILAVIVPALLVAGTANAAEIYNKDGNKVDLYGKAVGLHYFSKGNGENSYGGNGDMTYARLGFKGETQINSDLTGYGQWEYNFQGNNSEGADAQTGNKTRLAFAGLKYADVGSFDYGRNYGVVYDALGYTDMLPEFGGDTAYSDDFFVGRVGGVATYRNSNFFGLVDGLNFAVQYLGKNERDTARRSNGDGVGGSISYEYEGFGIVGAYGAADRTNLQEAQPLGNGKKAEQWATGLKYDANNIYLAANYGETRNATPITNKFTNTSGFANKTQDVLLVAQYQFDFGLRPSIAYTKSKAKDVEGIGDVDLVNYFEVGATYYFNKNMSTYVDYIINQIDSDNKLGVGSDDTVAVGIVYQF'}

#dict = {'proOmpC':'MKVKVLSLLVPALLVAGAANAAEVYNKDGNKLDLYGKVDGLHYFSDNKDVDGDQTYMRLGFKGETQVTDQLTGYGQWEYQIQGNSAENENNSWTRVAFAGLKFQDVGSFDYGRNYGVVYDVTSWTDVLPEFGGDTYGSDNFMQQRGNGFATYRNTDFFGLVDGLNFAVQYQGKNGNPSGEGFTSGVTNNGRDALRQNGDGVGGSITYDYEGFGIGGAISSSKRTDAQNTAAYIGNGDRAETYTGGLKYDANNIYLAAQYTQTYNATRVGSLGWANKAQNFEAVAQYQFDFGLRPSLAYLQSKGKNLGRGYDDEDILKYVDVGATYYFNKNMSTYVDYKINLLDDNQFTRDAGINTDNIVALGLVYQF'}

seqlen = len(dict['proOmpA'])
plotl = np.zeros(seqlen)

#i=0
#for u, aa in dict.iteritems():
#    plt.figure(1) 
#    i=i+1               
#    plt.plot(i,aa_props.local_flexibility(aa))
#    accessible_surface_area
#    hydrophilicity
#    volume
#            plt.figure(1) 
#        i=i+1               
#        plt.plot(i,aa_props.local_flexibility(aa))
#plt.show()


aas = {'F','W','Y'}
plt.figure(1)
for motif in aas:
    counter=0
    for u, seq in dict.iteritems():
        i=-1 #start search at the beginning of the sequence=
        while True:
            try:
                i= seq.index(motif, i+1) #get the index of the next occurrence
                counter +=1
                plotl[i]+=1
                if motif == 'F': color = 'r'; ms = 5
                if motif == 'W': color = 'g'; ms = 10
                if motif == 'Y': color = 'b'; ms = 15
                
                plt.plot(i,1,c=color,markersize=ms,marker='o',alpha = 0.5)                
                print "%s has been found in %d position of %s, %s"%(motif, i+1, u, counter)  
            except ValueError:
                break #no more motifs found

#plt.figure(2)
#plt.plot(plotl,'ro'); plt.ylim([0.5, 1.5])

hitind = [j for j, e in enumerate(plotl) if e != 0]
hitdiff = np.diff(hitind)
plt.figure(2)
plt.plot (hitdiff,'go-')
plt.figure(3)
plt.hist(hitdiff)