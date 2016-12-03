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

pOA='MKKTAIAIAVALAGFATVAQAAPKDNTWYTGAKLGWSQYHDTGFINNNGPTHENQLGAGAFGGYQVNPYVGFEMGYDWLGRMPYKGSVENGAYKAQGVQLTAKLGYPITDDLDIYTRLGGMVWRADTKSNVYGKNHDTGVSPVFAGGVEYAITPEIATRLEYQWTNNIGDAHTIGTRPDNGMLSLGVSYRFGQGEAAPVVAPAPAPAPEVQTKHFTLKSDVLFNFNKATLKPEGQAALDQLYSQLSNLDPKDGSVVVLGYTDRIGSDAYNQGLSERRAQSVVDYLISKGIPADKISARGMGESNPVTGNTCDNVKQRAALIDCLAPDRRVEIEVKGIKDVVTQPQA'

#dict={'proOmpF': 'MMKRNILAVIVPALLVAGTANAAEIYNKDGNKVDLYGKAVGLHYFSKGNGENSYGGNGDMTYARLGFKGETQINSDLTGYGQWEYNFQGNNSEGADAQTGNKTRLAFAGLKYADVGSFDYGRNYGVVYDALGYTDMLPEFGGDTAYSDDFFVGRVGGVATYRNSNFFGLVDGLNFAVQYLGKNERDTARRSNGDGVGGSISYEYEGFGIVGAYGAADRTNLQEAQPLGNGKKAEQWATGLKYDANNIYLAANYGETRNATPITNKFTNTSGFANKTQDVLLVAQYQFDFGLRPSIAYTKSKAKDVEGIGDVDLVNYFEVGATYYFNKNMSTYVDYIINQIDSDNKLGVGSDDTVAVGIVYQF'}

#dict = {'proOmpC':'MKVKVLSLLVPALLVAGAANAAEVYNKDGNKLDLYGKVDGLHYFSDNKDVDGDQTYMRLGFKGETQVTDQLTGYGQWEYQIQGNSAENENNSWTRVAFAGLKFQDVGSFDYGRNYGVVYDVTSWTDVLPEFGGDTYGSDNFMQQRGNGFATYRNTDFFGLVDGLNFAVQYQGKNGNPSGEGFTSGVTNNGRDALRQNGDGVGGSITYDYEGFGIGGAISSSKRTDAQNTAAYIGNGDRAETYTGGLKYDANNIYLAAQYTQTYNATRVGSLGWANKAQNFEAVAQYQFDFGLRPSLAYLQSKGKNLGRGYDDEDILKYVDVGATYYFNKNMSTYVDYKINLLDDNQFTRDAGINTDNIVALGLVYQF'}

locFlex =[]
volume = []
area = []
polarity = []
hydrophilicity = []
hydropathy = []
for index, letter in enumerate(pOA):
    locFlex.append(aa_props.local_flexibility(letter))
    volume.append(aa_props.volume(letter))
    area.append(aa_props.accessible_surface_area(letter))
    polarity.append(aa_props.polarity(letter))
    hydrophilicity.append(aa_props.hydrophilicity(letter))
    hydropathy.append(aa_props.hydropathy(letter))

#plt.figure()
#plt.plot(volume,area,'g.'); plt.xlabel('Volume'); plt.ylabel('area')
#for index,letter in enumerate (pOA):
#    plt.text(volume[index],area[index],letter)
    
plt.figure()
plt.plot(locFlex,volume,'r.'); plt.xlabel('LocFlex'); plt.ylabel('Volume')
for index,letter in enumerate (pOA):
    plt.text(locFlex[index],volume[index],letter)
    
#plt.figure()
#plt.plot(locFlex,area,'b.'); plt.xlabel('LocFlex'); plt.ylabel('Area')
#for index,letter in enumerate (pOA):
#    plt.text(locFlex[index],area[index],letter)

plt.show()
plt.plot(locFlex,polarity,'b.'); plt.xlabel('LocFlex'); plt.ylabel('polarity')
for index,letter in enumerate (pOA):
    plt.text(locFlex[index],polarity[index],letter)
    
plt.show()
plt.plot(volume,polarity,'b.'); plt.xlabel('volume'); plt.ylabel('polarity')
for index,letter in enumerate (pOA):
    plt.text(volume[index],polarity[index],letter)

plt.show()
plt.plot(volume,hydrophilicity,'b.'); plt.xlabel('Volume'); plt.ylabel('hydrophilicity')
for index,letter in enumerate (pOA):
    plt.text(volume[index],hydrophilicity[index],letter)
    
#plt.show()
#plt.plot(volume,polarity,'b.'); plt.xlabel('volume'); plt.ylabel('polarity')
#plt.show()
#plt.plot(volume,hydrophilicity,'b.'); plt.xlabel('volume'); plt.ylabel('hydrophilicity')
#plt.show()
#plt.plot(polarity,hydrophilicity,'b.'); plt.xlabel('polarity'); plt.ylabel('hydrophilicity')
#plt.show()
#plt.plot(hydropathy,hydrophilicity,'b.'); plt.xlabel('hydropathy'); plt.ylabel('hydrophilicity')
#plt.show()
#plt.plot(hydropathy,locFlex,'b.'); plt.xlabel('hydropathy'); plt.ylabel('locFlex')
#plt.show()