import sys
from astroquery.vizier import Vizier
from astropy import units as u
from astropy.coordinates import SkyCoord
import numpy as np
from math import *
import matplotlib as mpl
import matplotlib.pyplot as plt
plt.rc('font', family='serif')
mpl.rcParams.update({'font.size': 12})
mpl.rcParams.update({'legend.labelspacing':0.25, 'legend.fontsize': 12})
mpl.rcParams.update({'errorbar.capsize': 4})


Vizier.ROW_LIMIT = -1

catalogs = Vizier.get_catalogs('J/ApJS/247/66')

print (catalogs['J/ApJS/247/66/table6'])

#print (len(catalogs['J/ApJS/247/66/table6']))
#sys.exit(0)

first  = catalogs['J/ApJS/247/66/table4']
second = catalogs['J/ApJS/247/66/table5'] 
third  = catalogs['J/ApJS/247/66/table6']

inc = []
cos_inc = []
psep = []

psep1e4 = []
cos_i1e4 = []

sf = open ('cos_inc.txt', 'w')

for i in range (0, len(catalogs['J/ApJS/247/66/table6'])):

    if i % 100 == 0:
        print (i)
    #print ('Coordinates are: ', first['RA_ICRS'][i], first['DE_ICRS'][i])
    f = SkyCoord(ra=first['RA_ICRS'][i]*u.degree, dec=first['DE_ICRS'][i]*u.degree, frame='icrs')
    fg = f.transform_to('galactic')
    fg_l = fg.l
    fg_b = fg.b
    #print (fg_l)
    s = SkyCoord(ra=second['RA_ICRS'][i]*u.degree, dec=second['DE_ICRS'][i]*u.degree, frame='icrs')
    sg = s.transform_to('galactic')
    sg_l = sg.l
    sg_b = sg.b
    diff_l = sg_l - fg_l
    diff_b = sg_b - fg_b
    #print ('Coordinate differences are: ', diff_l.degree, diff_b.degree)
    inc_v = atan2 (diff_b.degree, diff_l.degree)
    inc.append (inc_v / pi * 180.0)
    cos_inc.append (cos(inc_v))
    psep.append (third['PSep'][i])
    sf.write (str(cos_inc[-1])+'\t' + str(psep[-1])+'\n')

    if psep[-1] > 1e4:
        psep1e4.append (psep[-1])
        cos_i1e4.append (cos_inc[-1])


#psep1e4 = [x for i in range (0, len(psep)) if psep[i] > 1e4: x = psep[i]]


plt.hist (cos_inc, 20)
plt.xlabel(r'$\cos i$')
plt.ylabel('Number of pairs')
plt.savefig('inc.pdf')
plt.show()

plt.scatter(cos_inc, psep)
plt.yscale('log')
plt.xlabel(r'$\cos i$')
plt.ylabel(r'Physical separation (au)')
plt.savefig('psep.pdf')
plt.show()

plt.hist(cos_i1e4, 20)
plt.xlabel('r$\cos i$')
plt.ylabel('Number of pairs')
plt.savefig('inc_1e4.pdf')
plt.show()
