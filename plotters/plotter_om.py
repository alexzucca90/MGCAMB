import os, sys
here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.normpath(os.path.join(here,'python/')))
import planckStyle as s
from pylab import *

#g = s.getSinglePlotter(plot_data='plot_data/')


import GetDistPlots



import planckStyle

g = planckStyle.getSinglePlotter(chain_dir = './chains', ratio=1)

roots = ['test_planck','test']
g.settings.solid_contour_palefactor = 0.8
g.plot_1d(roots, 'omegam', ls=['-','--','-.','-'], colors=['black','green','navy','darkred'])#, lims=[0, 0.16, 0, 1])
first = g.add_legend(['Planck','Planck+Euclid'],legend_loc='upper right',fontsize=6)

g.export('results_plots/omegam.pdf')






