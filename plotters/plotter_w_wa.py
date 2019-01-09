import os, sys
here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.normpath(os.path.join(here,'python/')))
import planckStyle as s
from pylab import *

#g = s.getSinglePlotter(plot_data='plot_data/')


import GetDistPlots



import planckStyle





g = planckStyle.getSinglePlotter(chain_dir = './chains', ratio=.7)

roots = ['EUmock_10bin_bfpk15_w0wa','pkTT+BSH_w0wa','pkTT+BSH_EUmock_10bin_bfpk15_w0wa']
g.settings.solid_contour_palefactor = 0.8
g.plot_2d(roots, 'w', 'wa', filled=True, colors=['navy','green','darkred','green'])#, lims=[0, 0.16, 0, 1.6])
g.add_legend(['Euclid WL', 'Planck+BSH','Planck+BSH+Euclid'], legend_loc='upper right', colored_text=True);
g.add_x_marker(-9.3111847E-01)
g.add_y_marker(-4.4578191E-01)

#plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
plt.show()
g.export('results_plots/w_wa.pdf')








