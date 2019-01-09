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
g.plot_2d(roots, 'omegam', 'H0', filled=True, colors=['navy','green','darkred','green'])#, lims=[0, 0.16, 0, 1.6])
g.add_legend(['Euclid WL', 'Planck+BSH','Planck+BSH+Euclid'], legend_loc='upper right', colored_text=True);
g.add_x_marker(0.31)
g.add_y_marker(68.40)

#plt.xticks([0.3, 0.31,0.32,0.33,0.34])
#plt.yticks([0,0.0004])
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
plt.show()
g.export('results_plots/om_H0.pdf')








