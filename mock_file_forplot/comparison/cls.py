# Configure Matplotlib options
from setup_matplotlib import *
from matplotlib.ticker import MaxNLocator
import pyfits



#load theoretical spectra
cls = np.loadtxt("comparison_cls_matteo.dat")
ste = np.loadtxt("stecls.dat")
win = np.loadtxt("wincls.dat")
# Create the plot


width = 10.
fig = plt.figure(figsize=(cm2inch(width), cm2inch(width*6/8.)))
# this should be changed for making a panel of multiple figures
ax = fig.add_subplot(111)
ax.set_yscale('log')
ax.set_xscale('log')

# power spectra
plt.plot(cls[:,0],cls[:,1], "k", lw=1, color='red', label=r"1-1")
plt.plot(cls[:,0],cls[:,2], "k", lw=1, color='green', label=r"2-2")
plt.plot(cls[:,0],cls[:,3], "k", lw=1, color='navy', label=r"1-2")

plt.plot(ste[:,0],ste[:,1], "k", lw=1, ls='--', color='red')
plt.plot(ste[:,0],ste[:,2], "k", lw=1, ls='--', color='green')
plt.plot(ste[:,0],ste[:,3], "k", lw=1, ls='--', color='navy')

plt.plot(win[:,0],win[:,1], "k", lw=1, ls=':', color='red')
plt.plot(win[:,0],win[:,2], "k", lw=1, ls=':', color='green')
plt.plot(win[:,0],win[:,3], "k", lw=1, ls=':', color='navy')

# legend
plt.legend(frameon=False,loc='lower right',ncol=1)

# axes labels
plt.xlabel(r"$\ell$"); plt.ylabel(r"$\ell(\ell+1) C_\ell/(2\pi)$")
ax.yaxis.labelpad = 10*width/17.; ax.xaxis.labelpad = 10*width/17. # distance of axis label to tick labels

    
# axes limits
#    plt.ylim([-1.e-5, 1.e5]); plt.xlim([1.8,2000]);
plt.xlim([10,3000]);
plt.ylim([1.e-7,5.e-4]);
# reduce white space around figure
plt.subplots_adjust(left=0.01, right=0.99, top=0.99, bottom=0.01)

# save to pdf with right bounding box
plt.savefig("cls_matteo.pdf", bbox_inches='tight', pad_inches=0.02)

plt.close()
