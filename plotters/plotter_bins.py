# Configure Matplotlib options
from setup_matplotlib import *
from matplotlib.ticker import MaxNLocator

# Load data

#load theoretical spectra
total = np.loadtxt("./mock_file_forplot/EUmock_10bin_bfpk15_window_total.dat") # l, C_l
bins = np.loadtxt("./mock_output/WL_mock/EUmock_10bin_bfpk15_window.dat") # l, C_l

# Create the plot


for width in [10.]:
    fig = plt.figure(figsize=(cm2inch(width), cm2inch(width*6/8.)))
    # this should be changed for making a panel of multiple figures
    ax = fig.add_subplot(111)
#    ax.set_xscale('log')
#    ax.set_yscale('log')

# power spectra
    plt.plot(total[:,0],total[:,1], "k", lw=1, ls='-',color='black',label="n(z)")
    plt.plot(bins[:,0],bins[:,1], "k", lw=1, ls='-',color='blue')
    plt.plot(bins[:,0],bins[:,2], "k", lw=1, ls='-',color='red')
    plt.plot(bins[:,0],bins[:,3], "k", lw=1, ls='-',color='green')
    plt.plot(bins[:,0],bins[:,4], "k", lw=1, ls='-',color='navy')
    plt.plot(bins[:,0],bins[:,5], "k", lw=1, ls='-',color='purple')
    plt.plot(bins[:,0],bins[:,6], "k", lw=1, ls='--',color='blue')
    plt.plot(bins[:,0],bins[:,7], "k", lw=1, ls='--',color='red')
    plt.plot(bins[:,0],bins[:,8], "k", lw=1, ls='--',color='green')
    plt.plot(bins[:,0],bins[:,9], "k", lw=1, ls='--',color='navy')
    plt.plot(bins[:,0],bins[:,10], "k", lw=1, ls='--',color='purple')



    # x axis
    plt.hlines(0, 0, 3300)

#    plt.title(r'variable $h^{cut}$, $r_{T/S}=0.05$, $\beta_1=0.48$, $\beta_4=0.94$', fontsize=10)

    # legend
    plt.legend(frameon=False,loc='upper right',ncol=1)

    # axes labels
    plt.xlabel(r"$z$"); plt.ylabel(r"$n(z)$")
    ax.yaxis.labelpad = 10*width/17.; ax.xaxis.labelpad = 10*width/17. # distance of axis label to tick labels

    
    # axes limits
    plt.xlim([0,2.5]); 
#plt.ylim([02,1000]);

    # reduce white space around figure
    plt.subplots_adjust(left=0.01, right=0.99, top=0.99, bottom=0.01)

    # set vertical y axis ticklables
#    for ticklabel in ax.yaxis.get_ticklabels():
#        ticklabel.set_rotation("vertical")

    # save to pdf with right bounding box
    plt.savefig("bins.pdf", bbox_inches='tight', pad_inches=0.02)
