import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import colorsys

from mpl_toolkits.axes_grid1.inset_locator import inset_axes

matplotlib.rc('figure', figsize=(3.15, 3.15), dpi=300)
matplotlib.rc('figure.subplot', left=0.06, bottom=0.05,
              right=0.98, top=0.96, wspace=0.15, hspace=0.15)
matplotlib.rc('font', family='Times New Roman', size=10)
#matplotlib.rc('axes', labelsize=10, labelpad=2)
matplotlib.rc('xtick', labelsize=8)
matplotlib.rc('ytick', labelsize=8)
matplotlib.rc('lines', linewidth=1, markersize=1.0)
matplotlib.rc('image', cmap='jet')
matplotlib.rc('xtick', direction='in')
matplotlib.rc('ytick', direction='in')

def plot_reconstruct_waveform (ax, srcDir, coef=0.011, diff=0.):

    axin = inset_axes(ax, width="40%", height="60%")
    axin.tick_params(labelleft=False)

    data = np.loadtxt (srcDir + 'res/field_50.dat')
    time_field, Ef_field, Af_field = data[:,0], data[:,1], data[:,2]
    axin.plot (time_field, 0.0001 * Ef_field, lw=1)

    data = np.loadtxt (srcDir + 'field_101.dat')
    time_field, Ef_field, Af_field = data[:,0], data[:,1], data[:,2]
    axin.plot (time_field, Ef_field, lw=1)
    ax.plot (time_field, np.abs (Ef_field), lw=1.5)

    data = np.loadtxt (srcDir + 'reconstruct.dat')
    ax.plot (data[:,0], coef*(10**data[:,1]-diff), 'rx', markersize=4)




from matplotlib.gridspec import GridSpec
#fig = plt.figure (constrained_layout=True)
fig = plt.figure ()
#gs = GridSpec (2, 1, figure=fig, hspace=0.1, wspace=0.1)
#gs.update(left=0.1,right=0.9,top=0.985,bottom=0.04,wspace=0.05,hspace=0.05)

""" ax1 = plt.subplot (gs[0,0])
ax2 = plt.subplot (gs[1,0]) """
ax1 = plt.subplot (211)
ax2 = plt.subplot (212)

plot_reconstruct_waveform (ax1, '../THz-reconstruct-scheme/', coef=3.2e-2, diff=5.8e-5)
plot_reconstruct_waveform (ax2, '../THz-reconstruct-params/', coef=0.011)

plt.show ()