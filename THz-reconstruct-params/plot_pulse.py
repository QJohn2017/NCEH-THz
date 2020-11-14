import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import colorsys

matplotlib.rc('figure', figsize=(3.15, 2.17), dpi=300)
matplotlib.rc('figure.subplot', left=0.15, bottom=0.15,
              right=0.96, top=0.96, wspace=0.15, hspace=0)
matplotlib.rc('font', family='Times New Roman', size=10)
matplotlib.rc('axes', labelsize=10, labelpad=2)
matplotlib.rc('xtick', labelsize=8)
matplotlib.rc('ytick', labelsize=8)
matplotlib.rc('lines', linewidth=1, markersize=1.0)
matplotlib.rc('image', cmap='jet')
#matplotlib.rc('text', usetex=True)
#matplotlib.rc('text.latex', unicode=True, preamble=['\usepackage{amsmath}'])

# Ef(t)
def plot_Ef (ax, coef=1):
    ax.plot(time_field, coef*Ef_field)
    #ax.tick_params(labelbottom=False)
    ax.grid()

from matplotlib.gridspec import GridSpec
fig = plt.figure (constrained_layout=True)
gs = GridSpec (1, 3, figure=fig)
ax1 = fig.add_subplot (gs[0,:])

dirData = './'

# fs pulse
data = np.loadtxt (dirData+'field_102.dat')
time_field, Ef_field, Af_field = data[:,0], data[:,1], data[:,2]
plot_Ef (ax1, coef=0.0001)

# THz field
data = np.loadtxt (dirData+'field_101.dat')
time_field, Ef_field, Af_field = data[:,0], data[:,1], data[:,2]
plot_Ef (ax1)


""" 
data = np.loadtxt (dirData+'field_2.dat')
time_field, Ef_field, Af_field = data[:,0], data[:,1], data[:,2]
plot_Ef (ax1, coef=1e3)
 """

plt.show()
