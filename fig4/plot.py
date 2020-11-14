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
matplotlib.rc('lines', linewidth=1.0, markersize=4.0)
matplotlib.rc('image', cmap='jet')
#matplotlib.rc('text', usetex=True)
#matplotlib.rc('text.latex', unicode=True, preamble=['\\usepackage{amsmath}'])


# compare gaussian-enveloped pulse and monochromatic laser
#data = np.loadtxt ('res.dat')
#plt.plot (data[:,0], data[:,1], data[:,0], data[:,2])

data = np.loadtxt ('res_mono.dat')
plt.plot (data[:,0], data[:,1], label='Phi_cw')
data = np.loadtxt ('res_gaus.dat')
plt.plot (data[:,0], data[:,1], '--', label='Phi_gauss')

data = np.loadtxt ('res_gaus_env.dat')
plt.plot (data[:,0], data[:,1], label='e*Phi_gaus')

plt.legend ()
#plt.ylim ([-0.05, 0.15])
plt.ylim ([0, 7])
plt.xlabel ('phi_t')
plt.ylabel ('value')

plt.show ()