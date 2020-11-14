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
matplotlib.rc('xtick', direction='in')
#matplotlib.rc('text', usetex=True)
#matplotlib.rc('text.latex', unicode=True, preamble=['\usepackage{amsmath}'])

# Ef(t)
def plot_Ef (ax):
    ax.plot(time_field, Ef_field)
    ax.set_ylabel('Ef')
    ax.tick_params(labelbottom=False)
    ax.grid()
# d(t)
def plot_dipole (ax):
    ax.plot(time, dreal)
    ax.set_xlabel('Time (a.u.)')
    ax.set_ylabel(r'$d(t)$')
    ax.grid()

# d(w)
def plot_hhg (ax, ylim0=-2.5, ylim1=0., lw=1):
    order_max = 100
    ax.plot(order, spec, lw=lw)
    ax.set_xlim (0, order_max)
    ax.set_xlim (60, 90)
    ax.set_ylim (ylim0, ylim1)
    [ax.axvline (x=i, lw=0.4, c='grey', linestyle='--') for i in range (order_max) [1:-1:2]]
    [ax.axvline (x=i, lw=0.2, c='grey', linestyle='--', dashes=[6,10]) for i in range (order_max) [2:-1:2]]
    ax.axvline(x=order_cutoff, lw=2, c='k', linestyle='--')
    ax.axvline(x=Ip/w0, lw=2, c='k')

    from scipy.signal import argrelextrema
    idx_valid = np.nonzero (spec > -6)
    spec_valid = spec[idx_valid]
    order_valid = order[idx_valid]
    idx_local_max = argrelextrema (spec_valid, np.greater)[0].tolist()
    idx_odd_order = [np.argmin (np.abs (order_valid - i))
                     for i in range (2, order_max, 2)]
    #ax.plot (order_valid[idx_odd_order], spec_valid[idx_odd_order])
    ax.scatter (order_valid[idx_odd_order],
                spec_valid[idx_odd_order],
                c='r', s=8)
    

def calc_spec (dipole):
    nt = np.size (dipole)
    w = 2. * np.pi * np.fft.fftfreq (nt, dt) [:nt//2]
    spec = np.fft.fft (2 * dipole) [:nt//2]
    #spec = np.log10(w*w*w*w*np.abs(spec)**2)
    spec = np.log10 (np.abs (spec)) # not squre
    order = w / w0
    return order, spec




def wrap_padding (t_sel, n_period):
    idx1 = np.argmin (np.abs (time-t_sel))
    idx2 = np.argmin (np.abs (time-(t_sel+2.*np.pi/w0)))
    time_seg = time[idx1:idx2]
    dreal_seg = dreal[idx1:idx2]
    n_seg = len (time_seg)
    dreal_pad = np.pad (dreal_seg, (n_seg*n_period,), 'wrap')
    return dreal_pad

from matplotlib.gridspec import GridSpec
fig = plt.figure (constrained_layout=True)
gs = GridSpec (2, 3, figure=fig)
ax1 = fig.add_subplot (gs[0,:])
#ax2 = fig.add_subplot (gs[1,:], sharex=ax1)
ax2 = fig.add_subplot (gs[1,:])


dirData = './res/phi-0.5/'

# parameters
with open (dirData+'para_0.dat', 'r') as f:
    l = np.array (f.readline().split(' '), float)
    Ip, w0, E0 = l[0], l[1], l[2]
    l = np.array (f.readline().split(' '), float)
    dt = l[2]
E_cutoff = Ip + 3.17*(E0**2/(4.*w0**2))
order_cutoff = E_cutoff / w0
data = np.loadtxt (dirData+'data_0.dat')
time, Ef, Af, alphaf, dreal, dimag = data[:,0], data[:,1], data[:,2], data[:,3], data[:,4], data[:,5]

# for arbitrary field
order, spec = calc_spec (dreal)
plot_hhg (ax1, ylim0=-4.5, ylim1=-2)

# parameters
with open (dirData+'para_1.dat', 'r') as f:
    l = np.array (f.readline().split(' '), float)
    Ip, w0, E0 = l[0], l[1], l[2]
    l = np.array (f.readline().split(' '), float)
    dt = l[2]
E_cutoff = Ip + 3.17*(E0**2/(4.*w0**2))
order_cutoff = E_cutoff / w0
data = np.loadtxt (dirData+'data_1.dat')
time, Ef, Af, alphaf, dreal, dimag = data[:,0], data[:,1], data[:,2], data[:,3], data[:,4], data[:,5]

# for arbitrary field
order, spec = calc_spec (dreal)
plot_hhg (ax1, ylim0=-4.5, ylim1=-2, lw=0.6)


dirData = './res/phi-1.0/'

# parameters
with open (dirData+'para_0.dat', 'r') as f:
    l = np.array (f.readline().split(' '), float)
    Ip, w0, E0 = l[0], l[1], l[2]
    l = np.array (f.readline().split(' '), float)
    dt = l[2]
E_cutoff = Ip + 3.17*(E0**2/(4.*w0**2))
order_cutoff = E_cutoff / w0
data = np.loadtxt (dirData+'data_0.dat')
time, Ef, Af, alphaf, dreal, dimag = data[:,0], data[:,1], data[:,2], data[:,3], data[:,4], data[:,5]

# for arbitrary field
order, spec = calc_spec (dreal)
plot_hhg (ax2, ylim0=-4.5, ylim1=-2)

# parameters
with open (dirData+'para_1.dat', 'r') as f:
    l = np.array (f.readline().split(' '), float)
    Ip, w0, E0 = l[0], l[1], l[2]
    l = np.array (f.readline().split(' '), float)
    dt = l[2]
E_cutoff = Ip + 3.17*(E0**2/(4.*w0**2))
order_cutoff = E_cutoff / w0
data = np.loadtxt (dirData+'data_1.dat')
time, Ef, Af, alphaf, dreal, dimag = data[:,0], data[:,1], data[:,2], data[:,3], data[:,4], data[:,5]

# for arbitrary field
order, spec = calc_spec (dreal)
plot_hhg (ax2, ylim0=-4.5, ylim1=-2, lw=0.6)

ax1.label_outer ()

plt.show()
