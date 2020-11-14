import numpy as np
import matplotlib.pyplot as plt

# Ef(t)
def plot_Ef (ax):
    ax.plot(time_field, Af_field)
    ax.set_ylabel('Af')
    ax.tick_params(labelbottom=False)
    ax.grid()
# d(t)
def plot_dipole (ax):
    ax.plot(time, dreal)
    ax.set_xlabel('Time (a.u.)')
    ax.set_ylabel(r'$d(t)$')
    ax.grid()

# d(w)
def plot_hhg (ax):
    order_max = 100
    ax.plot(order, spec, '.-')
    ax.set_xlim (0, order_max)
    #ax.set_ylim (-2.5, 0.)
    [ax.axvline (x=i, lw=0.5, c='r') for i in range (order_max) [1:-1:2]]
    ax.axvline(x=order_cutoff, lw=2, c='k')
    ax.axvline(x=Ip/w0, lw=2, c='k')
    ax.set_xlabel(r'Order ($\omega/\omega_c$)')
    ax.set_ylabel('Harmonics')

    from scipy.signal import argrelextrema
    idx_valid = np.nonzero (spec > -6)
    spec_valid = spec[idx_valid]
    order_valid = order[idx_valid]
    idx_local_max = argrelextrema (spec_valid, np.greater)[0].tolist()
    idx_odd_order = [np.argmin (np.abs (order_valid - i))
                     for i in range (1, order_max, 2)]
    import functools
    import operator
    threshold = 1
    idx_near_odd_order = functools.reduce (
        operator.add,
        [[j + i for j in range(-threshold, threshold+1)]
         for i in idx_odd_order])
    idx_peaks = (list (set (idx_local_max) & set (idx_near_odd_order)))
    idx_peaks.sort()
    ax.scatter (order_valid[idx_peaks],
                spec_valid[idx_peaks],
                c='r')
    ax.plot (order_valid[idx_peaks], spec_valid[idx_peaks])

def calc_spec (dipole):
    nt = np.size (dipole)
    w = 2. * np.pi * np.fft.fftfreq (nt, dt) [:nt//2]
    spec = np.fft.fft (2 * dipole) [:nt//2]
    #spec = np.log10(w*w*w*w*np.abs(spec)**2)
    spec = np.log10 (np.abs (spec)) # not squre
    order = w / w0
    return order, spec


dirData = './'
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

data = np.loadtxt (dirData+'field_0.dat')
time_field, Ef_field, Af_field = data[:,0], data[:,1], data[:,2]

def wrap_padding (t_sel, n_period):
    idx1 = np.argmin (np.abs (time-t_sel))
    idx2 = np.argmin (np.abs (time-(t_sel+2.*np.pi/w0)))
    time_seg = time[idx1:idx2]
    dreal_seg = dreal[idx1:idx2]
    n_seg = len (time_seg)
    dreal_pad = np.pad (dreal_seg, (n_seg*n_period,), 'wrap')
    return dreal_pad

# for arbitrary field
# order, spec = calc_spec (dreal)
# for monochromatic field
time_selected = 1735.59
order, spec = calc_spec (wrap_padding (time_selected, 100))

from matplotlib.gridspec import GridSpec
fig = plt.figure (constrained_layout=True)
gs = GridSpec (3, 3, figure=fig)
ax1 = fig.add_subplot (gs[0,:])
ax2 = fig.add_subplot (gs[1,:], sharex=ax1)
ax3 = fig.add_subplot (gs[2,:])

plot_Ef (ax1)
plot_dipole (ax2)
plot_hhg (ax3)

plt.show()
