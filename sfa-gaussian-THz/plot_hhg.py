import numpy as np
import matplotlib.pyplot as plt

# Ef(t)
def plot_Ef (ax):
    ax.plot(time, Af)
    ax.plot(time, Ef)
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
def plot_hhg (ax):
    order_max = 100
    ax.plot(order, spec, '.-')
    ax.set_xlim (0, order_max)
    ax.set_ylim (-4.5, -1.5)
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


def process_data (dirData='./', id=-1):
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

    if id > -1:
        filename_para = dirData + 'para_%d.dat' % id
        filename_data = dirData + 'data_%d.dat' % id
    else:
        filename_para = dirData + 'para.dat'
        filename_data = dirData + 'data.dat'

    with open (filename_para, 'r') as f:
        l = np.array (f.readline().split(' '), float)
        Ip, w0, E0 = l[0], l[1], l[2]
        l = np.array (f.readline().split(' '), float)
        #dt = l[2]
        dt, tc = l[2], l[3]

        E_cutoff = Ip + 3.17*(E0**2/(4.*w0**2))
        order_cutoff = E_cutoff / w0

        data = np.loadtxt (filename_data)
        time, Ef, Af, alphaf, dreal, dimag = data[:,0], data[:,1], data[:,2], data[:,3], data[:,4], data[:,5]

    # for arbitrary field
    order, spec = calc_spec (dreal)
    # # for monochromatic field
    # time_selected = 1735.59
    # order, spec = calc_spec (wrap_padding (time_selected, 100))        

    return time, Af, Ef, dreal, order, spec, Ip, w0, E0, tc, order_cutoff
    #return time, Af, Ef, dreal, order, spec, Ip, w0, E0, order_cutoff

###### figures ######

from matplotlib.gridspec import GridSpec
fig = plt.figure (constrained_layout=True)
gs = GridSpec (3, 3, figure=fig)
ax1 = fig.add_subplot (gs[0,:])
ax2 = fig.add_subplot (gs[1,:], sharex=ax1)
ax3 = fig.add_subplot (gs[2,:])


time, Af, Ef, dreal, order, spec, Ip, w0, E0, tc, order_cutoff = process_data ('./test-data/numeric_vs_analytic/full/')
plot_Ef (ax1)
plot_dipole (ax2)
plot_hhg (ax3)

time, Af, Ef, dreal, order, spec, Ip, w0, E0, tc, order_cutoff = process_data ('./test-data/using_faddeeva_eps_to_1st_order/')
plot_hhg (ax3)

""" 
time, Af, Ef, dreal, order, spec, Ip, w0, E0, tc, order_cutoff = process_data ('./test-data/numeric_without_THz/')
plot_hhg (ax3)
 """

""" 
time, Af, Ef, dreal, order, spec, Ip, w0, E0, tc, order_cutoff = process_data ('./test-data/using_faddeeva_eps_to_1st_order_large_z_limit_1/')
plot_hhg (ax3)
 """


time, Af, Ef, dreal, order, spec, Ip, w0, E0, tc, order_cutoff = process_data ('./test-data/using_faddeeva_eps_to_1st_order_large_z_limit_2_polynomial_highest_order_exp1/')
plot_hhg (ax3)


""" time, Af, Ef, dreal, order, spec, Ip, w0, E0, tc, order_cutoff = process_data ('./test-data/f_approx/')
plot_hhg (ax3)

time, Af, Ef, dreal, order, spec, Ip, w0, E0, tc, order_cutoff = process_data ('./test-data/derived-formula/')
plot_hhg (ax3) """

#ax3.set_xlim ([70, 84])
#ax3.set_ylim ([-4.5, -2])

plt.show()



# select a certain order harmonic to check the dependence on some parameter
# 1. generate data file from given order harmonic
""" with open ("even-harmonic-vs-tc.dat", 'w') as f:
    for i in range (100):
        time, Af, Ef, dreal, order, spec, Ip, w0, E0, tc, order_cutoff = process_data ('./test-data/scan_test/', i)
        # select the 78th even-order harmonic for comparison
        idx = np.argmin (np.abs (order - 78))
        f.write ("%f %lf" % (tc, spec[idx]))
        f.write ('\n')
 """
# 2. plot the data, comparing with the original THz electric field
""" data = np.loadtxt ("even-harmonic-vs-tc.dat")
#plt.plot (data[:,0], 10**data[:,1], np.abs (data[:,2]))
THz_electric_field = np.sin ((1e-4)*data[:,0])
#plt.plot (data[:,0], np.abs (electric_field))
plt.plot (data[:,0], 10**data[:,1], data[:,0], 3.7e-4*np.abs (THz_electric_field))
plt.show () """