import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import colorsys

matplotlib.rc('figure', figsize=(3.15, 3.15), dpi=300)
matplotlib.rc('figure.subplot', left=0.15, bottom=0.15,
              right=0.96, top=0.96, wspace=0.15, hspace=0)
matplotlib.rc('font', family='Times New Roman', size=10)
matplotlib.rc('axes', labelsize=10, labelpad=2)
matplotlib.rc('xtick', labelsize=8)
matplotlib.rc('ytick', labelsize=8)
matplotlib.rc('lines', linewidth=1, markersize=1.0)
matplotlib.rc('image', cmap='jet')
matplotlib.rc('xtick', direction='in')
matplotlib.rc('ytick', direction='in')
#matplotlib.rc('text', usetex=True)
#matplotlib.rc('text.latex', unicode=True, preamble=['\usepackage{amsmath}'])

order_reconstruct = 6
delay_0, delay_1, n_delay = 0., 6283., 100

# d(w)
def plot_hhg (ax, ylim0=-2.5, ylim1=0., lw=1):
    order_max = 100
    ax.plot(order, spec, lw=lw)
    #ax.set_xlim (0, order_max)
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

    idx_result = np.argmin (np.abs (order_valid - order_reconstruct))
    spec_valid[idx_result]

    #idx_local_max = argrelextrema (spec_valid, np.greater)[0].tolist()
    #idx_odd_order = [np.argmin (np.abs (order_valid - i))
    #                 for i in range (2, order_max, 2)]
    #ax.plot (order_valid[idx_odd_order], spec_valid[idx_odd_order])
    #ax.scatter (order_valid[idx_odd_order],
    #            spec_valid[idx_odd_order],
    #            c='r', s=8)
    

def calc_spec (dipole):
    global dt, w0
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
""" 
from matplotlib.gridspec import GridSpec
fig = plt.figure (constrained_layout=True)
gs = GridSpec (2, 3, figure=fig)
ax1 = fig.add_subplot (gs[0,:])
ax2 = fig.add_subplot (gs[1,:], sharex=ax1)
#ax2 = fig.add_subplot (gs[1,:])
 """

dirData = './res/'

#data = np.loadtxt (dirData+'field_0.dat')
#time_field, Ef_field, Af_field = data[:,0], data[:,1], data[:,2]
# plot_Ef (ax1)
# plot_dipole (ax2)

def process_data_of_reconstruction ():
    global dt, w0
    result = np.empty ([n_delay, 2])

    d_delay = (delay_1 - delay_0) / (n_delay - 1.)

    for i_pts in range (n_delay):
        print ("processing at " + str(i_pts) + "...")
        # parameters
        with open (dirData+'para_%d.dat' % i_pts, 'r') as f:
            l = np.array (f.readline().split(' '), float)
            Ip, w0, E0 = l[0], l[1], l[2]
            l = np.array (f.readline().split(' '), float)
            dt = l[2]
        E_cutoff = Ip + 3.17*(E0**2/(4.*w0**2))
        order_cutoff = E_cutoff / w0
        data = np.loadtxt (dirData+'data_%d.dat' % i_pts)
        time, Ef, Af, alphaf, dreal, dimag = data[:,0], data[:,1], data[:,2], data[:,3], data[:,4], data[:,5]

        # for arbitrary field
        order, spec = calc_spec (dreal)
        idx_valid = np.nonzero (spec > -6)
        spec_valid = spec[idx_valid]
        order_valid = order[idx_valid]

        idx_result = np.argmin (np.abs (order_valid - order_reconstruct))
        result[i_pts, 1] = spec_valid[idx_result]
        result[i_pts, 0] = delay_0 + i_pts * d_delay

    np.savetxt ('reconstruct.dat', result)


def generate_data_for_full_spectrum_vs_delay ():
    global dt, w0

    d_delay = (delay_1 - delay_0) / (n_delay - 1.)

    with open (dirData+'para_0.dat', 'r') as f:
        l = np.array (f.readline().split(' '), float)
        Ip, w0, E0 = l[0], l[1], l[2]
        l = np.array (f.readline().split(' '), float)
        dt = l[2]
    E_cutoff = Ip + 3.17*(E0**2/(4.*w0**2))
    order_cutoff = E_cutoff / w0
    data = np.loadtxt (dirData+'data_0.dat')
    time, Ef, Af, alphaf, dreal, dimag = data[:,0], data[:,1], data[:,2], data[:,3], data[:,4], data[:,5]

    result = np.empty ([len (time)//2, n_delay])

    order, spec = calc_spec (dreal)
    result[:,0] = spec

    for i_pts in range (1, n_delay):
        print ("processing at " + str(i_pts) + "...")
        # parameters
        with open (dirData+'para_%d.dat' % i_pts, 'r') as f:
            l = np.array (f.readline().split(' '), float)
            Ip, w0, E0 = l[0], l[1], l[2]
            l = np.array (f.readline().split(' '), float)
            dt = l[2]
        E_cutoff = Ip + 3.17*(E0**2/(4.*w0**2))
        order_cutoff = E_cutoff / w0
        data = np.loadtxt (dirData+'data_%d.dat' % i_pts)
        time, Ef, Af, alphaf, dreal, dimag = data[:,0], data[:,1], data[:,2], data[:,3], data[:,4], data[:,5]

        # for arbitrary field
        order, spec = calc_spec (dreal)
        result[:,i_pts] = spec

    np.savetxt ("spec2d_vs_delay.dat", result)
    np.savetxt ("order.dat", order)

def generate_data_for_full_spectrum_vs_delay_o_THz ():
    global dt, w0
    dirData = './'
    with open (dirData+'para_0.dat', 'r') as f:
        l = np.array (f.readline().split(' '), float)
        Ip, w0, E0 = l[0], l[1], l[2]
        l = np.array (f.readline().split(' '), float)
        dt = l[2]
    E_cutoff = Ip + 3.17*(E0**2/(4.*w0**2))
    order_cutoff = E_cutoff / w0
    data = np.loadtxt (dirData+'data_0.dat')
    time, Ef, Af, alphaf, dreal, dimag = data[:,0], data[:,1], data[:,2], data[:,3], data[:,4], data[:,5]

    result = np.empty ([len (time)//2, n_delay])

    order, spec = calc_spec (dreal)
    for i_pts in range (n_delay):
        result[:,i_pts] = spec
    np.savetxt ("spec2d_vs_delay_o_THz.dat", result)    
""" 
    idx_valid = np.nonzero (spec > -6)
    spec_valid = spec[idx_valid]
    order_valid = order[idx_valid]
    idx_result = np.argmin (np.abs (order_valid - 78))
    print (idx_result)
 """



def plot_reconstruct_waveform ():
    data = np.loadtxt ('field_101.dat')
    time_field, Ef_field, Af_field = data[:,0], data[:,1], data[:,2]
    plt.plot (time_field, np.abs (Ef_field), lw=1.5)
    data = np.loadtxt ('reconstruct.dat')
    plt.plot (data[:,0], 0.011*(10**data[:,1]), 'rx', markersize=4)
    plt.show ()

def plot_full_spectrum_vs_delay ():
    data_o_THz = np.loadtxt ('spec2d_vs_delay_o_THz.dat')
    data = np.loadtxt ('spec2d_vs_delay.dat')
    
    #plt.imshow ( (10**data), origin='lower', aspect=1000., extent=(0., 94274., 0, 100), vmin=0., vmax=0.01 )    
    plt.imshow ( (10**data - 10**data_o_THz), origin='lower', aspect=1000., extent=(delay_0, delay_1, 0, 100),  vmin=-0.00, vmax=0.0007 )
    plt.colorbar (orientation='horizontal')
    #plt.plot (data[:,0])
    #plt.plot (10**data[6422,:]-10**data_o_THz[6422,:])
    plt.show ()


def plot_colorbar ():
    data_o_THz = np.loadtxt ('spec2d_vs_delay_o_THz.dat')
    order = np.loadtxt ('order.dat')
    data = np.loadtxt ('spec2d_vs_delay.dat')
    data_field = np.loadtxt ('field_101.dat')
    data_reconstruct = np.loadtxt ('reconstruct.dat')
    
    t0, t1 = 0., 60000. # 94274

    idx_t1 = np.argmin (np.abs (data_reconstruct[:,0] - t1))
    t1 = data_reconstruct[idx_t1,0]
    idx_range = range(0, idx_t1)

    aspect = 350

    #plt.imshow ( data[:,idx_range], aspect=aspect, origin='lower', extent=(0., t1, 0, 100), vmin=-4, vmax=-2)
    #plt.colorbar ()
    #plt.show ()


    plt.imshow ( (10**data[:,idx_range] - 10**data_o_THz[:,idx_range]), aspect=aspect, origin='lower', extent=(0., t1, 0, 100),  vmin=-5e-4, vmax=5e-4,
                cmap=plt.get_cmap('bwr') )
    plt.colorbar ()
    plt.show ()


def plot_figure ():
    from matplotlib.gridspec import GridSpec
    #fig = plt.figure (constrained_layout=True)
    fig = plt.figure ()
    gs = GridSpec (5, 5, figure=fig, hspace=0.1, wspace=0.1)
    gs.update(left=0.1,right=0.9,top=0.985,bottom=0.04,wspace=0.05,hspace=0.05)

    ax1 = plt.subplot (gs[0,1:])
    ax2 = plt.subplot (gs[1:3,0])
    ax3 = plt.subplot (gs[1:3,1:], sharex=ax1, sharey=ax2)
    ax4 = plt.subplot (gs[3:5,1:], sharex=ax1)

    data_o_THz = np.loadtxt ('spec2d_vs_delay_o_THz.dat')
    order = np.loadtxt ('order.dat')
    data = np.loadtxt ('spec2d_vs_delay.dat')
    data_field = np.loadtxt ('field_101.dat')
    data_reconstruct = np.loadtxt ('reconstruct.dat')
    
    t0, t1 = 0., 60000. # 94274

    idx_t1 = np.argmin (np.abs (data_reconstruct[:,0] - t1))
    t1 = data_reconstruct[idx_t1,0]
    idx_range = range(0, idx_t1)

    ax2.plot (data[:,0], order, lw=0.2, c='k')
    ax2.set_xlim (0, -5)
    ax2.xaxis.tick_top()
    aspect = 350
    im = ax3.imshow ( data[:,idx_range], aspect=aspect, origin='lower', extent=(0., t1, 0, 100), vmin=-4, vmax=-2)
    ax3.axes.tick_params (axis='both', which='both', bottom=False, top=False, right=False, left=False, labelleft=False, labelbottom=False)
    # ax2.label_outer ()
    
    im = ax4.imshow ( (10**data[:,idx_range] - 10**data_o_THz[:,idx_range]), aspect=aspect, origin='lower', extent=(0., t1, 0, 100),  vmin=-5e-4, vmax=5e-4,
                        cmap=plt.get_cmap('bwr') )
    ax4.axes.tick_params (axis='both', which='both', top=True, right=True,  labelleft=False)
    # fig.colorbar (im, cax=ax4)

    idx_t0 = np.argmin (np.abs (data_field[:,0] - 0))
    if data_field[idx_t0, 0] < 0:
        idx_t0 = idx_t0 + 1
    idx_t1 = np.argmin (np.abs (data_field[:,0] - t1))
    if data_field[idx_t1, 0] > t1:
        idx_t1 = idx_t1 - 1
    
    ax1.plot (data_field[idx_t0:idx_t1,0], data_field[idx_t0:idx_t1,1], c='k', label='THz field')
    ax1.grid (dashes=(12,12), lw=0.2)
    #ax1.plot (data_reconstruct[idx_range,0], 0.1 * (10**data[6422,idx_range]-10**data_o_THz[6422,idx_range]), label='slice')
    #ax1.legend ()
    #print (np.shape (data[:,0]))
    
    with open (dirData+'para_0.dat', 'r') as f:
        l = np.array (f.readline().split(' '), float)
        Ip, w0, E0 = l[0], l[1], l[2]
        l = np.array (f.readline().split(' '), float)
        dt = l[2]
    E_cutoff = Ip + 3.17*(E0**2/(4.*w0**2))
    order_cutoff = E_cutoff / w0 
    ax2.axhline(y=order_cutoff, lw=0.5, c='k', linestyle='--')
    ax4.axhline(y=order_cutoff, xmin=0, xmax=0.2, lw=0.2, c='k', linestyle='--', dashes=(20, 20))

    order_E_2 = (Ip + 2.4*(E0**2/(4.*w0**2))) / w0
    ax4.axhline(y=order_E_2, xmin=0, xmax=0.2, lw=0.2, c='k', linestyle='--', dashes=(20, 20))

    order_E_Ip = (Ip) / w0
    ax4.axhline(y=order_E_Ip, xmin=0, xmax=0.2, lw=0.2, c='k', linestyle='--', dashes=(20, 20))    

    ax4.axhline(y=2, xmin=0, xmax=0.2, lw=0.2, c='k', linestyle='--', dashes=(20, 20))
    #ax4.grid (dashes=(12,12), lw=0.2)  

    plt.show ()




#process_data_of_reconstruction ()
plot_reconstruct_waveform ()

# generate_data_for_full_spectrum_vs_delay ()
#generate_data_for_full_spectrum_vs_delay_o_THz ()

#plot_full_spectrum_vs_delay ()
#plot_colorbar ()
#plot_figure ()


#plot_hhg (ax1, ylim0=-4.5, ylim1=-2)

#ax1.label_outer ()

#plt.show()
