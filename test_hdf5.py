import numpy as np
import h5py
import matplotlib.pyplot as plt

f = h5py.File('./all_waveforms.hdf5', 'r')
grp_train_x = f['train_x']
grp_train_y = f['train_y']

n_rows = len(list(grp_train_x.keys()))

# dset_train_x = f['train_x']
# dset_train_y = f['train_y']
# n_rows = len(dset_train_x)


while True:
    r = np.random.randint(0, n_rows)
    waveform = grp_train_x[str(r)]
    label = grp_train_y[str(r)]

    # waveform = dset_train_x[r,:,:]
    # label = dset_train_y[r,:,:]

    for attr in list(waveform.attrs):
        print(f'{attr}: {waveform.attrs[attr]}')

    p_sample = waveform.attrs['p_sample']
    s_sample = waveform.attrs['s_sample']

    xs_3000 = np.linspace(0,3000,3000)

    fig, axes = plt.subplots(4,1,figsize = (8,5))
    axes[0].plot(waveform[0], c='k', linewidth=0.5)
    axes[1].plot(waveform[1], c='k', linewidth=0.5)
    axes[2].plot(waveform[2], c='k', linewidth=0.5)
    axes[3].plot(xs_3000, label[0], c='g', linestyle='--')
    axes[3].plot(xs_3000, label[1], c='purple', linestyle='--')
    axes[3].plot(xs_3000, label[2], c='k', alpha=0.5, linestyle='--')

    fig.suptitle(f"nc{waveform.attrs['event_id']} -- {waveform.attrs['network']}.{waveform.attrs['station']}.{waveform.attrs['location']}.{waveform.attrs['channel']} -- n_samples: {len(waveform[0])}")

    axes[0].axvline(int(len(waveform[0])/2), c='orange', alpha=0.3, linestyle='--')
    axes[1].axvline(int(len(waveform[0])/2), c='orange', alpha=0.3, linestyle='--')
    axes[2].axvline(int(len(waveform[0])/2), c='orange', alpha=0.3, linestyle='--')

    for ax in axes:
        ax.margins(x=0)
        ax.axvline(p_sample, c='b')
        ax.axvline(s_sample, c='r')

    plt.tight_layout()
    plt.show()