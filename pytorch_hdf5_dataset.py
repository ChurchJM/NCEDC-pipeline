import h5py
import numpy as np
from torch.utils.data import DataLoader, Dataset
from scipy.stats import truncnorm
import matplotlib.pyplot as plt
from time import perf_counter

class HDF5Dataset(Dataset):
    def __init__(self, hdf5_path, single_dset, labels_on_fly):
        self.hdf5_file = h5py.File(hdf5_path, "r")
        self.labels_on_fly = labels_on_fly
        self.single_dset = single_dset
        if self.labels_on_fly:
            xs = np.linspace(-50, 50, num=101)
            self.trunc_norm = truncnorm.pdf(xs, loc=0, scale=10, a=-3, b=3)

    def __getitem__(self, index):
        if self.single_dset:
            x = self.hdf5_file['train_x'][index,:,:]
            y = self.hdf5_file['train_y'][index,:,:]
            return x,y
        if self.labels_on_fly:
            x = self.hdf5_file['train_x'][str(index)]
            
            p_sample = x.attrs['p_sample']
            s_sample = x.attrs['s_sample']
            arr_p = np.zeros(3000)
            arr_p[p_sample-50:p_sample+51] = self.trunc_norm
            arr_p = arr_p/np.max(self.trunc_norm)
            arr_s = np.zeros(3000)
            arr_s[s_sample-50:s_sample+51] = self.trunc_norm
            arr_s = arr_s/np.max(self.trunc_norm)
            arr_noise = 1-(arr_p + arr_s)
            y = np.asarray([arr_p, arr_s, arr_noise])
            
            return x[:],y
        else:
            x = self.hdf5_file['x_all'][str(index)][:]
            y = self.hdf5_file['y_all'][str(index)][:]
            return x,y

    def __len__(self):
        if self.single_dset:
            return len(self.hdf5_file['train_x'])
        else:
            return len(list(self.hdf5_file['train_x'].keys()))

options_grid = [
    ('./all_waveforms.hdf5', False, False),
    ('./all_waveforms.hdf5', False, True),
    ('./all_waveforms_gzip.hdf5', False, False),
    ('./all_waveforms_gzip.hdf5', False, True),
    ('./all_waveforms_no_attrs.hdf5', True, False),
]

for option in options_grid:
    dataloader = DataLoader(HDF5Dataset(hdf5_path=option[0], single_dset=option[1], labels_on_fly=option[2]), batch_size=32, num_workers=0, shuffle=True)
    n_batches = len(dataloader)

    start = perf_counter()
    for i, (data, target) in enumerate(dataloader):
        print(f'\rLoaded {i+1} batches of {n_batches}', end='', flush=True)
    end = perf_counter()
    print(f'\n{option}. Seconds elapsed: {end-start}')

# for i, (data, target) in enumerate(dataloader):
#     print(data.shape)
#     if i%100 == 0:
#         waveform = data[0,:,:]
#         label = target[0,:,:]

#         # for attr in list(waveform.attrs):
#         #     print(f'{attr}: {waveform.attrs[attr]}')

#         # p_sample = waveform.attrs['p_sample']
#         # s_sample = waveform.attrs['s_sample']

#         xs_3000 = np.linspace(0,3000,3000)

#         fig, axes = plt.subplots(4,1,figsize = (8,5))
#         axes[0].plot(waveform[0], c='k', linewidth=0.5)
#         axes[1].plot(waveform[1], c='k', linewidth=0.5)
#         axes[2].plot(waveform[2], c='k', linewidth=0.5)
#         axes[3].plot(xs_3000, label[0], c='g', linestyle='--')
#         axes[3].plot(xs_3000, label[1], c='purple', linestyle='--')
#         axes[3].plot(xs_3000, label[2], c='k', alpha=0.5, linestyle='--')

#         # fig.suptitle(f"nc{waveform.attrs['event_id']} -- {waveform.attrs['network']}.{waveform.attrs['station']}.{waveform.attrs['location']}.{waveform.attrs['channel']} -- n_samples: {len(waveform[0])}")

#         axes[0].axvline(int(len(waveform[0])/2), c='orange', alpha=0.3, linestyle='--')
#         axes[1].axvline(int(len(waveform[0])/2), c='orange', alpha=0.3, linestyle='--')
#         axes[2].axvline(int(len(waveform[0])/2), c='orange', alpha=0.3, linestyle='--')

#         for ax in axes:
#             ax.margins(x=0)
#             # ax.axvline(p_sample, c='b')
#             # ax.axvline(s_sample, c='r')

#         plt.tight_layout()
#         plt.show()