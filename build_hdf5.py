import pandas as pd
import sqlite3
from sklearn.model_selection import train_test_split
import h5py
from obspy import read
import numpy as np
from scipy.stats import truncnorm
from obspy import UTCDateTime
import os

np.random.seed(42)

def has_zero_component(row):
    stream = read(row['waveform_path'])
    for i in range(3):
        if np.std(stream[i].data) < 0.1: # Std and not mean because there are completely flat/constant waveforms with non-zero means.
            return True
    return False

# Load phase table (3-component) into pandas dataframe.
con = sqlite3.connect('./dbs/NCEDC.db')
df = pd.read_sql_query('select * from phase', con)
con.close()

if not os.path.isfile('./dataframes/df_complete.csv'):
    # Only use complete waveforms for now.
    df_complete = df[df['complete']==1]
    print(len(df_complete))

    # Waveforms with P/S lag < 0.35 s cause an issue where P + S probability is > 1 when computing truncated Gaussian labels.
    # P/S lags > 26 s are too large to fit in 30 s waveforms.
    df_complete['p_s_diff_sec'] = df_complete.apply(lambda x: UTCDateTime(x['S'])-UTCDateTime(x['P']), axis=1)
    df_complete = df_complete[(df_complete['p_s_diff_sec'] >= 0.35)&(df_complete['p_s_diff_sec'] <= 26)] # ToDo: Can we address the first constraint in a different way? (removes 6,605 total waveforms)

    # Find rows where one or more components returned from NCEDC are zeros (or some other completely flat, constant value).
    df_complete['has_zero_component'] = df_complete.apply(lambda x: has_zero_component(x), axis=1)
    df_complete = df_complete[~df_complete['has_zero_component']] # Filters out 4034 rows.

    df_complete.reset_index(drop=True, inplace=True)

    df_complete.to_csv('./dataframes/df_complete.csv', index=False)
else:
    print('Filtered dataframe already exists. Loading exisitng file.')
    df_complete = pd.read_csv('./dataframes/df_complete.csv')

print(f'Total complete rows: {len(df_complete)}.')
# print(df_complete.head())

df_complete['network.station'] = df_complete.apply(lambda x: f'{x["network"]}.{x["station"]}', axis=1)
station_counts = df_complete['network.station'].value_counts()
df_filtered = df_complete[df_complete['network.station'].isin(station_counts[station_counts > 9].index)]
print(f'Total rows where network.station occurrences >= 10: {len(df_filtered)}.')

if not os.path.isfile('./dataframes/df_train.csv'):
    # Stratified train/validation/test split.
    # The current sizes produce a dataset on the order of PhaseNet, but could be increased as we have additional data.
    # ToDo: sklearn seed?
    df_train_val, df_test = train_test_split(df_filtered, train_size=711000, test_size=79000, stratify=df_filtered["network.station"], random_state=42)
    df_train, df_val = train_test_split(df_train_val, train_size=632000, test_size=79000, stratify=df_train_val["network.station"], random_state=42)
    df_train.to_csv('./dataframes/df_train.csv', index=False)
    df_val.to_csv('./dataframes/df_val.csv', index=False)
    df_test.to_csv('./dataframes/df_test.csv', index=False)
else:
    print('Train, validation, and test split dataframes already exist. Loading existing files.')
    df_train = pd.read_csv('./dataframes/df_train.csv')
    df_val = pd.read_csv('./dataframes/df_val.csv')
    df_test = pd.read_csv('./dataframes/df_test.csv')

# Check the counts.
df_counts = pd.concat([df_train['network.station'].value_counts(), df_val['network.station'].value_counts(), df_test['network.station'].value_counts()],
                      axis=1)
df_counts.columns = ['train_count', 'val_count', 'test_count']
df_counts['val_count'].fillna(0, inplace=True)
df_counts['test_count'].fillna(0, inplace=True)
df_counts['val_count'] = df_counts['val_count'].astype(int)
df_counts['test_count'] = df_counts['test_count'].astype(int)
print(df_counts.head(20))

hdf5_file_path = './all_waveforms_gzip.hdf5' # ToDo: change this to a command line arg, or put in config.
if not os.path.isfile(hdf5_file_path):
    # Add ALL filtered waveforms to an HDF5 file, using links to store train/val/test splits.
    hdf5_file = h5py.File(hdf5_file_path, 'w')
    grp_x_all = hdf5_file.create_group('x_all')
    grp_y_all = hdf5_file.create_group('y_all')

    xs = np.linspace(-50, 50, num=101)
    xs_3000 = np.linspace(0,3000,3000)
    trunc_norm = truncnorm.pdf(xs, loc=0, scale=10, a=-3, b=3)

    # zero_component_ids = []

    # Save the dataframes out at this point? "Pandas processing/filtering takes around an hour, so these intermediate results are saved out..."
    # Also check for existence of hdf5 groups. Basically allow this script to be re-run in any state, only producing the missing assets.

    for index, row in df_complete.iterrows():
        id = row['id']
        stream = read(row['waveform_path'])
        
        p_dt = UTCDateTime(row['P'])
        s_dt = UTCDateTime(row['S'])
        p_s_offset = (s_dt-p_dt)
        mid_dt = p_dt + p_s_offset/2
        max_offset = mid_dt + 15 - s_dt - 2
        max_offset_samples = round(max_offset*100) # Resampling hasn't occurred yet, so need to hard code 100 Hz here. Not particularly safe.
        if max_offset_samples <= 0:
            print('OFFSET ERROR. Offset will be zero.', max_offset, p_dt, s_dt)
            random_offset_samples = 0
        else:
            random_offset_samples = np.random.randint(-max_offset_samples,max_offset_samples)

        start_idx = 1500 + random_offset_samples

        components = []
        for i in range(3):
            stream[i].resample(100) # Resample waveform to 100 Hz, as expected by PhaseNet.
            component = stream[i].data # Convert trace to numpy array.
            if len(component) < 6000: # Should hopefully never be true since we're using resample().
                print('Shorter than 6000: ', len(stream[i]))
                component = np.pad(component, (0,6000-len(component)), 'constant')

            # Trim waveforms to 30 s (with random offset), as expected by PhaseNet.
            component = component[start_idx:start_idx+3000]

            # Normalize the waveforms as specified by PhaseNet authors.
            component = (component-np.mean(component))/np.std(component)
            components.append(component)

        # ToDo: Add gzip option to config.
        dset_x = grp_x_all.create_dataset(str(index), data=np.asarray(components), compression='gzip')
        for column in df_complete.columns:
            if column in ['id', 'instrument_match', 'complete', 'waveform_path', 'phase_file', 'network.station', 'has_zero_component', 'p_s_diff_sec']:
                continue
            dset_x.attrs.create(column, str(row[column]))
        dset_x.attrs.create('db_id', id)
        starttime_30s = stream[0].stats.starttime + 15 + random_offset_samples/stream[0].stats.sampling_rate
        endtime_30s = stream[0].stats.endtime - 15 + random_offset_samples/stream[0].stats.sampling_rate
        dset_x.attrs.create('starttime', starttime_30s.isoformat())
        dset_x.attrs.create('endtime', endtime_30s.isoformat())
        
        p_sample = round((p_dt - starttime_30s) * stream[0].stats.sampling_rate)
        s_sample = round((s_dt - starttime_30s) * stream[0].stats.sampling_rate)

        dset_x.attrs.create('p_sample', p_sample)
        dset_x.attrs.create('s_sample', s_sample)

        arr_p = np.zeros(3000)
        arr_p[p_sample-50:p_sample+51] = trunc_norm
        arr_p = arr_p/np.max(trunc_norm)

        arr_s = np.zeros(3000)
        arr_s[s_sample-50:s_sample+51] = trunc_norm
        arr_s = arr_s/np.max(trunc_norm)

        arr_noise = 1-(arr_p + arr_s) # ToDo: check if arr_noise is ever negative!
        
        if np.any(arr_noise < 0):
            print('NEGATIVE VALUE IN NOISE DISTRIBUTION!')
            continue

        dset_y = grp_y_all.create_dataset(str(index), data=np.asarray([arr_p, arr_s, arr_noise]), compression='gzip')

        print(f'\rProcessed {index}', end='', flush=True)

    hdf5_file.close()
else:
    print(f'HDF5 file {hdf5_file_path} already exists. Loading existing file.')

# Before adding train/validation/test subgroups, ucompressed file size is 178.69 gb. Grew to 178.77 gb after adding links. To 178.75 gb after removin first train,val,test groups. 178.9?
# 178.85 gb
# For gzip: 73.16 gb -> 73.23 gb
hdf5_file = h5py.File(hdf5_file_path, 'a')

# print(list(hdf5_file.keys()))
# del hdf5_file['train']
# del hdf5_file['validation']
# del hdf5_file['test']
# print(list(hdf5_file.keys()))

if not 'train_x' in hdf5_file:
    grp_train_x = hdf5_file.create_group('train_x')
    grp_train_y = hdf5_file.create_group('train_y')
    for index, row in df_train.iterrows():
        db_id_train = row['id']
        id_x_all = df_complete.index[df_complete['id'] == db_id_train][0]
        grp_train_x[str(index)] = hdf5_file['x_all'][str(id_x_all)]
        grp_train_y[str(index)] = hdf5_file['y_all'][str(id_x_all)]
        print(f'\r(Train) processed {index}', end='', flush=True)
    print('\n(Train) finished.')

if not 'validation_x' in hdf5_file:
    grp_val_x = hdf5_file.create_group('validation_x')
    grp_val_y = hdf5_file.create_group('validation_y')
    for index, row in df_val.iterrows():
        db_id_val = row['id']
        id_x_all = df_complete.index[df_complete['id'] == db_id_val][0]
        grp_val_x[str(index)] = hdf5_file['x_all'][str(id_x_all)]
        grp_val_y[str(index)] = hdf5_file['y_all'][str(id_x_all)]
        print(f'\r(Validation) processed {index}', end='', flush=True)
    print('\n(Validation) finished.')

if not 'test_x' in hdf5_file:
    grp_test_x = hdf5_file.create_group('test_x')
    grp_test_y = hdf5_file.create_group('test_y')
    for index, row in df_test.iterrows():
        db_id_test = row['id']
        id_x_all = df_complete.index[df_complete['id'] == db_id_test][0]
        grp_test_x[str(index)] = hdf5_file['x_all'][str(id_x_all)]
        grp_test_y[str(index)] = hdf5_file['y_all'][str(id_x_all)]
        print(f'\r(Test) processed {index}', end='', flush=True)
    print('\n(Test) finished.')

# Some Validation Checks...
print(len(hdf5_file['train_x']), len(hdf5_file['validation_x']), len(hdf5_file['test_x']))

test_dset = hdf5_file['train_x'][str(0)]
for attr in list(test_dset.attrs):
    print(f'{attr}: {test_dset.attrs[attr]}')

hdf5_file.close()


# OPEN HDF5 file here and start processing the sub groups. (APPEND MODE!!!)

# print(f'Num zero component waveforms: {len(zero_component_ids)}')
# print(zero_component_ids[:100])

# Create second groups to hold stratified subsets using links. check diff in filesize.
    
# OR create hdf5 file with only the training data itself, if you want a smaller, less flexible file.


# check difference in size with compression on/off

# PyTorch dataloader.

# row_test = df_train.iloc[0]
# print(row_test)
# print(df_complete[df_complete['id'] == row_test['id']])
# row_id = df_complete.index[df_complete['id'] == row_test['id']][0]
# print(row_id)

# grp_all = hdf5_file['x_all']
# test_dset = grp_all[str(row_id)]

# for attr in list(test_dset.attrs):
#     print(f'{attr}: {test_dset.attrs[attr]}')