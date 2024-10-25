from obspy.clients.fdsn import Client
from obspy import UTCDateTime, Stream
import time
import os
import sqlite3
import os
import sys
from threading import Thread
from queue import Queue
import queue
import obspy

# Queues used for thread coordination.
# ToDo: add more extensive comment, including about the DownloadWorker class.
batch_queue = Queue()
sql_queue = Queue()    # For sqlite writes.
retry_queue = Queue()  # For sqlite writes.

class DownloadWorker(Thread):
    def __init__(self, queue):
        Thread.__init__(self)
        self.queue = queue
        self.log_path = ''
    
    def run(self):
        while True:
            try:
                file, batch_start, batch_end = self.queue.get()
                try:
                    process_lines(file, self.log_path, batch_start, batch_end)
                finally:
                    self.queue.task_done()
            except queue.Empty:
                pass

def log_message(message, log_path):
    with(open(log_path, 'a')) as logfile:
        logfile.write(f'{message}\n\n')

# Specific years or months to download may be specified via command line arguments.
# If ANY args exceed 4 characters, we treat all args as months (e.g., match a single file only).
timespans_to_process = []
for i in range(1, len(sys.argv)):
    timespans_to_process.append(sys.argv[i].strip())
print(f'Processing timespans: {timespans_to_process}')

# If command line args are years, we only want to match the first 4 characters of the file name.
# If args are months, we need to match the first 7.
endIdx = 4
maxArgLen = max(len(i) for i in timespans_to_process)
if maxArgLen > 4:
    endIdx = 7
    print('Command line arguments will be processed as months.')
else:
    print('Command line arguments will be processed as years.')

"""
Create two sqlite databases to associate waveform miniseed files with phase data, including p and s pick times.
A separate database is created for complete E/N/Z sets and Z "singletons" so that the two datasets can be distributed separately. 
"""
create_main_db = not os.path.isfile('./dbs/NCEDC.db')
create_singleton_db = not os.path.isfile('./dbs/NCEDC_s.db')
con_main = sqlite3.connect('./dbs/NCEDC.db')
con_single = sqlite3.connect('./dbs/NCEDC_s.db')

# Phase_file column is so that all waveforms corresponding to a single file can be deleted in the event that phase file must be re-processed.
cur_main = con_main.cursor()
if create_main_db:
    cur_main.execute('CREATE TABLE phase(id INTEGER PRIMARY KEY AUTOINCREMENT, event_id, network, station, location, channel, P, S, ' \
                     'instrument_match, complete, waveform_path, phase_file)')

cur_single = con_single.cursor()
if create_singleton_db:
    cur_single.execute('CREATE TABLE phase(id INTEGER PRIMARY KEY AUTOINCREMENT, event_id, network, station, location, channel, P, S, ' \
                       'instrument_match, complete, waveform_path, phase_file)')

"""
Create obspy client to download waveforms from NCEDC dataselect web service (https://service.ncedc.org/fdsnws/dataselect/1/).
Rather than querying the NCEDC station web service (https://service.ncedc.org/fdsnws/station/1/) for available instruments at each station, simply request everything of interest.
See https://ds.iris.edu/ds/nodes/dmc/data/formats/seed-channel-naming/ for channel naming conventions.
"""
client = Client("NCEDC")
bands = ['D', 'E', 'S', 'H', 'B']
instruments = ['H', 'L', 'N', 'P'] # Removed G, M

"""
Traces returned from NCEDC are labeled only by network, station, and start/end times. Based on these values, this function finds the matching
phase data, including p and s pick times, associated a trace.
When matching incomplete traces (e.g., those shorter than expected), a tolerance is allowed for start time to account for the missing data.
Any phase data that cannot be matched with a trace is logged to a retry file for further investigation.
"""
def get_matching_phase_data(trace, phase_data, matched_phase_data, log_path):
    tr_network = trace.stats.network
    tr_station = trace.stats.station
    tr_start = trace.stats.starttime
    tr_expected_length = 60 * trace.stats.sampling_rate

    match_idxs = []
    for j, data in enumerate(phase_data):
        if data['network'] == tr_network and data['station'] == tr_station and abs(tr_start-data['start']) < 0.1:
            match_idxs.append(j)
    # If no matches are found, try again with a start time tolerance, assuming all missing data is from the beginning of the trace (e.g., the maximum possible start time offset).
    if len(match_idxs) == 0:
        for j, data in enumerate(phase_data):
            if data['network'] == tr_network and data['station'] == tr_station and abs(tr_start-data['start']) < abs(tr_expected_length- len(trace))/trace.stats.sampling_rate + 0.1:
                match_idxs.append(j)
        
        if len(match_idxs) == 0:
            message = 'Zero phase data matches found for trace:\n'
            message += str(trace)
            log_message(message, log_path)
            return 0
    
    if len(match_idxs) > 1: 
        message = 'Multiple phase data matches found for trace:\n'
        message += str(trace) + '\n'
        message += 'Matches:\n'
        for match_idx in match_idxs:
            message += str(phase_data[match_idx]) + '\n'
        log_message(message, log_path)
        return 0
    else:
        matched_phase_data[match_idxs[0]] = 1
        return phase_data[match_idxs[0]]

def write_Z_stream(stream, phase_data, matched_phase_data, log_path):
    message = "Attempting to write all-Z stream:\n"
    message += "(Note that some traces may not be written if they cannot be matched with phase data.)\n"
    message += stream.__str__(extended=True)
    log_message(message, log_path)

    for trace in stream:
        complete = len(trace) == 60 * trace.stats.sampling_rate
        phase_row = get_matching_phase_data(trace, phase_data, matched_phase_data, log_path)
        if phase_row == 0: # Zero or multiple phase data matches.
            continue
        singleton_stream = Stream(trace)
        waveform_path = f'{dir_path_singleton}/{phase_row["event_id"]}.{trace.stats.network}.{trace.stats.station}' \
                        f'.{trace.stats.location}.{trace.stats.channel}.mseed'
        singleton_stream.write(waveform_path, format='MSEED')
        sql_queue.put({'event_id': phase_row['event_id'], 'network': trace.stats.network, 'station': trace.stats.station, 
                                'location': trace.stats.location, 'channel': trace.stats.channel, 'P': phase_row['p_text'], 'S': phase_row['s_text'], 
                                'inst_match': phase_row['inst_match'], 'complete': int(complete), 'file': file, 'waveform_path': waveform_path, 'n_components': 1})

def write_ENZ_streams(streams, phase_data, matched_phase_data, log_path):
    for stream in streams:
        complete = True
        longest_trace = stream[0]
        for trace in stream:
            if len(trace) != 60 * trace.stats.sampling_rate:
                complete = False
            if len(trace) > len(longest_trace):
                longest_trace = trace
        phase_row = get_matching_phase_data(longest_trace, phase_data, matched_phase_data, log_path)
        if phase_row == 0: # Zero or multiple phase data matches.
            continue

        message = "Writing E/N/Z set:\n"
        message += stream.__str__(extended=True)
        log_message(message, log_path)
        
        waveform_path = f'{dir_path_3comp}/{phase_row["event_id"]}.{longest_trace.stats.network}.{longest_trace.stats.station}' \
                        f'.{longest_trace.stats.location}.{longest_trace.stats.channel[:2]}.mseed'
        stream.write(waveform_path, format='MSEED')
        sql_queue.put({'event_id': phase_row['event_id'], 'network': trace.stats.network, 'station': trace.stats.station, 
                          'location': trace.stats.location, 'channel': trace.stats.channel[:2], 'P': phase_row['p_text'], 'S': phase_row['s_text'], 
                          'inst_match': phase_row['inst_match'], 'complete': int(complete), 'file': file, 'waveform_path': waveform_path, 'n_components': 3})

# Lines from parsed hypoinverse files.
lines = []

"""
Process one batch of lines from a parsed hypoinverse file.
All waveforms associated with the batch are downloaded from the NCEDC dataselect web service via a single bulk request to reduce overhead.
Testing so far indicates that a batch size of 100 is safe (requests will succeed before timing out).
"""
def process_lines(file, log_path, batch_start, batch_end):
    unconsumed_traces = Stream()

    batch = lines[batch_start:batch_end]

    # Tracks any phase data that cannot be associated with returned trace(s) from NCEDC. A zero value indicates unmatched data.
    matched_phase_data = [0] * len(batch)

    log_message(f'Processing lines. File: {file}, line ids (n={len(batch)}): {batch_start}-{batch_end}.', log_path)

    bulk = ''
    phase_data = []
    for line_id, line in enumerate(batch):
        split = line.split('|')
        event_id = split[0]
        network = split[1]
        station = split[2]
        p_text = split[3]
        p_dt = UTCDateTime(p_text)
        s_text = split[4]
        s_dt = UTCDateTime(s_text)
        instrument_match = int(split[5])
        dt_mid = p_dt + (s_dt - p_dt)/2

        """
        Retrieve 60 second waveforms centered on the midpoint between the p and s picks.
        30 second randomly-positioned waveforms containing both picks are required to train PhaseNet, but retrieving 60 second waveforms
        with consistent centering potentially increases the applications of the produced dataset.
        As mentioned previously, data from all bands and instruments of interest are requested rather than determining a station's available instruments.
        """
        dt_start = dt_mid - 30
        dt_end = dt_mid + 30

        for band in bands:
            for instrument in instruments:
                bulk += f'{network} {station} * {band}{instrument}* {dt_start} {dt_end}\n'

        # We need a dictionary of the phase data used to construct the bulk query in order to associate the returned waveforms with p and s pick times, etc.
        phase_data.append({'event_id': event_id, 'network': network, 'station': station, 'start': dt_start, 
                           'p_text': p_text, 's_text': s_text, 'p_dt': p_dt, 's_dt': s_dt, 'inst_match': instrument_match, 'line_id': line_id})

    # Bulk request to NCEDC dataselect web service.
    try:
        stream = client.get_waveforms_bulk(bulk)
    except obspy.clients.fdsn.header.FDSNNoDataException:
        message = 'No data available for batch:\n'
        message += '\n'.join(lines)
        log_message(message, log_path)
        return


    message = 'Downloaded stream.\n'
    message += stream.__str__(extended=True)
    log_message(message, log_path)

    incomplete_traces = Stream()
    for trace in stream:
        tr_expected_length = 60 * trace.stats.sampling_rate
        if len(trace)/tr_expected_length < 0.8:
            incomplete_traces.append(trace)
            stream.remove(trace)
    
    if len(incomplete_traces) > 0:
        message = "The following traces were removed because they were < 80% expected length:\n"
        message += incomplete_traces.__str__(extended=True)
        log_message(message, log_path)

    n_traces = len(stream)
    i = 0
    while i<n_traces:
        """
        Get all traces from the next "chunk."
        A chunk consists of all traces returned from a single station and instrument.
        """
        first_trace = stream[i]
        cur_stream = Stream(traces=[first_trace])
        while True:
            i += 1
            if i == n_traces:
                break

            next_trace = stream[i]
            network_match = (cur_stream[0].stats.network == next_trace.stats.network)
            station_match = (cur_stream[0].stats.station == next_trace.stats.station)
            instrument_match = (cur_stream[0].stats.channel[:2] == next_trace.stats.channel[:2])

            # ToDo: add note about assumption.
            if (not network_match) or (not station_match) or (not instrument_match):
                break

            cur_stream.append(next_trace)
        
        """
        Segment the current chunk into E/N/Z sets, Z singletons, and "leftovers" (e.g., E/N traces without a matching Z trace).
        """
        allZ = True
        for trace in cur_stream:
            if trace.stats.channel[-1] != 'Z' and trace.stats.channel[-1] != '3':
                allZ = False
                break

        if allZ:
                write_Z_stream(cur_stream, phase_data, matched_phase_data, log_path)
        else:
            # This is O(n^2), but total computation time is still much smaller than waiting for the data to come back from NCEDC.
            # So not worth the effort to optimize it.
            consumed_trace_idxs = [0] * len(cur_stream)
            streams_3comp = []
            # This nested loop = O(n^2)
            idx_E = -1
            idx_NZ = []
            for j, tracej in enumerate(cur_stream):
                if tracej.stats.channel[-1] == 'E' or tracej.stats.channel[-1] == '1':
                    idx_E = j
                    expected_trace_length = 60 * tracej.stats.sampling_rate
                    length_error_E = expected_trace_length - len(tracej)
                    for k, tracek in enumerate(cur_stream):
                        length_error_cur = expected_trace_length - len(tracek)
                        length_error_max = max(length_error_E, length_error_cur)
                        # Already know everything in the chunk has same network, station, and instrument.
                        # So just need to check location and start time here.
                        # length_error_max/sampling rate = seconds, 0.1 to account for misalignment...
                        # discuss tolerance.. same as previously described.
                        if (j != k
                            and tracek.stats.location == tracej.stats.location 
                            and abs(tracek.stats.starttime - tracej.stats.starttime) < length_error_max/tracej.stats.sampling_rate + 0.1):
                            idx_NZ.append(k)
                    if len(idx_NZ)==2:
                        traces = [cur_stream[idx_E], cur_stream[idx_NZ[0]], cur_stream[idx_NZ[1]]]
                        streams_3comp.append(Stream(traces=traces))
                        consumed_trace_idxs[idx_E] = 1
                        consumed_trace_idxs[idx_NZ[0]] = 1
                        consumed_trace_idxs[idx_NZ[1]] = 1
                    else:
                        log_stream = Stream()
                        log_stream.append(cur_stream[idx_E])
                        for trace_id in idx_NZ:
                            log_stream.append(cur_stream[trace_id])
                
                        message = 'Incomplete set, or set with too many components:\n'
                        message += log_stream.__str__(extended=True)
                        log_message(message, log_path)
                    idx_NZ = []
                    
            # Write E/N/Z sets
            if len(streams_3comp) > 0:
                write_ENZ_streams(streams_3comp, phase_data, matched_phase_data, log_path)

            # Write extra, "trailing" Z singletons
            stream_z = Stream()
            for j, trace in enumerate(cur_stream):
                if consumed_trace_idxs[j] == 0 and (trace.stats.channel[-1] == 'Z' or trace.stats.channel[-1] == '3'):
                    consumed_trace_idxs[j] = 1
                    stream_z.append(trace)
            if len(stream_z) > 0:
                write_Z_stream(stream_z, phase_data, matched_phase_data, log_path)
            
            for j, val in enumerate(consumed_trace_idxs):
                if val == 0:
                    unconsumed_traces.append(cur_stream[j])
            
    """
    Log traces that could not be associated with an E/N/Z set.
    This should not contain any Z traces, since these are written as singletons.
    """
    if len(unconsumed_traces) > 0:
        message = 'Unconsumed traces:\n'
        message += unconsumed_traces.__str__(extended=True)
        log_message(message, log_path)

    """
    Log any phase data which could not be matched to waveforms.
    These can be retried later, perhaps with different window start and end times.
    """
    for i in range(len(matched_phase_data)):
        if matched_phase_data[i] == 0:
            retry_queue.put(phase_data[i]['line_id'])

# --------------- PROGRAM START ------------------------

phases_root = './phases'

# Create logging directories if they don't exist.
dir_log = './logs'
dir_download_log = './logs/download'
dir_retry_log = './logs/retry'
if not os.path.exists(dir_log):
    os.makedirs(dir_log)
if not os.path.exists(dir_download_log):
    os.makedirs(dir_download_log)
if not os.path.exists(dir_retry_log):
    os.makedirs(dir_retry_log)

# ToDo: create dbs, waveforms folders...

# Create worker threads. Efficiency might be higher with > 10 threads, but this has not been tested.
workers = []
for x in range(10):
    worker = DownloadWorker(batch_queue)
    worker.daemon = True
    worker.start()
    workers.append(worker)

for root, dirs, files in os.walk(phases_root):
    if len(dirs) > 0:
        continue
    for file in files:
        if not file.endswith('parsed.txt'):
            continue
        if file[0:endIdx] not in timespans_to_process: # ToDo: default behavior if no years or months specified?
            continue
        if UTCDateTime(int(file[:4]), int(file[5:7]), 1) < UTCDateTime(1984, 3, 1):
            print(f'Data is unavailable for {file}. NCSN began digital recordings in March 1984.\n')
            continue
        
        print(f'Processing {file}...')
        
        start = time.perf_counter()

        log_path = f'./logs/{file}_{UTCDateTime.now()}.txt'
        retry_path = f'./logs/retry/{file}_{UTCDateTime.now()}.txt'

        """
        Update log paths of all worker threads.
        We want a separate file for each thread and parsed phase file, so that the log files don't grow to hundreds of thousands of lines.
        """
        now_time = UTCDateTime.now()
        for i, worker in enumerate(workers):
            worker.log_path = f'./logs/download/thread{i}_{file}_{now_time}.txt'

        # Create separate waveform folder for each month, so that we dont have millions of miniseeds in one directory.
        dir_path_3comp = f'./waveforms/3comp/{file[0:7]}'
        dir_path_singleton = f'./waveforms/singleton/{file[0:7]}'
        if not os.path.exists(dir_path_3comp):
            os.makedirs(dir_path_3comp)
        if not os.path.exists(dir_path_singleton):
            os.makedirs(dir_path_singleton)

        with open(os.path.join(root,file), 'r') as f:
            lines = f.readlines()
            if len(lines) == 0:
                continue

        batch_start = 0
        batch_size = 100
        more=True
        while more:
            batch_end = batch_start+batch_size
            if batch_end > len(lines):
                batch_end = len(lines)
                more = False
            if (batch_end - batch_start) == 0:
                break
            batch_queue.put((file, batch_start, batch_end))
            batch_start = batch_end

        batch_queue.join()

        main_commit = False
        singleton_commit = False
        while True:
            try:
                item = sql_queue.get(timeout=5)
                query = f'insert into phase (event_id, network, station, location, channel, P, S, instrument_match, complete, waveform_path, phase_file) ' \
                        f'values("{item["event_id"]}", "{item["network"]}", "{item["station"]}", "{item["location"]}", "{item["channel"]}", ' \
                        f'"{item["P"]}", "{item["S"]}", {item["inst_match"]}, {item["complete"]}, "{item["waveform_path"]}", "{item["file"]}")'
                if item['n_components'] == 3:
                    cur_main.execute(query)
                    main_commit = True
                else:
                    cur_single.execute(query)
                    singleton_commit = True
            except queue.Empty:
                break
        
        if main_commit:
            con_main.commit()
        if singleton_commit:
            con_single.commit()

        retry_lines = []
        while True:
            try:
                retry_line_id = retry_queue.get(timeout=5)
                retry_lines.append(lines[retry_line_id])
            except queue.Empty:
                break
        if len(retry_lines) > 0:
            with(open(retry_path, 'a')) as retry_file:
                retry_file.writelines(retry_lines)

        end = time.perf_counter()
        print(f'\n{file} processed in {end-start} seconds.\n')