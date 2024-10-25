import os
from obspy import UTCDateTime

def process_phase_file(dir, file):
    path = os.path.join(dir, file)
    lines = []
    event_header_line_ids = []
    with open(path, 'r') as f:
        lines = f.readlines()
    
    # Segment the file into separate events.
    next_line_is_header = True
    for i, line in enumerate(lines):
        if next_line_is_header:
            event_header_line_ids.append(i)
            next_line_is_header = False
        else:
            if line[0] == ' ': # This seems to work, but is there a better way to check for the end of an event block?
                next_line_is_header = True
    event_header_line_ids.append(len(lines))

    """ 
    Each iteration of the following loop parses a single event.
    p_picks : All recorded p picks for the event.
    s_picks : All recorded s picks for the event.
    pick_pairs : All pairs of p and s picks with matching event_id, network, and station. 
                We first look for picks recorded by the same instrument, but failing that we accept different instruments.
    Files are in Hypoinverse format (.phase). See the following documentation https://ncedc.org/ftp/pub/doc/cat5/ncsn.phase.txt.
    """ 
    pick_pairs = []
    for i in range(len(event_header_line_ids)-1):
        both = False
        pairmatch = False
        
        p_picks = []
        s_picks = []
        already_found = set()
        event_lines = lines[event_header_line_ids[i]:event_header_line_ids[i+1]]

        header_line = event_lines[0]
        event_id = header_line[138:146].strip()
        
        for line in event_lines[1:]:
            station = line[0:5].strip()
            network = line[5:9].strip()
            channel = line[9:13].strip()
            
            """
            If 'P' is present in column 15, the line represents a p pick. Likewise for 'S' in column 48.
            A single line can represent both a p and s pick.
            Duplicate (event_id, network, station) pick sets are rejected.
            The values for p and s seconds represent offsets from a shared datetime (down to the minute). These values can
            exceed 6000 (e.g., 60.00 sec) or be negative.
            """
            p = (line[14] == 'P')
            if p:
                p_sec = line[29:34].strip()
                factor = 1
                if '-' in p_sec:
                    factor = -1
                    p_sec = p_sec.replace('-', '')

                p_sec = '0'*(4-len(p_sec)) + p_sec # Pad value to four digits.
                p_sec = factor * float(f'{p_sec[0:2]}.{p_sec[2:4]}')

                p_dt = UTCDateTime(int(line[17:21].strip()), int(line[21:23].strip()), int(line[23:25].strip()), 
                                    int(line[25:27].strip()), int(line[27:29].strip()))
                p_dt = p_dt + p_sec # Add OR subtract seconds offset.
                p_dt_text = p_dt.isoformat()

            s = (line[47] == 'S')
            if s:
                s_sec = line[41:46].strip()
                factor = 1
                if '-' in s_sec:
                    factor = -1
                    s_sec = s_sec.replace('-', '')

                s_sec = '0'*(4-len(s_sec)) + s_sec # Pad value to four digits.
                s_sec = factor * float(f'{s_sec[0:2]}.{s_sec[2:4]}')

                s_dt = UTCDateTime(int(line[17:21].strip()), int(line[21:23].strip()), int(line[23:25].strip()), 
                                    int(line[25:27].strip()), int(line[27:29].strip()))
                s_dt = s_dt + s_sec # Add OR subtract seconds offset.
                s_dt_text = s_dt.isoformat()

            # If a line represents both a p and s pick, we can use it as a pick set.
            if p and s:
                if (event_id, network, station) not in already_found:
                    already_found.add((event_id, network, station))
                    pick_pairs.append([event_id, network, station, p_dt_text, s_dt_text, '1']) # Value of 1 (e.g., True) indicates instrument match.
                    both = True # ToDo: remove
            elif p:
                p_picks.append([network, station, channel, p_dt_text])
            elif s:
                s_picks.append([network, station, channel, s_dt_text])

        # Try to find a matching p pick for each s pick. S picks seem to be more rare, so this is most efficient.
        for s_pick in s_picks:
            if (event_id, s_pick[0], s_pick[1]) in already_found:
                continue
            p_match_found = False
            for p_pick in p_picks:
                if s_pick[0]==p_pick[0] and s_pick[1]==p_pick[1] and s_pick[2][:2]==p_pick[2][:2]:
                    already_found.add((event_id, s_pick[0], s_pick[1]))
                    pick_pairs.append([event_id, p_pick[0], p_pick[1], p_pick[3], s_pick[3], '1'])
                    p_match_found = True
                    pairmatch = True
                    break
            # Try to find matching p pick recorded by different instrument.
            if not p_match_found:
                for p_pick in p_picks:
                    if s_pick[0]==p_pick[0] and s_pick[1]==p_pick[1]:
                        already_found.add((event_id, s_pick[0], s_pick[1]))
                        pick_pairs.append([event_id, p_pick[0], p_pick[1], p_pick[3], s_pick[3], '0']) # Value of 0 (e.g., False) indicates different instrument.
                        pairmatch = True
                        break
    
            if both and pairmatch:
                print(f'BOTH! {file}: {event_id}')

    # Write pick sets to file.
    out_path = os.path.join(root, file.replace('phase', 'parsed.txt'))
    with open(out_path, 'w') as outfile:
        for pick_pair in pick_pairs:
            outfile.write('|'.join(pick_pair) + '\n')
           
phases_root = './phases'
processed_files = 0

for root, dirs, files in os.walk(phases_root):
    if len(dirs) > 0:
        continue
    for file in files:
        if file.endswith('.txt'):
            continue
        print(f'Processing file: {os.path.join(root, file)}. Processed so far: {processed_files}.')
        process_phase_file(root, file)
        processed_files += 1

print(f'Processed {processed_files} files.')