import sqlite3
import numpy as np
from obspy import read, UTCDateTime
import matplotlib.pyplot as plt

con = sqlite3.connect('./dbs/NCEDC.db')
cur = con.cursor()
res = cur.execute("select count(*) from phase")
num_rows = res.fetchone()[0]

# rs = [1712, 5272, 6433, 6592, 6673, 6822, 6922, 7042, 7043, 7170, 7171, 7277, 7329, 7358, 7500, 7561, 7562, 8970, 11000, 11386, 12134, 12795, 
#       23694, 32097, 37759, 39023, 49461, 62179, 88306, 88639, 104171, 124865, 125218, 126031, 128325, 128326, 128784, 128785, 128792, 129092, 
#       129093, 138978, 138979, 139074, 139076, 139078, 139080, 139081, 139329, 139330, 139331, 139457, 139570, 142829, 143721, 149207, 149605, 
#       149606, 149722, 149830, 151891, 152419, 154211, 154594, 154669, 154670, 155010, 163475, 164107, 177320, 180059, 183339, 183340, 183936, 
#       184423, 188510, 188842, 193738, 194735, 204783, 208389, 209314, 213258, 213353, 216015, 217257, 217376, 224588, 231011, 231120, 231121, 
#       231535, 231902, 233037, 236287, 237589, 237845, 237895, 241173, 261508, 262133, 379765, 379766, 381688, 384386, 384639, 385275, 385276, 
#       399127, 399691, 399848, 400160, 400161, 432882, 433127, 433339, 433462, 433599, 433600, 433601, 474692, 474693, 474694, 474695, 474696]

# rs = [1344774, 1344775, 1344776] # ToDo: KEEP 1344775 AS AN EXAMPLE OF WHY <0.1 is needed rather than ==0.0

while True:
# for r in rs:
    r = np.random.randint(0, num_rows)
    print(r)
    # r=1712
    # r=5272

    query = cur.execute(f"select * from phase where id ={r}")
    res = query.fetchone()

    event_id = res[1]
    network = res[2]
    station = res[3]
    location = res[4]
    channel = res[5]
    p_dt = UTCDateTime(res[6])
    s_dt = UTCDateTime(res[7])
    waveform_path = res[10]

    st = read(waveform_path)

    p_sample = round((p_dt - st[0].stats.starttime) * st[0].stats.sampling_rate)
    s_sample = round((s_dt - st[0].stats.starttime) * st[0].stats.sampling_rate)

    fig, axes = plt.subplots(3,1,figsize = (6,5))
    axes[0].plot(st[0], c='k', linewidth=0.5)
    axes[1].plot(st[1], c='k', linewidth=0.5)
    axes[2].plot(st[2], c='k', linewidth=0.5)

    fig.suptitle(f'nc{event_id} -- {network}.{station}.{location}.{channel} -- n_samples: {len(st[0])}')#{\nround(np.mean(st[0].data),2)}/{round(np.std(st[0].data),2)}, {round(np.mean(st[1].data),2)}/{round(np.std(st[1].data),2)}, {round(np.mean(st[2].data),2)}/{round(np.std(st[2].data),2)}')

    for ax in axes:
        ax.axvline(p_sample, c='b')
        ax.axvline(s_sample, c='r')
        ax.axvline(int(len(st[0])/2), c='orange', alpha=0.5, linestyle='--')
        ax.margins(x=0)

    plt.tight_layout()
    plt.show()
