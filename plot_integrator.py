import matplotlib.pyplot as plt
import csv

time = []
temp = []
pres = []
qc = []
qr = []
qv = []
qi = []
qh = []

with open('output/test.txt','r') as csvfile:
    plots = csv.reader(csvfile, delimiter=',')
    for row in plots:
        time.append(float(row[0]))
        temp.append(float(row[1]))
        pres.append(float(row[2]))
        qc.append(float(row[3]))
        qr.append(float(row[4]))
        qv.append(float(row[5]))
        qi.append(float(row[6]))
        qh.append(float(row[7]))

fig, ax1 = plt.subplots()
ax1.set_xscale('log')
ax1.set_xlabel(r'$q_{species}$')
ax1.set_ylabel('Time / s')
ax1.plot(qc,time, label=r'$q_c$')
ax1.plot(qr,time, label=r'$q_r$')
ax1.plot(qv,time, label=r'$q_v$')
ax1.plot(qi,time, label=r'$q_i$')
ax1.plot(qh,time, label=r'$q_h$')
ax1.tick_params(axis='y')
ax1.set_xlim(1.0e-8,1)
ax1.legend()

ax2 = ax1.twiny()
color = 'tab:red'
ax2.set_xlabel(r'Temperature / K', color=color)
ax2.plot(temp,time, label=r'T', color=color)
ax2.tick_params(axis='x', labelcolor=color)
ax2.legend()
plt.savefig('output/test.png', bbox_inches='tight')
