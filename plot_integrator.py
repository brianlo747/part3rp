import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import csv

run_data_dir = "output/"
run_data_list = ['01.txt','02.txt','03.txt','04.txt','05.txt','06.txt','07.txt','08.txt','09.txt','10.txt','11.txt','12.txt']
run_dfs = []

for run_data in run_data_list:
    dataframe = pd.read_csv(run_data_dir + run_data)
    dataframe.columns = ["height","temp","pres","qc","qr","qv","qi","qh","density"]
    run_dfs.append(dataframe)


plt.rcParams.update({'font.size': 16})


for i,df in enumerate(run_dfs):
    plt.close('all')
    fig1, ax1 = plt.subplots()
    ax1.set_xscale('log')
    ax1.set_xlabel(r'$q_{species}$ / $kgkg^{-1}$')
    ax1.set_ylabel('Height / m')
    ax1.plot(df["qc"],df["height"], label=r'$q_c$')
    ax1.plot(df["qr"],df["height"], label=r'$q_r$')
    ax1.plot(df["qv"],df["height"], label=r'$q_v$')
    ax1.plot(df["qi"],df["height"], label=r'$q_i$')
    ax1.plot(df["qh"],df["height"], label=r'$q_h$')
    ax1.tick_params(axis='y')
    ax1.set_xlim(1.0e-6,1.0e-1)
    ax1.set_ylim(0,15000, 200)
    if i==0 or i==6:
        ax1.legend(loc='upper right', prop={'size': 15}, bbox_to_anchor=(1.25, 1.03))
    ax1.tick_params(axis='both', which='major', pad=10)
    plt.savefig('output/' + run_data_list[i][0:2] + '_species.pdf', bbox_inches='tight')


plt.close('all')
fig1, ax1 = plt.subplots()
ax1.set_xlabel(r'Temperature / K', color='black')
ax1.set_ylabel('Height / m')

default_temp_precip = run_dfs[0]["temp"]
dalr_line = 298 - run_dfs[0]["height"] * 0.0098

ax1.plot(dalr_line,run_dfs[0]["height"], color='black', label=r'Dry Adiabatic Lapse Rate', linestyle='--')
ax1.plot(default_temp_precip,run_dfs[0]["height"], color='black', label=r'a')
ax1.plot(run_dfs[1]["temp"],run_dfs[1]["height"], linewidth=0.5, label=r'b')
ax1.plot(run_dfs[2]["temp"],run_dfs[2]["height"], linewidth=0.5, label=r'c')
ax1.plot(run_dfs[3]["temp"],run_dfs[3]["height"], linewidth=0.5, label=r'd')
ax1.plot(run_dfs[4]["temp"],run_dfs[4]["height"], linewidth=0.5, label=r'e')
ax1.plot(run_dfs[5]["temp"],run_dfs[5]["height"], linewidth=0.5, label=r'f')

ax1.set_xlim(230,298, 10)
ax1.set_ylim(0,15000, 200)
ax1.tick_params(axis='both', which='major', pad=10)
ax1.legend(loc="upper right", fontsize=9)
plt.savefig('output/temperature_absolute_precip.pdf', bbox_inches='tight')

plt.close('all')
fig1, ax1 = plt.subplots()
ax1.set_xlabel(r'Temperature / K', color='black')
ax1.set_ylabel('Height / m')

default_temp_no_precip = run_dfs[6]["temp"]

ax1.plot(dalr_line,run_dfs[0]["height"], color='black', label=r'Dry Adiabatic Lapse Rate', linestyle='--')
ax1.plot(default_temp_no_precip,run_dfs[6]["height"], color='black', label=r'a')
ax1.plot(run_dfs[7]["temp"],run_dfs[7]["height"], linewidth=0.5, label=r'b')
ax1.plot(run_dfs[8]["temp"],run_dfs[8]["height"], linewidth=0.5, label=r'c')
ax1.plot(run_dfs[9]["temp"],run_dfs[9]["height"], linewidth=0.5, label=r'd')
ax1.plot(run_dfs[10]["temp"],run_dfs[10]["height"], linewidth=0.5, label=r'e')
ax1.plot(run_dfs[11]["temp"],run_dfs[11]["height"], linewidth=0.5, label=r'f')

# #ax2.plot(T_dalr,height, '--', color='black')
ax1.set_xlim(230,298, 10)
ax1.set_ylim(0,15000, 200)
ax1.tick_params(axis='both', which='major', pad=10)
ax1.legend(loc="upper right", fontsize=9)
plt.savefig('output/temperature_absolute_no_precip.pdf', bbox_inches='tight')



#########
plt.close('all')
fig1, ax1 = plt.subplots()
ax1.set_xlabel(r'Temperature / K', color='black')
ax1.set_ylabel('Height / m')

default_temp_precip = run_dfs[0]["temp"]
dalr_line = 298 - run_dfs[0]["height"] * 0.0098

ax1.plot(dalr_line,run_dfs[0]["height"], color='black', label=r'Dry Adiabatic Lapse Rate', linestyle='--')
ax1.plot(default_temp_precip,run_dfs[0]["height"], color='black', label=r'a')
ax1.plot(run_dfs[1]["temp"],run_dfs[1]["height"], linewidth=0.5, label=r'b')
ax1.plot(run_dfs[2]["temp"],run_dfs[2]["height"], linewidth=0.5, label=r'c')
ax1.plot(run_dfs[3]["temp"],run_dfs[3]["height"], linewidth=0.5, label=r'd')
ax1.plot(run_dfs[4]["temp"],run_dfs[4]["height"], linewidth=0.5, label=r'e')
ax1.plot(run_dfs[5]["temp"],run_dfs[5]["height"], linewidth=0.5, label=r'f')

ax1.set_xlim(230,298, 10)
ax1.set_ylim(0,15000, 200)
ax1.tick_params(axis='both', which='major', pad=10)
ax1.legend(loc="upper right", fontsize=9)
plt.savefig('output/temperature_absolute_precip.pdf', bbox_inches='tight')

plt.close('all')
fig1, ax1 = plt.subplots()
ax1.set_xlabel(r'Temperature / K', color='black')
ax1.set_ylabel('Height / m')

default_temp_no_precip = run_dfs[6]["temp"]

ax1.plot(dalr_line,run_dfs[0]["height"], color='black', label=r'Dry Adiabatic Lapse Rate', linestyle='--')
ax1.plot(default_temp_no_precip,run_dfs[6]["height"], color='black', label=r'a')
ax1.plot(run_dfs[7]["temp"],run_dfs[7]["height"], linewidth=0.5, label=r'b')
ax1.plot(run_dfs[8]["temp"],run_dfs[8]["height"], linewidth=0.5, label=r'c')
ax1.plot(run_dfs[9]["temp"],run_dfs[9]["height"], linewidth=0.5, label=r'd')
ax1.plot(run_dfs[10]["temp"],run_dfs[10]["height"], linewidth=0.5, label=r'e')
ax1.plot(run_dfs[11]["temp"],run_dfs[11]["height"], linewidth=0.5, label=r'f')

# #ax2.plot(T_dalr,height, '--', color='black')
ax1.set_xlim(230,298, 10)
ax1.set_ylim(0,15000, 200)
ax1.tick_params(axis='both', which='major', pad=10)
ax1.legend(loc="upper right", fontsize=9)
plt.savefig('output/temperature_absolute_no_precip.pdf', bbox_inches='tight')


#########
min_idx_p = min([len(a_run) for a_run in run_dfs[0:6]]) - 1

default_temp_precip = run_dfs[0]["temp"][0:min_idx_p]

plt.close('all')
fig1, ax1 = plt.subplots()

ax1.set_xlabel(r'Temperature Difference / K', color='black')
ax1.set_ylabel('Height / m')

ax1.plot(default_temp_precip - default_temp_precip,run_dfs[0]["height"][0:min_idx_p], color='black', label=r'a')
ax1.plot(run_dfs[1]["temp"][0:min_idx_p] - default_temp_precip,run_dfs[1]["height"][0:min_idx_p], linewidth=0.5, label=r'b')
ax1.plot(run_dfs[2]["temp"][0:min_idx_p] - default_temp_precip,run_dfs[2]["height"][0:min_idx_p], linewidth=0.5, label=r'c')
ax1.plot(run_dfs[3]["temp"][0:min_idx_p] - default_temp_precip,run_dfs[3]["height"][0:min_idx_p], linewidth=0.5, label=r'd')
ax1.plot(run_dfs[4]["temp"][0:min_idx_p] - default_temp_precip,run_dfs[4]["height"][0:min_idx_p], linewidth=0.5, label=r'e')
ax1.plot(run_dfs[5]["temp"][0:min_idx_p] - default_temp_precip,run_dfs[5]["height"][0:min_idx_p], linewidth=0.5, label=r'f')
ax1.legend(loc='upper right', prop={'size': 15})
ax1.set_xlim(-3, 3, 0.5)
ax1.set_ylim(0,15000, 200)
ax1.tick_params(axis='both', which='major', pad=10)
plt.savefig('output/temperature_relative_precip.pdf', bbox_inches='tight')


#########
min_idx_np = min([len(a_run) for a_run in run_dfs[6:12]]) - 1

default_temp_no_precip = run_dfs[6]["temp"][0:min_idx_np]

plt.close('all')
fig1, ax1 = plt.subplots()

ax1.set_xlabel(r'Temperature Difference / K', color='black')
ax1.set_ylabel('Height / m')

ax1.plot(default_temp_no_precip - default_temp_no_precip,run_dfs[6]["height"][0:min_idx_np], color='black', label=r'a')
ax1.plot(run_dfs[7]["temp"][0:min_idx_np] - default_temp_no_precip,run_dfs[7]["height"][0:min_idx_np], linewidth=0.5, label=r'b')
ax1.plot(run_dfs[8]["temp"][0:min_idx_np] - default_temp_no_precip,run_dfs[8]["height"][0:min_idx_np], linewidth=0.5, label=r'c')
ax1.plot(run_dfs[9]["temp"][0:min_idx_np] - default_temp_no_precip,run_dfs[9]["height"][0:min_idx_np], linewidth=0.5, label=r'd')
ax1.plot(run_dfs[10]["temp"][0:min_idx_np] - default_temp_no_precip,run_dfs[10]["height"][0:min_idx_np], linewidth=0.5, label=r'e')
ax1.plot(run_dfs[11]["temp"][0:min_idx_np] - default_temp_no_precip,run_dfs[11]["height"][0:min_idx_np], linewidth=0.5, label=r'f')
ax1.legend(loc='upper right', prop={'size': 15})
ax1.set_xlim(-3, 3, 0.5)
ax1.set_ylim(0,15000, 200)
ax1.tick_params(axis='both', which='major', pad=10)
plt.savefig('output/temperature_relative_no_precip.pdf', bbox_inches='tight')

print(max(run_dfs[7]["temp"][0:5000] - default_temp_no_precip))
print(min(run_dfs[8]["temp"][0:5000] - default_temp_no_precip))

##########



#############################



all_qs = ["qc","qr","qv","qi","qh"]

for q in all_qs:
    plt.close('all')
    fig1, ax1 = plt.subplots()

    default_q_precip = run_dfs[0][q][0:min_idx_p]

    ax1.set_xlabel(r'$q_{' + q[1] + '}/q_{' + q[1] + ',ref}$', color='black')
    ax1.set_ylabel('Height / m')

    ax1.plot(default_q_precip / default_q_precip,run_dfs[0]["height"][0:min_idx_p], color='black', label=r'a')
    ax1.plot(run_dfs[1][q][0:min_idx_p] / default_q_precip,run_dfs[1]["height"][0:min_idx_p], linewidth=0.5, label=r'b')
    ax1.plot(run_dfs[2][q][0:min_idx_p] / default_q_precip,run_dfs[2]["height"][0:min_idx_p], linewidth=0.5, label=r'c')
    ax1.plot(run_dfs[3][q][0:min_idx_p] / default_q_precip,run_dfs[3]["height"][0:min_idx_p], linewidth=0.5, label=r'd')
    ax1.plot(run_dfs[4][q][0:min_idx_p] / default_q_precip,run_dfs[4]["height"][0:min_idx_p], linewidth=0.5, label=r'e')
    if q=="qi":
        print((run_dfs[3][q][0:min_idx_np] / default_q_precip)[5000])
        print((run_dfs[4][q][0:min_idx_np] / default_q_precip)[5000])
    if (q != "qi" or q != "qh"):
        ax1.plot(run_dfs[5][q][0:min_idx_p] / default_q_precip,run_dfs[5]["height"][0:min_idx_p], linewidth=0.5, label=r'f')

    ax1.set_xscale('log')
    ax1.set_xlim(1e-4, 1e+4)
    ax1.set_ylim(0,15000, 200)
    ax1.tick_params(axis='both', which='major', pad=10)
    ax1.legend()
    plt.savefig('output/' + q + '_relative_precip.pdf', bbox_inches='tight')

for q in all_qs:

    plt.close('all')
    fig1, ax1 = plt.subplots()

    default_q_no_precip = run_dfs[6][q][0:min_idx_np]

    ax1.set_xlabel(r'$q_{' + q[1] + '}/q_{' + q[1] + ',ref}$', color='black')
    ax1.set_ylabel('Height / m')

    ax1.plot(default_q_no_precip / default_q_no_precip,run_dfs[6]["height"][0:min_idx_np], color='black', label=r'a')
    ax1.plot(run_dfs[7][q][0:min_idx_np] / default_q_no_precip,run_dfs[7]["height"][0:min_idx_np], linewidth=0.5, label=r'b')
    ax1.plot(run_dfs[8][q][0:min_idx_np] / default_q_no_precip,run_dfs[8]["height"][0:min_idx_np], linewidth=0.5, label=r'c')
    ax1.plot(run_dfs[9][q][0:min_idx_np] / default_q_no_precip,run_dfs[9]["height"][0:min_idx_np], linewidth=0.5, label=r'd')
    ax1.plot(run_dfs[10][q][0:min_idx_np] / default_q_no_precip,run_dfs[10]["height"][0:min_idx_np], linewidth=0.5, label=r'e')
    if (q != "qi" or q != "qh"):
        ax1.plot(run_dfs[11][q][0:min_idx_np] / default_q_no_precip,run_dfs[11]["height"][0:min_idx_np], linewidth=0.5, label=r'f')

    ax1.set_xscale('log')
    ax1.set_xlim(1e-4, 1e+4)
    ax1.set_ylim(0,15000, 200)
    ax1.tick_params(axis='both', which='major', pad=10)
    ax1.legend()
    plt.savefig('output/' + q + '_relative_no_precip.pdf', bbox_inches='tight')




######################
min_idx_p = min([len(a_run) for a_run in run_dfs[0:6]]) - 1

default_dens_precip = run_dfs[0]["density"][0:min_idx_p]

plt.close('all')
fig1, ax1 = plt.subplots()

ax1.set_xlabel(r'Density Difference / $kgm^{-3}$', color='black')
ax1.set_ylabel('Height / m')

ax1.plot(default_dens_precip - default_dens_precip,run_dfs[0]["height"][0:min_idx_p], color='black', label=r'a')
ax1.plot(run_dfs[1]["density"][0:min_idx_p] - default_dens_precip,run_dfs[1]["height"][0:min_idx_p], linewidth=0.5, label=r'b')
ax1.plot(run_dfs[2]["density"][0:min_idx_p] - default_dens_precip,run_dfs[2]["height"][0:min_idx_p], linewidth=0.5, label=r'c')
ax1.plot(run_dfs[3]["density"][0:min_idx_p] - default_dens_precip,run_dfs[3]["height"][0:min_idx_p], linewidth=0.5, label=r'd')
ax1.plot(run_dfs[4]["density"][0:min_idx_p] - default_dens_precip,run_dfs[4]["height"][0:min_idx_p], linewidth=0.5, label=r'e')
ax1.plot(run_dfs[5]["density"][0:min_idx_p] - default_dens_precip,run_dfs[5]["height"][0:min_idx_p], linewidth=0.5, label=r'f')

ax1.set_xlim(-0.005, 0.005, 0.001)
ax1.set_ylim(0,15000, 200)
ax1.tick_params(axis='both', which='major', pad=10)
ax1.legend()
plt.savefig('output/density_relative_precip.pdf', bbox_inches='tight')


#########
min_idx_np = min([len(a_run) for a_run in run_dfs[6:12]]) - 1

default_dens_no_precip = run_dfs[6]["density"][0:min_idx_np]

plt.close('all')
fig1, ax1 = plt.subplots()

ax1.set_xlabel(r'Density Difference / $kgm^{-3}$', color='black')
ax1.set_ylabel('Height / m')

ax1.plot(default_dens_no_precip - default_dens_no_precip,run_dfs[6]["height"][0:min_idx_np], color='black', label=r'a')
ax1.plot(run_dfs[7]["density"][0:min_idx_np] - default_dens_no_precip,run_dfs[7]["height"][0:min_idx_np], linewidth=0.5, label=r'b')
ax1.plot(run_dfs[8]["density"][0:min_idx_np] - default_dens_no_precip,run_dfs[8]["height"][0:min_idx_np], linewidth=0.5, label=r'c')
ax1.plot(run_dfs[9]["density"][0:min_idx_np] - default_dens_no_precip,run_dfs[9]["height"][0:min_idx_np], linewidth=0.5, label=r'd')
ax1.plot(run_dfs[10]["density"][0:min_idx_np] - default_dens_no_precip,run_dfs[10]["height"][0:min_idx_np], linewidth=0.5, label=r'e')
ax1.plot(run_dfs[11]["density"][0:min_idx_np] - default_dens_no_precip,run_dfs[11]["height"][0:min_idx_np], linewidth=0.5, label=r'f')

ax1.set_xlim(-0.005, 0.005, 0.001)
ax1.set_ylim(0,15000, 200)
ax1.tick_params(axis='both', which='major', pad=10)
ax1.legend()
plt.savefig('output/density_relative_no_precip.pdf', bbox_inches='tight')

#
# plt.clf()
# fig1, ax1 = plt.subplots()
# ax1.set_xlabel(r'Density / $kgm^{-3}$', color='black')
# ax1.set_ylabel('Height / m')
# ax1.plot(density,height, color='black')
# ax1.set_xlim(0,1.4, 0.2)
# ax1.set_ylim(0,15000, 200)
# ax1.tick_params(axis='both', which='major', pad=10)
# plt.savefig('output/' + run_number + '_density.pdf', bbox_inches='tight')
