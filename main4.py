import json

from cobra.io import read_sbml_model
from metabolism import Metabolism
from other import *
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

###### INITIALIZE METABOLISM ######

sbmlMdl = read_sbml_model("code2v2/files/e_coli_core.xml")
rcnTrgs = readJson("code2v2/files/rcnTrgs.json") 
ecoli = Metabolism(sbmlMdl, rcnTrgs)
# print([i._id for i in ecoli.trgMtbs])
 
optFlxInts = ecoli.calOptFlxInts()
ecoli.setOptFlxInts(optFlxInts)

trgRcnWghs = readJson("code2v2/files/trgRcnWghs.json") 
ecoli.setWghs(trgRcnWghs)

mtbCncs = readJson("code2v2/files/mtbCncs2.json")
ecoli.setCncs(mtbCncs)

trgRcnGec0s = readJson("code2v2/files/rcnGec0s.json")
ecoli.setTrgRcnGec0s(trgRcnGec0s)

trgRcnHgr0s = readJson("code2v2/files/rcnHgr0s3.json")
ecoli.setTrgRcnHgr0s(trgRcnHgr0s)

rcnGec0Errs = readJson("code2v2/files/rcnGec0Errs.json") 
ecoli.setTrgRcnGec0Errs(rcnGec0Errs)

rcnHgr0Errs = {key: value * 0.1 for key, value in readJson("code2v2/files/rcnHgr0s3.json").items()}
ecoli.setTrgRcnHgr0Errs(rcnHgr0Errs)

# dr = pd.DataFrame({'Id':[i._id for i in ecoli.trgRcns],
#                 'Name':[i._nm for i in ecoli.trgRcns],
#                 'Weight':[round(i.wgh, 3) for i in ecoli.trgRcns],
#                 'dg':[round(i.gec0, 2) for i in ecoli.trgRcns],
#                 'dp':[round(i.hgr0, 2) for i in ecoli.trgRcns],})
# dm = pd.DataFrame({'Id':[i._id for i in ecoli.trgMtbs],
#                 'Name':[i._nm for i in ecoli.trgMtbs],
#                 'Concentration (mM)':[round(i.cnc, 3) for i in ecoli.trgMtbs],})
# print(dr)
# print(dm.sort_values(by='Id').to_latex(index=False))

# #gecx1, gecy1, gecx2, gecy2 
# plotGecs(ecoli)
# # print(len(ecoli.trgRcns))

fsss = pd.read_csv("code2v2/files/fluxSample10k.csv", index_col=False)

###### RUN GENETIC ALGORITHM ON SAMPLES ######

fldPth = 'code2v2/samples4'

# minDstThr = 7.28 # from previous analysis

# for i in range(352, 10000):
    
#     smp = fsss.iloc[i]
#     data = [ecoli.genAlg(smp)]

#     if data[0][1] <= minDstThr: 

#         print(f'WARNING: IMINENT SAVING OF {fldPth}/sample{i}.json')
#         for tmp in range(10, 200, 12): data.append(ecoli.genAlg(smp, tmp))
#         for tmp in range(200, 400, 4): data.append(ecoli.genAlg(smp, tmp))
#         for tmp in range(400, 800, 12): data.append(ecoli.genAlg(smp, tmp))
        
#         with open(f"{fldPth}/sample{i}.json", "w") as json_file: 
#             json.dump(data, json_file)

###### ANALYSIS ######

smps1 = [readJson(f'{fldPth}/{filename}') for filename in os.listdir(fldPth)]
minDst0 = [smp[0][1] for smp in smps1] # minimal distances at T0

thrPer = 100
thrVal = np.percentile(minDst0, thrPer)
# h2oCnc = 38850
h2oCnc = 41000

# plotDstsHist(minDst0, thrPer, thrVal)

smps2 = [smp for smp in smps1 if smp[0][1] <= thrVal]

avgDevAvgDsts = [avgDevAvgDst(smp) for smp in smps2]
# plotAvgDevAvgDstsHist(avgDevAvgDsts, nrm=False)

refCnc0, refKin0 = calRefs0(smps2) # references
refOsm0 = sum(refCnc0) - h2oCnc # 41000 from [h2o]

smps3 = [smp[1:] for smp in smps2] # no need for T0 anymore

tmps = [smpT[0] for smpT in smps3[0]]
cncDevTs = [cncDevT(smp, refOsm0, h2oCnc) for smp in smps3] 
avgCncDevTs = avg(cncDevTs)
# plotCncDevs(tmps, avgCncDevTs, refOsm0, n=len(smps3))




kinStrAvg = []
kinStrAvgj = []
kinStr37Avg = []
kinStr11Avg = []
cncSumAvg = []

for sample in smps3:

    tmp = []
    minDst = []
    cncSum = []
    cncStr = []
    
    kinStr = []
    kinStrj = []
    kinStr37 = []
    kinStr11 = []

    for entry in sample:
        tmp.append(entry[0])
        minDst.append(entry[1])
        cncSum.append(sum(entry[2])-h2oCnc)
        
        kinStrT = []
        for rcn in range(len(entry[3])):
            kinStrT.append((entry[3][rcn]+entry[4][rcn])*ecoli.trgRcns[rcn].wgh/refKin0[rcn])
        
        kinStr37.append((entry[3][37]+entry[4][37])*ecoli.trgRcns[37].wgh/refKin0[37])
        kinStr11.append((entry[3][11]+entry[4][11])*ecoli.trgRcns[11].wgh/refKin0[11])

        kinStr.append(sum(kinStrT))
        kinStrj.append(kinStrT)

        cncStrT = []
        for rcn in range(len(entry[3])):
            cncStrT.append((entry[3][rcn]+entry[4][rcn])/refKin0[rcn])
        cncStr.append(sum(cncStrT))
    
    kinStrAvg.append(kinStr)
    kinStrAvgj.append(kinStrj)
    kinStr37Avg.append(kinStr37)
    kinStr11Avg.append(kinStr11)
    cncSumAvg.append(cncSum)

# print(len(kinStrAvg), len(kinStrAvg[0]), len(kinStrAvgj), len(kinStrAvgj[0]))


# afasd = [[sum(x) / len(kinStrAvg) for x in zip(*i)][1:] for i in kinStrAvgj]


glbKinStrAvg = [sum(x) / len(kinStrAvg) for x in zip(*kinStrAvg)][1:]
KinStr37Avg = [sum(x) / len(kinStrAvg) for x in zip(*kinStr37Avg)][1:]
KinStr11Avg = [sum(x) / len(kinStrAvg) for x in zip(*kinStr11Avg)][1:]
_max = max(glbKinStrAvg)
_min = 0
glbKinStrAvg = [(i-_min)/(_max-_min) for i in glbKinStrAvg]
KinStr37Avg = [(i-_min)/(_max-_min) for i in KinStr37Avg]
KinStr11Avg = [(i-_min)/(_max-_min) for i in KinStr11Avg]


# for i in range(1584):
#     print(i)
#     plt.plot(tmp[1:], kinStrAvg[i][1:])
#     plt.plot(tmp[1:], cncSumAvg[i][1:])

#     # plt.plot(tmp[1:], KinStr37Avg)
#     # plt.plot(tmp[1:], KinStr11Avg)
#     # plt.xlabel('Temperature (K)')
#     # plt.ylabel('Normalized stress')
#     plt.show()

glbKinStrAvgj = [[sum(x)/len(i) for x in zip(*i)] for i in zip(*kinStrAvgj)][1:]

# for i in range(len(glbKinStrAvgj[0])):
#     if glbKinStrAvgj[0][i] > 10: print(i)

# print(ecoli.trgRcns[11].__dict__)
print(len(smps3))

plt.rcParams.update({'font.size': 11})
fig, ax1 = plt.subplots(figsize=(8,6))

ax1.set_xlabel('Temperature (K)')
ax1.set_ylabel('Normalized Stress', color='black')
ax1.plot(tmp[1:], glbKinStrAvg, color='tab:blue', label='Global')
ax1.plot(tmp[1:], KinStr11Avg, color='tab:red', label='ACONTb')
ax1.plot(tmp[1:], KinStr37Avg, color='tab:green', label='GAPD')
ax1.set_xlim(-10, 810)
ax1.legend(loc='lower left', bbox_to_anchor=(0, 0.1))

tmps = range(10, 800, 10)
gecs1, gecLow1, gecHgh1 = plotMtbGec(ecoli, tmps, 'ACONTb')
gecs2, gecLow2, gecHgh2 = plotMtbGec(ecoli, tmps, 'GAPD')
ax2 = ax1.twinx() 
ax2.set_yscale('log')
ax2.set_ylim(0.001, 1000)
ax2.set_ylabel('Equilibrium Constant', color='black') 
ax2.plot(tmps, gecs1, linestyle='dashed', color='tab:red', label='ACONTb')
ax2.fill_between(tmps, gecLow1, gecHgh1, color='tab:red', alpha=0.3)
# ax2.plot(tmps, gecHgh1, color='tab:red')
# ax2.plot(tmps, gecLow1, color='tab:red')
ax2.plot(tmps, gecs2, linestyle='dashed', color='tab:green', label='GAPD')
ax2.fill_between(tmps, gecLow2, gecHgh2, color='tab:green', alpha=0.3)
# ax2.plot(tmps, gecLow2, color='tab:green')
# ax2.plot(tmps, gecHgh2, color='tab:green')
ax2.spines['right'].set_linestyle((0,(4,4)))
ax2.legend(loc='lower right', bbox_to_anchor=(1, 0.1))

ax3 = ax1.twiny()
ax3.set_xlabel('Temperature (ºC)') 
ax3.set_xlim(-283.15, 536.85)

ax1.xaxis.set_major_locator(AutoLocator())
ax1.xaxis.set_minor_locator(AutoMinorLocator())
ax3.xaxis.set_major_locator(AutoLocator())
ax3.xaxis.set_minor_locator(AutoMinorLocator())
ax1.xaxis.grid(which='both', linestyle='--', linewidth=0.5)

#fig.suptitle(f'Stress analysis for {fldPth}')
fig.tight_layout()  # otherwise the right y-label is slightly clipped
plt.show()



indexs = [2, 8, 9, 12]
labels = ['$s_1$', '$s_2$', '$s_3$', '$s_4$']

plt.rcParams.update({'font.size': 11})
fig, axs = plt.subplots(1, 2, figsize=(12,6))

for i in range(len(indexs)):
    axs[0].plot(tmp[1:], kinStrAvg[indexs[i]][1:], label=labels[i])
    axs[1].plot(tmp[1:], cncSumAvg[indexs[i]][1:], label=labels[i])

for i in range(len(axs)):

    axs[i].set_xlabel('Temperature (K)')
    axs[i].set_xlim(-10, 810)

    ax_ = axs[i].twiny()
    ax_.set_xlabel('Temperature (ºC)') 
    ax_.set_xlim(-283.15, 536.85)

    axs[i].xaxis.set_major_locator(AutoLocator())
    axs[i].xaxis.set_minor_locator(AutoMinorLocator())
    ax_.xaxis.set_major_locator(AutoLocator())
    ax_.xaxis.set_minor_locator(AutoMinorLocator())
    axs[i].xaxis.grid(which='both', linestyle='--', linewidth=0.5)

axs[0].set_ylabel('Kinetic Parameter Stress', color='black')
axs[1].set_ylabel('Metabolite Concentration Stress', color='black')

axs[0].legend()#loc='center left', bbox_to_anchor=(1, 0.5))
axs[1].legend()#loc='center left', bbox_to_anchor=(1, 0.5))
fig.tight_layout()

plt.subplots_adjust(left=0.06, right=0.98, wspace=0.2)
plt.show()