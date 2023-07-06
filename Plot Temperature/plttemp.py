import pandas as pd
from matplotlib import pyplot as plt
import numpy as np

data1 = pd.read_excel('suhu.xlsx', sheet_name = 0)
data2 = pd.read_excel('suhu.xlsx', sheet_name = 1)

R0 = pd.DataFrame(data1, columns= ['Rad'])
T_bub = pd.DataFrame(data1, columns= ['Tmax'])
Mat = pd.DataFrame(data2, columns= ['Mat'])
T_melt = pd.DataFrame(data2, columns= ['Tmelt'])

plt.plot(R0,T_bub,'ko', label='$T_{max}$ Gelembung')
plt.ylabel('Suhu (Kelvin)', fontsize = 15)
plt.xlabel('Jari-jari awal ($\mu$m)', fontsize = 15)

plt.axhline(y=T_melt.iloc[0,0], color='gold', linestyle='-', label='$T_{melt}$ ' + str(Mat.iloc[0,0]))
plt.axhline(y=T_melt.iloc[1,0], color='tab:orange', linestyle='-', label='$T_{melt}$ ' + str(Mat.iloc[1,0]))
plt.axhline(y=T_melt.iloc[2,0], color='tab:red', linestyle='-', label='$T_{melt}$ ' + str(Mat.iloc[2,0]))
plt.axhline(y=T_melt.iloc[3,0], color='tab:blue', linestyle='-', label='$T_{melt}$ ' + str(Mat.iloc[3,0]))
plt.legend(loc="upper right")

xlab=[10,20,30]
plt.xticks(ticks=xlab, labels=xlab)


plt.show()