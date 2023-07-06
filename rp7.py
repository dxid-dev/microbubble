import numpy as np
from matplotlib import pyplot as plt
from scipy.integrate import odeint

def equation(r_init, t):
    R, dRdt = r_init
    d2Rdt2 = (P_g-P_0-P_a*np.cos(omega*t)-(2*sigma)/R-(4*mu*dRdt)/R+((2*sigma)/R_0+P_0-P_g)*(R_0/R)**(3*k)*(1-(((3*k)/c)*dRdt)))/(R*rho)-(3*dRdt**2)/(2*R)
    return dRdt, d2Rdt2

t = np.arange(0, 5e-5, 2.5e-10)                     # Jangka Waktu (s)
k = 5/3                                             # Eksponen dari hubungan barotropik
rho = 998.2                                         # Densitas Air (kg/m3)
mu = 8.9e-4                                         # Viskositas Air (Pa.s)
sigma = 0.0725                                      # Koefisien tegangan permukaan (N/m)
P_g = 2337                                          # Tekanan Gas (Pa)
P_0 = 1e5                                           # Tekanan Awal (Pa)
P_a = 1.155e5                                       # Tekanan Ultrasonik (Pa)
f = 26.5e3                                          # Frekuensi Ultrasonik (Hz)
omega = 2*np.pi*f                                   # Frekuensi Sudut (rad/s)
c = 1480                                            # Cepat Rambat Gelombang Akustik di Air (m/s)
T_0 = 300                                           # T_b air (Kelvin)
R_0 = 10e-6                                         # Jari-jari awal (m)
u_0 = 0                                             # Inisiasi dR/dT

R_s = odeint(equation, [R_0, u_0], t)

R = R_s[:,0]*(10**6)                                # Jari-jari (µm)
R_max = max(R)                                      # Jari-jari maksimum (µm)
mt = t*10**6                                        # Waktu (µs)
V = R_s[:,1]                                        # Kecepatan (µm/µs) 
P_in = [((P_0-P_g+((2*sigma)/R_0))*(R_0/x)**(3*k))/1000 for x in R_s[:,0]]
P_ex = (P_a*np.cos(omega*t))
T_b = [T_0*(R_max/x)**(1*(k-1)) for x in R]

RR0 = [x / (R_0*(10**6)) for x in R]

fig, axs = plt.subplots(2, 2)
fig.suptitle("Dinamika Gelembung dengan Persamaan Rayleigh–Plesset\n $R_0$ = " + str(R_0*10**6) + " $\mu$m, $P_a$ = " + str(P_a/1000) + " kPa, $f_{ultrasound}$ = " + str(f/1000) + " kHz")
fig.tight_layout()

# R vs t
axs[0, 0].plot(mt, RR0, linewidth = 1.5)
axs[0, 0].set(ylabel='Jari-jari gelembung - $R/R_0$', xlabel='Waktu ($\mu$s)')
axs[0, 0].yaxis.label.set_color('tab:blue')
axs[0, 0].tick_params(colors='tab:blue', axis='y')
axs[0, 0].grid()

# P external vs t
ax002 = axs[0, 0].twinx()
ax002.plot(mt, P_ex/1000, '--', color = "black", linewidth = 1)
ax002.spines['left'].set_color('tab:blue')
ax002.spines['left'].set_lw(1.5)
ax002.set(ylabel='Tekanan ultrasound (kPa)', xlabel='Waktu ($\mu$s)')

# dR/dt vs t
axs[0, 1].plot(mt, V, color = "black", linewidth = 1.5, label="$\dot{R}$")
axs[0, 1].set(ylabel='Kecepatan (m/s)', xlabel='Waktu ($\mu$s)')
axs[0, 1].legend(prop={'size':15})
axs[0, 1].grid()

# P internal vs t
axs[1, 0].plot(mt, P_in, color = "tab:green", linewidth = 1.5)
axs[1, 0].set(ylabel='Tekanan internal gelembung (kPa)', xlabel='Waktu ($\mu$s)')
axs[1, 0].yaxis.label.set_color('tab:green')
axs[1, 0].tick_params(colors='tab:green', axis='y')
axs[1, 0].grid()

# T vs t
axs[1, 1].plot(mt, T_b, color = "tab:orange", linewidth = 1.5)
axs[1, 1].set(ylabel='Suhu (Kelvin)', xlabel='Waktu ($\mu$s)')
axs[1, 1].yaxis.label.set_color('tab:orange')
axs[1, 1].tick_params(colors='tab:orange', axis='y')
axs[1, 1].grid()

# print('Pressure, Temprature:')
# print(max(P_in))
# print(max(T_b))

plt.show()