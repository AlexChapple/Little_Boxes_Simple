import matplotlib.pyplot as plt 
import numpy as np 
import cmath as cm
import matplotlib

import colours

# import data 
z_list = np.loadtxt("sigma_z.txt")
L_list = np.loadtxt("sigma_L.txt")
R_list = np.loadtxt("sigma_R.txt")
g2_list = np.loadtxt("g2.txt")

time_list = z_list[:,0]
z_list = z_list[:,1]
L_list = L_list[:,1]
R_list = R_list[:,1]
g2_list = g2_list[:,1]

z_list_analytical = []
L_list_analytical = []
R_list_analytical = []
g2_analytical_list = []

def sigma_z(t):

    Y = np.sqrt(2) * 2.3

    delta = (1/4) * cm.sqrt(1 - 8*Y**2)

    a = -1 * (1/(1 + Y**2))
    b = Y**2 * np.exp(-3*t/4) 
    c = np.cosh(delta*t) + (3/(4*delta))*np.sinh(delta*t)

    return a * (1 + b*c)

def sigma_LR(t, mode):

    p = 1

    if mode == "R":
        p = -1

    Y = np.sqrt(2) * 2.3
    delta = (1/4) * cm.sqrt(1 - 8*Y**2)

    a = p * (1.0j / np.sqrt(2)) * (Y / (1 + Y**2))
    b = np.exp(-3*t/4) * (np.cosh(delta*t) + (3/(4*delta))*np.sinh(delta*t))
    c = p * 1.0j * np.sqrt(2) * Y * np.exp(-3*t/4) * (1/(4*delta)) * np.sinh(delta*t)

    return np.imag(a * (1 - b) + c)

def g2_analytical(t):

    Y = np.sqrt(2) * 2.3
    delta = (1/4) * cm.sqrt(1 - 8*Y**2)

    a = np.exp(-3*t/4) 
    b = np.cosh(delta * t) + ((3/(4 * delta)) * np.sinh(delta * t))

    return 1 - a*b


### Prettify plots 
matplotlib.rcParams.update({'font.size': 24})

# fig1 = plt.figure(1) 
# fig1.set_size_inches(18.5, 10.5)
# ax = plt.axes()
# ax.spines['top'].set_visible(False)
# ax.spines['right'].set_visible(False)
# ax.spines['left'].set_color(colours.spanish_gray)
# ax.spines['bottom'].set_color(colours.spanish_gray)
# ax.tick_params(axis='x', colors=colours.spanish_gray)
# ax.tick_params(axis='y', colors=colours.spanish_gray)
# ax.xaxis.label.set_color(colours.spanish_gray)
# ax.yaxis.label.set_color(colours.spanish_gray)
# ax.set_facecolor(colours.cultured2)
# fig1.set_alpha(0)
# fig1.set_facecolor("none")
# plt.xlabel("Time (s)")
# plt.ylabel("$\langle \sigma_{z} (t) \\rangle$")
# plt.grid()
# plt.plot(time_list, z_list, label="numerical", c=colours.greek_red, linewidth=3)
# plt.plot(time_list, [sigma_z(t) for t in time_list], label="analytical", c=colours.greek_blue, linestyle="dashdot", linewidth=3)
# plt.legend()
# plt.savefig("sigma_z.pdf", dpi=600)


# fig2 = plt.figure(2) 
# fig2.set_size_inches(18.5, 10.5)
# ax2 = plt.axes()
# ax2.spines['top'].set_visible(False)
# ax2.spines['right'].set_visible(False)
# ax2.spines['left'].set_color(colours.spanish_gray)
# ax2.spines['bottom'].set_color(colours.spanish_gray)
# ax2.tick_params(axis='x', colors=colours.spanish_gray)
# ax2.tick_params(axis='y', colors=colours.spanish_gray)
# ax2.xaxis.label.set_color(colours.spanish_gray)
# ax2.yaxis.label.set_color(colours.spanish_gray)
# ax2.set_facecolor(colours.cultured2)
# fig2.set_alpha(0)
# fig2.set_facecolor("none")
# plt.xlabel("Time (s)")
# plt.ylabel("$\langle \sigma_{-} (t) \\rangle$")
# plt.grid()
# plt.plot(time_list, L_list, label="numerical", c=colours.greek_red, linewidth=3)
# plt.plot(time_list, [sigma_LR(t, "L") for t in time_list], label="analytical", c=colours.greek_blue, linestyle="dashdot", linewidth=3)
# plt.legend()
# plt.savefig("sigma_L.pdf", dpi=600)

# fig3 = plt.figure(3) 
# fig3.set_size_inches(18.5, 10.5)
# ax3 = plt.axes()
# ax3.spines['top'].set_visible(False)
# ax3.spines['right'].set_visible(False)
# ax3.spines['left'].set_color(colours.spanish_gray)
# ax3.spines['bottom'].set_color(colours.spanish_gray)
# ax3.tick_params(axis='x', colors=colours.spanish_gray)
# ax3.tick_params(axis='y', colors=colours.spanish_gray)
# ax3.xaxis.label.set_color(colours.spanish_gray)
# ax3.yaxis.label.set_color(colours.spanish_gray)
# ax3.set_facecolor(colours.cultured2)
# fig3.set_alpha(0)
# fig3.set_facecolor("none")
# plt.xlabel("Time (s)")
# plt.ylabel("$\langle \sigma_{+} (t) \\rangle$")
# plt.grid()
# plt.plot(time_list, R_list, label="numerical", c=colours.greek_red, linewidth=3)
# plt.plot(time_list, [sigma_LR(t, "R") for t in time_list], label="analytical", c=colours.greek_blue, linestyle="dashdot", linewidth=3)
# plt.legend()
# plt.savefig("sigma_R.pdf", dpi=600)

def modulo_func(z):

    a = z.real 
    b = z.imag 

    return np.sqrt(a**2 + b**2)

g2_list = [modulo_func(z) for z in g2_list]

fig4 = plt.figure(4)
plt.plot(time_list, g2_list)
plt.plot(time_list, [g2_analytical(t) for t in time_list])
plt.show()