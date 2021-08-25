import matplotlib.pyplot as plt 
import numpy as np 
import cmath as cm

# import data 
z_list = np.loadtxt("sigma_z.txt")
L_list = np.loadtxt("sigma_L.txt")
R_list = np.loadtxt("sigma_R.txt")

time_list = z_list[:,0]
z_list = z_list[:,1]
L_list = L_list[:,1]
R_list = R_list[:,1]

z_list_analytical = []
L_list_analytical = []
R_list_analytical = []

def sigma_z(t):

    Y = np.sqrt(2) * 2.3

    delta = (1/4) * cm.sqrt(1 - 8*Y**2)

    a = -1 * (1/(1 + Y**2))
    b = Y**2 * np.exp(-3*t/4) 
    c = np.cosh(delta*t) + (3/(4*delta))*np.sinh(delta*t)

    return a * (1 + b*c)

def sigma_LR(t, mode):

    p = 1

    # if mode == "R":
    #     p = -1

    Y = np.sqrt(2) * 2.3
    delta = (1/4) * cm.sqrt(1 - 8*Y**2)

    a = p * (1.0j / np.sqrt(2)) * (Y / (1 + Y**2))
    b = np.exp(-3*t/4) * (np.cosh(delta*t) + (3/(4*delta))*np.sinh(delta*t))
    c = p * 1.0j * np.sqrt(2) * Y * np.exp(-3*t/4) * (1/(4*delta)) * np.sinh(delta*t)

    return np.imag(a * (1 - b) + c)

plt.figure(1)
plt.plot(time_list, z_list)
plt.plot(time_list, [sigma_z(t) for t in time_list])
plt.xlabel("Time")
plt.ylabel("Probability")
plt.show()

plt.figure(2)
plt.plot(time_list, L_list)
plt.plot(time_list, [sigma_LR(t, "L") for t in time_list])
plt.xlabel("Time")
plt.ylabel("Probability")
plt.show()

plt.figure(3)
plt.plot(time_list, R_list)
plt.plot(time_list, [sigma_LR(t, "R") for t in time_list])
plt.xlabel("Time")
plt.ylabel("Probability")
plt.show()