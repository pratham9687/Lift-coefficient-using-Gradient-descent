import pandas as pd
import math
import matplotlib.pyplot as plt

m = 1.275
A = 2.0*0.15

set_num = 0
with open("log_number.txt", "r") as f:
    set_num = int(f.readline())

if(set_num<6):
    set_num = 6

test_cl_and_cl0_flag = bool(int(input("Enter 1 for testing existing CL(AOA) and Cl0 or 0 for finding new: ")))


def AOA(V_t, Vv, pitch_rad, roll_rad):
    # print(roll_rad, Vv, V_t)
    try:
        aoa = math.atan(math.tan(pitch_rad)/math.cos(roll_rad)) - math.asin(Vv/(V_t*math.cos(roll_rad)))
    except(ValueError):
        aoa = math.atan(math.tan(pitch_rad)/math.cos(roll_rad)) - math.asin(1.0)
    return aoa

def F_lr(acc_x, acc_z, m, AOA):
    return m*math.cos(AOA)*acc_z + m*math.sin(AOA)*acc_x

def gama(roh, v_t, A):
    return 0.5*roh*v_t*v_t*A

def dcost_by_dCl0(flr, fld, gama):
    return -2*(flr - fld)*gama

def dcost_by_dCl(flr, fld, gama, AOA):
    return -2*(flr - fld)*gama*AOA

L = 1e-5

flight_data = pd.read_csv(f"Controller_data{set_num}.csv")

alpha_dataset = [0]
F_lr_dataset = [0]
F_ld_dataset = [0]
C_dataset = [0]
total_velo = flight_data['estimated_total_velocity']
vertical_velo = flight_data['vertical_velocity']
pitch = flight_data['actual_pitch']
roll = flight_data['actual_roll']
acc_x = flight_data['a_x']
acc_z = flight_data['a_z']

if test_cl_and_cl0_flag:
    cl = 0.0
    cl0 = 0.0
    with open("coefficients.txt", 'r') as f:
       _ = f.readline().split(",")
    cl = float(_[0])
    cl0 = float(_[1])
else:
    cl = 2.0
    cl0 = 2.0

epochi = 1000
global_epochi = 10
counter = 0
total_cl = 0
total_cl0 = 0

for i in range(80, 497):
    print(i)
    alpha = AOA(total_velo[i], vertical_velo[i], pitch[i]*3.141/180.0, roll[i]*3.141/180.0)
    flr = F_lr(acc_x[i], acc_z[i], m, alpha)
    Gama = gama(1.23, total_velo[i], A)
    fld = cl*Gama*alpha + cl0*Gama
    
    alpha_dataset.append(alpha)
    F_lr_dataset.append(flr)
    F_ld_dataset.append(fld)
    C_dataset.append(math.pow(flr - fld, 2.0))

    total_cl += cl
    total_cl0 += cl0
    counter += 1

    if test_cl_and_cl0_flag == False:
        for _ in range(epochi):
        # while abs(flr - fld) >= 0.2:
            fld = cl*Gama*alpha + cl0*Gama
            cl -= L*dcost_by_dCl(flr,fld,Gama,alpha)
            cl0 -= L*dcost_by_dCl0(flr,fld,Gama)

cl = total_cl/counter
cl0 = total_cl0/counter

if test_cl_and_cl0_flag == False:
    with open("coefficients.txt", 'w') as f:
        f.write(f"{cl},{cl0}")

print(cl, cl0)
fig, axes = plt.subplots(3, 1, figsize=(8, 6), sharex=True)  # 2 rows, 1 column
axes[0].plot(alpha_dataset, label="AOA")
axes[0].legend()
axes[0].grid()
axes[1].plot(F_lr_dataset, label="F_actual")
axes[1].plot(F_ld_dataset, label=f"F_derived, for cl={cl} and cl0={cl0}")
axes[1].legend()
axes[1].grid()
axes[2].plot(C_dataset, label="cost")
axes[2].legend()
axes[2].grid()
plt.show()