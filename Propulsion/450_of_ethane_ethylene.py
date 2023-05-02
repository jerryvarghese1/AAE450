# -*- coding: utf-8 -*-
"""
Created on Tue Jan 24 17:40:03 2023

@author: john2530
"""

##450 ethane ethylene yuh

import numpy as np
import matplotlib.pyplot as plt
from pyfluids import Fluid, FluidsList, Input
import warnings
from CEA_Wrap import Fuel, Oxidizer, RocketProblem, HPProblem, ThermoInterface
from tabulate import tabulate
ethane = Fuel("C2H6", temp = 233)
ethylene = Fuel("C2H4", temp =233)
nitrousoxide = Oxidizer("N2O", temp=233)
pc = 101.16 #psi
E = 200

r = np.linspace(4,250,250)
e = np.linspace(4,250,250)
ispv1 = []
ivac_data = []
ivac_fro = []
temp_equ = []
temp_fro = []

# =============================================================================
# define problem with inital E guess of 40
# =============================================================================
j = np.linspace(4,250,250)
problem = RocketProblem(materials=[ethane, ethylene, nitrousoxide], \
                        pressure = pc, sup = E, analysis_type="frozen nfz=2")
for i in j:
    problem.set_o_f(i)                  #iterate through o:f values
    solution = problem.run()  
    cs = solution.cstar
    gam = solution.gamma
    pe = solution.p*100000
    cf = (2*gam**2/(gam-1)*(2/(gam+1))**((gam+1)/(gam-1))*(1-(pe/(pc*6894.76))**((gam-1)/gam)))**.5+(E*pe/(pc*6894.76))
    ivac_fro.append(cf*cs/9.81)
    temp_fro.append(solution.t)
    
max_ratio = j[np.argmax(ivac_fro)]         #find o:f yeilding highest isp-vacuum
print("ideal ratio o:f (frozen at throat) ", max_ratio)

    
problem = RocketProblem(materials=[ethane, ethylene, nitrousoxide], \
                        pressure = pc, sup = E)
for i in r:
    problem.set_o_f(i)                  #iterate through o:f values
    results = problem.run()  
    ispv1.append(results.ivac)
    temp_equ.append(results.t)

max_ratio = r[np.argmax(ispv1)]         #find o:f yeilding highest isp-vacuum
print("ideal ratio o:f (equ) ", max_ratio)

    
problem.set_o_f(max_ratio)              #set o:f to ideal
for i in e:
    problem.set_sup(i)                  #iterate through expansion ratios
    output = problem.run()
    
    ivac_data.append(output.ivac)
    
max_E = e[np.argmax(ivac_data)]
print('E at max isp', max_E)
problem.set_sup(max_ratio)
results = problem.run()

# =============================================================================
# find most common combustion products
# =============================================================================
product_tuples = results.prod_e.items()
n=0
for key, value in sorted(product_tuples, key=lambda value: value[1], reverse=True):
    print("{:20} : {:8.4%}".format(key, value))
    if n==5:
        break
    n+=1

# =============================================================================
# #PLOTS
# =============================================================================
#frozen flow ideal o:f based on optimizating chamber temperature
#equilibrium flow ideal o:f ratio based on optimizing ivac
plt.figure(1)
plt.rcParams['figure.facecolor']='white'
plt.plot(r,ispv1, color = 'xkcd:dark purple', label = 'equilibrium')
plt.plot(j, ivac_fro, color = 'xkcd:blood red', label = 'frozen at throat')
plt.grid()
plt.xlabel('o:f ratio')
plt.ylabel('vacuum isp (s)')
plt.title('o:f vs. ivac at E = 200')
plt.fill_between(r, ispv1, ivac_fro, color='xkcd:light grey')
plt.legend()

plt.figure(2)
plt.plot(e,ivac_data)
plt.grid()
plt.xlabel('expansion ratio')
plt.ylabel('vacuum isp (s)')
plt.title('E vs. ivac at ideal o:f')
plt.rcParams['figure.facecolor']='white'

plt.figure(20)
plt.plot(r, temp_equ, color='xkcd:dark blue', label='equilibrium')
plt.plot(j, temp_fro, color='xkcd:emerald', label='frozen at throat')
plt.legend()
plt.grid()
plt.xlabel('o:f')
plt.ylabel('exit temp (K)')
plt.title('exit temperature vs. o:f')
plt.fill_between(r, temp_equ, temp_fro, color='xkcd:light grey')
plt.show()


# =============================================================================
# mp prop estimates from Jerry's DV
# =============================================================================
max_ratio = 7.5

mpayload = np.linspace(0,94.3,20)     #kg (initial estimate from Mathew)
uranus_multi_moon = 2530                #m/s (moon tour DV, dependent on phasing)
isp_estimate = 310                      #(s) estimate for ethylene/ethane/N2O
mi = np.linspace(1000,919.17,20)           #kg estimate
drymass = [mpayload+mi]

#MR calc for uranus insertion burn
mprop_6 = (np.exp((uranus_multi_moon)/(9.81*isp_estimate)))*(mi+mpayload)-mi-mpayload
mr = np.exp((uranus_multi_moon)/(9.81*isp_estimate))
print(mr)
# =============================================================================
# calcularing mass of fuel and oxidizer, sizing at ambient
# =============================================================================
mw_n2o = 44.013     #g.mol
mw_c2h4 = 28.054    #g/mol
mw_c2h6 = 30.07     #g/mol

fuel_3 = (mprop_6/(max_ratio+1))
ox_3 = (mprop_6*max_ratio/(max_ratio+1))

# =============================================================================
# tank sizing
# =============================================================================
# =============================================================================
# c2h4_p = 1.142 #kg/m3 at 75 F, 1 bara
# c2h6_p =1.223 #kg/m3 at 75 F, 1 bara
# n2o_p = 1.141  #kg/m3 at 75 F, 1 bara
# =============================================================================

#c2h4_p = 1.142 #kg/m3 at 75 F, 1 bara
c2h4_p = 51.251 #kg/m3 30 bar guess
#c2h4_p35 = 66.8 #kg/m3
c2h4_p35 = 272.62 #kg/m3 at 75 F, 75 bar
#c2h4_p = 91.5 #kg/m3 at 75 F, 50 bara
c2h6_p = 51.251 #kg/m3 at 75 F, 30 bara
#c2h6_p35 = 66.8 #kg/m3 at 75 F, 35 bar
c2h6_p35 = 371.46 #kg/m3 at 75 F, 75 bar
#c2h6_p = 343.32 #kg/m3 at 75 F, 50 bara
#c2h6_p =1.223 #kg/m3 at 75 F, 1 bara
n2o_p = 65.851 #kg/m3 at 75 F, 30 bara
#n2o_p35 = 80.7
n2o_p35 = 795.78 #kg/m3 at 75 F, 75 bar
# n2o_p = 1.141  #kg/m3 at 75 F, 1 bara
#n2o_p = 143 #kg/m3 at 75 F, 50 bara

bulk_density = (c2h6_p35*c2h4_p35*2)/(c2h6_p35*2+c2h4_p35)
vfuel30 = fuel_3/((c2h4_p + c2h6_p)/2)
vox30 = ox_3/n2o_p

vfuel35 = fuel_3/((c2h4_p35 + c2h6_p35)/2)
vox35 = ox_3/n2o_p35

falconh_payload = 13.1*(5.2/2)*np.pi
starship_payload = 1100 

# =============================================================================
# print outputs
# =============================================================================
plt.figure(3)

plt.suptitle('pressurized to 30 bar ethane/ethylene')

plt.subplot(2,1,1)
plt.plot(drymass[0], vfuel30, label='fuel')
plt.plot(drymass[0], vox30, label='oxidizer')
plt.grid()
plt.xlabel('dry mass (kg)')
plt.ylabel('volume (m^3)')
plt.legend()

plt.subplot(2,1,2)
plt.plot(drymass[0], fuel_3, label='fuel')
plt.plot(drymass[0], ox_3, label='oxidizer')
plt.grid()
plt.xlabel('dry mass (kg)')
plt.ylabel('mass of fuel and ox (m)')
plt.legend()

plt.figure(4)

plt.suptitle('pressurized to 75 bar ethane/ethylene')

plt.subplot(2,1,1)
plt.plot(drymass[0], vfuel35, label='fuel')
plt.plot(drymass[0], vox35, label='oxidizer')
plt.grid()
plt.xlabel('dry mass (kg)')
plt.ylabel('volume (m^3)')
plt.legend()

plt.subplot(2,1,2)
plt.plot(drymass[0], fuel_3, label='fuel')
plt.plot(drymass[0], ox_3, label='oxidizer')
plt.grid()
plt.xlabel('dry mass (kg)')
plt.ylabel('mass of fuel and ox (m)')
plt.legend()

print('\nFalcon Heavy payload volume (m3): ', falconh_payload)
print('starship payload volume (m3)    : ', starship_payload)


# =============================================================================
# # tank mass calcs
# =============================================================================
# tank diameter must be less than starship diameter, use that to constrain 
# tanks must be less than starship fairing length

P_pa = 7500000 #pa (tank pressure)
SS_304_failure_stress = 515000000
SS_304_yield_stress = 250000000
fs_ss_failure = 1.1
copv_opp_press = 434000000#5000*6894.76
fs_copv = 1.1

#D = np.linspace(.2,2,50)
D=1

#calculating tank thinkness
t_ss_failure = fs_ss_failure*P_pa*(D/2)/SS_304_failure_stress
t_ss_yield = P_pa*(D/2)/SS_304_yield_stress
t_copv = fs_copv*P_pa*(D/2)/copv_opp_press

plt.figure(5)
plt.plot(D, t_ss_failure, label='304 sized at failure stress')
plt.plot(D, t_ss_yield, label = '304 sized at yeild')
plt.plot(D, t_copv, label = 'COPV')
plt.grid()
plt.legend()
plt.title('tank thickness vs. tank diameter')
plt.ylabel('thinkness (m)')
plt.xlabel('tank diameter (m)')
plt.show()

SS_304_rho = 7930 #kg/m3
copv_rho = 1750 #kg/m3

#D = .75

Df = .993
Do = 1
#calculating tank mass
#fuel and ox tank lengths
L_f = (vfuel35 - (4/3)*np.pi*(Df/2)**3)/(np.pi*(Df/2)**2)
L_o = (vox35 - (4/3)*np.pi*(Do/2)**3)/(np.pi*(Do/2)**2)

# SS at failure mass 
v_ssfail_f = (4*np.pi*(Df/2)**2+np.pi*Df*L_f)*t_ss_failure
v_ssfail_o = (4*np.pi*(Do/2)**2+np.pi*Do*L_o)*t_ss_failure
m_ssfail_f = v_ssfail_f * SS_304_rho
m_ssfail_o = v_ssfail_o * SS_304_rho

# SS at yeild mass
v_ssyeild_f = (4*np.pi*(Df/2)**2+np.pi*Df*L_f)*t_ss_yield
v_ssyeild_o = (4*np.pi*(Do/2)**2+np.pi*Do*L_o)*t_ss_yield
m_ssyeild_f = v_ssyeild_f * SS_304_rho
m_ssyeild_o = v_ssyeild_o * SS_304_rho

# copv mass
m_copv_f = P_pa*vfuel35/183000
m_copv_o = P_pa*vox35/183000

plt.figure(6)
plt.suptitle('propellant mass vs. tank mass')

plt.subplot(2,1,1)
plt.plot(fuel_3, m_ssfail_f, label='SS at failure')
plt.plot(fuel_3, m_ssyeild_f, label='SS at yield')
plt.plot(fuel_3, m_copv_f, label='copv')
plt.grid()
plt.title('fuel tank mass')
plt.xlabel('mass of fuel')
plt.ylabel('mass (kg)')
plt.legend(loc='upper left')

plt.subplot(2,1,2)
plt.plot(ox_3, m_ssfail_o, label='SS at failure')
plt.plot(ox_3, m_ssyeild_o, label='SS at yield')
plt.plot(ox_3, m_copv_o, label='copv')
plt.grid()
plt.title('ox tank mass')
plt.xlabel('mass of ox')
plt.ylabel('mass (kg)')
plt.legend(loc='upper left')

plt.subplots_adjust(left=.1,bottom=1.1,right=.9,top=2.5,wspace=.3,hspace=.4)

# =============================================================================
# #input diameter, get mass
# =============================================================================

#D = input('tank diameter (less than 20: ')
D = 1

t_ss_failure = fs_ss_failure*P_pa*(D/2)/SS_304_failure_stress
t_ss_yield = P_pa*(D/2)/SS_304_yield_stress

SS_304_rho = 7930 #kg/m3
copv_rho = 1750 #kg/m3

#calculating tank mass
#fuel and ox tank lengths
L_f = (vfuel35[19] - (4/3)*np.pi*(Df/2)**3)/(np.pi*(Df/2)**2)
L_o = (vox35[19] - (4/3)*np.pi*(Do/2)**3)/(np.pi*(Do/2)**2)

# SS at failure mass 
v_ssfail_f = (4*np.pi*(Df/2)**2+np.pi*Df*L_f)*t_ss_failure
v_ssfail_o = (4*np.pi*(Do/2)**2+np.pi*Do*L_o)*t_ss_failure
m_ssfail_f = v_ssfail_f * SS_304_rho
m_ssfail_o = v_ssfail_o * SS_304_rho

# SS at yeild mass
v_ssyeild_f = (4*np.pi*(Df/2)**2+np.pi*Df*L_f)*t_ss_yield
v_ssyeild_o = (4*np.pi*(Do/2)**2+np.pi*Do*L_o)*t_ss_yield
m_ssyeild_f = v_ssyeild_f * SS_304_rho
m_ssyeild_o = v_ssyeild_o * SS_304_rho

# copv mass
m_copv_f = P_pa*vfuel35[19]/183000
m_copv_o = P_pa*vox35[19]/183000

#new MR calc
MRnewcopv = (mpayload[19]+mi[19]+m_copv_f+m_copv_o+fuel_3[19]+ox_3[19])/(mpayload[19]+m_copv_f+m_copv_o+mi[19])
MRnewSSf = (mpayload[19]+mi[19]+ m_ssfail_o+ m_ssfail_f +fuel_3[19]+ox_3[19])/(mpayload[19]+m_ssfail_o+ m_ssfail_f+mi[19])
MRnewSSy = (mpayload[19]+mi[19]+ m_ssyeild_o+ m_ssyeild_f +fuel_3[19]+ox_3[19])/(mpayload[19]+m_ssyeild_o+ m_ssyeild_f+mi[19])

#new dV
copvdV = 9.81*310*np.log(MRnewcopv)
ssfdV = 9.81*310*np.log(MRnewSSf)
ssydV = 9.81*310*np.log(MRnewSSy)

print('Ethane tank data, (dry mass = 2000kg) D = ')
table = [['material', 'mass ox', 'mass fuel', 'mass ox tank', 'mass fuel tank', 'new MR', 'new dV'],\
         ['SS at failure', "{:.2f}".format(fuel_3[19]), "{:.2f}".format(ox_3[19]), "{:.2f}".format(m_ssfail_o), "{:.2f}".format(m_ssfail_f), "{:.2f}".format(MRnewSSf), "{:.2f}".format(ssfdV)],\
         ['SS at yield', "{:.2f}".format(fuel_3[19]), "{:.2f}".format(ox_3[19]), "{:.2f}".format(m_ssyeild_o), "{:.2f}".format(m_ssyeild_f), "{:.2f}".format(MRnewSSy), "{:.2f}".format(ssydV)],\
         ['COPV', "{:.2f}".format(fuel_3[19]), "{:.2f}".format(ox_3[19]), "{:.2f}".format(m_copv_o), "{:.2f}".format(m_copv_f), "{:.2f}".format(MRnewcopv), "{:.2f}".format(copvdV)]]
print(tabulate(table))

# =============================================================================
# phase change investigation
# =============================================================================
# =============================================================================
# ethylene_hvap = [[155, 156, 161, 167, 169.4, 175, 239, 258, 267],\
#                  [14.4, 14.4, 14.3, 14.1,13.544, 14, 13.6, 13.7, 14]]
#  
# plt.figure(7)
# plt.scatter(ethylene_hvap[0], ethylene_hvap[1])
# a,b = np.polyfit(ethylene_hvap[0], ethylene_hvap[1], 1)
# plt.plot(ethylene_hvap[0], a*ethylene_hvap[0]+b)
# plt.grid()
# plt.show()
# 
# #ethylene vapor pressure at -40C
# pvap_ethylene = P_pa/(np.exp((-13.6/.2964)-((1/297)-(1/233.5))))
# =============================================================================

# =============================================================================
# vapor state
# =============================================================================

#find liquid properties

t = np.linspace(233,298,100)
density_c2h4_array = []
density_c2h6_array = []
for i in t:
    ethane = Fluid(FluidsList.Ethane).with_state(Input.pressure(2.758e6), Input.temperature(233))
    ethylene = Fluid(FluidsList.Ethylene).with_state(Input.pressure(2.758e6), Input.temperature(233))

    rho_ethane = ethane.density
    rho_ethylene = ethylene.density
    
    density_c2h4_array.append(rho_ethylene)
    density_c2h6_array.append(rho_ethane)
    
plt.plot(t, density_c2h4_array)
plt.plot(t, density_c2h6_array)

ethane = Fluid(FluidsList.Ethane).with_state(Input.density(49.6), Input.temperature(233))

ethane.pressure

#mixture = Mixture([FluidsList.Ethane, FluidsList.Ethylene], [50,50]).with_state(Input.pressure(7.5e6), Input.temperature(297))

#mixture.density

from CoolProp.CoolProp import PhaseSI, PropsSI, get_global_param_string
PropsSI("D", "T", 298, "P", 7500000, "Ethane"), "kg/m^3"

PropsSI("P", "T", 233, "D", 368.69377862033497, "Ethane")


PropsSI("P", "T", 233, "D", 368.69377862033497, "Ethylene")
PropsSI("D", "T", 298, "P", 7500000, "N2O"), "kg/m^3"

PropsSI("P", "T", 233, "D", 787.8596935766616, "N2O")


# =============================================================================
# newton's solver 
# =============================================================================
mp = mprop_6[19]

def newton(guess):
    ########### Edit Value Below Depending on Estimates
    Isp = 300 # Specific Impulse [s], temporary value -> find max expansion ratio based on throat area and fairing area
    ########### Edit Value Below Depending on Estimates
    deltaV = 2.53 * 1000 # Delta V requirement [m/s]
    ############ Edit Value Below Depending on Estimates
    m_pl = 94.3 # Payload Mass [kg]; from Matthew's science payload
    ############ Edit Value Below Depending on Estimates
    m_i = 919.17+m_copv_o+m_copv_f # Inert Mass [kg]; Assume constant (Tank sizes change)
    ############

    g = 9.81 # Gravitational Constant [m/s^2]
    MR = np.exp(deltaV / (g * Isp)) # Mass Ratio, constant

    F = np.zeros(100) # Initialize vector for F
    F_prime = np.zeros(100) # Initialize vector for F_prime
    
    lamb = np.zeros(100) # Initalize vector for Propellant Mass Fraction
    lamb[0] = guess
    
    m_p = np.zeros(100) # Initalize vector for Propellant Mass
    m_p[0] = lamb[0] * m_i / (1 - lamb[0]) # Initial Propellant Mass [kg]
    
    # Looping
    index = 0
    error = 1 # Arbitrary Value to Start Loop
    while error >  0.0001:
        F[index] = m_p[index] - (lamb[index] * (MR - 1) * m_pl) / (1 - MR * (1 - lamb[index]))
        F_prime[index] = 1 + (m_i * m_pl * ((1 - MR) ** 2)) / ((m_p[index] + m_i * (1 - MR)) ** 2)
        m_p[index + 1] = m_p[index] - (F[index] / F_prime[index])
        error = np.absolute((m_p[index + 1] - m_p[index]) / (m_p[index]))
        lamb[index + 1] = m_p[index + 1] / (m_p[index + 1] + m_i) # Propellant Mass Fraction
        index = index + 1

    return(guess, index, m_p, lamb, F)

def results(plot):   
    # Ignore only runtime warnings
    warnings.filterwarnings("ignore", category=RuntimeWarning)

    # Looping to display outputs of all guesses
    guess = np.linspace(.1, 1, 100) # Create empty array for guess values

    for g in guess:
        initial_guess, index, m_p, lamb, F = newton(g)
        if not np.isnan(m_p[index]):
            if (m_p[index] != 0):
                print("----------")
                print("Initial Lambda Guess:", initial_guess)
                print("Iterations:", index)
                print("Propellant Mass:", m_p[index])
                print("Lambda:", lamb[index])
                print("----------")
                mp = m_p[index]
            
                if plot == 1:
                    plt.figure()
                    for i in range(index):
                        m_p_prev, F_prev = m_p[i], F[i]
                        plt.plot(m_p_prev, F_prev, marker='o', linestyle='-', color='blue')
                        plt.text(m_p_prev + 0.1, F_prev + 0.1, str(i), fontsize=12, color='red')
                        plt.xlabel('Propellant Mass [kg]') 
                        plt.ylabel('F')
                        plt.title(f"F vs. Propellant Mass [kg] for an Initial Lambda Guess of {initial_guess:.6f}");
                    plt.show()
    return(mp)
                    
mp = results(0) #outputs new prop mass

#solve for new copv mass
def newmi(mp, c2h4_p35, c2h6_p35, max_ratio, mi, mpayload):
    #new propellant volume
    fuel_3 = mp/(1+max_ratio)
    ox_3 = mp/(1+(1/max_ratio))
    vfuel = fuel_3/((c2h4_p35 + c2h6_p35)/2)
    vox = ox_3/n2o_p35

    # new copv mass
    m_copv_f = P_pa*vfuel/183000
    m_copv_o = P_pa*vox/183000

    
    return(m_copv_f, m_copv_o, vfuel, vox)

tankmass = newmi(mp, c2h4_p35, c2h6_p35, max_ratio, mi, mpayload)

m_copv_f = tankmass[0]
m_copv_o = tankmass[1]