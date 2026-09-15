# -*- coding: utf-8 -*-
"""
Created on Tue Sep 24 11:18:25 2024

@author: Bipul Biswas EE24MTECH12003
"""
import numpy as np
import matplotlib.pyplot as plt
import math

#...................Constant.....................
q = 1.60217e-19
ef = 0
kt = 0.026*q
vt = 1E-10
v_max = vt+1
ni=8.6e15

#...................Data.....................
Nc_data = {'Si':2.8e25, 'Ge':1.04e25,'GaAs':4.7e23 }  # m^-3
Nv_data = {'Si':1.04e25, 'Ge':6.0e24,'GaAs':7.0e24}  # m^-3
Eg_Data = {'Si':1.12, 'Ge':0.66,'GaAs':1.43 }  # eV
#EA_Data= {'Si':4.05,  'Ge':4.0}  # eV
epsilon_Data = {'Si':11.7,'Ge':16.0,'GaAs':12.9} 

#..................Take Input....................

v_initial = float(input("Enter v_initial voltage: "))
v_final = float(input("Ente v_final voltage: "))
N_temp =0
Nodes = 0
Z = []
nd = []
nd_temp = []
na = []
rhos=[]
na_temp = []
#D_volt = []
d_volt = []
test1=[]
test2=[]
# Nc = 10e25
# Nv = 10e25
epsilon = []
delta = []
domain = int(input("Enter the number of Materials: "))
Length=[]  #just to capture total length
nodes=[]   #just to capture total nodes 
Acceptor_conc=[]
Donner_conc=[]
hole_conc = np.zeros(N_temp, dtype = float)
electron_conc = np.zeros(N_temp, dtype = float)

for d in range(domain):
    print("Materials: ", d+1)
    material=input("Enter Material(Si or Ge or GaAs)")
    len =float(input("Enter the length of the material in m^-3: "))
    Length.append(len)
    node = 500
    nodes.append(node)
    Na = float(input("Enter Na value in cm^-3: "))*1e6
    Acceptor_conc.append(Na)
    Nd = float(input("Enter Nd value in cm^-3: "))*1e6
    Donner_conc.append(Nd)
    
    
    ep = epsilon_Data[material] * 8.85e-12
    print("the epsilon of the Material epsilon_Data",ep) 
    
    Nc = Nc_data[material]
    print("the Conduction Band DOS of the Material",Nc) 
    
    Nv = Nv_data[material]
    print("the Valance Band DOS of the Material",Nv) 
    
    Eg = Eg_Data[material]
    print("the Energy Band of the Material",Eg) 
    
    

    Eg = Eg*q      #ev to V
    dels = len/node
    
    delta = delta + [dels]*node
    epsilon = epsilon + [ep]*node
    
    if Na>Nd:
        na = na + [Na]*node
        nd = nd + [Nd]*node
        v_guess = (-(1/q))*(Eg+kt*math.log(Na/Nv))
        z = (-(Eg -(kt*math.log(Nv/Na)))/q) 
        
    else:
        na = na + [Na]*node
        nd = nd + [Nd]*node
        v_guess = ((kt/q)*math.log((Nd/Nc)))
        z = ((kt*math.log((Nd/Nc)))/q)
        
    
    # if na[0]>nd[0]:
    #     v_guess = (-(1/q))*(Eg+kt*math.log(na[0]/Nv))
    #     z = ((Eg -(kt*math.log(Nv/na[0])))/q) 
    # elif nd[0]>na[0]:
    #     v_guess = ((kt/q)*math.log((nd[0]/Nc)))
    #     z=((kt*math.log((nd[0]/Nc)))/q)
        
    d_volt = d_volt + [v_guess]*node
    N_temp +=node
    na_temp = na_temp + [Na]*node
    nd_temp = nd_temp + [Nd]*node
    Z.append(z)
    #D_volt.append(d_volt)
total_nodes=sum(nodes) # sum of length of different materials
total_len=sum(Length)  # sum of total node of different materials
first_v = v_initial+Z[0]
last_v = v_final+Z[-1]
d_volt[0] = first_v
d_volt[N_temp-1] = last_v
hole_conc = np.zeros(N_temp, dtype = float)
electron_conc = np.zeros(N_temp, dtype = float)


    
    
while(v_max > vt):

      coef_matrix = np.identity(N_temp)
      sol = np.zeros(N_temp)

      for i in range(N_temp):
          if i > 0 and i < (N_temp-1):
              a = (epsilon[i] + epsilon[i+1])/delta[i]
              b = -1*((epsilon[i] + epsilon[i+1])/delta[i] + (epsilon[i-1] + epsilon[i])/delta[i-1])
              c = (epsilon[i-1] + epsilon[i])/delta[i-1]
              d = (delta[i] + delta[i-1])

              r = Nv * math.pow(math.e, (-q*d_volt[i] -Eg - ef)/kt)
              #test2.append(r)
              s = Nc * math.pow(math.e, (q*d_volt[i] + ef)/kt)
              #test1.append(s)
              
                       
              rho =( nd[i] - na[i] + r - s)*q
              
              rhos.append(rho)
              diff_rho = -1*(q*q/kt)*(r + s)

              f = a*d_volt[i+1] + b*d_volt[i] + c*d_volt[i-1] + d*rho
                 
              sol[i] = -f
              
              
              
              hole_conc[i] = r
              
              electron_conc[i] = s
              if nd[0]<na[0]:
                  hole_conc[0]=na[1]
                  electron_conc[0]=electron_conc[1]
                  
              elif nd[0]>na[0]:    
                  hole_conc[0]=hole_conc[1]
                  electron_conc[0]=nd[0]
             
              if nd[total_nodes-1]>na[total_nodes-1]:
                  hole_conc[total_nodes-1]=hole_conc[total_nodes-2]
                  electron_conc[total_nodes-1]=nd[total_nodes-1]
              
              elif nd[total_nodes-1]<na[total_nodes-1]:
                  hole_conc[total_nodes-1]=na[total_nodes-1]
                  electron_conc[total_nodes-1]=electron_conc[total_nodes-2]
              
                
              # making TriDiagonal  matrix
              for j in range(N_temp):
                  if j == (i-1):
                      coef_matrix[i][j] = c
                  if j == i:
                      coef_matrix[i][j] = b + d*diff_rho
                  if j == (i+1):
                      coef_matrix[i][j] = a

     
      A = coef_matrix
      B = sol
    # X = (A**-1)*B
      matrixInv = np.linalg.inv(A)
      X = np.matmul(matrixInv, B)

      print("A: ", A)
      print("B: ", B)
      print("X: ", X)

      max_temp = [0]*N_temp
      temp = []
      for k in range(N_temp):
          d_volt[k] = d_volt[k] + X[k]
          m = abs(X[k])
          temp.append(m)
          if max_temp<temp:
              max_temp=temp
      v_max =max(max_temp)
      print("Voltages: ", d_volt, "\n")

x_values = np.linspace(0, total_len,total_nodes)



electric_field = np.zeros(N_temp, dtype=float)
for i in range(1, N_temp - 1):
    electric_field[i] = -(d_volt[i + 1] - d_volt[i - 1]) / (2 * delta[i])

# Boundary Conditions
electric_field[0] = -(d_volt[1] - d_volt[0]) / delta[0]
electric_field[-1] = -(d_volt[-1] - d_volt[-2]) / delta[-1]



# Plot Potential Profile in a separate window

plt.figure(figsize=(10, 6))
plt.plot(x_values, d_volt, label='Potential', color='blue')
plt.title('Potential Profile')
plt.xlabel('Position (m)')
plt.ylabel('Voltage (V)')
plt.grid()
plt.legend()  
plt.show()

# Plot Electric Field in a separate window
plt.figure(figsize=(10, 6))
plt.plot(x_values, electric_field, label='Electric Field', color='orange')
plt.title('Electric Field Distribution')
plt.xlabel('Position (m)')
plt.ylabel('Electric Field (V/m)')
plt.grid()
plt.legend()
plt.show()

# Plot Electron and Hole Concentrations in the same window

plt.figure(figsize=(10, 6))
plt.plot(x_values, electron_conc, label='Electron Concentration', color='green')
plt.plot(x_values, hole_conc, label='Hole Concentration', color='red')
plt.title('Carrier Concentration Distribution')
plt.xlabel('Position (m)')
plt.ylabel('Concentration (m^-3)')
plt.grid()
plt.legend()
plt.show()

plt.figure(figsize=(10, 6))
plt.yscale('log')
plt.plot(x_values, electron_conc, label='Electron Concentration', color='green')
plt.plot(x_values, hole_conc, label='Hole Concentration', color='red')
plt.title('Carrier Concentration Distribution (Log Scale)')
plt.xlabel('Position (m)')
plt.ylabel('Concentration (m^-3)')
plt.grid()
plt.legend()
plt.show()

# Plot charge Concentrations in the same window
plt.figure(figsize=(10, 6))
plt.plot(x_values, rhos[3000:4000], label='Charge', color='magenta')
plt.title('Charge Distribution')
plt.xlabel('Position (m)')
plt.ylabel('Charge (columb/m^-3)')
plt.grid()
plt.legend()
plt.show()