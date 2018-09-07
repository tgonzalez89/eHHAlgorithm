#!/usr/bin/env python3

################################
###                          ###
### Hodgkin-Huxley algorithm ###
### Tomás González Aragón    ###
###                          ###
################################

import math
import copy
import time

##### CONSTANTS #####
DELTA = 0.05 # Duration of 1 simulation step (in miliseconds)
CONDUCTANCE = 0.04 # Conductance for neighbors' coupling
C_M = 1 # Capacitance
# Somatic conductances (mS/cm2)
G_NA_S = 150 # Na gate conductance (=90 in Schweighofer code, 70 in paper) 120 too little
G_KDR_S = 9.0 # K delayed rectifier gate conductance (alternative value: 18)
G_K_S = 5 # Voltage-dependent (fast) potassium
G_LS = 0.016 # Leak conductance (0.015)
# Dendritic conductances (mS/cm2)
G_K_CA = 35 # Potassium gate conductance (35)
G_CAH = 4.5 # High-threshold Ca gate conductance (4.5)
G_LD = 0.016 # Dendrite leak conductance (0.015)
G_H = 0.125 # H current gate conductance (1.5) (0.15 in SCHWEIGHOFER 2004)
# Axon hillock conductances (mS/cm2)
G_NA_A = 240 # Na gate conductance (according to literature: 100 to 200 times as big as somatic conductance)
G_NA_R = 0 # Na (resurgent) gate conductance
G_K_A = 20 # K voltage-dependent
G_LA = 0.016 # Leak conductance
# Cell morphology
P1 = 0.25 # Cell surface ratio soma/dendrite (0.2)
P2 = 0.15 # Cell surface ratio axon(hillock)/soma (0.1)
G_INT = 0.13 # Cell internal conductance (0.13)
# Reversal potentials
V_NA = 55 # Na reversal potential (55)
V_K = -75 # K reversal potential
V_CA = 120 # Ca reversal potential (120)
V_H = -43 # H current reversal potential
V_L = 10 # Leak current

##### DEND #####
class Dend:
    def __init__(self, V_dend, Hcurrent_q, Calcium_r, Potassium_s, I_CaH, Ca2Plus):
        self.V_dend = V_dend
        self.Hcurrent_q = Hcurrent_q
        self.Calcium_r = Calcium_r
        self.Potassium_s = Potassium_s
        self.I_CaH = I_CaH
        self.Ca2Plus = Ca2Plus

    def CalcNextState(self, V_soma, iApp, NeighVdend, Connectivity):
        self.CalcHcurrent_q()
        self.CalcCalcium_r()
        self.CalcPotassium_s()
        self.CalcCa2Plus()
        Ic = self.IcNeighbors(NeighVdend, Connectivity)
        self.CalcV_dendI_CaH(Ic, iApp, V_soma)

    def CalcHcurrent_q(self): #Hcurrent_q, V_dend
        q_inf = 1 / (1 + math.exp(0.25 * (80 + self.V_dend)))
        tau_q = 1 / (math.exp(-0.086 * self.V_dend - 14.6) + math.exp(0.070 * self.V_dend - 1.87))
        dq_dt = (q_inf - self.Hcurrent_q) / tau_q
        nextHcurrent_q = DELTA * dq_dt + self.Hcurrent_q
        self.Hcurrent_q = nextHcurrent_q

    def CalcCalcium_r(self): #Calcium_r, V_dend
        alpha_r = 1.7 / (1 + math.exp((5 - self.V_dend) / 13.9))
        beta_r = 0.02 * (8.5 + self.V_dend) / math.expm1(0.2 * (8.5 + self.V_dend))
        r_inf = alpha_r / (alpha_r + beta_r)
        tau_r = 5 / (alpha_r + beta_r)
        dr_dt = (r_inf - self.Calcium_r) / tau_r
        nextCalcium_r = DELTA * dr_dt + self.Calcium_r
        self.Calcium_r = nextCalcium_r

    def CalcPotassium_s(self): #Potassium_s, Ca2Plus
        alpha_s = 0.00002 * self.Ca2Plus
        if alpha_s < 0.01:
            alpha_s = 0.01
        beta_s = 0.015
        s_inf = alpha_s / (alpha_s + beta_s)
        tau_s = 1 / (alpha_s + beta_s)
        ds_dt = (s_inf - self.Potassium_s) / tau_s
        nextPotassium_s = DELTA * ds_dt + self.Potassium_s
        self.Potassium_s = nextPotassium_s

    def CalcCa2Plus(self): #Ca2Plus, I_CaH
        dCa_dt = -3 * self.I_CaH - 0.075 * self.Ca2Plus
        nextCa2Plus = DELTA * dCa_dt + self.Ca2Plus
        self.Ca2Plus = nextCa2Plus

    def IcNeighbors(self, NeighVdend, Connectivity): #V_dend
        F_acc = 0
        V_acc = 0
        for i, currNeighVdend in enumerate(NeighVdend):
            V = self.V_dend - currNeighVdend
            f = V * math.exp(-0.01 * V * V)
            F_acc += f * Connectivity[i]
            V_acc += V * Connectivity[i]
        I_c = 0.8 * F_acc + 0.2 * V_acc
        return I_c

    def CalcV_dendI_CaH(self, I_c, I_app, V_soma):
        # Dentritic currents
        # Soma-dendrite interaction current I_sd
        I_sd = (G_INT / (1 - P1)) * (self.V_dend - V_soma)
        # Inward high-threshold Ca current I_CaH
        I_CaH = G_CAH * self.Calcium_r * self.Calcium_r * (self.V_dend - V_CA)
        # Outward Ca-dependent K current I_K_Ca
        I_K_Ca = G_K_CA * self.Potassium_s * (self.V_dend - V_K)
        # Leakage current I_ld
        I_ld = G_LD * (self.V_dend - V_L)
        # Inward anomalous rectifier I_h
        I_h = G_H * self.Hcurrent_q * (self.V_dend - V_H)

        dVd_dt = (I_app - (I_c + I_CaH + I_sd + I_ld + I_K_Ca + I_h)) / C_M
        nextV_dend = DELTA * dVd_dt + self.V_dend
        nextI_CaH = I_CaH
        self.V_dend = nextV_dend
        self.I_CaH = nextI_CaH


##### SOMA #####
class Soma:
    def __init__(self, g_CaL, V_soma, Sodium_m, Sodium_h, Calcium_k, Calcium_l, Potassium_n, Potassium_p, Potassium_x_s):
        self.g_CaL = g_CaL
        self.V_soma = V_soma
        self.Sodium_m = Sodium_m
        self.Sodium_h = Sodium_h
        self.Calcium_k = Calcium_k
        self.Calcium_l = Calcium_l
        self.Potassium_n = Potassium_n
        self.Potassium_p = Potassium_p
        self.Potassium_x_s = Potassium_x_s

    def CalcNextState(self, V_dend, V_axon):
        self.CalcCalcium_kCalcium_l()
        self.CalcSodium_mSodium_h()
        self.CalcPotassium_nPotassium_p()
        self.CalcPotassium_x_s()
        self.CalcV_soma(V_dend, V_axon)

    def CalcCalcium_kCalcium_l(self): #Calcium_k, Calcium_l, V_soma
        k_inf = 1 / (1 + math.exp((self.V_soma + 61) / -4.2))
        l_inf = 1 / (1 + math.exp((self.V_soma + 85.5) / 8.5))
        tau_k = 1
        tau_l = 35 + 20 * math.exp((self.V_soma + 160) / 30) / (1 + math.exp((self.V_soma + 84) / 7.3))
        dk_dt = (k_inf - self.Calcium_k) / tau_k;
        dl_dt = (l_inf - self.Calcium_l) / tau_l;
        nextCalcium_k = DELTA * dk_dt + self.Calcium_k;
        nextCalcium_l = DELTA * dl_dt + self.Calcium_l;
        self.Calcium_k = nextCalcium_k
        self.Calcium_l = nextCalcium_l

    def CalcSodium_mSodium_h(self): #Sodium_m, Sodium_h, V_soma
        m_inf = 1 / (1 + math.exp((30 + self.V_soma) / -5.5))
        h_inf = 1 / (1 + math.exp((70 + self.V_soma) / 5.8))
        tau_h = 3 * math.exp((40 + self.V_soma) / -33)
        dh_dt = (h_inf - self.Sodium_h) / tau_h
        nextSodium_m = m_inf
        nextSodium_h = DELTA * dh_dt + self.Sodium_h
        self.Sodium_m = nextSodium_m
        self.Sodium_h = nextSodium_h

    def CalcPotassium_nPotassium_p(self): #Potassium_n, Potassium_p, V_soma
        n_inf = 1 / (1 + math.exp(-0.1 * (3 + self.V_soma)))
        p_inf = 1 / (1 + math.exp((51 + self.V_soma) / 12))
        tau_n = 5 + 47 * math.exp((50 + self.V_soma) / 900)
        tau_p = tau_n
        dn_dt = (n_inf - self.Potassium_n) / tau_n
        dp_dt = (p_inf - self.Potassium_p) / tau_p
        nextPotassium_n = DELTA * dn_dt + self.Potassium_n
        nextPotassium_p = DELTA * dp_dt + self.Potassium_p
        self.Potassium_n = nextPotassium_n
        self.Potassium_p = nextPotassium_p

    def CalcPotassium_x_s(self): #Potassium_x_s, V_soma
        alpha_x_s = -0.13 * (25 + self.V_soma) / math.expm1(-0.1 * (25 + self.V_soma))
        beta_x_s = 1.69 * math.exp(-0.0125 * (35 + self.V_soma))
        x_inf_s = alpha_x_s / (alpha_x_s + beta_x_s)
        tau_x_s = 1 / (alpha_x_s + beta_x_s)
        dx_dt_s = (x_inf_s - self.Potassium_x_s) / tau_x_s
        nextPotassium_x_s = DELTA * dx_dt_s + self.Potassium_x_s
        self.Potassium_x_s = nextPotassium_x_s

    def CalcV_soma(self, V_dend, V_axon):
        # Somatic currents
        # Dendrite-soma interaction current I_ds
        I_ds = (G_INT / P1) * (self.V_soma - V_dend)
        # Inward low-threshold Ca current I_CaL
        I_CaL = self.g_CaL * self.Calcium_k * self.Calcium_k * self.Calcium_k * self.Calcium_l * (self.V_soma - V_CA)
        # Inward Na current I_Na_s
        I_Na_s = G_NA_S * self.Sodium_m * self.Sodium_m * self.Sodium_m * self.Sodium_h * (self.V_soma - V_NA)
        # Leakage current I_ls
        I_ls = G_LS * (self.V_soma - V_L)
        # Outward delayed potassium current I_Kdr
        I_Kdr_s = G_KDR_S * self.Potassium_n * self.Potassium_n * self.Potassium_n * self.Potassium_n * (self.V_soma - V_K)
        # I_K_s
        I_K_s = G_K_S * self.Potassium_x_s * self.Potassium_x_s * self.Potassium_x_s * self.Potassium_x_s * (self.V_soma - V_K)
        # Axon-soma interaction current I_as
        I_as = (G_INT / (1 - P2)) * (self.V_soma - V_axon)

        dVs_dt = -(I_CaL + I_ds + I_as + I_Na_s + I_ls + I_Kdr_s + I_K_s) / C_M
        nextV_soma = DELTA * dVs_dt + self.V_soma;
        self.V_soma = nextV_soma


##### AXON #####
class Axon:
    def __init__(self, V_axon, Sodium_m_a, Sodium_h_a, Potassium_x_a):
        self.V_axon = V_axon
        self.Sodium_m_a = Sodium_m_a
        self.Sodium_h_a = Sodium_h_a
        self.Potassium_x_a = Potassium_x_a

    def CalcNextState(self, V_soma):
        self.CalcSodium_m_aSodium_h_a()
        self.CalcPotassium_x_a()
        self.CalcV_axon(V_soma)

    def CalcSodium_m_aSodium_h_a(self): #Sodium_h_a, V_axon
        m_inf_a = 1 / (1 + math.exp((30 + self.V_axon) / -5.5))
        h_inf_a = 1 / (1 + math.exp((60 + self.V_axon) / 5.8))
        tau_h_a = 1.5 * math.exp((40 + self.V_axon) / -33)
        dh_dt_a = (h_inf_a - self.Sodium_h_a) / tau_h_a
        nextSodium_m_a = m_inf_a
        nextSodium_h_a = DELTA * dh_dt_a + self.Sodium_h_a
        self.Sodium_m_a = nextSodium_m_a
        self.Sodium_h_a = nextSodium_h_a

    def CalcPotassium_x_a(self): #Potassium_x_a, V_axon
        alpha_x_a = -0.13 * (25 + self.V_axon) / math.expm1(-0.1 * (25 + self.V_axon))
        beta_x_a = 1.69 * math.exp(-0.0125 * (35 + self.V_axon))
        x_inf_a = alpha_x_a / (alpha_x_a + beta_x_a)
        tau_x_a = 1 / (alpha_x_a + beta_x_a)
        dx_dt_a = (x_inf_a - self.Potassium_x_a) / tau_x_a
        nextPotassium_x_a = DELTA * dx_dt_a + self.Potassium_x_a
        self.Potassium_x_a = nextPotassium_x_a

    def CalcV_axon(self, V_soma):
        # Axonal currents
        # Sodium
        I_Na_a = G_NA_A * self.Sodium_m_a * self.Sodium_m_a * self.Sodium_m_a * self.Sodium_h_a * (self.V_axon - V_NA)
        # Leak
        I_la = G_LA * (self.V_axon - V_L)
        # Soma-axon interaction current I_sa
        I_sa = (G_INT / P2) * (self.V_axon - V_soma)
        # Potassium (transient)
        I_K_a = G_K_A * self.Potassium_x_a * self.Potassium_x_a * self.Potassium_x_a * self.Potassium_x_a * (self.V_axon - V_K)

        dVa_dt = -(I_K_a + I_sa + I_la + I_Na_a) / C_M
        nextV_axon = DELTA * dVa_dt + self.V_axon
        self.V_axon = nextV_axon


##### NEURON #####
class Neuron:
    def __init__(self, dend, soma, axon):
        self.dend = dend
        self.soma = soma
        self.axon = axon

    def CalcNextState(self, iApp, NeighVdend, Connectivity):
        prevV_dend = self.dend.V_dend
        prevV_soma = self.soma.V_soma
        prevV_axon = self.axon.V_axon
        self.dend.CalcNextState(prevV_soma, iApp, NeighVdend, Connectivity)
        self.soma.CalcNextState(prevV_dend, prevV_axon)
        self.axon.CalcNextState(prevV_soma)


##### NETWORK #####
class Network:
    def __init__(self, size, initStateNeuron, ConnectivityMatrix):
        self.size = size
        self.ConnectivityMatrix = ConnectivityMatrix
        self.neurons = []
        for _ in range(size):
            self.neurons.append(copy.deepcopy(initStateNeuron))

    def CalcNextState(self, iAppArray):
        NeighVdend = []
        for i in range(self.size):
            NeighVdend.append(self.neurons[i].dend.V_dend)
        for i in range(self.size):
            self.neurons[i].CalcNextState(iAppArray[i], NeighVdend, self.ConnectivityMatrix[i])


################################################
##################### MAIN #####################
################################################

##### INITIAL STATE #####
# Initial dendritic parameters
V_dend = -60
Hcurrent_q = 0.0337836 # H current
Calcium_r = 0.0112788 # High-threshold calcium
Potassium_s = 0.0049291 # Calcium-dependent potassium
I_CaH = 0.5 # High-threshold calcium current
Ca2Plus = 3.7152 # Calcium concentration
# Initial somatic parameters
g_CaL = 0.68 # Default arbitrary value but it should be randomized per cell
V_soma = -60
Sodium_m = 1.0127807 # Sodium (artificial)
Sodium_h = 0.3596066
Calcium_k = 0.7423159 # Low-threshold calcium
Calcium_l = 0.0321349
Potassium_n = 0.2369847 # Potassium (delayed rectifier)
Potassium_p = 0.2369847
Potassium_x_s = 0.1 # Potassium (voltage-dependent)
# Initial axonal parameters
V_axon = -60
Sodium_m_a = 0.003596066 # Sodium (thalamocortical)
Sodium_h_a = 0.9
Potassium_x_a = 0.2369847 # Potassium (transient)

initStateDend = Dend(V_dend, Hcurrent_q, Calcium_r, Potassium_s, I_CaH, Ca2Plus)
initStateSoma = Soma(g_CaL, V_soma, Sodium_m, Sodium_h, Calcium_k, Calcium_l, Potassium_n, Potassium_p, Potassium_x_s)
initStateAxon = Axon(V_axon, Sodium_m_a, Sodium_h_a, Potassium_x_a)
initStateNeuron = Neuron(initStateDend, initStateSoma, initStateAxon)

##### SIMULATION #####
networkSize = 1
ConnectivityMatrix = [[CONDUCTANCE]*networkSize for _ in range(networkSize)]
neuronNetwork = Network(networkSize, initStateNeuron, ConnectivityMatrix)

#simTime = 1 # In miliseconds
simSteps = 20 #math.ceil(simTime/DELTA)
iAppArray = [[0]*networkSize for _ in range(simSteps)]
print("Simulating", simSteps, "steps for a network of size", networkSize)
out = open("InferiorOlive_Output.txt", "w")
print("SimStep SimTime Neuron# iApp V_axon")
out.write("SimStep SimTime Neuron# iApp V_axon\n")
start = time.time()
for i in range(simSteps):
    neuronNetwork.CalcNextState(iAppArray[i])
    for j in range(networkSize):
        print('{0:7d}'.format(i+1)+" "+'{0:7.2f}'.format((i+1)*DELTA)+" "+'{0:7d}'.format(j)+" "+'{0:4.2f}'.format(iAppArray[i][j])+" "+'{0:.8f}'.format(neuronNetwork.neurons[j].axon.V_axon))
        out.write('{0:7d}'.format(i+1)+" "+'{0:7.2f}'.format((i+1)*DELTA)+" "+'{0:7d}'.format(j)+" "+'{0:4.2f}'.format(iAppArray[i][j])+" "+'{0:.8f}'.format(neuronNetwork.neurons[j].axon.V_axon)+"\n")
end = time.time()
out.close()
print("Time elapsed:", 1000*(end-start), "ms,", 1000*(end-start)/simSteps, "ms per step")

