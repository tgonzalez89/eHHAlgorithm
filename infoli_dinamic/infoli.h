#include <stddef.h> //size_t

#ifndef INFOLI_H_INCLUDED
#define INFOLI_H_INCLUDED

// Cell properties
#define DELTA 0.05 // Duration of 1 simulation step (in miliseconds)
#define CONDUCTANCE 0.04 // Conductance for neighbors' coupling
#define C_M 1 // Capacitance
// Somatic conductances (mS/cm2)
#define G_NA_S 150 // Na gate conductance (=90 in Schweighofer code, 70 in paper) 120 too little
#define G_KDR_S 9.0 // K delayed rectifier gate conductance (alternative value: 18)
#define G_K_S 5 // Voltage-dependent (fast) potassium
#define G_LS 0.016 // Leak conductance (0.015)
// Dendritic conductances (mS/cm2)
#define G_K_CA 35 // Potassium gate conductance (35)
#define G_CAH 4.5 // High-threshold Ca gate conductance (4.5)
#define G_LD 0.016 // Dendrite leak conductance (0.015)
#define G_H 0.125 // H current gate conductance (1.5) (0.15 in SCHWEIGHOFER 2004)
// Axon hillock conductances (mS/cm2)
#define G_NA_A 240 // Na gate conductance (according to literature: 100 to 200 times as big as somatic conductance)
#define G_NA_R 0 // Na (resurgent) gate conductance
#define G_K_A 20 // K voltage-dependent
#define G_LA 0.016 // Leak conductance
// Cell morphology
#define P1 0.25 // Cell surface ratio soma/dendrite (0.2)
#define P2 0.15 // Cell surface ratio axon(hillock)/soma (0.1)
#define G_INT 0.13 // Cell internal conductance (0.13)
// Reversal potentials
#define V_NA 55 // Na reversal potential (55)
#define V_K -75 // K reversal potential
#define V_CA 120 // Ca reversal potential (120)
#define V_H -43 // H current reversal potential
#define V_L 10 // Leak current

typedef float fpn;

typedef struct Dend {
    fpn V_dend;
    fpn Hcurrent_q;
    fpn Calcium_r; // High-threshold calcium
    fpn Potassium_s; // Calcium-dependent potassium
    fpn I_CaH; // High-threshold calcium current
    fpn Ca2Plus; // Calcium concentration
} Dend;

void Dend_CalcNextState(Dend* dend, fpn V_soma, fpn I_app, fpn* NeighV_dend, fpn* Connectivity, size_t size);
void Dend_CalcHcurrent_q(fpn* Hcurrent_q, fpn V_dend);
void Dend_CalcCalcium_r(fpn* Calcium_r, fpn V_dend);
void Dend_CalcPotassium_s(fpn* Potassium_s, fpn Ca2Plus);
void Dend_CalcCa2Plus(fpn* Ca2Plus, fpn I_CaH);
fpn Dend_CalcI_cNeighbors(fpn V_dend, fpn* NeighV_dend, fpn* Connectivity, size_t size);
void Dend_CalcV_dendI_CaH(Dend* dend, fpn I_c, fpn I_app, fpn V_soma);

typedef struct Soma {
    fpn g_CaL;
    fpn V_soma;
    fpn Sodium_m; // Sodium (artificial)
    fpn Sodium_h;
    fpn Calcium_k; // Low-threshold calcium
    fpn Calcium_l;
    fpn Potassium_n; // Potassium (delayed rectifier)
    fpn Potassium_p;
    fpn Potassium_x_s; // Potassium (voltage-dependent)
} Soma;

void Soma_CalcNextState(Soma* soma, fpn V_dend, fpn V_axon);
void Soma_CalcCalcium_kCalcium_l(fpn* Calcium_k, fpn* Calcium_l, fpn V_soma);
void Soma_CalcSodium_mSodium_h(fpn* Sodium_m, fpn* Sodium_h, fpn V_soma);
void Soma_CalcPotassium_nPotassium_p(fpn* Potassium_n, fpn* Potassium_p, fpn V_soma);
void Soma_CalcPotassium_x_s(fpn* Potassium_x_s, fpn V_soma);
void Soma_CalcV_soma(Soma* soma, fpn V_dend, fpn V_axon);

typedef struct Axon {
    fpn V_axon;
    fpn Sodium_m_a; // Sodium (thalamocortical)
    fpn Sodium_h_a;
    fpn Potassium_x_a; // Potassium (transient)
} Axon;

void Axon_CalcNextState(Axon* axon, fpn V_soma);
void Axon_CalcSodium_m_aSodium_h_a(fpn* Sodium_m_a, fpn* Sodium_h_a, fpn V_axon);
void Axon_CalcPotassium_x_a(fpn* Potassium_x_a, fpn V_axon);
void Axon_CalcV_axon(Axon* axon, fpn V_soma);

typedef struct Neuron {
    Dend dend;
    Soma soma;
    Axon axon;
} Neuron;

void Neuron_CalcNextState(Neuron* neuron, fpn I_app, fpn* NeighV_dend, fpn* Connectivity, size_t size);
void Network_CalcNextState(Neuron* network, fpn* NeighV_dend, fpn* I_app_array, fpn** ConnectivityMatrix, size_t size);

#endif // INFOLI_H_INCLUDED

