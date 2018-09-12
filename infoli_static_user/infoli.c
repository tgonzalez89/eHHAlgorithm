#include <math.h>
#include <omp.h>
#include <stdio.h>
#include "infoli.h"

/// /// /// NETWORK /// /// ///

// Calculate the next state of the network
void Network_CalcNextState(size_t network_size, Neuron network[], fpn I_app_array[], fpn ConnectivityMatrix[][network_size]) {
    fpn NeighV_dend[network_size];
    for (size_t i = 0; i < network_size; i++) {
        NeighV_dend[i] = network[i].dend.V_dend;
    }
    // Calculate the next state for each neuron
    #pragma omp parallel for
    for (size_t i = 0; i < network_size; i++) {
    	#ifdef DEBUG
        printf("Thread number: %llu, i=%llu\n", (unsigned long long)omp_get_thread_num(), (unsigned long long)i);
        #endif
        fpn prevV_dend = network[i].dend.V_dend;
        fpn prevV_soma = network[i].soma.V_soma;
        fpn prevV_axon = network[i].axon.V_axon;
        network[i].dend = Dend_CalcNextState(network[i].dend, prevV_soma, I_app_array[i], network_size, NeighV_dend, ConnectivityMatrix[i]);
        network[i].soma = Soma_CalcNextState(network[i].soma, prevV_dend, prevV_axon);
        network[i].axon = Axon_CalcNextState(network[i].axon, prevV_soma);
    }
}

/// /// /// DEND /// /// ///

// Calculate the next state of the dentrite
Dend Dend_CalcNextState(Dend dend, fpn V_soma, fpn I_app, size_t network_size, fpn NeighV_dend[], fpn Connectivity[]) {
    fpn_pair temp;
    dend.Hcurrent_q = Dend_CalcHcurrent_q(dend.Hcurrent_q, dend.V_dend);
    dend.Calcium_r = Dend_CalcCalcium_r(dend.Calcium_r, dend.V_dend);
    dend.Potassium_s = Dend_CalcPotassium_s(dend.Potassium_s, dend.Ca2Plus);
    dend.Ca2Plus = Dend_CalcCa2Plus(dend.Ca2Plus, dend.I_CaH);
    fpn I_c = Dend_CalcI_cNeighbors(dend.V_dend, network_size, NeighV_dend, Connectivity);
    temp = Dend_CalcV_dendI_CaH(dend, I_c, I_app, V_soma);
    dend.V_dend = temp.first;
    dend.I_CaH = temp.second;
    return dend;
}

// Calculate the next state for the dentrite's Hcurrent_q
fpn Dend_CalcHcurrent_q(fpn Hcurrent_q, fpn V_dend) {
    fpn const div_four = 0.25;

    fpn q_inf = 1 / (1 + expf((V_dend + 80) * div_four));
    fpn tau_q = 1 / (expf(-0.086 * V_dend - 14.6) + expf(0.070 * V_dend - 1.87));
    fpn dq_dt = (q_inf - Hcurrent_q) / tau_q;
    fpn result = DELTA * dq_dt + Hcurrent_q;
    return result;
}

// Calculate the next state for the dentrite's Calcium_r
fpn Dend_CalcCalcium_r(fpn Calcium_r, fpn V_dend) {
    fpn const div_five = 0.2;
    fpn const div_thirteen = 1 / 13.9;

    fpn alpha_r = 1.7 / (1 + expf(-(V_dend - 5) * div_thirteen));
    fpn beta_r = 0.02 * (V_dend + 8.5) / (expf((V_dend + 8.5) * div_five) - 1);
    fpn r_inf = alpha_r / (alpha_r + beta_r);
    fpn tau_r = 5 / (alpha_r + beta_r);
    fpn dr_dt = (r_inf - Calcium_r) / tau_r;
    fpn result = DELTA * dr_dt + Calcium_r;
    return result;
}

// Calculate the next state for the dentrite's Potassium_s
fpn Dend_CalcPotassium_s(fpn Potassium_s, fpn Ca2Plus) {
    fpn alpha_s = ((0.00002 * Ca2Plus) < 0.01) ? (0.00002 * Ca2Plus) : 0.01;
    fpn beta_s = 0.015;
    fpn s_inf = alpha_s / (alpha_s + beta_s);
    fpn tau_s = 1 / (alpha_s + beta_s);
    fpn ds_dt = (s_inf - Potassium_s) / tau_s;
    fpn result = DELTA * ds_dt + Potassium_s;
    return result;
}

// Calculate the next state for the dentrite's Ca2Plus
fpn Dend_CalcCa2Plus(fpn Ca2Plus, fpn I_CaH) {
    fpn dCa_dt = -3 * I_CaH - 0.075 * Ca2Plus;
    fpn result = DELTA * dCa_dt + Ca2Plus;
    return result;
}

// Calculate I_c utilizing all of the neighbor neuron's V_dends
fpn Dend_CalcI_cNeighbors(fpn V_dend, size_t network_size, fpn NeighV_dend[], fpn Connectivity[]) {
    fpn const hundred = -1 / 100.0;

    fpn V, f;
    fpn V_acc = 0;
    fpn F_acc = 0;
    for (size_t i = 0; i < network_size; i++) {
        V = V_dend - NeighV_dend[i];
        f = V * expf(V * V * hundred);
        V_acc += V * Connectivity[i];
        F_acc += f * Connectivity[i];
    }
    fpn I_c = 0.8 * F_acc + 0.2 * V_acc;
    return I_c;
}

// Calculate the next state for the dentrite's V_dend and I_CaH
fpn_pair Dend_CalcV_dendI_CaH(Dend dend, fpn I_c, fpn I_app, fpn V_soma) {
    fpn_pair result;
    // Dentritic currents
    // Soma-dendrite interaction current I_sd
    fpn I_sd = (G_INT / (1 - P1)) * (dend.V_dend - V_soma);
    // Inward high-threshold Ca current I_CaH
    fpn I_CaH = G_CAH * dend.Calcium_r * dend.Calcium_r * (dend.V_dend - V_CA);
    // Outward Ca-dependent K current I_K_Ca
    fpn I_K_Ca = G_K_CA * dend.Potassium_s * (dend.V_dend - V_K);
    // Leakage current I_ld
    fpn I_ld = G_LD * (dend.V_dend - V_L);
    // Inward anomalous rectifier I_h
    fpn I_h = G_H * dend.Hcurrent_q * (dend.V_dend - V_H);

    fpn dVd_dt = (-(I_CaH + I_sd + I_ld + I_K_Ca + I_c + I_h) + I_app) / C_M;
    result.first = DELTA * dVd_dt + dend.V_dend;
    result.second = I_CaH;
    return result;
}

/// /// /// SOMA /// /// ///

// Calculate the next state of the soma
Soma Soma_CalcNextState(Soma soma, fpn V_dend, fpn V_axon) {
    fpn_pair temp;
    temp = Soma_CalcCalcium_kCalcium_l(soma.Calcium_k, soma.Calcium_l, soma.V_soma);
    soma.Calcium_k = temp.first;
    soma.Calcium_l = temp.second;
    temp = Soma_CalcSodium_mSodium_h(soma.Sodium_h, soma.V_soma);
    soma.Sodium_m = temp.first;
    soma.Sodium_h = temp.second;
    temp = Soma_CalcPotassium_nPotassium_p(soma.Potassium_n, soma.Potassium_p, soma.V_soma);
    soma.Potassium_n = temp.first;
    soma.Potassium_p = temp.second;
    soma.Potassium_x_s = Soma_CalcPotassium_x_s(soma.Potassium_x_s, soma.V_soma);
    soma.V_soma = Soma_CalcV_soma(soma, V_dend, V_axon);
    return soma;
}

// Calculate the next state for the soma's Calcium_k and Calcium_l
fpn_pair Soma_CalcCalcium_kCalcium_l(fpn Calcium_k, fpn Calcium_l, fpn V_soma) {
    fpn_pair result;
    fpn const four = 1 / 4.2;
    fpn const eight = 1 / 8.5;
    fpn const thirty = 1 / 30.0;
    fpn const seven = 1 / 7.3;

    fpn k_inf = 1 / (1 + expf(-1 * (V_soma + 61) * four));
    fpn l_inf = 1 / (1 + expf((V_soma + 85.5) * eight));
    fpn tau_k = 1;
    fpn tau_l = 20 * expf((V_soma + 160) * thirty) / (1 + expf((V_soma + 84) * seven)) + 35;
    fpn dk_dt = (k_inf - Calcium_k) / tau_k;
    fpn dl_dt = (l_inf - Calcium_l) / tau_l;
    result.first = DELTA * dk_dt + Calcium_k;
    result.second = DELTA * dl_dt + Calcium_l;
    return result;
}

// Calculate the next state for the soma's Sodium_m and Sodium_h
fpn_pair Soma_CalcSodium_mSodium_h(fpn Sodium_h, fpn V_soma) {
    fpn_pair result;
    fpn const five_five = 1 / 5.5;
    fpn const five_eight = -1 / 5.8;
    fpn const thirty_three = 1 / 33.0;

    fpn m_inf = 1 / (1 + expf((-30 - V_soma) * five_five));
    fpn h_inf = 1 / (1 + expf((-70 - V_soma) * five_eight));
    fpn tau_h = 3 * expf((-40 - V_soma) * thirty_three);
    fpn dh_dt = (h_inf - Sodium_h) / tau_h;
    result.first = m_inf;
    result.second = Sodium_h + DELTA * dh_dt;
    return result;
}

// Calculate the next state for the soma's Potassium_n and Potassium_p
fpn_pair Soma_CalcPotassium_nPotassium_p(fpn Potassium_n, fpn Potassium_p, fpn V_soma) {
    fpn_pair result;
    fpn const ten = 0.1;
    fpn const twelve = -1 / 12.0;
    fpn const nine_hundred = 1 / 900.0;

    fpn n_inf = 1 / (1 + expf(( -3 - V_soma) * ten));
    fpn p_inf = 1 / (1 + expf((-51 - V_soma) * twelve));
    fpn tau_n = 5 + 47 * expf(-(-50 - V_soma) * nine_hundred);
    fpn tau_p = tau_n;
    fpn dn_dt = (n_inf - Potassium_n) / tau_n;
    fpn dp_dt = (p_inf - Potassium_p) / tau_p;
    result.first = DELTA * dn_dt + Potassium_n;
    result.second = DELTA * dp_dt + Potassium_p;
    return result;
}

// Calculate the next state for the soma's Potassium_x_s
fpn Soma_CalcPotassium_x_s(fpn Potassium_x_s, fpn V_soma) {
    fpn const ten = 0.1;

    fpn alpha_x_s = 0.13 * (V_soma + 25) / (1 - expf(-(V_soma + 25) * ten));
    fpn beta_x_s = 1.69 * expf(-0.0125 * (V_soma + 35));
    fpn x_inf_s = alpha_x_s / (alpha_x_s + beta_x_s);
    fpn tau_x_s = 1 / (alpha_x_s + beta_x_s);
    fpn dx_dt_s = (x_inf_s - Potassium_x_s) / tau_x_s;
    fpn result = 0.05 * dx_dt_s + Potassium_x_s;
    return result;
}

// Calculate the next state for the soma's V_soma
fpn Soma_CalcV_soma(Soma soma, fpn V_dend, fpn V_axon) {
    // Somatic currents
    // Dendrite-soma interaction current I_ds
    fpn I_ds = (G_INT / P1) * (soma.V_soma - V_dend);
    // Inward low-threshold Ca current I_CaL
    fpn I_CaL = soma.g_CaL * soma.Calcium_k * soma.Calcium_k * soma.Calcium_k * soma.Calcium_l * (soma.V_soma - V_CA);
    // Inward Na current I_Na_s
    fpn I_Na_s = G_NA_S * soma.Sodium_m * soma.Sodium_m * soma.Sodium_m * soma.Sodium_h * (soma.V_soma - V_NA);
    // Leakage current I_ls
    fpn I_ls = G_LS * (soma.V_soma - V_L);
    // Outward delayed potassium current I_Kdr
    fpn I_Kdr_s = G_KDR_S * soma.Potassium_n * soma.Potassium_n * soma.Potassium_n * soma.Potassium_n * (soma.V_soma - V_K);
    // I_K_s
    fpn I_K_s = G_K_S * soma.Potassium_x_s * soma.Potassium_x_s * soma.Potassium_x_s * soma.Potassium_x_s * (soma.V_soma - V_K);
    // Axon-soma interaction current I_as
    fpn I_as = (G_INT / (1 - P2)) * (soma.V_soma - V_axon);

    fpn dVs_dt = -(I_CaL + I_ds + I_as + I_Na_s + I_ls + I_Kdr_s + I_K_s) / C_M;
    fpn result = DELTA * dVs_dt + soma.V_soma;
    return result;
}

/// /// /// AXON /// /// ///

// Calculate the next state of the axon
Axon Axon_CalcNextState(Axon axon, fpn V_soma) {
    fpn_pair temp;
    temp = Axon_CalcSodium_m_aSodium_h_a(axon.Sodium_h_a, axon.V_axon);
    axon.Sodium_m_a = temp.first;
    axon.Sodium_h_a = temp.second;
    axon.Potassium_x_a = Axon_CalcPotassium_x_a(axon.Potassium_x_a, axon.V_axon);
    axon.V_axon = Axon_CalcV_axon(axon, V_soma);
    return axon;
}

// Calculate the next state for the axon's Sodium_m_a and Sodium_h_a
fpn_pair Axon_CalcSodium_m_aSodium_h_a(fpn Sodium_h_a, fpn V_axon) {
    fpn_pair result;
    fpn const five_five = 1 / 5.5;
    fpn const five_eight = -1 / 5.8;
    fpn const thirty_three = 1 / 33.0;

    fpn m_inf_a = 1 / (1 + expf((-30 - V_axon) * five_five));
    fpn h_inf_a = 1 / (1 + expf((-60 - V_axon) * five_eight));
    fpn tau_h_a = 1.5 * expf((-40 - V_axon) * thirty_three);
    fpn dh_dt_a = (h_inf_a - Sodium_h_a) / tau_h_a;
    result.first = m_inf_a;
    result.second = Sodium_h_a + DELTA * dh_dt_a;
    return result;
}

// Calculate the next state for the axon's Potassium_x_a
fpn Axon_CalcPotassium_x_a(fpn Potassium_x_a, fpn V_axon) {
    fpn const ten = 0.1;

    fpn alpha_x_a = 0.13 * (V_axon + 25) / (1 - expf(-(V_axon + 25) * ten));
    fpn beta_x_a = 1.69 * expf(-0.0125 * (V_axon + 35));
    fpn x_inf_a = alpha_x_a / (alpha_x_a + beta_x_a);
    fpn tau_x_a = 1 / (alpha_x_a + beta_x_a);
    fpn dx_dt_a = (x_inf_a - Potassium_x_a) / tau_x_a;
    fpn result = 0.05 * dx_dt_a + Potassium_x_a;
    return result;
}

// Calculate the next state for the axon's V_axon
fpn Axon_CalcV_axon(Axon axon, fpn V_soma) {
    //Axonal currents
    // Sodium
    fpn I_Na_a = G_NA_A * axon.Sodium_m_a * axon.Sodium_m_a * axon.Sodium_m_a * axon.Sodium_h_a * (axon.V_axon - V_NA);
    // Leak
    fpn I_la = G_LA * (axon.V_axon - V_L);
    // Soma-axon interaction current
    fpn I_sa = (G_INT / P2) * (axon.V_axon - V_soma);
    // Potassium (transient)
    fpn I_K_a = G_K_A * axon.Potassium_x_a * axon.Potassium_x_a * axon.Potassium_x_a * axon.Potassium_x_a * (axon.V_axon - V_K);

    fpn dVa_dt = -(I_K_a + I_sa + I_la + I_Na_a) / C_M;
    fpn result = DELTA * dVa_dt + axon.V_axon;
    return result;
}
