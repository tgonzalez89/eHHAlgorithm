#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <sys/time.h>
#include "infoli.h"

#ifndef SIMULATION_STEPS
#define SIMULATION_STEPS 20
#endif

fpn ConnectivityMatrix[NETWORK_SIZE][NETWORK_SIZE];
fpn I_app_array[SIMULATION_STEPS][NETWORK_SIZE];
Neuron network[NETWORK_SIZE];

// Create a network with some initial values
void createInitialStateNetwork(Neuron network[], size_t size) {
    Dend dend;
    dend.V_dend = -60;
    dend.Hcurrent_q = 0.0337836;
    dend.Calcium_r = 0.0112788;
    dend.Potassium_s = 0.0049291;
    dend.I_CaH = 0.5;
    dend.Ca2Plus = 3.7152;

    Soma soma;
    soma.g_CaL = 0.68; // Default arbitrary value but it should be randomized per cell
    soma.V_soma = -60;
    soma.Sodium_m = 1.0127807;
    soma.Sodium_h = 0.3596066;
    soma.Calcium_k = 0.742315;
    soma.Calcium_l = 0.0321349;
    soma.Potassium_n = 0.2369847;
    soma.Potassium_p = 0.2369847;
    soma.Potassium_x_s = 0.1;

    Axon axon;
    axon.V_axon = -60;
    axon.Sodium_m_a = 0.003596066;
    axon.Sodium_h_a = 0.9;
    axon.Potassium_x_a = 0.2369847;

    Neuron neuron;
    neuron.dend = dend;
    neuron.soma = soma;
    neuron.axon = axon;

    for (size_t i = 0; i < size; i++) {
        network[i] = neuron;
    }
}

// Main
int main(int argc, char* argv[]) {
    // Initialize network
    createInitialStateNetwork(network, NETWORK_SIZE);
    // Initialize connectivity matrix
    for (size_t i = 0; i < NETWORK_SIZE; i++) {
        for (size_t j = 0; j < NETWORK_SIZE; j++) {
            ConnectivityMatrix[i][j] = CONDUCTANCE;
        }
    }
    /// TO DO: Read ConnectivityMatrix from input file
    // Initialize I_app array
    for (size_t i = 0; i < SIMULATION_STEPS; i++) {
        for (size_t j = 0; j < NETWORK_SIZE; j++) {
            I_app_array[i][j] = 0;
        }
    }
    /// TO DO: Read I_app_array from input file

    // Simulation
    #ifndef NO_OUTPUT
    FILE* fp;
    fp = fopen("InferiorOlive_Output.txt", "w");
    #endif
    printf("Simulating %u steps for a network of size %u.\n", SIMULATION_STEPS, NETWORK_SIZE);
    #ifndef NO_OUTPUT
    fprintf(fp, "Simulating %u steps for a network of size %u.\n", SIMULATION_STEPS, NETWORK_SIZE);
    printf("Simulation Step | Simulation Time (ms) | Neuron # | Input (I_app) | Output (V_axon)\n");
    fprintf(fp, "Simulation Step | Simulation Time (ms) | Neuron # | Input (I_app) | Output (V_axon)\n");
    #endif
    #ifndef NO_BENCH
    clock_t begin = clock();
    struct timeval lbegin;
    gettimeofday(&lbegin, NULL);
    #endif
    for (size_t i = 0; i < SIMULATION_STEPS; i++) {
        Network_CalcNextState(network, I_app_array[i], ConnectivityMatrix);
        #ifndef NO_OUTPUT
        for (size_t j = 0; j < NETWORK_SIZE; j++) {
            printf("%15u | %20.2f | %8u | %13.3f | %15.8f\n", i+1, (fpn)(i+1) * DELTA, j, I_app_array[i][j], network[j].axon.V_axon);
            fprintf(fp, "%15u | %20.2f | %8u | %13.3f | %15.8f\n", i+1, (fpn)(i+1) * DELTA, j, I_app_array[i][j], network[j].axon.V_axon);
        }
        #endif
    }
    #ifndef NO_BENCH
    clock_t end = clock();
    struct timeval lend;
    gettimeofday(&lend, NULL);
    float ms_passed = 1000 * (end - begin) / CLOCKS_PER_SEC;
    printf("CPU time elapsed:  %.3f ms, %.3f ms per step.\n", ms_passed, ms_passed / SIMULATION_STEPS);
    float lms_passed = ((lend.tv_sec - lbegin.tv_sec)*1000000ul + (lend.tv_usec - lbegin.tv_usec))*0.001;
    printf("Wall time elapsed: %.3f ms, %.3f ms per step.\n", lms_passed, lms_passed / SIMULATION_STEPS);
    #endif
    #ifndef NO_OUTPUT
    #ifndef NO_BENCH
    fprintf(fp, "CPU time elapsed: %.3f ms, %.3f ms per step.\n", ms_passed, ms_passed / SIMULATION_STEPS);
    fprintf(fp, "Wall time elapsed: %.3f ms, %.3f ms per step.\n", lms_passed, lms_passed / SIMULATION_STEPS);
    #endif
    fclose(fp);
    #endif

    return 0;
}
