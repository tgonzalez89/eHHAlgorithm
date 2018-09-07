#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <errno.h>
#include <time.h>
#include "infoli.h"

// Parse unsigned long
int str2ulong(const char* str, unsigned long* val) {
    char* endptr;
    errno = 0;
    *val = strtoul(str, &endptr, 0);
    if (endptr == str) {
        return -1; // Empty string
    } else if (*endptr != '\0') {
        return -2; // Garbage
    } else if (*val == ULONG_MAX && errno == ERANGE) {
        return -3; // Overflow
    } else if (errno != 0) {
        return 0; // Unidentified error
    }
    return 1;
}

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

// Help
void help() {
    printf("Usage: infoli_new <network_size> <simulation_steps>\n");
}

// Main
int main(int argc, char* argv[]) {
    // Get arguments: network size and simulation steps
    if (argc < 3) {
        printf("ERROR: Wrong number of arguments. Got %d, expected 2.\n", argc - 1);
        help();
        return 0;
    }
    unsigned long size, steps;
    if(str2ulong(argv[1], &size) <= 0) {
        printf("ERROR: Cannot parse argument <network_size> with value '%s'.\n", argv[1]);
        return -2;
    }
    if(str2ulong(argv[2], &steps) <= 0) {
        printf("ERROR: Cannot parse argument <simulation_steps> with value '%s'.\n", argv[2]);
        return -2;
    }

    // Create and initialize network
    Neuron* network = (Neuron*)malloc(size * sizeof(Neuron));
    createInitialStateNetwork(network, size);

    // Create and initialize connectivity matrix
    fpn** ConnectivityMatrix = (fpn**)malloc(size * sizeof(fpn*));
    for (size_t i = 0; i < size; i++) {
        ConnectivityMatrix[i] = (fpn*)malloc(size * sizeof(fpn));
    }
    for (size_t i = 0; i < size; i++) {
        for (size_t j = 0; j < size; j++) {
            ConnectivityMatrix[i][j] = CONDUCTANCE;
        }
    }

    // Create and initialize I_app array
    fpn** I_app_array = (fpn**)malloc(steps * sizeof(fpn*));
    for (size_t i = 0; i < steps; i++) {
        I_app_array[i] = (fpn*)malloc(size * sizeof(fpn));
    }
    for (size_t i = 0; i < steps; i++) {
        for (size_t j = 0; j < size; j++) {
            I_app_array[i][j] = 0;
        }
    }

    // Create NeighV_dend
    fpn* NeighV_dend = (fpn*)malloc(size * sizeof(fpn));

    // Simulation
    #ifndef NO_OUTPUT
    FILE* fp;
    fp = fopen("InferiorOlive_Output.txt", "w");
    #endif
    printf("Simulating %lu steps for a network of size %lu.\n", steps, size);
    #ifndef NO_OUTPUT
    fprintf(fp, "Simulating %lu steps for a network of size %lu.\n", steps, size);
    printf("Simulation Step | Simulation Time (ms) | Neuron # | Input (I_app) | Output (V_axon)\n");
    fprintf(fp, "Simulation Step | Simulation Time (ms) | Neuron # | Input (I_app) | Output (V_axon)\n");
    #endif
    clock_t begin = clock();
    for (size_t i = 0; i < steps; i++) {
        Network_CalcNextState(network, NeighV_dend, I_app_array[i], ConnectivityMatrix, size);
        #ifndef NO_OUTPUT
        for (size_t j = 0; j < size; j++) {
            printf("%15lu | %20.2f | %8lu | %13.3f | %15.8f\n", i+1, (fpn)(i+1) * DELTA, j, I_app_array[i][j], network[j].axon.V_axon);
            fprintf(fp, "%15lu | %20.2f | %8lu | %13.3f | %15.8f\n", i+1, (fpn)(i+1) * DELTA, j, I_app_array[i][j], network[j].axon.V_axon);
        }
        #endif
    }
    clock_t end = clock();
    float ms_passed = 1000 * (end - begin) / CLOCKS_PER_SEC;
    printf("Time elapsed: %.0f ms, %.3f ms per step.\n", ms_passed, ms_passed / steps);
    #ifndef NO_OUTPUT
    fprintf(fp, "Time elapsed: %.0f ms, %.3f ms per step.\n", ms_passed, ms_passed / steps);
    fclose(fp);
    #endif

    // Free stuff and return
    free(network);
    for (size_t i = 0; i < size; i++) {
        free(ConnectivityMatrix[i]);
    }
    free(ConnectivityMatrix);
    for (size_t i = 0; i < steps; i++) {
        free(I_app_array[i]);
    }
    free(I_app_array);
    free(NeighV_dend);
    return 0;
}

