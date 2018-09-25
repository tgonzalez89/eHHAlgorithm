#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include "infoli.h"

#ifndef MAX_NETWORK_SIZE
#define MAX_NETWORK_SIZE 10000
#endif

#ifndef MAX_SIMULATION_STEPS
#define MAX_SIMULATION_STEPS 10000
#endif

fpn ConnectivityMatrix[MAX_NETWORK_SIZE][MAX_NETWORK_SIZE];
fpn I_app_array[MAX_SIMULATION_STEPS][MAX_NETWORK_SIZE];
Neuron network[MAX_NETWORK_SIZE];

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

void print_help() {
	printf("USAGE:\n");
	printf("./infoli <network_size> <simulation_steps>\n");
	printf("1 <=   network_size   <= %llu\n", (unsigned long long)MAX_NETWORK_SIZE);
	printf("1 <= simulation_steps <= %llu\n", (unsigned long long)MAX_SIMULATION_STEPS);
}

// Main
int main(int argc, char* argv[]) {
	if (argc != 3) {
		printf("ERROR: Wrong number of arguments.\n");
		print_help();
		return -1;
	}

	int tmp = atoi(argv[1]);
	size_t network_size = tmp >= 0 ? tmp : 0;
	tmp = atoi(argv[2]);
	size_t simulation_steps = tmp >= 0 ? tmp : 0;
	if (network_size < 1 || network_size > MAX_NETWORK_SIZE) {
		printf("ERROR: network_size (%llu) outside of range.\n", (unsigned long long)network_size);
		print_help();
		return -1;
	}
	if (simulation_steps < 1 || simulation_steps > MAX_SIMULATION_STEPS) {
		printf("ERROR: simulation_steps (%llu) outside of range.\n", (unsigned long long)simulation_steps);
		print_help();
		return -1;
	}

    // Initialize network
    createInitialStateNetwork(network, network_size);
    // Initialize connectivity matrix
    /// TO DO: Option to read ConnectivityMatrix from input file
    for (size_t i = 0; i < network_size; i++) {
        for (size_t j = 0; j < network_size; j++) {
            ConnectivityMatrix[i][j] = CONDUCTANCE;
        }
    }
    // Initialize I_app array
    /// TO DO: Option to read I_app_array from input file
    for (size_t i = 0; i < simulation_steps; i++) {
        for (size_t j = 0; j < network_size; j++) {
            I_app_array[i][j] = 0;
        }
    }

    // Simulation
    #ifdef OUTPUT_RESULTS
    FILE* fp;
    fp = fopen("InferiorOlive_Output.txt", "w");
    fprintf(fp, "Simulation Step | Simulation Time (ms) | Neuron # | Input (I_app) | Output (V_axon)\n");
    #endif
    printf("Simulating %llu steps for a network of size %llu.\n", (unsigned long long)simulation_steps, (unsigned long long)network_size);
    struct timeval lbegin;
    gettimeofday(&lbegin, NULL);

    for (size_t i = 0; i < simulation_steps; i++) {
        Network_CalcNextState(network_size, network, I_app_array[i], ConnectivityMatrix);
        #ifdef OUTPUT_RESULTS
        for (size_t j = 0; j < network_size; j++) {
            fprintf(fp, "%15llu | %20.2f | %8llu | %13.3f | %15.8f\n", (unsigned long long)i+1, (float)(i+1) * DELTA, (unsigned long long)j, I_app_array[i][j], network[j].axon.V_axon);
        }
        #endif
        #ifdef STEP_CHECKPOINT
		struct timeval lend;
		gettimeofday(&lend, NULL);
		float ms_elapsed = ((lend.tv_sec - lbegin.tv_sec)*1000000ull + (lend.tv_usec - lbegin.tv_usec))*0.001f;
        printf("Done simulating step #%llu at %.3f ms.\n", (unsigned long long)i+1, ms_elapsed);
        #endif
    }

    struct timeval lend;
    gettimeofday(&lend, NULL);
    float us_elapsed = (lend.tv_sec - lbegin.tv_sec)*1000000ull + (lend.tv_usec - lbegin.tv_usec);
    printf("Time elapsed: %.3f ms, %.3f us per step.\n", us_elapsed/1000.0f, us_elapsed / simulation_steps);
    #ifdef OUTPUT_RESULTS
    fclose(fp);
    #endif

    return 0;
}

