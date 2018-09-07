#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/time.h>
#include <time.h>
#include <stdbool.h>
#include "infoli.h"

//typedef unsigned long long timestamp_t;

//static timestamp_t get_timestamp ()
//{
//    struct timeval now;
//    gettimeofday (&now, NULL);
//    return  now.tv_usec + (timestamp_t)now.tv_sec * 1000000;
//}

int main(int argc, char *argv[]){

    char *inFileName;
    char *outFileName = "InferiorOlive_Output.txt";
    FILE *pInFile;
    FILE *pOutFile;
    char *iAppBuf;
    const int iAppBufSize =  IAPP_MAX_CHARS*HW_CELLS+1;
    mod_prec iAppArray[IO_NETWORK_SIZE];
    int i, j, k, p, q, n;
    bool ini , new_matrix;
    int simSteps = 0;
    int simTime = 0;
    int inputFromFile = 0;
    int initSteps;
    returnState cellOut;
    cellState IniArray[IO_NETWORK_SIZE];
    cellCompParams cellCompParamsPtr;
    int seedvar;
    char temp[100];//warning: this buffer may overflow
    mod_prec iApp;
    mod_prec Connectivity_Matrix[CONN_MATRIX_SIZE];
    //timestamp_t t0, t1, secs;
    //double secs;


    printf("Inferior Olive Model (%d cell network)\n", IO_NETWORK_SIZE);

    //Process command line arguments
    if(argc == 1){
        inputFromFile = 0;
        printf("Warning: No input file has been specified. A one-pulse input will be used.\n");
    }else if(argc == 2){
        inputFromFile = 1;
        inFileName = argv[1];//comment out for a hardcoded name
        pInFile = fopen(inFileName,"r");
        if(pInFile==NULL){
            printf("Error: Couldn't open %s\n", inFileName);
            exit(EXIT_FAILURE);
        }
    }else{
        printf("Error: Too many arguments.\nUsage: ./InferiorOlive <Iapp_input_file> or ./InferiorOlive\n");
        exit(EXIT_FAILURE);
    }

    //Open output file
    pOutFile = fopen(outFileName,"w");
    if(pOutFile==NULL){
        printf("Error: Couldn't create %s\n", outFileName);
        exit(EXIT_FAILURE);
    }
    sprintf(temp, "#simSteps Time(ms) Input(Iapp) Output(V_dend) Output(V_soma) Output(V_axon)\n");
    fputs(temp, pOutFile);

    //Malloc for iAppBuffer holding iApp arrays, one 2D array (a single line in the file though) at the time
    printf("Malloc'ing memory...\n");
    printf("iAppBuf: %dB\n", iAppBufSize);
    iAppBuf = (char *)malloc(iAppBufSize);
    if(iAppBuf==NULL){
        printf("Error: Couldn't malloc for iAppBuf\n");
        exit(EXIT_FAILURE);
    }


    for(j=0;j<IO_NETWORK_SIZE;j++){


    	IniArray[j] = InitState();

        }

    //Initialize g_CaL
    seedvar = 1;
    for(j=0;j<IO_NETWORK_SIZE;j++){

            srand(seedvar++);   // use this for debugging, now there is difference

            IniArray[j].soma.g_CaL = 0.68;



    }


    //initialize connection Matrix
    for (j=0;j<CONN_MATRIX_SIZE; j++){
    	Connectivity_Matrix[j] = CONDUCTANCE;
    }
    j = 2;

    if(inputFromFile){
        simSteps = 0;
        //Read full lines until end of file. Every iteration (line) is one simulation step.
        while(ReadFileLine(iAppBuf, iAppBufSize, pInFile, iAppArray)){
            //Compute one sim step for all cells
            for(j=0;j<IO_NETWORK_SIZE;j++){

            	 //ComputeNetwork(ini,new_matrix, IniArray, iAppArray,IO_NETWORK_SIZE,TIME_MUX_FACTOR,Connectivity_Matrix,CONN_MATRIX_SIZE,cellOut.axonOut);
                    //Store results
                    sprintf(temp, "%d %.3f %.3f %.8f\n", simSteps, (float)simSteps/20000, iAppArray[j+k], cellOut.axonOut[0]);
                    fputs(temp, pOutFile);


            }
            simSteps++;
        }
    }else{
        simTime = SIMTIME; // in miliseconds
        simSteps = ceil(simTime/DELTA);

        for(i=0;i<simSteps;i++){
          
            if(i>20000-1 && i<20500-1){ iApp = 6;} // start @ 1 because skipping initial values
            else{ iApp = 0;}


           for(j=0;j<IO_NETWORK_SIZE;j++){

                	 iAppArray[j] = iApp;

           }
           sprintf(temp, "%d %.2f %.1f ", i+1, i*0.05,  iAppArray[0]); // start @ 1 because skipping initial values
           fputs(temp, pOutFile);

                    n = 0;

                    //Compute Network...
                    if (i==0){
                    	ini=1;
                    	new_matrix = "true";
                    }
                    else{
                    	ini=0;
                    	new_matrix = "false";
                    }
                    ComputeNetwork(IniArray, iAppArray, TIME_MUX_FACTOR, Connectivity_Matrix, cellOut.axonOut);

                   for(j=0;j<IO_NETWORK_SIZE;j++){
                	   sprintf(temp, "%d: %.8f ",j,cellOut.axonOut[j]);
                	   fputs(temp, pOutFile);
                   }



                   sprintf(temp, "\n");
               	   fputs(temp, pOutFile);
        }
    }

    //t1 = get_timestamp();
    //secs = (t1 - t0);// / 1000000;
    printf("%d ms of brain time in %d simulation steps\n", simTime, simSteps);
    //printf(" %lld usecs real time \n", secs);


    free(iAppBuf);
    fclose (pOutFile);
    if(inputFromFile){ fclose (pInFile);}

    return 0;
}

int ReadFileLine(char *iAppBuf, int iAppBufSize, FILE *pInFile, mod_prec *iAppArray){
    //FIXME: make this function more robust
    char *strNumber;
    int i = 0;
    //Get one line
    if(fgets(iAppBuf, iAppBufSize, pInFile)){
        //Convert the ASCII string of one element to a double precision floating point value
        strNumber = strtok(iAppBuf," ");
        i = 0;
        //printf("Line:\n");
        while ((strNumber != NULL) && (i<IO_NETWORK_SIZE)){
            iAppArray[i] = atof(strNumber);//atof() should change if using integers or fixed point
            //printf ("(%s) %0.2f ", strNumber, iAppArray[i]);
            strNumber = strtok(NULL, " ");
            i++;
        }
        //printf("i: %d\n", i);
        if(i<IO_NETWORK_SIZE){
            //BUG: if only one element is missing but the line ends in a space, the error is not detected
            printf("Error: Input line doesn't have enough elements, only %d\n", i);
            exit(EXIT_FAILURE);
        }
        return 1;//success
    }else{
        if(!feof(pInFile)){
        printf("Error: Reading from input file didn't finish successfully\n");
        exit(EXIT_FAILURE);
        }
        return 0;//end of file
    }
}

cellState InitState(){
    //int j, k;
    cellState initState;
    //Initial dendritic parameters
    initState.dend.V_dend = -60;
    initState.dend.Calcium_r = 0.0112788;// High-threshold calcium
    initState.dend.Potassium_s = 0.0049291;// Calcium-dependent potassium
    initState.dend.Hcurrent_q = 0.0337836;// H current
    initState.dend.Ca2Plus = 3.7152;// Calcium concentration
    initState.dend.I_CaH   = 0.5;// High-threshold calcium current
    //Initial somatic parameters
    initState.soma.g_CaL = 0.68; //default arbitrary value but it should be randomized per cell
    initState.soma.V_soma = -60;
    initState.soma.Sodium_m = 1.0127807;// Sodium (artificial)
    initState.soma.Sodium_h = 0.3596066;
    initState.soma.Potassium_n = 0.2369847;// Potassium (delayed rectifier)
    initState.soma.Potassium_p = 0.2369847;
    initState.soma.Potassium_x_s = 0.1;// Potassium (voltage-dependent)
    initState.soma.Calcium_k = 0.7423159;// Low-threshold calcium
    initState.soma.Calcium_l = 0.0321349;
    // Initial axonal parameters
    initState.axon.V_axon = -60;
    //sisaza: Sodium_m_a doesn't have a state, therefore this assignment doesn'thave any effect
    initState.axon.Sodium_m_a = 0.003596066;// Sodium (thalamocortical)
    initState.axon.Sodium_h_a = 0.9;
    initState.axon.Potassium_x_a = 0.2369847;// Potassium (transient)

    //Copy init sate to all cell states
  //  for(j=0;j<IO_NETWORK_DIM1;j++){
  //     for(k=0;k<IO_NETWORK_DIM2;k++){
  //      	cellCompParamsPtr[j][k].prevCellState = initState;
  //      }
  //  }

    return (initState);
}
