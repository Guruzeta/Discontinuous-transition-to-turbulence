#include <stdio.h>
#include <time.h>

#define NSTEPS 200000  // Number of time steps
#define DT 0.001     // Time step size
#define NUM_P_VALUES 100 // Number of p values
#define MIN_P 0.1      // Minimum p value
#define MAX_P 1.0      // Maximum p value
#define GAMMA 0.500       // Diffusion coefficient Gamma
#define qweu 0.6000      // Value of parameter q   

double rho[NUM_P_VALUES]; // Density array

// Function to initialize the density and height arrays
void initialize() {
    for (int i = 0; i < NUM_P_VALUES; i++) {
        rho[i] = 1; // Initial density is 1 forall p i.e. setting fully turbulent initial condition
    }
}


// Function to perform a single time step of the simulation
void step(FILE *file) {
    
    double rho_new[NUM_P_VALUES];

    for (int i = 0; i < NUM_P_VALUES; i++) { //loop to update the rho as a function of p array
        double p = MIN_P + (MAX_P - MIN_P) * i / (NUM_P_VALUES - 1); // Calculate current p value
        double A = 2 * p - 1;
        double B = p;
        rho_new[i] = rho[i] + DT * ((A) * rho[i] - B * rho[i] * rho[i] - qweu*rho[i] * rho[i] * rho[i] +((2-p)*qweu)*rho[i] * rho[i]);
    }

    for (int i = 0; i < NUM_P_VALUES; i++) {
        rho[i] = rho_new[i];
    }
    // rho[0]=1;

    // Output the density values for this time step to the file
    for (int i = 0; i < NUM_P_VALUES; i++) {
        fprintf(file, "%.3f ", rho[i]);
    }
    fprintf(file, "\n");
}

// Function to run the simulation for a range of p values
void run_simulation() {
    char newFilename[100]; // Adjust size as needed
    sprintf(newFilename, "tricritical_DP_qweu_%g_gamma_%g_DT_%g_NT%d_T_%g.txt", qweu, GAMMA,DT,NSTEPS,DT*NSTEPS);
    FILE *output_file = fopen(newFilename, "w");
    if (output_file == NULL) {
        printf("Error: Unable to open output file\n");
        return;
    }
    clock_t start_time = clock(); // Record start time
    initialize();
    for (int t = 0; t < NSTEPS; t++) {
        step(output_file);
    }
    // }

    clock_t end_time = clock(); // Record end time
    double total_time = (double)(end_time - start_time) / CLOCKS_PER_SEC; // Calculate total time in seconds
    fclose(output_file);
    printf("Simulation complete. Total execution time: %.2f seconds.\n", total_time);
}


int main() {
    // Initialize the density and height arrays
    // double alpha = GAMMA*DT/(DX*DX);
    // double beta = u*DT/(DX);

    initialize();

    // Run the simulation and output results to file
    run_simulation();

    printf("Simulation complete. Results written to output.txt.\n");
        /* code */
 }











