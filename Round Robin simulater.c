/*
***********************************************************************************************
*                RR simulater | Group 1 | salhah Ibrahim & Lujain Alharbi
**********************************************************************************************
*/

//needed libraries 
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <stdbool.h>

//the needed constants for thecalculations in the program
#define NUM_FILES 5
#define NUM_PROCESSES 50
#define ARRIVAL_RATE 0.5
#define MEAN_CPU_TIME 10.0
#define STD_CPU_TIME 1.0
#define MEAN_DISK_TIME 5.0
#define STD_DISK_TIME 0.5
#define TIME_QUANTUM 10
//data structrure for the process state
typedef enum {
    NEW, READY, RUNNING, WAITING
} process_state_t;

//the basic informations for the process 
typedef struct {
    int processID, 
    arrivalTime, 
    CPUTime, diskTime, 
    remainingCPUTime, 
    turnAroundTime, 
    waitingTime;
    process_state_t state;
} process_struct;

//structre for our devices(cpu and disk) flags
typedef struct {
    bool isBusy;
    process_struct *currentProcess;
} Device;

//Arrival queue for the processes
typedef struct {
    process_struct *arr;
    int front, rear, capacity, size;
} ArrivalQueue;

//structure for managing the scheduling on cpu
typedef struct {
    ArrivalQueue *arrivalQueue;
    Device *cpu;
    int timeQuantum;
} CPUScheduler;

//structure for managing the scheduling on disk
typedef struct {
    ArrivalQueue *diskQueue;
    Device *disk;
} DiskScheduler;

//the random numbers generation for input files (arrival time by poisson and cpu and disk time by log-normal)
int generatePoisson(double lambda) {
    double L = exp(-0.5), p = 1.0;
    int k = 0;
    do {
        k++;
        //p is random number between 0 and 1 using the rand() function and the RAND_MAX macro
        p *= ((double)rand() / RAND_MAX);
    } while (p >= L);
    return k - 1;
}
int generateLogNormal(double mean, double stdDev) {
    double standardNormal = 0.0;
    for (int i = 0; i <= 50; i++) {
        standardNormal += (double)rand() / RAND_MAX;
    }
    //round to the nearest integer
    return (int)round(mean + stdDev * standardNormal);
}
//creat the arrival queue entering the sysytem
ArrivalQueue *createArrivalQueue(int capacity) {
    ArrivalQueue *queue = (ArrivalQueue *)malloc(sizeof(ArrivalQueue));
    if (queue == NULL) {
        fprintf(stderr, "Memory allocation error for ArrivalQueue\n");
        exit(EXIT_FAILURE);
    }
    //memory allocation array for arrivel queue 
    queue->arr = (process_struct *)malloc(capacity * sizeof(process_struct));
    if (queue->arr == NULL) {
        fprintf(stderr, "Memory allocation error for ArrivalQueue array\n");
        exit(EXIT_FAILURE);
    }
    queue->capacity = capacity;
    queue->front = queue->size = 0;
    queue->rear = capacity - 1;

    return queue;
}
// function for enetering processes to the arrival queue
void enqueue(ArrivalQueue *queue, process_struct process) {
    if (queue->size == queue->capacity) {
        fprintf(stderr, "Enqueue operation failed. ArrivalQueue is full.\n");
        exit(EXIT_FAILURE);
    }
    queue->rear = (queue->rear + 1) % queue->capacity;
    queue->arr[queue->rear] = process;
    queue->size++;
}
//to give a message in case the proceses failed leaving the queue
process_struct dequeue(ArrivalQueue *queue) {
    if (queue->size == 0) {
        fprintf(stderr, "Dequeue operation failed. ArrivalQueue is empty.\n");
        exit(EXIT_FAILURE);
    }

    process_struct process = queue->arr[queue->front];
    queue->front = (queue->front + 1) % queue->capacity;
    queue->size--;

    return process;
}
//memory allocation for the scheduler and declares a pointer variable named createCPUScheduler for the CPUScheduler
CPUScheduler *createCPUScheduler(ArrivalQueue *arrivalQueue, int timeQuantum) {
    CPUScheduler *cpuScheduler = (CPUScheduler *)malloc(sizeof(CPUScheduler));
    //condition to check an error
    if (cpuScheduler == NULL) {
        fprintf(stderr, "Memory allocation error for CPUScheduler\n");
        exit(EXIT_FAILURE);
    }
    //initialization for the cpuScheduler
    cpuScheduler->arrivalQueue = arrivalQueue;
    cpuScheduler->cpu = (Device *)malloc(sizeof(Device));
    cpuScheduler->cpu->isBusy = false;
    cpuScheduler->cpu->currentProcess = NULL;
    cpuScheduler->timeQuantum = timeQuantum;

    return cpuScheduler;
}
//memory allocation for the scheduler and declares a pointer variable named createDiskScheduler for the DiskScheduler
DiskScheduler *createDiskScheduler() {
    DiskScheduler *diskScheduler = (DiskScheduler *)malloc(sizeof(DiskScheduler));
    //to check if an error occur
    if (diskScheduler == NULL) {
        fprintf(stderr, "Memory allocation error for DiskScheduler\n");
        exit(EXIT_FAILURE);
    }
    //initialization for the diskScheduler
    diskScheduler->diskQueue = createArrivalQueue(NUM_PROCESSES);
    diskScheduler->disk = (Device *)malloc(sizeof(Device));
    diskScheduler->disk->isBusy = false;
    diskScheduler->disk->currentProcess = NULL;

    return diskScheduler;
}
//the round Robin function 
void roundRobinScheduling(CPUScheduler *cpuScheduler, DiskScheduler *diskScheduler, process_struct processes[]) {
    int currentTime = 0;
    int remainingProcesses = NUM_PROCESSES;
    //the initialization for RR
    while (remainingProcesses > 0) {
        for (int i = 0; i < NUM_PROCESSES; i++) {
            if (processes[i].arrivalTime <= currentTime && processes[i].remainingCPUTime > 0) {
                int executionTime = (processes[i].remainingCPUTime > cpuScheduler->timeQuantum)
                                        ? cpuScheduler->timeQuantum
                                        : processes[i].remainingCPUTime;

                processes[i].remainingCPUTime -= executionTime;
                currentTime += executionTime;

                if (processes[i].remainingCPUTime == 0) {
                    processes[i].turnAroundTime = currentTime - processes[i].arrivalTime;
                    processes[i].waitingTime = processes[i].turnAroundTime - processes[i].CPUTime - processes[i].diskTime;
                    processes[i].state = WAITING;

                    enqueue(diskScheduler->diskQueue, processes[i]);
                    remainingProcesses--;
                }
            }
        }
    }
}
//function to creat an outout file where all the results of the scheduler will be stored in 
void printOutputToFile(const char *output, process_struct processes[], int numProcesses, int fileIndex, double averageWaitTime, double averageCompletionTime) {
    FILE *file = fopen(output, "a");
    if (file == NULL) {
        //check if file did not created
        fprintf(stderr, "Error opening file %s for writing.\n", output);
        exit(EXIT_FAILURE);
    }
    //print statement to the file
    fprintf(file, "\nROUND ROBIN SCHEDULING ON CPU AND DISK with TIME_QUANTUM 10 \n");
    fprintf(file, "\n-------------------------------------------------------------------------------------------------\n");
    fprintf(file, "\nScheduling for the input file %d:\n", fileIndex);
    fprintf(file, "ProcessID\tArrivalTime\tBurstTime\tCompletionTime\tWaitingTime\tCPUTime\tDiskTime\n");
    //loop for the printing process of the numbers 
    for (int i = 0; i < numProcesses; i++) {
        fprintf(file, "%d\t\t%d\t\t%d\t\t%d\t\t%d\t\t%d\t\t%d\n",
                processes[i].processID,
                processes[i].arrivalTime,
                processes[i].CPUTime + processes[i].diskTime,
                processes[i].turnAroundTime,
                processes[i].waitingTime,
                processes[i].CPUTime,
                processes[i].diskTime);
    }
    //printing the averages that will be calculated in the bellow function
    fprintf(file, "\nAverage Waiting Time: %.2f\n", averageWaitTime);
    fprintf(file, "Average Completion Time: %.2f\n", averageCompletionTime);
    fclose(file);
}
// function to calculate the averages Waiting Time and Completion Time for each input file reuslts
void calculateAverages(process_struct processes[], int numProcesses, double *averageWaitTime, double *averageCompletionTime) {
    for (int i = 0; i < numProcesses; i++) {
        *averageWaitTime += processes[i].waitingTime;
        *averageCompletionTime += processes[i].turnAroundTime;
    }
    *averageWaitTime /= numProcesses;
    *averageCompletionTime /= numProcesses;
}
//creating the input files that contains the number of processes sequantially and the arrival time generated by poisson and cpu, dik time by log-normal that been calculated above
void generateInputFile(const char *inputs, int fileIndex) {
    srand(time(NULL) + fileIndex);
    //open a file  in writing mode to write the numbers in 5 inout trace files
    FILE *file = fopen(inputs, "w");
    if (file == NULL) {
        fprintf(stderr, "Error opening file %s for writing.\n", inputs);
        exit(EXIT_FAILURE);
    }
    //for loop for writing in the input files
    for (int i = 0; i < NUM_PROCESSES; i++) {
        fprintf(file, "%d %d %d %d\n",
                i + 1,
                generatePoisson(ARRIVAL_RATE),
                generateLogNormal(MEAN_CPU_TIME, STD_CPU_TIME),
                generateLogNormal(MEAN_DISK_TIME, STD_DISK_TIME));
    }

    fclose(file);
}
// functin for generating the input And output files for reading and writing
void generateFiles(const char *output, int fileIndex) {
    double totalWaitingTime = 0;
    double totalCompletionTime = 0;

    process_struct processes[NUM_PROCESSES];
    ArrivalQueue *arrivalQueue = createArrivalQueue(NUM_PROCESSES);
    CPUScheduler *cpuScheduler = createCPUScheduler(arrivalQueue, TIME_QUANTUM);
    DiskScheduler *diskScheduler = createDiskScheduler();

    char inputs[20];
    snprintf(inputs, sizeof(inputs), "input%d.txt", fileIndex);

    generateInputFile(inputs, fileIndex);
    //open the input files on reading mode to read them for the outbut file calculations
    FILE *inputFile= fopen(inputs, "r");
    if (inputFile == NULL) {
        fprintf(stderr, "Error opening file %s for reading.\n", inputs);
        exit(EXIT_FAILURE);
    }

    for (int i = 0; i < NUM_PROCESSES; i++) {
        fscanf(inputFile, "%d %d %d %d",
               &processes[i].processID,
               &processes[i].arrivalTime,
               &processes[i].CPUTime,
               &processes[i].diskTime);

        processes[i].remainingCPUTime = processes[i].CPUTime;
        processes[i].state = NEW;
        enqueue(arrivalQueue, processes[i]);
    }

    fclose(inputFile);

    for (int i = 0; i < NUM_PROCESSES - 1; i++) {
        for (int j = 0; j < NUM_PROCESSES - i - 1; j++) {
            if (arrivalQueue->arr[j].arrivalTime > arrivalQueue->arr[j + 1].arrivalTime) {
                process_struct temp = arrivalQueue->arr[j];
                arrivalQueue->arr[j] = arrivalQueue->arr[j + 1];
                arrivalQueue->arr[j + 1] = temp;
            }
        }
    }
    // the assigned parameters to each function, the basic processing the code
    roundRobinScheduling(cpuScheduler, diskScheduler, arrivalQueue->arr);
    calculateAverages(arrivalQueue->arr, NUM_PROCESSES, &totalWaitingTime, &totalCompletionTime);
    printOutputToFile(output, arrivalQueue->arr, NUM_PROCESSES, fileIndex, totalWaitingTime, totalCompletionTime);

    //this part is to free the dynamically allocated memory to for the arrival queue, cpu scheduler, and disk scheduler for memory management and to prevent leaking
    free(arrivalQueue->arr);
    free(arrivalQueue);
    free(cpuScheduler->cpu);
    free(cpuScheduler);
    free(diskScheduler->diskQueue->arr);
    free(diskScheduler->diskQueue);
    free(diskScheduler->disk);
    free(diskScheduler);
}
// function to calculate the sum of averages been calculated in calculateAverages function for the 5 files every iteration 
void calculateOverallAverages(const char *output, int numFiles) {
    double sumTotalWaitingTime = 0.0;
    double sumTotalCompletionTime = 0.0;

    for (int fileIndex = 1; fileIndex <= numFiles; fileIndex++) {
        //open the output file for reading the averages
        FILE *file = fopen(output, "r");
        //error message if the file not open
        if (file == NULL) {
            fprintf(stderr, "Error opening file %s for reading.\n", output);
            exit(EXIT_FAILURE);
        }
        //initzilize the variables to 0
        double fileAverageWaitingTime = 0.0;
        double fileAverageCompletionTime = 0.0;

        // skip to the relevant section in the file, so it reads only the avreges in the output file
        for (int i = 0; i < (NUM_PROCESSES + 5) * (fileIndex - 1); i++) {
            char line[256];
            fgets(line, sizeof(line), file);
        }

        // read average waiting time and average completion time
        fscanf(file, "\nAverage Waiting Time: %lf", &fileAverageWaitingTime);
        fscanf(file, "\nAverage Completion Time: %lf", &fileAverageCompletionTime);

        // accumulate the times from each file
        sumTotalWaitingTime += fileAverageWaitingTime;
        sumTotalCompletionTime += fileAverageCompletionTime;

        fclose(file);
    }

    // Calculate the overall averages
    double overallAverageWaitingTime = sumTotalWaitingTime;
    double overallAverageCompletionTime = sumTotalCompletionTime;
    // open the output file in append mode
    FILE *file = fopen(output, "a");
    if (file == NULL) {
        // error check if file did not open
        fprintf(stderr, "Error opening file %s for writing.\n", output);
        exit(EXIT_FAILURE);
    }
    // print statements at the end of the output file after all the 5 input files been processed
    fprintf(file, "\n-------------------------------------------------------\n");
    fprintf(file, "\nSum Average Waiting Time for %d files: %.2f\n", numFiles, overallAverageWaitingTime);
    fprintf(file, "Sum Average Completion Time for %d files: %.2f\n", numFiles, overallAverageCompletionTime);
    fprintf(file, "\n--------------------------------------------------------\n");
    fclose(file);
}
//the main function to the outout
int main() {
    const char *output = "output.txt";

    for (int fileIndex = 1; fileIndex <= NUM_FILES; fileIndex++) {
        generateFiles(output, fileIndex);
    }
    calculateOverallAverages(output, NUM_FILES);
    return 0;
}


