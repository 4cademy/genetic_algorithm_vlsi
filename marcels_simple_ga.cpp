#include "Header.h"
#include <sys/time.h>
#include <cstdio>
#include <cstdlib>
#include <unistd.h>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <omp.h>

double getrnd() {
  return (double)rand() / (double)RAND_MAX;
}

void create_population(double** pop, const unsigned pop_size, const unsigned dim) {
    #pragma omp parallel for
    for (unsigned i=0; i<pop_size; i++) {
        for (unsigned j=0; j<dim; j++) {
            pop[i][j] = getrnd() * 200 - 100;
        }
    }
}

// SELECTION
void compute_fitness(double** pop, double* fitness, double &min_fitness, double &max_fitness, double &convergence, double &best_fitness, const unsigned pop_size, Benchmarks* fp) {
    fitness[0] = fp->compute(pop[0]);
    min_fitness = fitness[0];
    max_fitness = fitness[0];
    // compute new fitness and find min and max fitness
    #pragma omp parallel for reduction(min:min_fitness) reduction(max:max_fitness)
    for (unsigned i=1; i<pop_size; i++) {
        fitness[i] = fp->compute(pop[i]);
        if (fitness[i] < min_fitness) {
            min_fitness = fitness[i];
        } else if (fitness[i] > max_fitness) {
            max_fitness = fitness[i];
        }
    }
    // check if min_fitness is better than current best_fitness
    if (min_fitness < best_fitness) {
        best_fitness = min_fitness;
    }
    convergence = 1-min_fitness/max_fitness;
}

// CROSSOVER

void crossover_uniform(double** mating_list, double** pop, const unsigned pop_size, const unsigned dim) {
    for (unsigned  i=0; i<pop_size; i++) {
        for (unsigned  j=0; j<dim; j++) {
            if (getrnd() < 0.5) {       // choose gene of parent A
                pop[i][j] = mating_list[2*i][j];
            } else {                    // choose gene of parent B
                pop[i][j] = mating_list[2*i+1][j];
            }
        }
    }
}

// MUTATION

void mutation_random_resetting(double** pop, const unsigned pop_size, const unsigned dim, double mutation_rate) {
    for (unsigned i=0; i<pop_size; i++) {
        for (unsigned j=0; j<dim; j++) {
            if (getrnd() < mutation_rate) {
                pop[i][j] = getrnd() * 200 - 100;
            }
        }
    }
}

double genetic_algorithm(Benchmarks*  fp, int maxevals) {
    unsigned tries = 1;
    const unsigned pop_size=10'000;
    const unsigned dim=1000;
    const double convergence_threshold = 0.1;

    auto* best_fitnesses = new double[tries];

    for (unsigned t = 0; t < tries; t++) {
        //create population
        auto** pop = new double*[pop_size];
        for (unsigned i = 0; i < pop_size; i++) {
            pop[i] = new double[dim];
        }

        // create fitness array
        auto* fitness = new double[pop_size];
        double best_fitness;

        // create mating list
        auto** mating_list = new double*[2 * pop_size];
        for (unsigned i = 0; i < 2 * pop_size; i++) {
            mating_list[i] = new double[dim];
        }

        // create fitness metrics
        double min_fitness;
        double max_fitness;
        double convergence;
        double total_fitness;

        auto* offset = new double[pop_size];
        double roulette_random;

        // create and open export file
        std::ofstream fitness_export;
        fitness_export.open ("fitness.csv");

        fp->nextRun();

        // create INITIAL POPULATION
        create_population(pop, pop_size, dim);

        // compute INITIAL FITNESS
        best_fitness = fp->compute(pop[0]);
        compute_fitness(pop, fitness, min_fitness, max_fitness, convergence, best_fitness, pop_size, fp);

        // write initial fitness metrics to export file
        string minf = std::to_string(min_fitness);
        std::replace(minf.begin(), minf.end(), '.', ',');
        string maxf = std::to_string(max_fitness);
        std::replace(maxf.begin(), maxf.end(), '.', ',');
        string conv = std::to_string(convergence);
        std::replace(conv.begin(), conv.end(), '.', ',');

        fitness_export << "0" << ";" << minf << ";" << maxf << ";" << conv << "\n";

        // run actual EVOLUTION
        int evals = 0;
        int convergence_counter = 0;

        printf("Gen: %u\n", 0);
        while (evals < maxevals && convergence_counter < 10) {
            evals++;

            // SELECTION (roulette wheel selection)
            // shift fitness to all positive values and invert it so the minimal value has the highest fitness
            // printf("Shifted fitness:\n");
            for (unsigned i=0; i<pop_size; i++) {
                // fitness[i] = (max_fitness + abs(min_fitness)) - (fitness[i]+abs(min_fitness));
                fitness[i] = max_fitness - fitness[i];
                // printf("%e\n", fitness[i]);
            }

            // calculate total shifted fitness
            total_fitness = 0;
            for (unsigned i=0; i<pop_size; i++) {
                total_fitness += fitness[i];
            }
            offset[0] = fitness[0]/total_fitness;
            // printf("Roulette Offsets:\n");
            for (unsigned i=1; i<pop_size; i++) {
                offset[i] = offset[i-1] + (fitness[i]/total_fitness);
                // printf("%f\n", offset[i]);
            }

            // do roulette selection
            for (unsigned i=0; i < 2*pop_size; i++) {
                roulette_random = getrnd();
                for (unsigned j=0; j<pop_size; j++) {
                    if (roulette_random < offset[j]) {
                        for (unsigned k=0; k < dim; k++) {
                            mating_list[i][k] = pop[j][k];          // ToDo prevent mating with itself
                        }
                        break;
                    }
                }
            }

            // CROSSOVER
            crossover_uniform(mating_list, pop, pop_size, dim);


            // MUTATION
            double mutation_rate = 0.0001;
            // mutation_random_resetting(pop, pop_size, dim, mutation_rate);

            // compute new FITNESS values
            compute_fitness(pop, fitness, min_fitness, max_fitness, convergence, best_fitness, pop_size, fp);

            if (convergence < convergence_threshold) {
                convergence_counter++;
            } else {
                convergence_counter = 0;
            }

            printf("Min Fitness: %e\n", min_fitness);
            printf("Max Fitness: %e\n", max_fitness);
            printf("Convergence: %f\n", convergence);

            // write next fitness metrics to export file
            minf = std::to_string(min_fitness);
            std::replace(minf.begin(), minf.end(), '.', ',');
            maxf = std::to_string(max_fitness);
            std::replace(maxf.begin(), maxf.end(), '.', ',');
            conv = std::to_string(convergence);
            std::replace(conv.begin(), conv.end(), '.', ',');

            fitness_export << std::to_string(evals) << ";" << minf << ";" << maxf << ";" << conv << "\n";


            printf("---------------------\n");
            printf("Gen: %u\n", evals);
        }

        // free allocated memory
        for (unsigned i = 0; i < pop_size; i++) {
            delete[] pop[i];
        }
        delete[] pop;
        for (unsigned i = 0; i < 2 * pop_size; i++) {
            delete[] mating_list[i];
        }
        delete[] mating_list;
        delete[] fitness;
        delete[] offset;

        // close export file
        fitness_export.close();

        // store best fitness value of current try
        best_fitnesses[t] = best_fitness;
    }

    // return global best fitness value over all trys
    double global_best_fitness = best_fitnesses[0];
    for (unsigned i=0; i<tries; i++)
    {
        printf("Best fitness %i: %e\n", i, best_fitnesses[i]);
    }
    for (unsigned i=0; i<tries; i++) {
        if (best_fitnesses[i] < global_best_fitness) {
            global_best_fitness = best_fitnesses[i];
        }
    }
    delete[] best_fitnesses;

    return global_best_fitness;
}

int main(){
  /*  Test the basic benchmark function */
  Benchmarks* fp=NULL;
  unsigned funToRun[] = {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15};
  // unsigned funToRun[] = {1};
  // unsigned funToRun[] = {15};
  unsigned funNum = 1;
  double best_fitness;

  vector<double> runTimeVec;
  struct timeval start, end;
  long seconds, useconds;
  double mtime;
  unsigned maxevals = 1'000;

  srand(time(NULL));

  for (unsigned i=0; i<funNum; i++){
    fp = generateFuncObj(funToRun[i]);
    gettimeofday(&start, NULL);
    best_fitness = genetic_algorithm(fp, maxevals);
    gettimeofday(&end, NULL);

    printf("F %d value = %1.20E\n", fp->getID(), best_fitness);
    seconds  = end.tv_sec  - start.tv_sec;
    useconds = end.tv_usec - start.tv_usec;

    mtime = (((seconds) * 1000 + useconds/1000.0) + 0.5)/1000;

    runTimeVec.push_back(mtime);
    printf ( "F %d, Running Time = %f s\n\n", fp->getID(), mtime);

    delete fp;
  }

  // for (unsigned i=0; i<runTimeVec.size(); i++){
  // 	printf ( "%f\n", runTimeVec[i] );
  // }

  return 0;
}

// create new object of class with default setting
Benchmarks* generateFuncObj(int funcID){
  Benchmarks *fp;
  // run each of specified function in "configure.ini"
  if (funcID==1){
    fp = new F1();
  }else if (funcID==2){
    fp = new F2();
  }else if (funcID==3){
    fp = new F3();
  }else if (funcID==4){
    fp = new F4();
  }else if (funcID==5){
    fp = new F5();
  }else if (funcID==6){
    fp = new F6();
  }else if (funcID==7){
    fp = new F7();
  }else if (funcID==8){
    fp = new F8();
  }else if (funcID==9){
    fp = new F9();
  }else if (funcID==10){
    fp = new F10();
  }else if (funcID==11){
    fp = new F11();
  }else if (funcID==12){
    fp = new F12();
  }else if (funcID==13){
    fp = new F13();
  }else if (funcID==14){
    fp = new F14();
  }else if (funcID==15){
    fp = new F15();
  }else{
    cerr<<"Fail to locate Specified Function Index"<<endl;
    exit(-1);
  }
  return fp;
}
