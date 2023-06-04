#include "Header.h"
#include <sys/time.h>
#include <cstdio>
#include <cstdlib>
#include <unistd.h>
#include <algorithm>

double getrnd() {
  return (double)rand() / (double)RAND_MAX;
}

double genetic_algorithm(Benchmarks*  fp, int maxevals) {
    const unsigned pop_size=100;
    const unsigned dim=1000;
    double pop[pop_size][dim];
    double fitness[pop_size], best_fitness, best_individual[dim];
    // variables for selection
    double mating_list[2*pop_size][dim];
    double min_fitness;
    double max_fitness;
    double total_fitness;
    double offset[pop_size];
    double total_offset;
    double roulette_random;

    fp->nextRun();

    // create INITIAL POPULATION
    for (unsigned i=0; i<pop_size; i++) {
        for (unsigned j = 0; j < dim; j++) {
          pop[i][j] = getrnd() * 200 - 100;
        }
    }

    // compute INITIAL FITNESS
    fitness[0] = fp->compute(pop[0]);
    best_fitness = fitness[0];
    // std::copy(best_individual, best_individual+dim, pop[0]);
    for (unsigned i=1; i<pop_size; i++) {
        fitness[i] = fp->compute(pop[i]);
        if (fitness[i] < best_fitness) {
          best_fitness = fitness[i];
          // std::copy(best_individual, best_individual+dim, pop[i]);
        }
    }

    // run actual EVOLUTION
    for (unsigned evals = 1; evals < maxevals; evals++) {
        printf("Evals: %u\n", evals);

        // SELECTION (roulette wheel selection)
        // find minimal (best) fitness and maximal (worst) fitness
        min_fitness = fitness[0];
        max_fitness = fitness[0];
        for (unsigned i=0; i<pop_size; i++) {
            if (fitness[i] < min_fitness) {
                min_fitness = fitness[i];
            } else if (fitness[i] > min_fitness) {
                max_fitness = fitness[i];
            }
        }
        // shift fitness to all positive values and invert it so the minimal value has the highest fitness
        for (unsigned i=0; i<pop_size; i++) {
            // fitness[i] = (max_fitness + abs(min_fitness)) - (fitness[i]+abs(min_fitness));
            fitness[i] = max_fitness - fitness[i];
        }

        // calculate total fitness
        total_fitness = 0;
        for (unsigned i=0; i<pop_size; i++) {
          total_fitness += fitness[i];
        }
        offset[0] = fitness[0]/total_fitness;
        for (unsigned i=1; i<pop_size; i++) {
            offset[i] = offset[i-1] + (fitness[i]/total_fitness);
        }

        // do roulette selection
        for (unsigned i=0; i < 2*pop_size; i++) {
            roulette_random = getrnd();
            for (unsigned j=0; j<pop_size; j++) {
                if (roulette_random < offset[j]) {
                    std::copy(mating_list[i],mating_list[i]+dim,pop[j]);
                    break;
                }
            }
        }

        // CROSSOVER (uniform crossover)
        for (unsigned  i=1; i<pop_size; i++) {
            for (unsigned  j=0; j<dim; j++) {
                if (getrnd() < 0.5) {       // choose gene of parent A
                    pop[i][j] = mating_list[2*i][j];
                } else {                    // choose gene of parent B
                    pop[i][j] = mating_list[2*i+1][j];
                }
            }
        }

        // MUTATION
        //
        //

        // compute FITNESS of new generation
        for (unsigned i=0; i<pop_size; i++) {
            fitness[i] = fp->compute(pop[i]);
            if (fitness[i] < best_fitness) {
                best_fitness = fitness[i];
                // std::copy(best_individual, best_individual+dim, pop[i]);
            }
        }
        printf("%e\n", best_fitness);
    }

    return best_fitness;
}

int main(){
  /*  Test the basic benchmark function */
  Benchmarks* fp=NULL;
  unsigned funToRun[] = {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15};
  // unsigned funToRun[] = {1};
  // unsigned funToRun[] = {15};
  unsigned funNum = 3;
  double best_fitness;

  vector<double> runTimeVec;
  struct timeval start, end;
  long seconds, useconds;
  double mtime;
  unsigned maxevals = 1000;

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
