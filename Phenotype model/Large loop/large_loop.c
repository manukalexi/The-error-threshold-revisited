#include <stdio.h>
#include <string.h>
#include <ViennaRNA/fold.h>
#include <ViennaRNA/utils/basic.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>

#define L (25) //sequence length
#define N (100000) //size of the population
#define sigma (10) //masters fitness (reproductive superiority)
#define GEN (10000000) //number of generations
#define REP (50) //number of parallel simulations
#define n (0.2) //fitness landscape geometry factor
#define K (7.0) //the maximum error score, reducing fitness to the baseline level 

int m; //number of master sequences in the population
int seed_number; //seed for the random number generator
double mu; //mutation rate
char *master= "CCCCAAAAAGUGGCAUGGUAGGGGG\0"; //RNA sequence of the master of length L containing a large loop in its secondary structure and the active region in a previously defined position
int a = ((L-1)/2)-9; //start base of loop (with neck) -> counting from 0!
int b = ((L-1)/2)+8; //end base of loop (with neck) -> counting from 0!
int c = ((L-1)/2)-2; //start base of active region -> counting from 0!
int d = ((L-1)/2)+1; //end base of active region -> counting from 0!
char PM[N][L]; //population matrix
double fitness[N]; //fitness array 

/******* RANDOM NUMBER GENERATOR **********/
#define AA 471
#define BB 1586
#define CC 6988
#define DD 9689
#define M 16383
#define RIMAX 2147483648.0        /* = 2^31 */
#define RandomInteger (++nd, ra[nd & M] = ra[(nd-AA) & M] ^ ra[(nd-BB) & M] ^ ra[(nd-CC) & M] ^ ra[(nd-DD) & M])
void seed(long seed);
static long ra[M+1], nd;

void seed(long seed){
 int  i;

 if(seed<0) { puts("SEED error."); exit(1); }
 ra[0]= (long) fmod(16807.0*(double)seed, 2147483647.0);
 for(i=1; i<=M; i++){
  ra[i] = (long)fmod( 16807.0 * (double) ra[i-1], 2147483647.0);
 }
}

//random number between 0 and num-1
long randl(long num){  
 return(RandomInteger % num);
}

//random number between 0 and 1
double randd(void){ 
 return((double) RandomInteger / RIMAX);
}

/********* END OF RANDOM NUMBER GENERATOR ********/

//filling up the population matrix with master sequences
void init(){
    int i, j;
    for(i=0; i<N; i++){
        for(j=0; j<L; j++){
            PM[i][j] = master[j];
        }
    }
    for(i=0; i<N; i++){
        fitness[i] = sigma;
    }
    m = N;
} 

//selecting an individual for replication with the broken stick method (the chance for replication is based on its fitness proportionate to the total fitness)
int broken_stick(){
    double sumfit=0.0, presfit=0.0;
    double ranfit;
    int i;

    for(i=0; i<N; i++){
        sumfit = sumfit + fitness[i];
    }

    ranfit = randd() * sumfit;

    for(i=0; i<N; i++){
        presfit = presfit + fitness[i];
        if(presfit >= ranfit)
            break;
    }
    return(i);
} 

int type(int i){
    int j;
    int sum = 0;
    int loop_size = 5; //not a specific loop size, only represents junk
    char seq[L+1];

    for(j = 0; j < L; j++){
        seq[j] = PM[i][j];
    }

    seq[L] = '\0'; // null-terminating the sequence to treat it as a string

    char *structure = (char *)vrna_alloc(sizeof(char) * (L + 1));
    float mfe = vrna_fold(seq, structure);

    char *functions[] = {"(..............)","(...............)","(................)","(.................)","(..................)"};

    for(int size = 18; size >= 14; size--){
        int loop_start, loop_end;

        // defining the neccessary loop boundaries (position of the start and end base) based on the size of the loop
        switch(size){
            case 14:
                loop_start = a + 1;
                loop_end = b - 1;
                break;
            case 15:
                loop_start = a;
                loop_end = b - 1;
                break;
            case 16:
                loop_start = a;
                loop_end = b;
                break;
            case 17:
                loop_start = a - 1;
                loop_end = b;
                break;
            case 18:
                loop_start = a - 1;
                loop_end = b + 1;
                break;
            default:
                break;
        }

        // checking if the folded structure matches the expected loop structure of the given size
        for(j = loop_start; j <= loop_end; j++){
            if(structure[j] == functions[size - 14][j -loop_start]){ 
                sum++;
            }
        }

        if(sum == size + 2){ // loop bases (size) + neck brackets (2) 
            loop_size = size; // if the structure matches the loop structure exactly then we treath the structure as a [size] length loop
            free(structure);
            return loop_size;
        }else{
            sum = 0; // reseting the match counter for the next size
        }
    }

    free(structure);
    return loop_size;
}

double fit(int i){
    int j;
    double f = 1.0; //the fitness of the sequence, default is 1
    int loop_size; // identifying it with function 'type'
    double fit1 = 5.0; //the fitness component based on the secondary structure of the sequence -> 5.0 is default, represents junk
    double fit2 = 0.0; //the fitness component based on the the number of erroneous bases in the active region
    char seq[L+1];

    for(j = 0; j < L; j++){
        seq[j] = PM[i][j];
    }

    seq[L] = '\0';

    char *structure = (char *)vrna_alloc(sizeof(char) * (L + 1));
    float mfe = vrna_fold(seq, structure);

    loop_size = type(i);

    //if there is no acceptable loop in the defined positions
    if(loop_size == 5){
        f=1; //baseline fitness -> junk!
        free(structure);
        return f;
    }

    //the difference in loop size results in different error values (master: loop 8 -> no error score)  
    switch(loop_size){
        case 18:
            fit1 = 4;
            break;
        case 17:
            fit1 = 2;
            break;
        case 16:
            fit1 = 0;
            break;
        case 15:
            fit1 = 1;
            break;
        case 14:
            fit1 = 2;
            break;
    }

    // checking for errors in the active region
    for(j = c; j <= d; j++){
        if(seq[j] != master[j]){ 
            fit2 += 1;
        }
    }

    //no correct secondary structure is present
    if(fit1 == 5){ 
        f=1; //baseline fitness -> junk!
        free(structure);
        return f;
    }

    if(fit1 !=5){ // then it means that the structure is something meaningful
        if(fit2 == 4){ //there is no correct base in the active region
        f = 1; //baseline fitness -> junk!
        free(structure);
        return f;
        }else{ //there is correct base in the active region
        f = (sigma-1) * pow((1-pow((((double)(fit1+fit2))/((double)K)),n)),(1/n)) + 1;
        }
    }

    if(f < 1){
        f = 1;
    }

    free(structure);
    return f;
    
}

void replicate_with_error(int i){
    int j;
    double p; //probability
    char newgenotype[L+1]; //the genotype of the new sequence of length L + \0 null terminator
    double fitj, fitk; // fitness of sequence j and k
    int k;

    for(j = 0; j < L; j++){
    p = randd(); 
        if(p < mu){
            int r = randl(3);
            //mutations
            switch(PM[i][j]){
                case 'A':
                    newgenotype[j] = (r == 0) ? 'G' : ((r == 1) ? 'U' : 'C');
                    break;
                case 'G':
                    newgenotype[j] = (r == 0) ? 'U' : ((r == 1) ? 'C' : 'A');
                    break;
                case 'U':
                    newgenotype[j] = (r == 0) ? 'C' : ((r == 1) ? 'A' : 'G');
                    break;
                case 'C':
                    newgenotype[j] = (r == 0) ? 'A' : ((r == 1) ? 'G' : 'U');
                    break;
            }
        }else{
            newgenotype[j] = PM[i][j]; // no mutation with probability 1-mu
        }
    }   

    newgenotype[L] = '\0';

    //the new copy will overwrite a random member of the population (k)
    k = randl(N); 
    while (k == i){
        k = randl(N);
    }
    fitk = fit(k);

    for(j = 0; j < L; j++){
        PM[k][j] = newgenotype[j]; 
    }
    
    fitj = fit(k);
    fitness[k] = fitj; 

    if(fitk == sigma){ // if the new sequence is a master
        m = (fitj == 10) ? m : m - 1; // updating the number of masters in the population
    }else if (fitj == 10){ // if the new sequence is not a master 
        m = (fitk != 10) ? m + 1 : m; // updating the number of masters in the population
    }
}

int main(int argc, char** argv){
    int i, ex, ism, coex = 0, ext = 0, nmu;
    double ave, aveave;
    char c[50];
    FILE *output;
    for(nmu = 0; nmu < 200; nmu = nmu+10){
        mu = 0.000 + 0.001 * nmu;
        sprintf(c, "L%d_n%.1lf_large_loop_REP1-%d_%lf.txt", L, n, REP, mu);
        output = fopen(c, "wt");
        coex = 0;
        ext = 0;
        aveave = 0.0;
        for (ism = 0; ism < REP; ism++) {
            seed_number = 168750 + ism;
            seed(seed_number);
            init();
            ave = 0.0;
            ex = 0;

            for(int t = 0; t < GEN; t++){
                i = broken_stick();
                replicate_with_error(i);

                if(m==0){ //the master is lost from the system
                    ex = 1;
                    break;
                }

                if(t >= 0.9 * GEN){
                    ave += (double)m;
                    if (m == 0) {
                        ex = 1;
                        break;
                    }
                }
            }

            if(ex == 1){
                ext++;
                fprintf(output, "# EXTINCTION\n");
            }else{
                coex++;
                aveave += ave / (0.1 * GEN) / N;
                fprintf(output, "%lf\t%d\n", ave / (0.1 * GEN) / N, seed_number);
            }
            fflush(output);
        }

        if(coex > 0){
            aveave /= (double)coex;
        }else{
            aveave = 0.0;
        }
        fprintf(output, "#%lf\t%lf\t%d\t%d\n", mu, aveave, coex, ext);
        fclose(output);
    }

    return 0;
}
