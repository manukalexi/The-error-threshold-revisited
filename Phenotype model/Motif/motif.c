#include <stdio.h>
#include <string.h>
#include <ViennaRNA/fold.h>
#include <ViennaRNA/utils/basic.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>

#define L (25) //genome length
#define N (100000) //size of the population
#define sigma (10) //masters fitness (reproductive superiority)
#define GEN (10000000) //number of generations
#define REP (50) //number of parallel simulations
#define n (3.0) //fitness landscape geometry factor
#define K (7.0) //the maximum error score, reducing fitness to the baseline level  

int m; //number of master sequences in the population
int seed_number; //seed for the random number generator
double mu; // mutation rate
char *master= "AUCUUACUGGCUACGUCUACUUAGG"; //RNA sequence of the master of length L containing a linear motif and the active region in the specified position
int c = 5; //start base of active region
int d = 12; //end base of active region
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

// filling up the population matrix with master sequences
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

double fit(int i)
{
    int j;
    int fit2 = 0; //the fitness component based on the the number of erroneous bases in the active region
    double f = 1.0;
    int count_bop = 0, count_bcp = 0, count_aop = 0, count_acp = 0;
    char seq[L+1];

    for(j = 0; j < L; j++){
        seq[j] = PM[i][j];
    }

    seq[L] = '\0';

    char *structure = (char *)vrna_alloc(sizeof(char) * (L + 1));
    float mfe = vrna_fold(seq, structure);

    //checking if the folded structure meets the requirements (the motif should be a linear region)
    for (j = 0; j < L; j++) {
        if (j >= c - 2 && j <= d + 2) {
            if (structure[j] != '.') {
            f = 1;
            free(structure);
            return f;
            }
        }
    }

    //checking if the folded structure meets the requirements (the motif does not participate in any secondary structure formation)
    for(j = 0; j < L; j++){
        if(j < c - 2){
            if(structure[j] == '('){
                count_bop+=1;
            }
            if(structure[j] == ')'){
                count_bcp+=1;
            }
        }
        if(j > d + 2){
            if(structure[j] == '('){
                count_aop+=1;
            }

            if(structure[j] == ')'){
                count_acp+=1;
            }
        }   
    }

    //checcking if all loops starting before the motif end before the motif 
    if(count_bop != count_bcp){
        f = 1;
        free(structure);
        return f;
    }

    //checking if all loops starting before the motif end after the motif
    if(count_aop != count_acp){
        f = 1;
        free(structure);
        return f;
        }
    
    //checking for erroneus bases in the active region
    for(j = c; j <= d; j++){
        if (seq[j] != master[j]) {
            fit2 += 1;   
        }
    }

    if(fit2<8){
        f = (sigma-1) * pow((1-pow((((double)(fit2))/((double)K)),n)),(1/n)) + 1;
    }else{
        f=1;
    }

    if(f < 1){
        f = 1;
    }

    free(structure);
    return f;
}


void replicate_with_error(int i)
{
    int j;
    double p; // probability
    char newgenotype[L+1]; // the genotype of the new sequence of length L + \0 null terminator
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

    // the new copy will overwrite a random member of the population (k)
    k = randl(N); 
    while (k == i){
        k = randl(N);
    }
    fitk = fit(k);

    for (j = 0; j < L; j++){
        PM[k][j] = newgenotype[j]; 
    }
    fitj = fit(k);
    fitness[k] = fitj; 

    if(fitk == sigma){ // if the new sequence is a master
        m = (fitj == 10) ? m : m - 1; // updating the number of masters in the population
    }else if(fitj == 10){ // if the new sequence is not a master 
        m = (fitk != 10) ? m + 1 : m; // updating the number of masters in the population
    }
}

int main(int argc, char** argv){
    int i, ex, ism, coex = 0, ext = 0, nmu;
    double ave, aveave;
    char c[50];
    FILE *output;
    for(nmu =  0; nmu < 200; nmu = nmu+1){
        mu = 0.000 + 0.001 * nmu;
        sprintf(c, "L%d_n%.1lf_seq_motif_REP1-%d_%lf.txt", L, n, REP, mu);
        output = fopen(c, "wt");
        coex = 0;
        ext = 0;
        aveave = 0.0;
        for(ism = 0; ism < REP; ism++){
            seed_number = 168750 + ism;
            seed(seed_number);
            init();
            ave = 0.0;
            ex = 0;

            for(int t = 0; t < GEN; t++){
                i = broken_stick();
                replicate_with_error(i);
                
                if(m==0){
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
                fprintf(output, "%lf\n", ave / (0.1 * GEN) / N);
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
