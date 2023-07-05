#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <math.h>
#include "rand.c"
#include <sys/time.h>
#include <time.h>

//########################################################################################
// Main parameters:
#define mu (0.01) // mutation rate {0.01 (default) 0.05,  0.1, 0.15} 
#define SEQ_NUM (100000) // population size
#define Delta (0.005) // fitness distance -> master fitnesses: 0.005 0.010 0.015 0.020 0.025 0.030 0.035 0.040 0.045 0.050
//#########################################################################################

// constant expressions: 
#define MAX_REPLICATIONS (10000000)    
#define N_REPLICATOR_TYPES (10) // number of competing replicator types
#define L (10) // sequence length
#define PRINT_TIME (2000) 
#define init_pop_sizes (SEQ_NUM/N_REPLICATOR_TYPES)
#define min_binding_loci (2) 
#define dupl_num (21-min_binding_loci) // number of double-stranded chemical species types (=possible binding energies = 2L-2 since BE=0 and BE=1 are not possible)
#define infinity (1.0/0.0)
#define sum_assoc_and_dissoc (1+dupl_num)
#define minimum_HD_between_masters (2)

// static parameters:
char output_filename[100] = "data_delta_0_005__N100000_seed_";
double minimum_master_fitness=0.005;
double full_compl_diss_prob=0.01; // dissociation probability of double strands with maximum binding energy (BEmax=20)
// BE     p_diss
// 0    1.00000000
// 1    0.79432823
// 2    0.63095734
// 3    0.50118723
// 4    0.39810717
// 5    0.31622777
// 6    0.25118864
// 7    0.19952623
// 8    0.15848932
// 9    0.12589254
// 10   0.10000000
// 11   0.07943282
// 12   0.06309573
// 13   0.05011872
// 14   0.03981072
// 15   0.03162278
// 16   0.02511886
// 17   0.01995262
// 18   0.01584893
// 19   0.01258925
// 20   0.01000000
 

// further globals:
int reaction_num; 
double fitness_array[30]; //N*3 vector for replication rates
double *sorted_fitness_array;
double baseline_fitness;
int molecular_count_array[30];
int unique_fitnesses=0;
double sum_simp_dup=0.0;
int sum_simp_dup_counter=0;
double add_rel_sum_of_all_masters=0.0;
int duplex_species[dupl_num]; // vector for concentrations of double-stranded chemical species
double c_d[dupl_num]; // duplex dissociation rates
double master_fitnesses[N_REPLICATOR_TYPES]; // Master (Hd0) fitness values for different types:
double selection_coefficients[]={1.0, 0.2, 0.05}; // master fitnesses are multiplied by these values in case of different mutant classes 

// basis equivalents: 0-> A   2-> U 
// basis equivalents: 1-> G   3-> C

int BE_lookup[4][4] = {
//   A	G  U  C -> 0,1,2,3
    {0, 0, 1, 0}, //A     
    {0, 0, 0, 2}, //G     
    {1, 0, 0, 0}, //U    
    {0, 2, 0, 0}, //C  
};

int sequence_matrix[N_REPLICATOR_TYPES][L]; //N*10 matrix for initial (master) sequences: rows->species, columns->sequences
int del_pos=-333; 
int replication_counter=0; // counts the total number of replication reaction-steps occured
int simplex_state=SEQ_NUM;
int duplex_state=0;
int master_HDs[45];
int simplex_lookup[SEQ_NUM]; // lookup table for simplex indexes 
int master_concentrations[N_REPLICATOR_TYPES]; // concentrations of the original information carrying molecules; 
// so we measure the maintainability of their diversity:
// (the fact that the Gillespie algorithm treats the associated Hamming-distance=0 strands as duplex species  
// and in this way reduces their concentrations at association, necessitates that we follow the master concentrations in an additional array)   
int survival_vector[]={0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
int surv_num=0;
int sum_of_surviving_indeces=0;

typedef struct _COMP
{
	int ind_fitness_index; 
	int master; // if master: it takes the value of the master index, if not master: it is -888  
	int seq[L]; // sequence
	int pair;  //  if x strand is unbinded: pair=-999, if x is binded, pair=y where y is the sequence No. (seq[y]) to which x strand is binded   
	int energy_level; // binding energy if it constitutes a duplex
	int deleted; // 0=no, 1=yes 
}COMP;

COMP comp[SEQ_NUM];

typedef struct _REACTION
{
	int reaction_fitness_index;
	int duplex_bindingE;
	double propensity;
}REACTION;

typedef struct _BE_FEATURES
{
	int energy;
	int binding_loci;
}BE_FEATURES;

BE_FEATURES BE_features;

// functions:
unsigned long int random_seed()
{
	struct timeval tv;
	gettimeofday(&tv,0);
	return (tv.tv_sec + tv.tv_usec);
}

int generate_HD_larger_than_x_seqs(int matr_row){               // function that checks the Hamming-distace of a sequence (i.e. a matrix row) 
	int HammingDist, condition_counter=0, former_row, z;	   // from former sequences (rows) in the sequence matrix
	for(former_row = 0; former_row < matr_row; former_row++) 
	{
		HammingDist=0;
		for(z=0;z<L;z++)
		{
			if(sequence_matrix[matr_row][z]!=sequence_matrix[former_row][z]) HammingDist++;
		}
		if(HammingDist<minimum_HD_between_masters) condition_counter++; 
		if(condition_counter>0) return(1); // return 1 if there is at least one sequence pair 
										   // between which the HD is smaller than "minimum_HD_between_masters" parameter
	}
	return(0); 
}

void HDs_between_masters(void){  //                      
	int Hamming_D, sp_1st_col, pair, z, counter=0;	  
	for(sp_1st_col = 0; sp_1st_col < N_REPLICATOR_TYPES-1; sp_1st_col++) 
	{
		for(pair = sp_1st_col+1; pair < N_REPLICATOR_TYPES; pair++) 
		{
			Hamming_D=0;
			for(z=0; z<L; z++) 
			{
				if(sequence_matrix[sp_1st_col][z]!=sequence_matrix[pair][z]) Hamming_D++;
			}
			master_HDs[counter]=Hamming_D;
			counter++;
		}
	}
}

void rand_fit_and_seq(void)
{	
	int jed, wf;
	int i, j, deco;
	double temp_dec, temp, array1[30];

	master_fitnesses[0]=minimum_master_fitness;
	master_fitnesses[1]=minimum_master_fitness+Delta;
	master_fitnesses[2]=minimum_master_fitness+2*Delta;
	master_fitnesses[3]=minimum_master_fitness+3*Delta;
	master_fitnesses[4]=minimum_master_fitness+4*Delta;
	master_fitnesses[5]=minimum_master_fitness+5*Delta;
	master_fitnesses[6]=minimum_master_fitness+6*Delta;
	master_fitnesses[7]=minimum_master_fitness+7*Delta;
	master_fitnesses[8]=minimum_master_fitness+8*Delta;
	master_fitnesses[9]=minimum_master_fitness+9*Delta;

	baseline_fitness=minimum_master_fitness*selection_coefficients[2];
	
	// all possible fitness values determined here: 
	for(wf = 0; wf < N_REPLICATOR_TYPES; wf++) fitness_array[wf]=master_fitnesses[wf]*selection_coefficients[0];
	for(wf = 0; wf < N_REPLICATOR_TYPES; wf++) fitness_array[10+wf]=master_fitnesses[wf]*selection_coefficients[1];
	for(wf = 0; wf < N_REPLICATOR_TYPES; wf++) fitness_array[20+wf]=master_fitnesses[wf]*selection_coefficients[2];

	// a vector of unique fitness values for replication Gillespie steps made by discarding the possible identical fitness values:
	for(i = 0; i < 30; i++) {
		for(j = i+1; j < 30; j++) {
			if(fitness_array[i] == fitness_array[j]) {
            /* Duplicate element found */
				break;
			}
		}
      /* If j is equal to size, it means we traversed whole
      array and didn't found a duplicate of array[i] */
	if(j == 30) {
		array1[unique_fitnesses++] = fitness_array[i];
		}
	}
	
	//sorting the array1 where only the distinct values are stored
	for ( i = 0; i < unique_fitnesses-1; i++) {
		for ( j = i+1; j < unique_fitnesses; j++) {
			if(array1[i]>array1[j]) {
				temp = array1[i];
				array1[i] = array1[j];
				array1[j] = temp;
			}
		}
	}
   
	// inverse the sorted vector so that it will be in a decreasing order:
	for(deco = 0; deco<unique_fitnesses/2; deco++){
	        temp_dec = array1[deco];
	        array1[deco] = array1[unique_fitnesses-deco-1];
	        array1[unique_fitnesses-deco-1] = temp_dec;
	}	  
   
	sorted_fitness_array = (double*) calloc(unique_fitnesses, sizeof(double));
	for ( i = 0; i < unique_fitnesses; ++i) {
		sorted_fitness_array[i]=array1[i];
	}

	// generate random sequences until they are dissimilar enough (Hamming-distance between all pairs > "minimum_HD_between_masters" parameter)
	for(wf = 0; wf < N_REPLICATOR_TYPES; wf++)
	{
		if(wf==0) for(jed=0;jed<L;jed++) sequence_matrix[wf][jed]=randl(4); // First row of the matrix is certainly not identical to any other row
		else
		{
			for(jed=0;jed<L;jed++) sequence_matrix[wf][jed]=randl(4);                            // In case of further rows, we check if the 
			while(generate_HD_larger_than_x_seqs(wf)==1) for(jed=0;jed<L;jed++) sequence_matrix[wf][jed]=randl(4); // recent row is identical to any former rows 	
		}	
	}	
}

void init_population_vectors(void){
	int s, dd, spp, iii, i, kk, ssr, cc;
	
	// concentrations of duplex Gillespie agents with different binding energies duplex_species[0]: binding energy=0
	for (s = 0; s < dupl_num; s++) duplex_species[s]=0;   
	
	// dissociation rates of duplexes according to the Boltzmann function, maximum binding energy: 2L (E=20) (starts from 2 as BE min. for min binding loci=2)
	for (dd = 0; dd < dupl_num; dd++) c_d[dd]=exp(((double)-(dd+min_binding_loci))/((2.0*(double)L)/(log(full_compl_diss_prob))*(-1.0))); // dd denotes binding energy here

	// initial master concentrations (measured independently from the Gillespie algorithm):
	for(spp = 0; spp < N_REPLICATOR_TYPES; spp++) master_concentrations[spp]=init_pop_sizes; 
		
	// Init properties of the sequences which are initially the same in the entire population:
	for(iii = 0; iii < SEQ_NUM; iii++)
	{
		comp[iii].pair=-999;
		comp[iii].deleted=0;
		comp[iii].energy_level=-11;
		simplex_lookup[iii]=iii;
	}
	
	// Init species-specific properties: 
	for(ssr = 1; ssr < N_REPLICATOR_TYPES+1; ssr++)
	{
		for(i = init_pop_sizes*(ssr-1); i < init_pop_sizes*ssr; i++) 
		{
			comp[i].master=ssr-1;
			comp[i].ind_fitness_index=N_REPLICATOR_TYPES-(ssr-1)-1; // master sp9 -> best fitness -> fitness_index=0
			for(kk=0;kk<L;kk++) comp[i].seq[kk]=sequence_matrix[ssr-1][kk];
		}
	}
	
	for(cc = 0; cc < unique_fitnesses; cc++) molecular_count_array[cc]=0;
	
	// (initial) population densities of simple stranded chemical species:
	for(iii = 0; iii < SEQ_NUM; iii++)
	{	
		for(cc = 0; cc < unique_fitnesses; cc++) 
		{
			if(comp[iii].ind_fitness_index==cc) 
			{
				molecular_count_array[cc]++; // (molecular counts for the Gillespie algorithm)
			}
		}
	}
}

void bindingEnergy2(int a,int b) {
	/*x: pointer to forward strand bases
	 * y: same for reverse bases
	 * x_max: pointer ro first non-vector position of x ptr. WARNING: denomening of x_max causes SEQFAULT!
	 * */
	BE_features.energy = 0;
	BE_features.binding_loci = 0;
	for(int *x= comp[a].seq, *y= comp[b].seq + L-1, *x_max= comp[a].seq + L; x != x_max; ++x, --y) 
	{
		BE_features.energy += BE_lookup[*x][*y];
		if(BE_lookup[*x][*y]>0) BE_features.binding_loci++;
	}
}

void replicate_and_mutate(int iii, int jjj){
	int bazis;//,x,temp;
	
	// complementary replication:
	for(bazis=0; bazis<L; bazis++) 
	{
		if(comp[iii].seq[L - 1 - bazis] % 2 == 0 ) comp[jjj].seq[bazis]=-comp[iii].seq[L -1-bazis]+2; // this "function" makes 2 from 0, and 0 from 2
		else comp[jjj].seq[bazis]=-comp[iii].seq[L -1 -bazis]+4; // this "function" makes 3 from 1, and 1 from
		if(randd() < mu) comp[jjj].seq[bazis] = (comp[jjj].seq[bazis] + randl(3) + 1) & 3;
	} 
}

void define_fitness(int ind)	
{
	int sqnc, sp, uf;
	int HD_vector[N_REPLICATOR_TYPES];
	double fitness_vector_to_choose[N_REPLICATOR_TYPES];
	double max;
  	
  	comp[ind].master=-888;
	for(sp=0;sp<N_REPLICATOR_TYPES;sp++)
	{
		HD_vector[sp]=0;
		for(sqnc=0;sqnc<L;sqnc++)
		{
			if(comp[ind].seq[sqnc]!=sequence_matrix[sp][sqnc]) HD_vector[sp]++; //  the number of loci that is different from the master in the given sequence
			if(HD_vector[sp]>2) break; 
		}	
		if(HD_vector[sp]>2) fitness_vector_to_choose[sp]=baseline_fitness;
		else 
		{	
			fitness_vector_to_choose[sp]=master_fitnesses[sp]*selection_coefficients[HD_vector[sp]];
			if(HD_vector[sp]==0) comp[ind].master=sp;
		}
	}	
	
	// choose the largest fitness value as the fitness of the new replica:
    max = fitness_vector_to_choose[0];
    for (int i = 1; i < N_REPLICATOR_TYPES; i++) 
    {     
        //Compare elements of array with max    
		if(fitness_vector_to_choose[i] > max) 
			max = fitness_vector_to_choose[i];    
	}

	for(uf = 0; uf < unique_fitnesses; uf++) 
	{
		if(max==sorted_fitness_array[uf]) 
		{
			comp[ind].ind_fitness_index=uf;
			break;
		}
	}		
}

void shift(int nn) // shift function: it is used when the element number is changing in the simplex_lookup vector
{
	int s;  
	for(s=nn;s<simplex_state;s++)
		simplex_lookup[s]=simplex_lookup[s+1];    
}

int main(int argc, char** argv)
{
	double sum, ran, prop, tau=0.0, t=0.0, sum_of_master_HDs=0.0, sum_of_all_masters; //tau: step of time;	
	int ik, jk, chbe;
	long random_integer_seed_variable;
	int randdiss, randpair, duplex_index;
	int i_sum, a_one, z_one, i, j, ind_sh;
	int event;			// reaction number selected
	int nulls, di, reps;
	int output_t = 0;
	int i_s = 0, j_s = 0;
	char str_container[100];
	char merged_strings[100];
	random_integer_seed_variable=random_seed();
	sprintf(str_container, "%li", random_integer_seed_variable);	
	// Insert the first string in the new string 
    while (output_filename[i_s] != '\0') { 
        merged_strings[j_s] = output_filename[i_s]; 
        i_s++; 
        j_s++; 
    } 
    // Insert the second string in the new string 
    i_s = 0; 
    while (str_container[i_s] != '\0') { 
        merged_strings[j_s] = str_container[i_s]; 
        i_s++; 
        j_s++; 
    } 
    merged_strings[j_s] = '\0'; 
	seed(random_integer_seed_variable);
	rand_fit_and_seq();
	init_population_vectors();
	HDs_between_masters();
	
	//############### # init_reaction_struct ############################################################################
	reaction_num=sum_assoc_and_dissoc+unique_fitnesses;
	REACTION reaction[reaction_num];
	
	// Init properties of the reactions (if a property remains -11, it means that it doesn't play a role in the reaction):
	for(nulls = 0; nulls < reaction_num; nulls++) reaction[nulls].duplex_bindingE=-11;

	// association of simplexes:
	reaction[0].propensity=(double)simplex_state*((double)simplex_state-1.0)/2.0; 
													   
	// dissociation of duplexes (monomolecular reactions):
	for(di = 0; di < dupl_num; di++)	
	{	
		reaction[1+di].propensity=0.0;
		reaction[1+di].duplex_bindingE=di+min_binding_loci;
	}

	// replication of simplexes (monomolecular reactions):
	for(reps = 0; reps < unique_fitnesses; reps++) 
	{
		reaction[sum_assoc_and_dissoc+reps].reaction_fitness_index = reps;
		reaction[sum_assoc_and_dissoc+reps].propensity = sorted_fitness_array[reps]*molecular_count_array[reps];
	}						
	//###############################################################################################################
	
	FILE *parameters;
	parameters=fopen("master_fit_and_seq.txt","a");
	fprintf(parameters, "seed %li: \n", random_integer_seed_variable);
	for(ik=0; ik<N_REPLICATOR_TYPES; ik++) {
		fprintf(parameters, "sp%d HD0 fitness: %f, sequence: ", ik, master_fitnesses[ik]);
		for(jk=0;jk<L;jk++) {
			fprintf(parameters, "%d ", sequence_matrix[ik][jk]);
			if(jk==L-1){
				fprintf(parameters, "\n");
			}
		}
	}
	fprintf(parameters, "Hamming-distances between master sequences: \n");
	for(chbe=0; chbe<45; chbe++) 
	{
		sum_of_master_HDs+=master_HDs[chbe];
		fprintf(parameters, "%d, ", master_HDs[chbe]);
	}
	fprintf(parameters, "mean Hamming-distances between masters: %f\n", (double)sum_of_master_HDs/(double)45.0);	
	fprintf(parameters, "\n");
	fclose(parameters);
	FILE *file = NULL;
	file = fopen(merged_strings, "w");
	while(replication_counter<MAX_REPLICATIONS){ // The main loop
		
		// select reaction:
		sum=0.0;		
		ran=0.0;
		for(i_sum=0; i_sum<reaction_num; i_sum++) sum+=reaction[i_sum].propensity; // sum all propensities		
		ran = randd()*sum;
		prop = 0.0;
		for(event=0; event<reaction_num; event++){
			prop += reaction[event].propensity;
			if(prop >= ran){
				break;
			}
		} // chosen reaction number: event
		
		// sample tau		
		tau=(1/sum)*log(1/randd()); 
		while(tau==infinity) tau=(1/sum)*log(1/randd()); // this is  because randd() generate (machine) null sometimes and therefore tau would be Inf
		
		// association of simplexes 
		if(event==0)  
		{
			a_one = randl(simplex_state); //comp[simplex_lookup[a_one]]: simplex_lookup[a_one]-th sequence: 
			z_one = randl(simplex_state); 
			while(a_one == z_one) z_one = randl(simplex_state);			
			bindingEnergy2(simplex_lookup[a_one], simplex_lookup[z_one]);	
			if(BE_features.binding_loci>=min_binding_loci)
			{
				comp[simplex_lookup[a_one]].pair=simplex_lookup[z_one]; //associate
				comp[simplex_lookup[z_one]].pair=simplex_lookup[a_one];
				comp[simplex_lookup[a_one]].energy_level=BE_features.energy;
				comp[simplex_lookup[z_one]].energy_level=comp[simplex_lookup[a_one]].energy_level;
				duplex_index=comp[simplex_lookup[a_one]].energy_level-min_binding_loci;
				duplex_species[duplex_index]++; // duplex species for Gillespie increases
				reaction[1+duplex_index].propensity=c_d[duplex_index]*duplex_species[duplex_index];  // update propensities related to duplex concentrations
				duplex_state+=2; // total concentration of duplexes increases 
				
				//##### update replication propensities:
				
				// simplex1:
				molecular_count_array[comp[simplex_lookup[a_one]].ind_fitness_index]--; 
				reaction[sum_assoc_and_dissoc+comp[simplex_lookup[a_one]].ind_fitness_index].propensity 
				 = sorted_fitness_array[comp[simplex_lookup[a_one]].ind_fitness_index]*molecular_count_array[comp[simplex_lookup[a_one]].ind_fitness_index];	
	
				// simplex2:
				molecular_count_array[comp[simplex_lookup[z_one]].ind_fitness_index]--; 
				reaction[sum_assoc_and_dissoc+comp[simplex_lookup[z_one]].ind_fitness_index].propensity 
				 = sorted_fitness_array[comp[simplex_lookup[z_one]].ind_fitness_index]*molecular_count_array[comp[simplex_lookup[z_one]].ind_fitness_index];					
					
				//##### update association propensities:
				if(a_one>z_one) {shift(a_one); simplex_state--; shift(z_one); simplex_state--;}
				else {shift(a_one); simplex_state--; shift(z_one-1); simplex_state--;}
				reaction[0].propensity = (double)simplex_state*((double)simplex_state-1.0)/2.0; // update association propensity
				t += tau; // time update
			}						
		} 
				
		// dissociation of duplexes:		
		else if(event>0 && event<sum_assoc_and_dissoc)   
		{   // indexing of duplex_species from 0 to 18: [0] denotes BE=2 , [18] denotes BE=20
			randdiss = randl(SEQ_NUM);// Choose a random individual from the pool of duplexes with the given energy level:
			while((comp[randdiss].pair==-999) || (comp[randdiss].energy_level!=reaction[event].duplex_bindingE)) randdiss = randl(SEQ_NUM);
			randpair=comp[randdiss].pair;
			comp[randpair].pair=-999;
			comp[randdiss].pair=-999; //dissociate
			duplex_index=reaction[event].duplex_bindingE-min_binding_loci;
			duplex_species[duplex_index]--; // duplex species for Gillespie decreases
			reaction[event].propensity=c_d[duplex_index]*duplex_species[duplex_index];  // update propensities related to duplex concentrations
			duplex_state-=2; // total concentration of duplexes decreases 
	
			//##### update replication propensities:
			
			// simplex1:
			molecular_count_array[comp[randpair].ind_fitness_index]++; 
			reaction[sum_assoc_and_dissoc+comp[randpair].ind_fitness_index].propensity 
			 = sorted_fitness_array[comp[randpair].ind_fitness_index]*molecular_count_array[comp[randpair].ind_fitness_index];	

			// simplex2:
			molecular_count_array[comp[randdiss].ind_fitness_index]++; 
			reaction[sum_assoc_and_dissoc+comp[randdiss].ind_fitness_index].propensity 
			 = sorted_fitness_array[comp[randdiss].ind_fitness_index]*molecular_count_array[comp[randdiss].ind_fitness_index];					
		
			//##### update association propensities:	
			simplex_lookup[simplex_state]=randpair;
			simplex_state++;
			simplex_lookup[simplex_state]=randdiss;
			simplex_state++;
			reaction[0].propensity = (double)simplex_state*((double)simplex_state-1.0)/2.0; // update association propensity
			t += tau; // time update
		} 

		// Replication of simplexes:
		else   
		{ 
			switch(del_pos) 
			{
				case -333: // population size = N
					i = randl(simplex_state); 
					while(comp[simplex_lookup[i]].ind_fitness_index!=reaction[event].reaction_fitness_index) i = randl(simplex_state); // choose individual for replication
					j = randl(SEQ_NUM); // Choose an individual to be overwrited (Moran-process)
					while(simplex_lookup[i] == j) j = randl(SEQ_NUM); // Choose another individual than i (we know already that there is no deleted sequence now: del_pos=-333)
					if(comp[j].pair==-999)  //the overwritten individual is simplex
					{  
						molecular_count_array[comp[j].ind_fitness_index]--; 
						reaction[sum_assoc_and_dissoc+comp[j].ind_fitness_index].propensity 
						 = sorted_fitness_array[comp[j].ind_fitness_index]*molecular_count_array[comp[j].ind_fitness_index];							
						if(comp[j].master>-1) master_concentrations[comp[j].master]--;
						replicate_and_mutate(simplex_lookup[i], j); // replication with overwriting 
						define_fitness(j);
						if(comp[j].master>-1) master_concentrations[comp[j].master]++;
						bindingEnergy2(simplex_lookup[i],j);
						if(BE_features.binding_loci>=min_binding_loci)
						{ 
							comp[j].pair=simplex_lookup[i]; // after replication the strands stay binded to each other 
							comp[simplex_lookup[i]].pair=j;
							comp[simplex_lookup[i]].energy_level=BE_features.energy;
							comp[j].energy_level=BE_features.energy;
							duplex_index=BE_features.energy-min_binding_loci;
							duplex_species[duplex_index]++; // duplex species for Gillespie increases 
							reaction[1+duplex_index].propensity=c_d[duplex_index]*duplex_species[duplex_index];  // update duplex propensities	
							duplex_state+=2; // total concentration of duplexes increases
							molecular_count_array[event-sum_assoc_and_dissoc]--; // replicating i takes duplex form, so simple stranded N Gillespie agents decrease	
							reaction[event].propensity = sorted_fitness_array[event-sum_assoc_and_dissoc]*molecular_count_array[event-sum_assoc_and_dissoc];							
							shift(i); 
							simplex_state--;
							for(ind_sh=0; ind_sh<simplex_state; ind_sh++)
							{
								if(simplex_lookup[ind_sh]==j) // we have to find the lookup table index of j for remove from it
								{
									shift(ind_sh);
									break;
								}  
							}
							simplex_state--;
							reaction[0].propensity = (double)simplex_state*((double)simplex_state-1.0)/2.0; // update association propensity
						}else{	 // overwritten Gillespie agent after replication and mutation stays simple stranded
							molecular_count_array[comp[j].ind_fitness_index]++; 
							reaction[sum_assoc_and_dissoc+comp[j].ind_fitness_index].propensity 
							 = sorted_fitness_array[comp[j].ind_fitness_index]*molecular_count_array[comp[j].ind_fitness_index];							
						}									  
						replication_counter++; // update the number of occured reactions
						t += tau; // time update
						if(output_t++>=PRINT_TIME) 
						{
							sum_of_all_masters=0.0;
							for(ik = 0; ik < N_REPLICATOR_TYPES ; ik++) 
							{
								sum_of_all_masters+=(double)master_concentrations[ik];
								fprintf(file, "%e ", ((double)master_concentrations[ik])/((double)SEQ_NUM)); // masters
							}
							fprintf(file, "%e %f %d %f\n", 
							((double)sum_of_all_masters)/((double)SEQ_NUM), // proportion of all masters relative to the total molecular concentration
							((double)(simplex_state-duplex_state))/((double)(simplex_state+duplex_state)),  //  simplex-duplex ratio
							replication_counter,  //  total number of replication reactions
							t); // reaction times (t+tau)
							fflush(file);
							output_t=0;
							if(replication_counter>=MAX_REPLICATIONS-10000) 
							{
								for(ik = 0; ik < N_REPLICATOR_TYPES ; ik++)
								{
									if(((double)master_concentrations[ik])/((double)SEQ_NUM)>=0.001) survival_vector[ik]++;
								}
								sum_simp_dup += ((double)(simplex_state-duplex_state))/((double)(simplex_state+duplex_state));
								add_rel_sum_of_all_masters+= ((double)sum_of_all_masters)/((double)SEQ_NUM);
								sum_simp_dup_counter++;								
							}
						}
					}
					else //the overwitten individual is in duplex state -> both duplex constituting strands undergo "decay" 	
					{
						duplex_index=comp[j].energy_level-min_binding_loci;
						duplex_species[duplex_index]--; // duplex species for Gillespie decreases 
						reaction[1+duplex_index].propensity=c_d[duplex_index]*duplex_species[duplex_index];  // update duplex propensities									
						del_pos=comp[j].pair;
						if(comp[del_pos].master>-1) master_concentrations[comp[del_pos].master]--;
						comp[del_pos].deleted=1; // "delete" pair of "j" 
						comp[del_pos].pair=-999;
						if(comp[j].master>-1) master_concentrations[comp[j].master]--;
						replicate_and_mutate(simplex_lookup[i], j); // replication 
						define_fitness(j);
						if(comp[j].master>-1) master_concentrations[comp[j].master]++;
						bindingEnergy2(simplex_lookup[i],j);
						if(BE_features.binding_loci>=min_binding_loci)
						{ 
							comp[j].pair=simplex_lookup[i]; // after replication the strands stay binded to each other 
							comp[simplex_lookup[i]].pair=j;
							comp[simplex_lookup[i]].energy_level=BE_features.energy;
							comp[j].energy_level=BE_features.energy;
							duplex_index=BE_features.energy-min_binding_loci;
							duplex_species[duplex_index]++; // duplex species for Gillespie increases 
							reaction[1+duplex_index].propensity=c_d[duplex_index]*duplex_species[duplex_index];  // update duplex propensities	
							molecular_count_array[event-sum_assoc_and_dissoc]--; // replicating i takes duplex form, so simple stranded N Gillespie agents decrease	
							reaction[event].propensity = sorted_fitness_array[event-sum_assoc_and_dissoc]*molecular_count_array[event-sum_assoc_and_dissoc];							
							shift(i); 
							simplex_state--;													 
						}else{   // overwritten Gillespie agent after replication and mutation stays simple stranded
							comp[j].pair=-999;	
							molecular_count_array[comp[j].ind_fitness_index]++; 
							reaction[sum_assoc_and_dissoc+comp[j].ind_fitness_index].propensity 
							 = sorted_fitness_array[comp[j].ind_fitness_index]*molecular_count_array[comp[j].ind_fitness_index];																										
							duplex_state-=2; // total concentration of duplexes decreases	
							simplex_lookup[simplex_state]=j;	
							simplex_state++;			
						}  
						reaction[0].propensity = (double)simplex_state*((double)simplex_state-1.0)/2.0; // update association propensity										   
						replication_counter++; // update the number of occured reactions
						t += tau; // time update
						if(output_t++>=PRINT_TIME) 
						{
							sum_of_all_masters=0.0;
							for(ik = 0; ik < N_REPLICATOR_TYPES ; ik++) 
							{
								sum_of_all_masters+=(double)master_concentrations[ik];
								fprintf(file, "%e ", ((double)master_concentrations[ik])/((double)SEQ_NUM-1)); // masters
							}
							fprintf(file, "%e %f %d %f\n", 
							((double)sum_of_all_masters)/((double)SEQ_NUM-1), // proportion of all masters relative to the total molecular concentration
							((double)(simplex_state-duplex_state))/((double)(simplex_state+duplex_state)),  //  simplex-duplex ratio
							replication_counter,  //  total number of replication reactions
							t); // reaction times (t+tau)
							fflush(file);
							output_t=0;
							if(replication_counter>=MAX_REPLICATIONS-10000) 
							{
								for(ik = 0; ik < N_REPLICATOR_TYPES ; ik++)
								{
									if(((double)master_concentrations[ik])/((double)SEQ_NUM-1)>=0.001) survival_vector[ik]++;
								}
								sum_simp_dup += ((double)(simplex_state-duplex_state))/((double)(simplex_state+duplex_state));
								add_rel_sum_of_all_masters+= ((double)sum_of_all_masters)/((double)SEQ_NUM-1);
								sum_simp_dup_counter++;
							}
						}
					} 
					break;  // end of case 1 in the switch
					
				default: // population size = N-1: replication with no overwriting, ie. "inserting" a new sequence to the array
					i = randl(simplex_state); 
					while(comp[simplex_lookup[i]].ind_fitness_index!=reaction[event].reaction_fitness_index) i = randl(simplex_state); // choose individual for replication
					comp[del_pos].deleted=0;  
					replicate_and_mutate(simplex_lookup[i], del_pos);
					define_fitness(del_pos);
					if(comp[del_pos].master>-1) master_concentrations[comp[del_pos].master]++;
					bindingEnergy2(simplex_lookup[i], del_pos);
					if(BE_features.binding_loci>=min_binding_loci)
					{
						comp[del_pos].pair=simplex_lookup[i]; // after replication the strands stay binded to each other 
						comp[simplex_lookup[i]].pair=del_pos;
						comp[simplex_lookup[i]].energy_level=BE_features.energy;
						comp[del_pos].energy_level=BE_features.energy;
						duplex_index=BE_features.energy-min_binding_loci;
						duplex_species[duplex_index]++; //replicating simple stranded Gillespie agent takes duplex form so duplex population increases
						molecular_count_array[event-sum_assoc_and_dissoc]--; // replicating i takes duplex form, so simple stranded N Gillespie agents decrease	
						reaction[event].propensity = sorted_fitness_array[event-sum_assoc_and_dissoc]*molecular_count_array[event-sum_assoc_and_dissoc];						
						reaction[1+duplex_index].propensity=c_d[duplex_index]*duplex_species[duplex_index];  // update duplex propensities	
						duplex_state+=2; // total concentration of duplexes increases	
						shift(i); 
						simplex_state--;			
					}else{ // the new replica immediately dissociates and stays simple stranded
						molecular_count_array[comp[del_pos].ind_fitness_index]++; 
						reaction[sum_assoc_and_dissoc+comp[del_pos].ind_fitness_index].propensity 
						 = sorted_fitness_array[comp[del_pos].ind_fitness_index]*molecular_count_array[comp[del_pos].ind_fitness_index];													
						simplex_lookup[simplex_state]=del_pos;	
						simplex_state++;						
					}	
					reaction[0].propensity = (double)simplex_state*((double)simplex_state-1.0)/2.0; // update association propensity					  							  
					del_pos=-333;
					replication_counter++; // update the number of occured reactions
					t += tau; // time update
					if(output_t++>=PRINT_TIME) 
					{
						sum_of_all_masters=0.0;
						for(ik = 0; ik < N_REPLICATOR_TYPES ; ik++) 
						{
							sum_of_all_masters+=(double)master_concentrations[ik];
							fprintf(file, "%e ", ((double)master_concentrations[ik])/((double)SEQ_NUM)); // masters
						}
						fprintf(file, "%e %f %d %f\n", 
						((double)sum_of_all_masters)/((double)SEQ_NUM), // proportion of all masters relative to the total molecular concentration
						((double)(simplex_state-duplex_state))/((double)(simplex_state+duplex_state)),  //  simplex-duplex ratio
						replication_counter,  //  total number of replication reactions
						t); // reaction times (t+tau)
						fflush(file);
						output_t=0;
						if(replication_counter>=MAX_REPLICATIONS-10000) 
						{
							for(ik = 0; ik < N_REPLICATOR_TYPES ; ik++)
							{
								if(((double)master_concentrations[ik])/((double)SEQ_NUM)>=0.001) survival_vector[ik]++;
							}
							sum_simp_dup += ((double)(simplex_state-duplex_state))/((double)(simplex_state+duplex_state));
							add_rel_sum_of_all_masters+= ((double)sum_of_all_masters)/((double)SEQ_NUM);
							sum_simp_dup_counter++;
						}
					}
			} // end of the del_pos switch 
		} // end of replication events block	
	} //end time
	free(sorted_fitness_array);
	fprintf(file, "\nsurvived types:\n");
	for(ik = 0; ik < N_REPLICATOR_TYPES ; ik++) 
	{
		if(survival_vector[ik]>1) 
		{
			surv_num++;
			sum_of_surviving_indeces+=ik+1;
			fprintf(file, "%d ", ik+1);
		}
	}
	fprintf(file, " (%f)", ((double)surv_num)/((double)10));
	fprintf(file, "\n");	
	fprintf(file, "\nsurvival rank additive score:\n");
	fprintf(file, "%d", sum_of_surviving_indeces);	
	fprintf(file, "\n");
	fprintf(file, "\nsurvival rank mean:\n");
	fprintf(file, "%f", ((double)sum_of_surviving_indeces)/((double)surv_num));
	fclose(file);
	FILE *surv = NULL;
	surv = fopen("survival_stat__delta_0_005__N100000.txt", "a");
	fprintf(surv, "%.1f, ", ((double)surv_num)/((double)10));
	fclose(surv);
	FILE *sim_dup_file = NULL;
	sim_dup_file = fopen("simplex_duplex_ratio_mean__delta_0_005__N100000.txt", "a");
	fprintf(sim_dup_file, "%f, ", ((double)sum_simp_dup)/((double)sum_simp_dup_counter));
	fclose(sim_dup_file);	
	FILE *sum_master_file = NULL;
	sum_master_file = fopen("mean_sum_of_all_masters__delta_0_005__N100000.txt", "a");
	fprintf(sum_master_file, "%f, ", ((double)add_rel_sum_of_all_masters)/((double)sum_simp_dup_counter));
	fclose(sum_master_file);
	return(0);
}//main end

// compile: gcc fitness_distance_0_005__N100000.c -lm -Wall -Ofast -o delta_0_005__N100000
// run: nohup ./delta_0_005__N100000 &           

