#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <math.h>
#include "rand.c"
#include <sys/time.h>
#include <time.h>

//########################################################################################
// screened parameters:
#define INFLOW_PROPENSITY (0.01) // constant resource inflow 
#define duplex_decay_factor (0.01) // {0.1, 0.5, 1}   scales the duplex decay rates relative to that of the simplexes
#define Delta (0.001) // fitness distances{0.001,  0.005,  0.01} 
	//######  master fitnesses:
	// delta=0.001: 0.005 0.006 0.007 0.008 0.009 0.010 0.011 0.012 0.013 0.014
	// delta=0.005: 0.005 0.010 0.015 0.020 0.025 0.030 0.035 0.040 0.045 0.050
	// delta=0.01:  0.005 0.015 0.025 0.035 0.045 0.055 0.065 0.075 0.085 0.095
//#########################################################################################

// constant parameters:
#define OUTFLOW_PROPENSITY (0.0000019) // constant resource and replicator outflow  
#define simplex_decay (0.0001)  
#define exponent_base (0.8) // 
#define SEQ_NUM (1500) // max population size, if the actual population size excesses this number it gives error
#define mu (0.01) // mutation rate 

// constant expressions: 
#define PRINT_TIME (2000) 
#define MAX_REPLICATIONS (10000000)   
#define inflow_volume (1) // the number of nucleotides that flow into the system at once
#define outflow_volume (1) // the number of nucleotides that flow out from the system at once WARNING: it cannot be larger than one, otherwise it could drop below 0
#define min_binding_loci (2) 
#define N_REPLICATOR_TYPES (10) // number of competing replicator types
#define L (10) // sequence length
#define INIT_POP_SIZE (200)
#define init_type_sizes (INIT_POP_SIZE/N_REPLICATOR_TYPES)
#define dupl_num (21-min_binding_loci) // number of double-stranded chemical species types (=possible binding energies = 2L-2 since BE=0 and BE=1 are not possible)
#define infinity (1.0/0.0)
#define minimum_HD_between_masters (2)

// static parameters:
char output_filename[100] = "data_delta_0_001__dup_dec_fact_0_01_seed_";
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
int	max_decay_react; 
int	first_dupl_diss;
int	sum_assoc_dissoc_decay_inflow;
double fitness_array[30]; //N*3 vector for replication rates
double *sorted_fitness_array;
double baseline_fitness;
int molecular_count_array[30];
int unique_fitnesses=0; // it will be calculated as the actual length of the sorted_fitness_array and molecular_count_array 
double sum_simp_dup=0.0;
int sum_simp_dup_counter=0;
double add_rel_sum_of_all_masters=0.0;
int duplex_species[dupl_num]; // vector for concentrations of double-stranded chemical species
double c_d[dupl_num]; // duplex dissociation rates
double c_dec[dupl_num]; // duplex decay rates
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
int replication_counter=0; // counts the total number of replication reaction-steps occured
int simplex_state=INIT_POP_SIZE;
int duplex_state=0;
int master_HDs[45];
int simplex_lookup[SEQ_NUM]; // lookup table for simplex indexes 
int master_concentrations[N_REPLICATOR_TYPES]; // concentrations of the original information carrying molecules; 
// so we measure the maintainability of their diversity:
// (the fact that the Gillespie algorithm treats the associated Hamming-distance=0 strands as duplex species  
// and in this way reduces their concentrations at association, necessitates that we follow the master concentrations in an additional array)   
int survival_vector[]={0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
int surv_num=0;
int resource_vector[]={100, 100, 100, 100}; // [0]: A(0), [1]: U(2),  [2]: G(1),  [3]: C(3),
int ps=INIT_POP_SIZE;// population size
int sum_pop=0;
int sum_pop_counter=0;

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

int generate_HD_larger_than_x_seqs(int matr_row){               // function that checks what is the Hamming-distace of a sequence (i.e. a matrix row) 
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
	int s, d_d, dd, spp, iii, i, i_ps, d_ps, kk, ssr, cc;
	
	// concentrations of duplex Gillespie agents with different binding energies duplex_species[0]: binding energy=0
	for (s = 0; s < dupl_num; s++) duplex_species[s]=0;   

	// decay rates of duplexes:
	for (d_d = 0; d_d < dupl_num; d_d++) c_dec[d_d]=((double)simplex_decay*(double)duplex_decay_factor)*pow((double)exponent_base, ((double)(d_d+min_binding_loci)));
	
	// dissociation rates of duplexes according to the Boltzmann function, maximum binding energy: 2L (E=20) (starts from 2 as BE min. for min binding loci=2)
	for (dd = 0; dd < dupl_num; dd++) c_d[dd]=exp(((double)-(dd+min_binding_loci))/((2.0*(double)L)/(log(full_compl_diss_prob))*(-1.0))); // dd denotes binding energy here

	// initial master concentrations (measured independently from the Gillespie algorithm):
	for(spp = 0; spp < N_REPLICATOR_TYPES; spp++) master_concentrations[spp]=init_type_sizes; 
		
				
	// Init properties of the sequences which are initially the same in the entire population:
	for(iii = 0; iii < SEQ_NUM; iii++)
	{
		comp[iii].pair=-999;
		comp[iii].energy_level=-11;
	}


    // Init properties of the initial population:
	for(i_ps = 0; i_ps < INIT_POP_SIZE; i_ps++) 
	{
		simplex_lookup[i_ps]=i_ps;
		comp[i_ps].deleted=0;
	}
	
	for(d_ps = 0; d_ps < SEQ_NUM-INIT_POP_SIZE; d_ps++) 
	{
		comp[INIT_POP_SIZE+d_ps].deleted=1;
	}
	
	for(ssr = 1; ssr < N_REPLICATOR_TYPES+1; ssr++)
	{
		for(i = init_type_sizes*(ssr-1); i < init_type_sizes*ssr; i++) 
		{
			comp[i].master=ssr-1; // init without mutant fitnesses (there are only master sequences in the beggining)
			comp[i].ind_fitness_index=N_REPLICATOR_TYPES-(ssr-1)-1; // master sp9 -> best fitness -> fitness_index=0
			for(kk=0;kk<L;kk++) comp[i].seq[kk]=sequence_matrix[ssr-1][kk];
		}
	}
	
	for(cc = 0; cc < unique_fitnesses; cc++) molecular_count_array[cc]=0;
	
	// (initial) population densities of simple stranded chemical species:
	for(iii = 0; iii < INIT_POP_SIZE; iii++)
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
	
	// complementer replication:
	for(bazis=0; bazis<L; bazis++) 
	{
		if(comp[iii].seq[L - 1 - bazis] % 2 == 0 ) comp[jjj].seq[bazis]=-comp[iii].seq[L -1-bazis]+2; // this "function" makes 2 from 0, and 0 from 2
		else comp[jjj].seq[bazis]=-comp[iii].seq[L -1 -bazis]+4; // this "function" makes 3 from 1, and 1 from
		if(randd() < mu) comp[jjj].seq[bazis] = (comp[jjj].seq[bazis] + randl(3) + 1) & 3;
	} 
}

int IsThereEnoughResource(int aa)	// Test whether there is enough nucleotide basis for COMPLEMENTARY replication
{
	int bazis,A=0,U=0,G=0,C=0; // resource_vector[0]: A(0), resource_vector[1]: U(2),  resource_vector[2]: G(1),  resource_vector[3]: C(3),
	for(bazis=0; bazis<L; bazis++)  
	{
		if(comp[aa].seq[bazis]==0) U++; // calculate resource demand in terms of the number of complementary nucleotides
		else if(comp[aa].seq[bazis]==1) C++;
		else if(comp[aa].seq[bazis]==2) A++;
		else if(comp[aa].seq[bazis]==3) G++;
	}
	if(resource_vector[0]>=A && resource_vector[1]>=U && resource_vector[2]>=G && resource_vector[3]>=C) 
	{
		resource_vector[0]-=A; // substract the nucleotides required for replication from their pool
		resource_vector[1]-=U;
		resource_vector[2]-=G;
		resource_vector[3]-=C;
		return(1); // there were enough resources for replication
	}
	else return(0); // there were not enough resources for replication
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
	int ik, jk, chbe, rv;
	long random_integer_seed_variable;
	int ran_dup, del_randpair, randdiss, randpair, duplex_index, ran_sim;
	int i_sum, a_one, z_one, i, j;
	int event;			// reaction number selected
	int nulls, di, reps;
	int output_t = 0;
	int i_s = 0, j_s = 0;
	char str_container[100];
	char merged_strings[100];
	random_integer_seed_variable=random_seed();
	//random_integer_seed_variable=1677156176;
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
	max_decay_react=9+dupl_num; 
	first_dupl_diss=10+dupl_num;
	sum_assoc_dissoc_decay_inflow=first_dupl_diss+dupl_num;
	reaction_num=sum_assoc_dissoc_decay_inflow+unique_fitnesses;	
	REACTION reaction[reaction_num];
		
	// Init properties of the reactions (if a property remains -11, it means that it doesn't play a role in the reaction):
	for(nulls = 0; nulls < reaction_num; nulls++) reaction[nulls].duplex_bindingE=-11;
	
	reaction[0].propensity=(double)simplex_state*((double)simplex_state-1.0)/2.0; // association
	reaction[1].propensity=(double)INFLOW_PROPENSITY; // inflow of A
	reaction[2].propensity=(double)INFLOW_PROPENSITY; // inflow of U
	reaction[3].propensity=(double)INFLOW_PROPENSITY; // inflow of G
	reaction[4].propensity=(double)INFLOW_PROPENSITY; // inflow of C

	reaction[5].propensity=(double)resource_vector[0]*(double)OUTFLOW_PROPENSITY; // outflow of A
	reaction[6].propensity=(double)resource_vector[1]*(double)OUTFLOW_PROPENSITY; // outflow of U
	reaction[7].propensity=(double)resource_vector[2]*(double)OUTFLOW_PROPENSITY; // outflow of G
	reaction[8].propensity=(double)resource_vector[3]*(double)OUTFLOW_PROPENSITY; // outflow of C

	// decay of duplexes (monomolecular reactions):
	for(di = 0; di < dupl_num; di++)	
	{	
		reaction[9+di].propensity=0.0;
		reaction[9+di].duplex_bindingE=di+min_binding_loci;
	}
		
	reaction[max_decay_react].propensity=(double)simplex_state*(simplex_decay+OUTFLOW_PROPENSITY); // simplex decay
													   
	// dissociation of duplexes (monomolecular reactions):
	for(di = 0; di < dupl_num; di++)	
	{	
		reaction[first_dupl_diss+di].propensity=0.0;
		reaction[first_dupl_diss+di].duplex_bindingE=di+min_binding_loci;
	}

	// replication of simplexes (monomolecular reactions):
	for(reps = 0; reps < unique_fitnesses; reps++) 
	{
		reaction[sum_assoc_dissoc_decay_inflow+reps].reaction_fitness_index = reps;
		reaction[sum_assoc_dissoc_decay_inflow+reps].propensity = (double)sorted_fitness_array[reps]*(double)molecular_count_array[reps];
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

	//######### The main loop
	while(replication_counter<MAX_REPLICATIONS){ 
		
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
				duplex_state+=2; // total concentration of duplexes increases 
				
				//##### update replication propensities:
				
				// simplex1:
				molecular_count_array[comp[simplex_lookup[a_one]].ind_fitness_index]--;
				reaction[sum_assoc_dissoc_decay_inflow+comp[simplex_lookup[a_one]].ind_fitness_index].propensity 
				 = (double)sorted_fitness_array[comp[simplex_lookup[a_one]].ind_fitness_index]*(double)molecular_count_array[comp[simplex_lookup[a_one]].ind_fitness_index];	
	
				// simplex2:
				molecular_count_array[comp[simplex_lookup[z_one]].ind_fitness_index]--; 
				reaction[sum_assoc_dissoc_decay_inflow+comp[simplex_lookup[z_one]].ind_fitness_index].propensity 
				 = (double)sorted_fitness_array[comp[simplex_lookup[z_one]].ind_fitness_index]*(double)molecular_count_array[comp[simplex_lookup[z_one]].ind_fitness_index];					
					
				//##### update association propensities:
				if(a_one>z_one) {shift(a_one); simplex_state--; shift(z_one); simplex_state--;}
				else {shift(a_one); simplex_state--; shift(z_one-1); simplex_state--;}
				reaction[0].propensity = (double)simplex_state*((double)simplex_state-1.0)/2.0; // update association propensity
				reaction[max_decay_react].propensity=(double)simplex_state*(simplex_decay+OUTFLOW_PROPENSITY); // update simplex decay
				reaction[9+duplex_index].propensity=(OUTFLOW_PROPENSITY+c_dec[duplex_index])*(double)duplex_species[duplex_index];  // update duplex decay propensities 
				reaction[first_dupl_diss+duplex_index].propensity=c_d[duplex_index]*(double)duplex_species[duplex_index];  // update duplex dissociation propensities
				t += tau; // time update
			}						
		} 

		// constant resource inflow:
		else if(event==1) 
		{
			resource_vector[0]+=inflow_volume; // inflow of A
			reaction[5].propensity=(double)resource_vector[0]*(double)OUTFLOW_PROPENSITY; // update the propensity for outflow of A
		    t += tau; // time update
		}
		else if(event==2) 
		{
			resource_vector[1]+=inflow_volume; // inflow of U
			reaction[6].propensity=(double)resource_vector[1]*(double)OUTFLOW_PROPENSITY; // outflow of U
		    t += tau; // time update
		}	
		else if(event==3) 
		{
			resource_vector[2]+=inflow_volume; // inflow of G
			reaction[7].propensity=(double)resource_vector[2]*(double)OUTFLOW_PROPENSITY; // outflow of G
		    t += tau; // time update
		}	
		else if(event==4) 
		{
			resource_vector[3]+=inflow_volume; // inflow of C
			reaction[8].propensity=(double)resource_vector[3]*(double)OUTFLOW_PROPENSITY; // outflow of C
		    t += tau; // time update
		}	
	
		// resource outflow:
		else if(event==5) 
		{
			resource_vector[0]-=outflow_volume; 
			reaction[5].propensity=(double)resource_vector[0]*(double)OUTFLOW_PROPENSITY; // outflow of A
		    t += tau; // time update
		}
		else if(event==6) 
		{
			resource_vector[1]-=outflow_volume; 
			reaction[6].propensity=(double)resource_vector[1]*(double)OUTFLOW_PROPENSITY; // outflow of U
		    t += tau; // time update
		}
		else if(event==7) 
		{
			resource_vector[2]-=outflow_volume; 
			reaction[7].propensity=(double)resource_vector[2]*(double)OUTFLOW_PROPENSITY; // outflow of G
		    t += tau; // time update
		}
		else if(event==8) 
		{
			resource_vector[3]-=outflow_volume; 
			reaction[8].propensity=(double)resource_vector[3]*(double)OUTFLOW_PROPENSITY; // outflow of C
		    t += tau; // time update
		}

		// duplex decay		
		else if(event>8 && event<max_decay_react)   
		{
			ran_dup = randl(SEQ_NUM);// Choose a random individual from the pool of duplexes with the given energy level:
			while((comp[ran_dup].pair==-999) || (comp[ran_dup].energy_level!=reaction[event].duplex_bindingE)) ran_dup = randl(SEQ_NUM);
			del_randpair=comp[ran_dup].pair;
			comp[ran_dup].pair=-999;
			comp[del_randpair].pair=-999; 
			comp[ran_dup].deleted=1;
			comp[del_randpair].deleted=1; //
			duplex_index=comp[ran_dup].energy_level-min_binding_loci;
			duplex_species[duplex_index]--; // duplex species for Gillespie decreases
			duplex_state-=2; // total concentration of duplexes decreases 
			if(comp[ran_dup].master>-1) master_concentrations[comp[ran_dup].master]--;
			if(comp[del_randpair].master>-1) master_concentrations[comp[del_randpair].master]--;
			ps-=2;
			reaction[9+duplex_index].propensity=(OUTFLOW_PROPENSITY+c_dec[duplex_index])*(double)duplex_species[duplex_index];  // update duplex decay propensities
			reaction[first_dupl_diss+duplex_index].propensity=c_d[duplex_index]*(double)duplex_species[duplex_index];  // update duplex dissociation propensities
			t += tau; // time update
		}

		// simplex decay (only one reaction channel including all simplex molecules)		
		else if(event==max_decay_react)   
		{   //  Choose a random simplex - if it is in the simplex lookup table, cannot be previously deleted 
			ran_sim = randl(simplex_state);
			comp[simplex_lookup[ran_sim]].pair=-999;
			comp[simplex_lookup[ran_sim]].deleted=1;
			if(comp[simplex_lookup[ran_sim]].master>-1) master_concentrations[comp[simplex_lookup[ran_sim]].master]--;
			//##### update replication propensities:
			molecular_count_array[comp[simplex_lookup[ran_sim]].ind_fitness_index]--; 				
			reaction[sum_assoc_dissoc_decay_inflow+comp[simplex_lookup[ran_sim]].ind_fitness_index].propensity 
			 = (double)sorted_fitness_array[comp[simplex_lookup[ran_sim]].ind_fitness_index]*(double)molecular_count_array[comp[simplex_lookup[ran_sim]].ind_fitness_index];					
			shift(ran_sim); 
			simplex_state--;
			ps--;	
			reaction[0].propensity=(double)simplex_state*((double)simplex_state-1.0)/2.0; // update simplex association propensity
			reaction[max_decay_react].propensity=(double)simplex_state*(simplex_decay+OUTFLOW_PROPENSITY); // update simplex decay
			t += tau; // time update
		}
				
		// dissociation of duplexes:		
		else if(event>max_decay_react && event<sum_assoc_dissoc_decay_inflow)  
		{   // indexing of duplex_species from 0 to 18: [0] denotes BE=2 , [18] denotes BE=20
			randdiss = randl(SEQ_NUM);// Choose a random individual from the pool of duplexes with the given energy level:
			while((comp[randdiss].pair==-999) || (comp[randdiss].energy_level!=reaction[event].duplex_bindingE)) randdiss = randl(SEQ_NUM);
			randpair=comp[randdiss].pair;
			comp[randpair].pair=-999;
			comp[randdiss].pair=-999; //dissociate
			duplex_index=reaction[event].duplex_bindingE-min_binding_loci;
			duplex_species[duplex_index]--; // duplex species for Gillespie decreases
			duplex_state-=2; // total concentration of duplexes decreases 
	
			//##### update replication propensities:
			
			// simplex1:
			molecular_count_array[comp[randpair].ind_fitness_index]++; 
			reaction[sum_assoc_dissoc_decay_inflow+comp[randpair].ind_fitness_index].propensity 
			 = (double)sorted_fitness_array[comp[randpair].ind_fitness_index]*(double)molecular_count_array[comp[randpair].ind_fitness_index];	

			// simplex2:
			molecular_count_array[comp[randdiss].ind_fitness_index]++; 
			reaction[sum_assoc_dissoc_decay_inflow+comp[randdiss].ind_fitness_index].propensity 
			 = (double)sorted_fitness_array[comp[randdiss].ind_fitness_index]*(double)molecular_count_array[comp[randdiss].ind_fitness_index];					
		
			//##### update association propensities:	
			simplex_lookup[simplex_state]=randpair;
			simplex_state++;
			simplex_lookup[simplex_state]=randdiss;
			simplex_state++;
			reaction[0].propensity=(double)simplex_state*((double)simplex_state-1.0)/2.0; // update simplex association propensity
			reaction[max_decay_react].propensity=(double)simplex_state*(simplex_decay+OUTFLOW_PROPENSITY); // update simplex decay
			reaction[9+duplex_index].propensity=(OUTFLOW_PROPENSITY+c_dec[duplex_index])*(double)duplex_species[duplex_index];  // update duplex decay propensities
			reaction[first_dupl_diss+duplex_index].propensity=c_d[duplex_index]*(double)duplex_species[duplex_index];  // update duplex dissociation propensities
			t += tau; // time update
		} 

		// Replication of simplexes:
		else   
		{ 			
			i = randl(simplex_state); 
			while(comp[simplex_lookup[i]].ind_fitness_index!=reaction[event].reaction_fitness_index) i = randl(simplex_state); // choose individual for replication
			if(IsThereEnoughResource(simplex_lookup[i])==1)
			{
				j = randl(SEQ_NUM); 
				while(comp[j].deleted==0) j = randl(SEQ_NUM); 
				replicate_and_mutate(simplex_lookup[i], j); // replication with overwriting 
				define_fitness(j);
				comp[j].deleted=0; 
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
					duplex_state+=2; // total concentration of duplexes increases
					reaction[9+duplex_index].propensity=(OUTFLOW_PROPENSITY+c_dec[duplex_index])*(double)duplex_species[duplex_index];  // update duplex decay propensities
					reaction[first_dupl_diss+duplex_index].propensity=c_d[duplex_index]*(double)duplex_species[duplex_index];  // update duplex dissociation propensities
					molecular_count_array[event-sum_assoc_dissoc_decay_inflow]--; // replicating i takes duplex form, so simple stranded N Gillespie agents decrease	
					reaction[event].propensity 
					= (double)sorted_fitness_array[event-sum_assoc_dissoc_decay_inflow]*(double)molecular_count_array[event-sum_assoc_dissoc_decay_inflow];							
					shift(i); 
					simplex_state--;
				}else{	 // newly formed Gillespie agent after replication and mutation stays simple stranded				
					molecular_count_array[comp[j].ind_fitness_index]++;
					reaction[sum_assoc_dissoc_decay_inflow+comp[j].ind_fitness_index].propensity 
					 = (double)sorted_fitness_array[comp[j].ind_fitness_index]*(double)molecular_count_array[comp[j].ind_fitness_index];						
					simplex_lookup[simplex_state]=j;
					simplex_state++;							
				}	
				reaction[0].propensity=(double)simplex_state*((double)simplex_state-1.0)/2.0; // update simplex association propensity
				reaction[max_decay_react].propensity=(double)simplex_state*(simplex_decay+OUTFLOW_PROPENSITY); // update simplex decay

				// update resource outflow propensities since with replication their concentration decreased:
				reaction[5].propensity=(double)resource_vector[0]*(double)OUTFLOW_PROPENSITY; // outflow of A
				reaction[6].propensity=(double)resource_vector[1]*(double)OUTFLOW_PROPENSITY; // outflow of U
				reaction[7].propensity=(double)resource_vector[2]*(double)OUTFLOW_PROPENSITY; // outflow of G
				reaction[8].propensity=(double)resource_vector[3]*(double)OUTFLOW_PROPENSITY; // outflow of C
				ps++;										  
				replication_counter++; // update the number of occured reactions
				t += tau; // time update
				if(output_t++>=PRINT_TIME) 
				{
					sum_of_all_masters=0.0;
					for(ik = 0; ik < N_REPLICATOR_TYPES ; ik++) 
					{
						sum_of_all_masters+=(double)master_concentrations[ik];
						fprintf(file, "%e ", ((double)master_concentrations[ik])/((double)ps)); // masters
					}
					fprintf(file, "%e %f %d %f %d %d %d %d %d\n", 
					((double)sum_of_all_masters)/((double)ps), // proportion of all masters relative to the total molecular concentration
					((double)(simplex_state-duplex_state))/((double)(simplex_state+duplex_state)),  //  simplex-duplex ratio
					replication_counter,  //  total number of replication reactions
					t,
					ps,
					resource_vector[0], resource_vector[1], resource_vector[2], resource_vector[3]); 
					fflush(file);
					output_t=0;
					if(replication_counter>=MAX_REPLICATIONS-10000) 
					{
						for(ik = 0; ik < N_REPLICATOR_TYPES ; ik++)
						{
							if(((double)master_concentrations[ik])/((double)ps)>=0.001) survival_vector[ik]++;
						}
						sum_simp_dup += ((double)(simplex_state-duplex_state))/((double)(simplex_state+duplex_state));
						add_rel_sum_of_all_masters+= ((double)sum_of_all_masters)/((double)ps);
						sum_simp_dup_counter++;	
						sum_pop += ps;						
					}
				}
			}
		} // end of replication events block	
		if (ps < 1) {	fprintf(file, "\nERROR. Population size is too small. Aborting.\n"); fclose(file); exit(0); }
		if (ps > SEQ_NUM-1) {	fprintf(file, "\nERROR. Population size is too large. Aborting.\n"); fclose(file); exit(0); }
		for(rv = 0; rv < 4 ; rv++) 
		{
			if(resource_vector[rv] < 0) {	fprintf(file, "\nERROR. Resource vector %d element smaller than 0. Aborting.\n", rv); fclose(file); exit(0); }
		}		
	} //end time
	free(sorted_fitness_array);
	
	// statistics:
	fprintf(file, "\nsurvived types:\n");
	for(ik = 0; ik < N_REPLICATOR_TYPES ; ik++) 
	{
		if(survival_vector[ik]>1) 
		{
			surv_num++;
			fprintf(file, "%d ", ik+1);
		}
	}
	fprintf(file, "\n%f", ((double)surv_num)/((double)10));
	fclose(file);
	FILE *surv = NULL;
	surv = fopen("survival_stat__delta_0_001__dupl_dec_fact_0_01.txt", "a");
	fprintf(surv, "%.1f, ", ((double)surv_num)/((double)10));
	fclose(surv);
	FILE *sim_dup_file = NULL;
	sim_dup_file = fopen("simplex_duplex_ratio_mean__delta_0_001__dupl_dec_fact_0_01.txt", "a");
	fprintf(sim_dup_file, "%f, ", ((double)sum_simp_dup)/((double)sum_simp_dup_counter));
	fclose(sim_dup_file);	
	FILE *sum_master_file = NULL;
	sum_master_file = fopen("mean_sum_of_all_masters__delta_0_001__dupl_dec_fact_0_01.txt", "a");
	fprintf(sum_master_file, "%f, ", ((double)add_rel_sum_of_all_masters)/((double)sum_simp_dup_counter));
	fclose(sum_master_file);
	FILE *popul_file = NULL;
	popul_file = fopen("pop_size_mean.txt", "a");
	fprintf(popul_file, "%f, ", ((double)sum_pop)/((double)sum_simp_dup_counter));
	fclose(popul_file);
	return(0);
}//main end

// compile: gcc OPEN_SYSTEM___delta_0_001____dupl_dec_fact_0_01___INFLOW_0_01.c -lm -Wall -Ofast -o delta_0_001__dupl_dec_fact_0_01
// run: nohup ./delta_0_001__dupl_dec_fact_0_01 &            

