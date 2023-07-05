#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <math.h>
#include "rand.c"
#include <sys/time.h>
#include <time.h>

//########################################################################################
// Main parameters:
#define mu (0.01) // mutation rate 
#define SEQ_NUM (10000) // population size
#define Delta (0.01) // fitness distance (now it is irrelevant, since all master sequences have the same fitness)

// constant expressions: 
#define MAX_REPLICATIONS (2000000)   
#define N_REPLICATOR_TYPES (10) // number of competing replicator types
#define L (10) // sequence length
#define PRINT_TIME (2000) 
#define init_pop_sizes (SEQ_NUM/N_REPLICATOR_TYPES)
#define min_binding_loci (2) 
#define dupl_num (21-min_binding_loci) // number of double-stranded chemical species types (=possible binding energies = 2L-2 since BE=0 and BE=1 are not possible)
#define infinity (1.0/0.0)
#define sum_assoc_and_dissoc (1+dupl_num)
#define minimum_HD_between_masters (4)
#define unique_fitnesses (4)

// static parameters:
char output_filename[100] = "data_GC_gradient__seed_";
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
double fitness_array[unique_fitnesses]; //constant fitnesses for HD=0, HD=1, HD=2, HD=3 (and larger-> baseline) 
double molecular_count_array[unique_fitnesses]; // for Gillespie, variable molecular concentration for HD=0, HD=1, HD=2, HD=3 (and larger-> baseline) 
double baseline_fitness;

double sum_simp_dup=0.0;
int sum_simp_dup_counter=0;
double add_rel_sum_of_all_masters=0.0;
int duplex_species[dupl_num]; // vector for concentrations of double-stranded chemical species
double c_d[dupl_num]; // duplex dissociation rates
double master_fitnesses; // Master (Hd0) fitness equal values for different types:
double selection_coefficients[]={1.0, 0.2, 0.05}; // master fitnesses are multiplied by these values in case of different mutant classes 

int BE_lookup[4][4] = {
//   A	G  U  C -> 0,1,2,3
    {0, 0, 1, 0}, //A     
    {0, 0, 0, 2}, //G     
    {1, 0, 0, 0}, //U    
    {0, 2, 0, 0}, //C  
};

int sequence_matrix[N_REPLICATOR_TYPES][L]; //N*10 matrix for initial (master) sequences: rows->species, columns->sequences
int complementary_sequence_matrix[N_REPLICATOR_TYPES][L]; //N*10 matrix for complementary master sequences: rows->species, columns->sequences

int del_pos=-333; 
int replication_counter=0; // counts the total number of replication reaction-steps occured
int simplex_state=SEQ_NUM;
int duplex_state=0;
int master_HDs[55];
int master_BEs[55];
int master_complement_HDs[100];
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
	int master; // if master: it takes the value of the master index, if not master: it is -888  -> this is sequence specific
	int ind_H_Dist; 
	int seq[L]; // sequence
	int pair;  //  if x strand is unbinded: pair=-999, if x is binded, pair=y where y is the sequence No. (seq[y]) to which x strand is binded   
	int energy_level; // binding energy if it constitutes a duplex
	int deleted; // 0=no, 1=yes 
}COMP;

COMP comp[SEQ_NUM];

typedef struct _REACTION
{
	int reaction_Hd;
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

void shuffle(int *array, int n) {
    int i, j, tmp;
 
    for (i = n - 1; i > 0; i--) {
        j = randl(i + 1);
        tmp = array[j];
        array[j] = array[i];
        array[i] = tmp;
   }
}

int generate_HD_larger_3_seqs(int matr_row){               // function that checks what is the Hamming-distace of a sequence (i.e. a matrix row) 
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

int generate_master_complementary_HD_larger_3_seqs(void){               // function that checks what is the Hamming-distace of a sequence (i.e. a matrix row) 
	int condition_counter=0, z;	   // from former sequences (rows) in the sequence matrix
	for(z=0;z<100;z++)
	{
		if(master_complement_HDs[z]<4) condition_counter++;
	}
	if(condition_counter>0) return(1); // return 1 if there is at least one sequence pair 
									   // between which the HD is smaller than "minimum_HD_between_masters" parameter
	return(0); 
}

void HDs_between_masters(void){  //                      
	int Hamming_D, sp_1st_col, pair, z, counter=0;	  
	for(sp_1st_col = 0; sp_1st_col < N_REPLICATOR_TYPES; sp_1st_col++) 
	{
		for(pair = sp_1st_col; pair < N_REPLICATOR_TYPES; pair++) 
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

void HDs_between_masters_and_complementary_strands(void){  //                      
	int Hamming_D, sp_1st_col, pair, z, counter=0;	  
	for(sp_1st_col = 0; sp_1st_col < N_REPLICATOR_TYPES; sp_1st_col++) 
	{
		for(pair = 0; pair < N_REPLICATOR_TYPES; pair++) 
		{
			Hamming_D=0;
			for(z=0; z<L; z++) 
			{
				if(sequence_matrix[sp_1st_col][z]!=complementary_sequence_matrix[pair][z]) Hamming_D++;
			}
			master_complement_HDs[counter]=Hamming_D;
			counter++;
		}
	}
}

void BEs_between_masters(void){  //                      
	int Binding_E, sp_1st_col, pair, counter=0;	  
	for(sp_1st_col = 0; sp_1st_col < N_REPLICATOR_TYPES; sp_1st_col++) 
	{
		for(pair = sp_1st_col; pair < N_REPLICATOR_TYPES; pair++) 
		{
			Binding_E=0;
			for(int *x= sequence_matrix[sp_1st_col], *y= sequence_matrix[pair] + L-1, *x_max= sequence_matrix[sp_1st_col] + L; x != x_max; ++x, --y) Binding_E += BE_lookup[*x][*y];
			master_BEs[counter]=Binding_E;
			counter++;
		}
	}
}

void generate_complementary_master_sequences(void){
	int wf, bazis;

	for(wf = 0; wf < N_REPLICATOR_TYPES; wf++)
	{
		for(bazis=0;bazis<L;bazis++) 
		{
			if(sequence_matrix[wf][L - 1 - bazis] % 2 == 0 ) complementary_sequence_matrix[wf][bazis]=-sequence_matrix[wf][L - 1 - bazis]+2; // this "function" makes 2 from 0, and 0 from 2
			else complementary_sequence_matrix[wf][bazis]=-sequence_matrix[wf][L - 1 - bazis]+4; // this "function" makes 3 from 1, and 1 from			
		}
	}
}

void rand_fit_and_seq(void)
{	
	int wf;
	int jed;
	int rm;

	master_fitnesses=minimum_master_fitness+0*Delta;
	baseline_fitness=minimum_master_fitness*selection_coefficients[2];
	
	// all possible fitness values determined here: 
	fitness_array[0]=master_fitnesses*selection_coefficients[0];
	fitness_array[1]=master_fitnesses*selection_coefficients[1];
	fitness_array[2]=master_fitnesses*selection_coefficients[2];
	fitness_array[3]=baseline_fitness;
	
	// generate prepared but random sequences according to a decreasing GC content gradient:	
	for(wf = 0; wf < N_REPLICATOR_TYPES; wf++)
	{	
		if(wf==0) // First row of the matrix is (the first master) has 100% GC content and it is certainly not identical to any other row
		{
			for(jed=0;jed<L;jed++) 
			{                       // A G  U  C -> 0,1,2,3
				if(randd() < 0.5) sequence_matrix[wf][jed]=1; 
				else sequence_matrix[wf][jed]=3;			
			}
			shuffle (sequence_matrix[wf], L);
		}
		else 
		{
			for(jed=0;jed<L-wf;jed++) 
			{                       // A G  U  C -> 0,1,2,3
				if(randd() < 0.5) sequence_matrix[wf][jed]=1; 
				else sequence_matrix[wf][jed]=3;								
			}
			for(rm=L-wf;rm<L;rm++) 
			{                       // A G  U  C -> 0,1,2,3
				if(randd() < 0.5) sequence_matrix[wf][rm]=0;
				else sequence_matrix[wf][rm]=2;							
			}			
			shuffle (sequence_matrix[wf], L);
			while(generate_HD_larger_3_seqs(wf)==1) shuffle (sequence_matrix[wf], L);						
		}
	}
	generate_complementary_master_sequences();
	HDs_between_masters_and_complementary_strands();
	
	while(generate_master_complementary_HD_larger_3_seqs()==1)  	
	{
		// generate prepared but random sequences according to a decreasing GC content gradient:	
		for(wf = 0; wf < N_REPLICATOR_TYPES; wf++)
		{	
			if(wf==0) // First row of the matrix is (the first master) has 100% GC content and it is certainly not identical to any other row
			{
				for(jed=0;jed<L;jed++) 
				{                       // A G  U  C -> 0,1,2,3
					if(randd() < 0.5) sequence_matrix[wf][jed]=1; 
					else sequence_matrix[wf][jed]=3;			
				}
				shuffle (sequence_matrix[wf], L);
			}
			else 
			{
				for(jed=0;jed<L-wf;jed++) 
				{                       // A G  U  C -> 0,1,2,3
					if(randd() < 0.5) sequence_matrix[wf][jed]=1; 
					else sequence_matrix[wf][jed]=3;								
				}
				for(rm=L-wf;rm<L;rm++) 
				{                       // A G  U  C -> 0,1,2,3
					if(randd() < 0.5) sequence_matrix[wf][rm]=0;
					else sequence_matrix[wf][rm]=2;							
				}			
				shuffle (sequence_matrix[wf], L);
				while(generate_HD_larger_3_seqs(wf)==1) shuffle (sequence_matrix[wf], L);						
			}
		}
		generate_complementary_master_sequences();
		HDs_between_masters_and_complementary_strands();		
	}				
}

void init_population_vectors(void){
	int s, dd, spp, iii, i, kk, ssr;
	
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
		comp[iii].ind_H_Dist=0;
		comp[iii].deleted=0;
		comp[iii].energy_level=-11;
		simplex_lookup[iii]=iii;
	}
	
	// Init species-specific properties: 
	for(ssr = 1; ssr < N_REPLICATOR_TYPES+1; ssr++)
	{
		for(i = init_pop_sizes*(ssr-1); i < init_pop_sizes*ssr; i++) 
		{
			comp[i].master=ssr-1; // sequence specific feature
			for(kk=0;kk<L;kk++) comp[i].seq[kk]=sequence_matrix[ssr-1][kk];
		}
	}

	molecular_count_array[0]=SEQ_NUM;  // for Gillespie, which treats all master replications one reaction type
	molecular_count_array[1]=0; // mutant class HD=1
	molecular_count_array[2]=0; // mutant class HD=2
	molecular_count_array[3]=0; // mutant class HD>2
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
	int sqnc, sp;
	int HD_vector[N_REPLICATOR_TYPES];
	int min;
  	
  	comp[ind].master=-888;
	for(sp=0;sp<N_REPLICATOR_TYPES;sp++)
	{
		HD_vector[sp]=0;
		for(sqnc=0;sqnc<L;sqnc++)
		{
			if(comp[ind].seq[sqnc]!=sequence_matrix[sp][sqnc]) HD_vector[sp]++; //  the number of loci that is different from the master in the given sequence
			if(HD_vector[sp]>2) break; 
		}	
		if(HD_vector[sp]==0) comp[ind].master=sp;
	}	
	
	// choose the smallest Hamming Distance of the new replica:
    min = HD_vector[0];
    for (int i = 1; i < N_REPLICATOR_TYPES; i++) 
    {     
        //Compare elements of array with max    
		if(HD_vector[i] < min) 
			min = HD_vector[i];    
	}
	comp[ind].ind_H_Dist=min;		
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
	random_integer_seed_variable=1674356367;
	//random_integer_seed_variable=random_seed();
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
	BEs_between_masters();
	
	//############### # init_reaction_struct ############################################################################
	reaction_num=sum_assoc_and_dissoc+unique_fitnesses;
	REACTION reaction[reaction_num];
	
	// Init properties of the reactions (if a property remains -11, it means that it doesn't play a role in the reaction):
	for(nulls = 0; nulls < reaction_num; nulls++) 
	{
		reaction[nulls].duplex_bindingE=-11;
		reaction[nulls].reaction_Hd=-11;
	}

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
		reaction[sum_assoc_and_dissoc+reps].reaction_Hd = reps;
		reaction[sum_assoc_and_dissoc+reps].propensity = fitness_array[reps]*molecular_count_array[reps];
	}						
	//###############################################################################################################
	
	FILE *parameters;
	parameters=fopen("master_fit_and_seq.txt","a");
	fprintf(parameters, "seed %li: \n", random_integer_seed_variable);
	for(ik=0; ik<N_REPLICATOR_TYPES; ik++) {
		fprintf(parameters, "sp%d HD0 fitness: %f, sequence: ", ik, master_fitnesses);
		for(jk=0;jk<L;jk++) {
			fprintf(parameters, "%d ", sequence_matrix[ik][jk]);
			if(jk==L-1){
				fprintf(parameters, "\n");
			}
		}
	}

	fprintf(parameters, "Hamming-distances between master sequences: \n");	
	fprintf(parameters, "sp0 HDs: ");
	fprintf(parameters, "%d ", master_HDs[0]); // 0-0
	fprintf(parameters, "%d ", master_HDs[1]); // 0-1
	fprintf(parameters, "%d ", master_HDs[2]); // 0-2
	fprintf(parameters, "%d ", master_HDs[3]); // 0-3
	fprintf(parameters, "%d ", master_HDs[4]); // 0-4
	fprintf(parameters, "%d ", master_HDs[5]); // 0-5
	fprintf(parameters, "%d ", master_HDs[6]); // 0-6
	fprintf(parameters, "%d ", master_HDs[7]); // 0-7
	fprintf(parameters, "%d ", master_HDs[8]); // 0-8
	fprintf(parameters, "%d ", master_HDs[9]); // 0-9
	fprintf(parameters, "\n");
	
	fprintf(parameters, "sp1 HDs: ");
	fprintf(parameters, "%d ", master_HDs[1]);  // 1-0
	fprintf(parameters, "%d ", master_HDs[10]); // 1-1
	fprintf(parameters, "%d ", master_HDs[11]); // 1-2
	fprintf(parameters, "%d ", master_HDs[12]); // 1-3
	fprintf(parameters, "%d ", master_HDs[13]); // 1-4
	fprintf(parameters, "%d ", master_HDs[14]); // 1-5
	fprintf(parameters, "%d ", master_HDs[15]); // 1-6
	fprintf(parameters, "%d ", master_HDs[16]); // 1-7
	fprintf(parameters, "%d ", master_HDs[17]); // 1-8
	fprintf(parameters, "%d ", master_HDs[18]);	// 1-9	
	fprintf(parameters, "\n");		

	fprintf(parameters, "sp2 HDs: ");
	fprintf(parameters, "%d ", master_HDs[2]);  // 2-0
	fprintf(parameters, "%d ", master_HDs[11]); // 2-1
	fprintf(parameters, "%d ", master_HDs[19]); // 2-2
	fprintf(parameters, "%d ", master_HDs[20]); // 2-3
	fprintf(parameters, "%d ", master_HDs[21]); // 2-4
	fprintf(parameters, "%d ", master_HDs[22]); // 2-5
	fprintf(parameters, "%d ", master_HDs[23]); // 2-6
	fprintf(parameters, "%d ", master_HDs[24]); // 2-7
	fprintf(parameters, "%d ", master_HDs[25]); // 2-8
	fprintf(parameters, "%d ", master_HDs[26]);	// 2-9	
	fprintf(parameters, "\n");	
	
	fprintf(parameters, "sp3 HDs: ");
	fprintf(parameters, "%d ", master_HDs[3]);  // 3-0
	fprintf(parameters, "%d ", master_HDs[12]); // 3-1
	fprintf(parameters, "%d ", master_HDs[20]); // 3-2
	fprintf(parameters, "%d ", master_HDs[27]); // 3-3
	fprintf(parameters, "%d ", master_HDs[28]); // 3-4
	fprintf(parameters, "%d ", master_HDs[29]); // 3-5
	fprintf(parameters, "%d ", master_HDs[30]); // 3-6
	fprintf(parameters, "%d ", master_HDs[31]); // 3-7
	fprintf(parameters, "%d ", master_HDs[32]); // 3-8
	fprintf(parameters, "%d ", master_HDs[33]);	// 3-9	
	fprintf(parameters, "\n");		
	
	fprintf(parameters, "sp4 HDs: ");
	fprintf(parameters, "%d ", master_HDs[4]);  // 4-0
	fprintf(parameters, "%d ", master_HDs[13]); // 4-1
	fprintf(parameters, "%d ", master_HDs[21]); // 4-2
	fprintf(parameters, "%d ", master_HDs[28]); // 4-3
	fprintf(parameters, "%d ", master_HDs[34]); // 4-4
	fprintf(parameters, "%d ", master_HDs[35]); // 4-5
	fprintf(parameters, "%d ", master_HDs[36]); // 4-6
	fprintf(parameters, "%d ", master_HDs[37]); // 4-7
	fprintf(parameters, "%d ", master_HDs[38]); // 4-8
	fprintf(parameters, "%d ", master_HDs[39]);	// 4-9	
	fprintf(parameters, "\n");	

	fprintf(parameters, "sp5 HDs: ");
	fprintf(parameters, "%d ", master_HDs[5]);  // 5-0
	fprintf(parameters, "%d ", master_HDs[14]); // 5-1
	fprintf(parameters, "%d ", master_HDs[22]); // 5-2
	fprintf(parameters, "%d ", master_HDs[29]); // 5-3
	fprintf(parameters, "%d ", master_HDs[35]); // 5-4
	fprintf(parameters, "%d ", master_HDs[40]); // 5-5
	fprintf(parameters, "%d ", master_HDs[41]); // 5-6
	fprintf(parameters, "%d ", master_HDs[42]); // 5-7
	fprintf(parameters, "%d ", master_HDs[43]); // 5-8
	fprintf(parameters, "%d ", master_HDs[44]);	// 5-9	
	fprintf(parameters, "\n");	

	fprintf(parameters, "sp6 HDs: ");
	fprintf(parameters, "%d ", master_HDs[6]);  // 6-0
	fprintf(parameters, "%d ", master_HDs[15]); // 6-1
	fprintf(parameters, "%d ", master_HDs[23]); // 6-2
	fprintf(parameters, "%d ", master_HDs[30]); // 6-3
	fprintf(parameters, "%d ", master_HDs[36]); // 6-4
	fprintf(parameters, "%d ", master_HDs[41]); // 6-5
	fprintf(parameters, "%d ", master_HDs[45]); // 6-6
	fprintf(parameters, "%d ", master_HDs[46]); // 6-7
	fprintf(parameters, "%d ", master_HDs[47]); // 6-8
	fprintf(parameters, "%d ", master_HDs[48]);	// 6-9	
	fprintf(parameters, "\n");	

	fprintf(parameters, "sp7 HDs: ");
	fprintf(parameters, "%d ", master_HDs[7]);  // 7-0
	fprintf(parameters, "%d ", master_HDs[16]); // 7-1
	fprintf(parameters, "%d ", master_HDs[24]); // 7-2
	fprintf(parameters, "%d ", master_HDs[31]); // 7-3
	fprintf(parameters, "%d ", master_HDs[37]); // 7-4
	fprintf(parameters, "%d ", master_HDs[42]); // 7-5
	fprintf(parameters, "%d ", master_HDs[46]); // 7-6
	fprintf(parameters, "%d ", master_HDs[49]); // 7-7
	fprintf(parameters, "%d ", master_HDs[50]); // 7-8
	fprintf(parameters, "%d ", master_HDs[51]);	// 7-9	
	fprintf(parameters, "\n");	

	fprintf(parameters, "sp8 HDs: ");
	fprintf(parameters, "%d ", master_HDs[8]);  // 8-0
	fprintf(parameters, "%d ", master_HDs[17]); // 8-1
	fprintf(parameters, "%d ", master_HDs[25]); // 8-2
	fprintf(parameters, "%d ", master_HDs[32]); // 8-3
	fprintf(parameters, "%d ", master_HDs[38]); // 8-4
	fprintf(parameters, "%d ", master_HDs[43]); // 8-5
	fprintf(parameters, "%d ", master_HDs[47]); // 8-6
	fprintf(parameters, "%d ", master_HDs[50]); // 8-7
	fprintf(parameters, "%d ", master_HDs[52]); // 8-8
	fprintf(parameters, "%d ", master_HDs[53]);	// 8-9	
	fprintf(parameters, "\n");	

	fprintf(parameters, "sp9 HDs: ");
	fprintf(parameters, "%d ", master_HDs[9]);  // 9-0
	fprintf(parameters, "%d ", master_HDs[18]); // 9-1
	fprintf(parameters, "%d ", master_HDs[26]); // 9-2
	fprintf(parameters, "%d ", master_HDs[33]); // 9-3
	fprintf(parameters, "%d ", master_HDs[39]); // 9-4
	fprintf(parameters, "%d ", master_HDs[44]); // 9-5
	fprintf(parameters, "%d ", master_HDs[48]); // 9-6
	fprintf(parameters, "%d ", master_HDs[51]); // 9-7
	fprintf(parameters, "%d ", master_HDs[53]); // 9-8
	fprintf(parameters, "%d ", master_HDs[54]);	// 9-9	
	fprintf(parameters, "\n");

	for(chbe=0; chbe<55; chbe++) 
	{
		sum_of_master_HDs+=master_HDs[chbe];
	}
	fprintf(parameters, "mean Hamming-distances between masters: %f\n", (double)sum_of_master_HDs/(double)45.0);	
	
	fprintf(parameters, "Binding energies between master sequences: \n");	
	fprintf(parameters, "sp0 BEs: ");
	fprintf(parameters, "%d ", master_BEs[0]); // 0-0
	fprintf(parameters, "%d ", master_BEs[1]); // 0-1
	fprintf(parameters, "%d ", master_BEs[2]); // 0-2
	fprintf(parameters, "%d ", master_BEs[3]); // 0-3
	fprintf(parameters, "%d ", master_BEs[4]); // 0-4
	fprintf(parameters, "%d ", master_BEs[5]); // 0-5
	fprintf(parameters, "%d ", master_BEs[6]); // 0-6
	fprintf(parameters, "%d ", master_BEs[7]); // 0-7
	fprintf(parameters, "%d ", master_BEs[8]); // 0-8
	fprintf(parameters, "%d ", master_BEs[9]); // 0-9
	fprintf(parameters, "  mean BE of sp0: %f", 
	(double)(master_BEs[0]+master_BEs[1]+master_BEs[2]+master_BEs[3]+master_BEs[4]+master_BEs[5]+master_BEs[6]+master_BEs[7]+master_BEs[8]+master_BEs[9])/(double)10.0);
	fprintf(parameters, "\n");
	
	fprintf(parameters, "sp1 BEs: ");
	fprintf(parameters, "%d ", master_BEs[1]);  // 1-0
	fprintf(parameters, "%d ", master_BEs[10]); // 1-1
	fprintf(parameters, "%d ", master_BEs[11]); // 1-2
	fprintf(parameters, "%d ", master_BEs[12]); // 1-3
	fprintf(parameters, "%d ", master_BEs[13]); // 1-4
	fprintf(parameters, "%d ", master_BEs[14]); // 1-5
	fprintf(parameters, "%d ", master_BEs[15]); // 1-6
	fprintf(parameters, "%d ", master_BEs[16]); // 1-7
	fprintf(parameters, "%d ", master_BEs[17]); // 1-8
	fprintf(parameters, "%d ", master_BEs[18]);	// 1-9	
	fprintf(parameters, "  mean BE of sp1: %f", 
	(double)(master_BEs[1]+master_BEs[10]+master_BEs[11]+master_BEs[12]+master_BEs[13]+master_BEs[14]+master_BEs[15]+master_BEs[16]+master_BEs[17]+master_BEs[18])/(double)10.0);
	fprintf(parameters, "\n");

	fprintf(parameters, "sp2 BEs: ");
	fprintf(parameters, "%d ", master_BEs[2]);  // 2-0
	fprintf(parameters, "%d ", master_BEs[11]); // 2-1
	fprintf(parameters, "%d ", master_BEs[19]); // 2-2
	fprintf(parameters, "%d ", master_BEs[20]); // 2-3
	fprintf(parameters, "%d ", master_BEs[21]); // 2-4
	fprintf(parameters, "%d ", master_BEs[22]); // 2-5
	fprintf(parameters, "%d ", master_BEs[23]); // 2-6
	fprintf(parameters, "%d ", master_BEs[24]); // 2-7
	fprintf(parameters, "%d ", master_BEs[25]); // 2-8
	fprintf(parameters, "%d ", master_BEs[26]);	// 2-9	
	fprintf(parameters, "  mean BE of sp2: %f", 
	(double)(master_BEs[2]+master_BEs[11]+master_BEs[19]+master_BEs[20]+master_BEs[21]+master_BEs[22]+master_BEs[23]+master_BEs[24]+master_BEs[25]+master_BEs[26])/(double)10.0);
	fprintf(parameters, "\n");
	
	fprintf(parameters, "sp3 BEs: ");
	fprintf(parameters, "%d ", master_BEs[3]);  // 3-0
	fprintf(parameters, "%d ", master_BEs[12]); // 3-1
	fprintf(parameters, "%d ", master_BEs[20]); // 3-2
	fprintf(parameters, "%d ", master_BEs[27]); // 3-3
	fprintf(parameters, "%d ", master_BEs[28]); // 3-4
	fprintf(parameters, "%d ", master_BEs[29]); // 3-5
	fprintf(parameters, "%d ", master_BEs[30]); // 3-6
	fprintf(parameters, "%d ", master_BEs[31]); // 3-7
	fprintf(parameters, "%d ", master_BEs[32]); // 3-8
	fprintf(parameters, "%d ", master_BEs[33]);	// 3-9	
	fprintf(parameters, "  mean BE of sp3: %f", 
	(double)(master_BEs[3]+master_BEs[12]+master_BEs[20]+master_BEs[27]+master_BEs[28]+master_BEs[29]+master_BEs[30]+master_BEs[31]+master_BEs[32]+master_BEs[33])/(double)10.0);
	fprintf(parameters, "\n");	
	
	fprintf(parameters, "sp4 BEs: ");
	fprintf(parameters, "%d ", master_BEs[4]);  // 4-0
	fprintf(parameters, "%d ", master_BEs[13]); // 4-1
	fprintf(parameters, "%d ", master_BEs[21]); // 4-2
	fprintf(parameters, "%d ", master_BEs[28]); // 4-3
	fprintf(parameters, "%d ", master_BEs[34]); // 4-4
	fprintf(parameters, "%d ", master_BEs[35]); // 4-5
	fprintf(parameters, "%d ", master_BEs[36]); // 4-6
	fprintf(parameters, "%d ", master_BEs[37]); // 4-7
	fprintf(parameters, "%d ", master_BEs[38]); // 4-8
	fprintf(parameters, "%d ", master_BEs[39]);	// 4-9	
	fprintf(parameters, "  mean BE of sp4: %f", 
	(double)(master_BEs[4]+master_BEs[13]+master_BEs[21]+master_BEs[28]+master_BEs[34]+master_BEs[35]+master_BEs[36]+master_BEs[37]+master_BEs[38]+master_BEs[39])/(double)10.0);
	fprintf(parameters, "\n");	

	fprintf(parameters, "sp5 BEs: ");
	fprintf(parameters, "%d ", master_BEs[5]);  // 5-0
	fprintf(parameters, "%d ", master_BEs[14]); // 5-1
	fprintf(parameters, "%d ", master_BEs[22]); // 5-2
	fprintf(parameters, "%d ", master_BEs[29]); // 5-3
	fprintf(parameters, "%d ", master_BEs[35]); // 5-4
	fprintf(parameters, "%d ", master_BEs[40]); // 5-5
	fprintf(parameters, "%d ", master_BEs[41]); // 5-6
	fprintf(parameters, "%d ", master_BEs[42]); // 5-7
	fprintf(parameters, "%d ", master_BEs[43]); // 5-8
	fprintf(parameters, "%d ", master_BEs[44]);	// 5-9	
	fprintf(parameters, "  mean BE of sp5: %f", 
	(double)(master_BEs[5]+master_BEs[14]+master_BEs[22]+master_BEs[29]+master_BEs[35]+master_BEs[40]+master_BEs[41]+master_BEs[42]+master_BEs[43]+master_BEs[44])/(double)10.0);
	fprintf(parameters, "\n");	

	fprintf(parameters, "sp6 BEs: ");
	fprintf(parameters, "%d ", master_BEs[6]);  // 6-0
	fprintf(parameters, "%d ", master_BEs[15]); // 6-1
	fprintf(parameters, "%d ", master_BEs[23]); // 6-2
	fprintf(parameters, "%d ", master_BEs[30]); // 6-3
	fprintf(parameters, "%d ", master_BEs[36]); // 6-4
	fprintf(parameters, "%d ", master_BEs[41]); // 6-5
	fprintf(parameters, "%d ", master_BEs[45]); // 6-6
	fprintf(parameters, "%d ", master_BEs[46]); // 6-7
	fprintf(parameters, "%d ", master_BEs[47]); // 6-8
	fprintf(parameters, "%d ", master_BEs[48]);	// 6-9	
	fprintf(parameters, "  mean BE of sp6: %f", 
	(double)(master_BEs[6]+master_BEs[15]+master_BEs[23]+master_BEs[30]+master_BEs[36]+master_BEs[41]+master_BEs[45]+master_BEs[46]+master_BEs[47]+master_BEs[48])/(double)10.0);
	fprintf(parameters, "\n");	

	fprintf(parameters, "sp7 BEs: ");
	fprintf(parameters, "%d ", master_BEs[7]);  // 7-0
	fprintf(parameters, "%d ", master_BEs[16]); // 7-1
	fprintf(parameters, "%d ", master_BEs[24]); // 7-2
	fprintf(parameters, "%d ", master_BEs[31]); // 7-3
	fprintf(parameters, "%d ", master_BEs[37]); // 7-4
	fprintf(parameters, "%d ", master_BEs[42]); // 7-5
	fprintf(parameters, "%d ", master_BEs[46]); // 7-6
	fprintf(parameters, "%d ", master_BEs[49]); // 7-7
	fprintf(parameters, "%d ", master_BEs[50]); // 7-8
	fprintf(parameters, "%d ", master_BEs[51]);	// 7-9	
	fprintf(parameters, "  mean BE of sp7: %f", 
	(double)(master_BEs[7]+master_BEs[16]+master_BEs[24]+master_BEs[31]+master_BEs[37]+master_BEs[42]+master_BEs[46]+master_BEs[49]+master_BEs[50]+master_BEs[51])/(double)10.0);
	fprintf(parameters, "\n");	

	fprintf(parameters, "sp8 BEs: ");
	fprintf(parameters, "%d ", master_BEs[8]);  // 8-0
	fprintf(parameters, "%d ", master_BEs[17]); // 8-1
	fprintf(parameters, "%d ", master_BEs[25]); // 8-2
	fprintf(parameters, "%d ", master_BEs[32]); // 8-3
	fprintf(parameters, "%d ", master_BEs[38]); // 8-4
	fprintf(parameters, "%d ", master_BEs[43]); // 8-5
	fprintf(parameters, "%d ", master_BEs[47]); // 8-6
	fprintf(parameters, "%d ", master_BEs[50]); // 8-7
	fprintf(parameters, "%d ", master_BEs[52]); // 8-8
	fprintf(parameters, "%d ", master_BEs[53]);	// 8-9	
	fprintf(parameters, "  mean BE of sp8: %f", 
	(double)(master_BEs[8]+master_BEs[17]+master_BEs[25]+master_BEs[32]+master_BEs[38]+master_BEs[43]+master_BEs[47]+master_BEs[50]+master_BEs[52]+master_BEs[53])/(double)10.0);
	fprintf(parameters, "\n");	

	fprintf(parameters, "sp9 BEs: ");
	fprintf(parameters, "%d ", master_BEs[9]);  // 9-0
	fprintf(parameters, "%d ", master_BEs[18]); // 9-1
	fprintf(parameters, "%d ", master_BEs[26]); // 9-2
	fprintf(parameters, "%d ", master_BEs[33]); // 9-3
	fprintf(parameters, "%d ", master_BEs[39]); // 9-4
	fprintf(parameters, "%d ", master_BEs[44]); // 9-5
	fprintf(parameters, "%d ", master_BEs[48]); // 9-6
	fprintf(parameters, "%d ", master_BEs[51]); // 9-7
	fprintf(parameters, "%d ", master_BEs[53]); // 9-8
	fprintf(parameters, "%d ", master_BEs[54]);	// 9-9	
	fprintf(parameters, "  mean BE of sp9: %f", 
	(double)(master_BEs[9]+master_BEs[18]+master_BEs[26]+master_BEs[33]+master_BEs[39]+master_BEs[44]+master_BEs[48]+master_BEs[51]+master_BEs[53]+master_BEs[54])/(double)10.0);
	fprintf(parameters, "\n");
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
				molecular_count_array[comp[simplex_lookup[a_one]].ind_H_Dist]--; 
				reaction[sum_assoc_and_dissoc+comp[simplex_lookup[a_one]].ind_H_Dist].propensity 
				 = fitness_array[comp[simplex_lookup[a_one]].ind_H_Dist]*molecular_count_array[comp[simplex_lookup[a_one]].ind_H_Dist];	
	
				// simplex2:
				molecular_count_array[comp[simplex_lookup[z_one]].ind_H_Dist]--; 
				reaction[sum_assoc_and_dissoc+comp[simplex_lookup[z_one]].ind_H_Dist].propensity 
				 = fitness_array[comp[simplex_lookup[z_one]].ind_H_Dist]*molecular_count_array[comp[simplex_lookup[z_one]].ind_H_Dist];					
					
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
			molecular_count_array[comp[randpair].ind_H_Dist]++; 
			reaction[sum_assoc_and_dissoc+comp[randpair].ind_H_Dist].propensity 
			 = fitness_array[comp[randpair].ind_H_Dist]*molecular_count_array[comp[randpair].ind_H_Dist];	

			// simplex2:
			molecular_count_array[comp[randdiss].ind_H_Dist]++;
			reaction[sum_assoc_and_dissoc+comp[randdiss].ind_H_Dist].propensity 
			 = fitness_array[comp[randdiss].ind_H_Dist]*molecular_count_array[comp[randdiss].ind_H_Dist];					
		
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
					while(comp[simplex_lookup[i]].ind_H_Dist!=reaction[event].reaction_Hd) i = randl(simplex_state); // choose individual for replication
					j = randl(SEQ_NUM); // Choose an individual to be overwrited (Moran-process)
					while(simplex_lookup[i] == j) j = randl(SEQ_NUM); // Choose another individual than i (we know already that there is no deleted sequence now: del_pos=-333)
					if(comp[j].pair==-999)  //the overwritten individual is simplex
					{  
						molecular_count_array[comp[j].ind_H_Dist]--; 
						reaction[sum_assoc_and_dissoc+comp[j].ind_H_Dist].propensity 
						 = fitness_array[comp[j].ind_H_Dist]*molecular_count_array[comp[j].ind_H_Dist];							
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
							reaction[event].propensity = fitness_array[event-sum_assoc_and_dissoc]*molecular_count_array[event-sum_assoc_and_dissoc];							
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
							molecular_count_array[comp[j].ind_H_Dist]++; 
							reaction[sum_assoc_and_dissoc+comp[j].ind_H_Dist].propensity 
							 = fitness_array[comp[j].ind_H_Dist]*molecular_count_array[comp[j].ind_H_Dist];							
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
							reaction[event].propensity = fitness_array[event-sum_assoc_and_dissoc]*molecular_count_array[event-sum_assoc_and_dissoc];							
							shift(i); 
							simplex_state--;													 
						}else{   // overwritten Gillespie agent after replication and mutation stays simple stranded
							comp[j].pair=-999;	
							molecular_count_array[comp[j].ind_H_Dist]++; // 
							reaction[sum_assoc_and_dissoc+comp[j].ind_H_Dist].propensity 
							 = fitness_array[comp[j].ind_H_Dist]*molecular_count_array[comp[j].ind_H_Dist];																										
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
					while(comp[simplex_lookup[i]].ind_H_Dist!=reaction[event].reaction_Hd) i = randl(simplex_state); // choose individual for replication
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
						reaction[event].propensity = fitness_array[event-sum_assoc_and_dissoc]*molecular_count_array[event-sum_assoc_and_dissoc];						
						reaction[1+duplex_index].propensity=c_d[duplex_index]*duplex_species[duplex_index];  // update duplex propensities	
						duplex_state+=2; // total concentration of duplexes increases	
						shift(i); 
						simplex_state--;			
					}else{ // the new replica immediately dissociates and stays simple stranded
						molecular_count_array[comp[del_pos].ind_H_Dist]++; 
						reaction[sum_assoc_and_dissoc+comp[del_pos].ind_H_Dist].propensity 
						 = fitness_array[comp[del_pos].ind_H_Dist]*molecular_count_array[comp[del_pos].ind_H_Dist];													
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
	//free(sorted_fitness_array);
	fprintf(file, "\nsurvived types:\n");
	for(ik = 0; ik < N_REPLICATOR_TYPES ; ik++) 
	{
		if(survival_vector[ik]>1) 
		{
			surv_num++;
			sum_of_surviving_indeces+=ik+1,
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
	surv = fopen("survival_stat__delta_0_005__N10000.txt", "a");
	fprintf(surv, "%.1f, ", ((double)surv_num)/((double)10));
	fclose(surv);
	FILE *sim_dup_file = NULL;
	sim_dup_file = fopen("simplex_duplex_ratio_mean__delta_0_005__N10000.txt", "a");
	fprintf(sim_dup_file, "%f, ", ((double)sum_simp_dup)/((double)sum_simp_dup_counter));
	fclose(sim_dup_file);	
	FILE *sum_master_file = NULL;
	sum_master_file = fopen("mean_sum_of_all_masters__delta_0_005__N10000.txt", "a");
	fprintf(sum_master_file, "%f, ", ((double)add_rel_sum_of_all_masters)/((double)sum_simp_dup_counter));
	fclose(sum_master_file);
	return(0);
}//main end

// compile: gcc equal_fit_code__grad_GC_BUT_RANDOM_sequences_shuffled.c -lm -Wall -Ofast -o GC_gradient__N10000
// run: nohup ./GC_gradient__N10000 &  


