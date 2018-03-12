

#ifndef MAIN_C
#define MAIN_C

#include "../heads/flowsort.h"

/* Uncomment the line below to use PROM_II. Otherwise PROM I is the default PROMETHEE used. */


/* display the results */
#define PRINT_STUFFS

/** read_data
 *
 * Read instances given as arguments and return the right form to be used by the application.
 *
 * @param : the text file (see ./data/README for more.) depicting the data for the instance.
 * @return : a usable form of data for PROMETHEE : a vector of matrices of size M. Each matrix is of size (K,N) and represents the preference values given by the current expert l in 1..M to the project i in 1..N on the criteria i in 1..K. Such matrix is denoted E_l, and the vector of all such matrices is denoted E.
 *
*/
data read_data(char** argv)
{
	int i=0,j=0,k=0;
	float buff=0;
	data E;
	FILE* fp;

	fp=fopen(argv[4],"r");

	if( fscanf(fp,"%d %d %d", &M, &K, &N) < 0 ){LEAVE_FAIL(__LINE__, __FILE__);}

	/* assign values for thresholds and grades of level criteria */
	THRESHOLDS=malloc(K*sizeof(float*));
	if(THRESHOLDS == NULL){ /* memory allocation failure */ PRINT_MEM_FAIL(__LINE__, __FILE__); }

	if(STOCH)
	{
		STDs=malloc(K*sizeof(float*));
		if(STDs == NULL){ /* memory allocation failure */ PRINT_MEM_FAIL(__LINE__, __FILE__); }
	}

	for(i=0;i<K;i++)
	{/* get the thresholds for each criterion */
		THRESHOLDS[i]=calloc(3,sizeof(float));/* -q, q and p for each criteria */
		if(THRESHOLDS[i] == NULL){ /* memory allocation failure */ PRINT_MEM_FAIL(__LINE__, __FILE__); }
		if( fscanf(fp,"%f %f", &THRESHOLDS[i][1], &THRESHOLDS[i][2]) < 0 ){LEAVE_FAIL(__LINE__, __FILE__);}
		THRESHOLDS[i][0]= -(THRESHOLDS[i][1]);

		if(STOCH)
		{
			STDs[i]=calloc(N,sizeof(float));
			if(STDs[i] == NULL){ /* memory allocation failure */ PRINT_MEM_FAIL(__LINE__, __FILE__); }
		}
	}

	LEV_CRIT_GRADS=calloc(4,sizeof(float));
	if(LEV_CRIT_GRADS == NULL){ /* memory allocation failure */ PRINT_MEM_FAIL(__LINE__, __FILE__); }

	LEV_CRIT_GRADS[0]=NO_PREF; LEV_CRIT_GRADS[1]=INDIFF; LEV_CRIT_GRADS[2]=WEAK_PREF; LEV_CRIT_GRADS[3]=PREF;

	E=malloc(M*sizeof(expert));
	if(E == NULL){ /* memory allocation failure */ PRINT_MEM_FAIL(__LINE__, __FILE__); }


	if(COMPLETE_PREORDER)
	{
		/* allocate memory for the ranking and "S" when complete preorder PROM II is used */
		RANKS=malloc(M*sizeof(int*));
		if(RANKS == NULL){ /* memory allocation failure */ PRINT_MEM_FAIL(__LINE__, __FILE__); }

		COMPLETE_PREORDER_S=calloc(N,sizeof(float));
		if(COMPLETE_PREORDER_S == NULL){ /* memory allocation failure */ PRINT_MEM_FAIL(__LINE__, __FILE__); }

		for(i=0;i<M;i++)
		{
			RANKS[i]=calloc(N,sizeof(int));
			if(RANKS[i] == NULL){ /* memory allocation failure */ PRINT_MEM_FAIL(__LINE__, __FILE__); }
		}
	}

	/* alloc memory for each expert */
	for(i=0;i<M;i++)
	{
		E[i].p_l=0.;

		E[i].e_ij=malloc(K*sizeof(float*));
		if(E[i].e_ij == NULL){ /* memory allocation failure */ PRINT_MEM_FAIL(__LINE__, __FILE__); }

		E[i].S_l=NULL;

		/* alloc memory for each criteria */
		for(j=0;j<K;j++)
		{
			E[i].e_ij[j]=calloc(N,sizeof(float));
			if(E[i].e_ij[j] == NULL){ /* memory allocation failure */ PRINT_MEM_FAIL(__LINE__, __FILE__); }
		}
	}

	/* alloc memory for the criterion weights */
	criterion_weights=calloc(K,sizeof(float));
	if(criterion_weights == NULL){ /* memory allocation failure */ PRINT_MEM_FAIL(__LINE__, __FILE__); }

	buff=0;

	/* get criterion weights */
	for(i=0;i<K;i++)
	{
		if( fscanf(fp,"%f", &buff) < 0 ){LEAVE_FAIL(__LINE__, __FILE__);};
		criterion_weights[i]=buff;
	}

	/* gather infos from experts */
	for(i=0;i<M;i++)
	{
		if( fscanf(fp,"%f", &buff) < 0 ){LEAVE_FAIL(__LINE__, __FILE__);};
		E[i].p_l=buff;

		/* for each criteria... */
		for(j=0;j<K;j++)
		{
			/* get the evaluation done by expert l to action k with respect to criterion j */
			for(k=0;k<N;k++)
			{
				if(STOCH)
				{/* If stochastic case : then all values are combined expected evalutations and their standard deviation (separated by a comma) */
					if( fscanf(fp,"%f,%f", &buff, &STDs[j][k]) < 0 ){LEAVE_FAIL(__LINE__, __FILE__);};
				}
				else
					if( fscanf(fp,"%f", &buff) < 0 ){LEAVE_FAIL(__LINE__, __FILE__);};
				E[i].e_ij[j][k]=buff;
			}
		}
	}
	
	fclose(fp);

	return E;
}


/** read_sorting_data
 *
 * Read the reference profiles
 *
 * @param : the text file (see ./data/README for more.) depicting the data for the instance.
 * @return : a usable form of data for PROMETHEE : a vector of matrices of size M. Each matrix is of size (K,N) and represents the preference values given by the current expert l in 1..M to the project i in 1..N on the criteria i in 1..K. Such matrix is denoted E_l, and the vector of all such matrices is denoted E.
 *
*/
void read_sorting_data(char** argv)
{
	int i=0,j=0;
	FILE* fp;

	fp=fopen(argv[5],"r");

	if( fscanf(fp,"%d", &N_CAT) < 0 ){LEAVE_FAIL(__LINE__, __FILE__);}

	/* assign values for the reference profiles */
	R=malloc(N_CAT*sizeof(float*));
	C=malloc(3*sizeof(ptr_cat*));
	act_to_cat=malloc(N*sizeof(int*));

	if(STOCH)
	{/* referecence profiles are also pair of : an expected value and its standard deviation */
		R_STDs=malloc(N_CAT*sizeof(float*));
		if(R_STDs == NULL){ /* memory allocation failure */ PRINT_MEM_FAIL(__LINE__, __FILE__); }
	}

	for(j=0; j<3; j++)
	{
		C[j]=malloc((N_CAT-1)*sizeof(ptr_cat*));
		for(i=0; i<N_CAT-1; i++)
			C[j][i]=NULL;/* at first a category is secured by assigning a NULL pointer */
	}

    for(i=0; i<N_CAT; i++)
    {
		act_to_cat[i]=calloc(3,sizeof(int));
		R[i]=calloc(K,sizeof(float));
		if(STOCH)
		{
			R_STDs[i]=calloc(K,sizeof(float));
			if(R_STDs[i] == NULL){ /* memory allocation failure */ PRINT_MEM_FAIL(__LINE__, __FILE__); }
			for(j=0; j<K; j++)  if( fscanf(fp,"%f,%f", &R[i][j],&R_STDs[i][j]) < 0 ){LEAVE_FAIL(__LINE__, __FILE__);}
		}
		else
			for(j=0; j<K; j++)  if( fscanf(fp,"%f", &R[i][j]) < 0 ){LEAVE_FAIL(__LINE__, __FILE__);}
	}

	fclose(fp);
}


int main(int argc, char** argv)
{
	data E;
	int l=0, i=0;

	float*** P_l=NULL;
	float** PI=NULL;
	float** PHI=NULL;

	/* used if COMPLETE_PREORDER */
	float** S;

	/* at least 4 arguments are needed : at 0 the executable, 1st. the shape of the preference function, 2nd. PROM I/II, 3rd. the numerical instance, */
	if( argc < 5 ){printf("Please add more parameters\n");return 0;};

	/* stochastic case */
	if( atoi(argv[2]) == 1 )
		STOCH=true;

	/* Assign the available shapes for the preference function */
	if(STOCH)
	{
		TBL_PREF_FUNC_EPF[0]=EPF_level_criterion;
		TBL_PREF_FUNC_EPF[1]=EPF_linear_criterion;
	}
	else
	{
		TBL_PREF_FUNC[0]=level_criterion;
		TBL_PREF_FUNC[1]=linear_criterion;
	}
	/* By default, use the level criterion type */
	TYPE_CRITERION=0;
	/* Change otherwise */
	if( atoi(argv[1]) == 1 )
		TYPE_CRITERION=1;

	if( atoi(argv[3]) == 1 )
		COMPLETE_PREORDER=true;

	E=read_data(argv);
	
	/* When a sorting is required, read all data as well. */
	if(argc > 5)
		read_sorting_data(argv);


#ifdef PRINT_STUFFS
print_data(E);/* Since data is of type "expert*" one can send directly the object that is a pointer, and no copies will be created. */
#endif

	for(l=0;l<M;l++)
	{

#ifdef PRINT_STUFFS
printf("Expert %d\n",l);
#endif

		P_l=malloc(K*sizeof(float**));
		if(P_l == NULL){ /* memory allocation failure */ PRINT_MEM_FAIL(__LINE__, __FILE__); }

		for(i=0;i<K;i++)
		{/* computes P_1_l,..., P_K_l and gathers result into P_l, wich means, for each expert, compute the multicriteria preference degree */
			if(STOCH)
				P_l[i]=EPF(E[l].e_ij[i], i, false);
			else
				P_l[i]=pref_func(E[l].e_ij[i], i, false);
		}

#ifdef PRINT_STUFFS
for(i=0;i<K;i++)
print_P_l(P_l[i]);
#endif

		/* Reminder : PI(a,b) = [ ( sum{j in 1 to K}(w[j]*P[j](a,b)) )/( sum{i in 1 to K} (w[j]) ) ] and is of size n*n.
			PI(a,b) ~ P_j(a,b) considering simultaneously all the criteria */
		PI=compute_pref_indices(P_l, false);

		/* dealloc memory provided for P_l[i] */
		for(i=0;i<K;i++)
            free_float_matrix(P_l[i],N);

		/* dealloc memory allocated to P_l, caution it has to be deallocated after deallocating all P_l[i] */
		free(P_l);
		P_l=NULL;

		PHI=compute_phi(PI, false);
		/* Here PI is not freed for it will be needed later */

#ifdef PRINT_STUFFS
print_PI(PI);
#endif

		/* apply PROMETHEE I or II ? */
		if(COMPLETE_PREORDER)
			PROM_2((E+l), PHI);
		else
			PROM_1((E+l), PHI);

#ifdef PRINT_STUFFS
print_PHI(PHI);
#endif

#ifdef PRINT_STUFFS
print_S_l(E[l].S_l);
#endif

		/* if COMPLETE_PREORDER, a ranking can be made */
		if(COMPLETE_PREORDER)
			rank_single(E[l].S_l, l);/* ranks the projects with respect to opinion of l and store the ranking into "RANKS" list. */

		write_prom_results(l, PI, PHI, E[l].S_l, argv);

		/* dealloc memory for PI */
		free_float_matrix(PI, N);

		/* dealloc PHI, 3 : pos,neg,net flows */
		free_float_matrix(PHI, 3);
	}

	/* PROMETHEE II */
	if(COMPLETE_PREORDER)
	{
		get_S_a(E);/* if the ranking is a complete preorder, the computing can be done based on the rankings. */

#ifdef PRINT_STUFFS
printf("complete preorder (PROM II) :\n");
print_ranks();
#endif
	}
	else
	{
		/* PROMETHEE I */
		S=aggregate_S_l(E);/* deallocating E[l].S_l in memory is done during the aggregation */

#ifdef PRINT_STUFFS
print_aggregated_S_l(S);
#endif
	}

	write_prom_global_results(argv, S);

	/* deallocate with respect to the type of preorder performed */
	if(COMPLETE_PREORDER)
		free_alloc_for_PROM_II();
	else
		free_float_matrix(S, N+1);/* N+1 as (N+1)th contains the results of the summation */

	/* Apply flowsort if appropriate */
	if(argc > 5)
	{
		for(l=0;l<M;l++)
		{
			P_l=malloc(K*sizeof(float**));
			if(P_l == NULL){ /* Checking memory allocation failures */ PRINT_MEM_FAIL(__LINE__, __FILE__); }

			for(i=0;i<K;i++)
			{/* computes P_1_l,..., P_K_l and... */
				if(STOCH)
					P_l[i]=EPF(E[l].e_ij[i], i, true);
				else
					P_l[i]=pref_func(E[l].e_ij[i], i, true);
			}

			/* gathers all of these results into a single P_l (but still for each expert separately) */
			PI=compute_pref_indices(P_l, true);

			/* Dealloc P_l[i] */
			for(i=0;i<K;i++)
				free_float_matrix(P_l[i],N);

			/* and P_l */
			free(P_l);
			P_l=NULL;

#ifdef PRINT_STUFFS
printf("Expert %d\n",l);
print_PI_FlowSort(PI);
#endif

			PHI=compute_phi(PI, true);
			flowsort(PHI);

#ifdef PRINT_STUFFS
print_PHI_FlowSort(PHI);
#endif

			/* write results obtained by flowsort */
			write_flowsort_results(PI, PHI, argv);

			/* dealloc PI and PHI */
			free_float_matrix(PI,(2*N));
			free_float_matrix(PHI,3);
		}
	}

	/* free remaining data ;P */
	free_data(E);

	return EXIT_SUCCESS;
}

#endif
