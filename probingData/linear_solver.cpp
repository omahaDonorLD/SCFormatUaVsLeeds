
#include "linear_solver.hpp"

double euclDistance(double *node1, double *node2)
{
	int i=0;
	double norm=0;
	for(i=0;i<dim;i++)	norm+=pow(node1[i]-node2[i],2);
    return sqrt(norm);
};


void solve_linear_model(sln* net, double range, double lb, double* soln)
{
	int i=0,j=0;

	/* Build it first */
	glp_prob *lp;
	int ia[nbr_grnds+1],ja[net->n_uavs+1];
	double ar[nbr_grnds+1], z=0.;
	
	lp = glp_create_prob();
	glp_set_prob_name(lp, "set cover");
	glp_set_obj_dir(lp, GLP_MIN);

	/* row names */
	char row_names[nbr_grnds][100];
	char col_names[net->n_uavs][100];
	char *buff;

	glp_add_rows(lp, nbr_grnds);
	for(i=1;i<=nbr_grnds;i++)
	{
		sprintf(buff,"r%d",i);
		strcpy(row_names[i],buff);
		glp_set_row_name(lp, i, row_names[i]);
		glp_set_row_bnds(lp, 1, GLP_LO, lb, 0.0);
	}

	glp_add_cols(lp, net->n_uavs);
	for(j=1;j<=net->n_uavs;j++)
	{
		sprintf(buff,"x%d",j);
		strcpy(col_names[j],buff);
		glp_set_col_name(lp, j, col_names[j]);
		glp_set_col_bnds(lp, j, GLP_DB, 0.0, 1.0);
		glp_set_obj_coef(lp, j, 1.0);
		glp_set_col_kind(lp, j, GLP_BV);
		//glp_set_col_kind(lp, j, GLP_BV);
	}

	int counter=1;
	for(i=1;i<=nbr_grnds;i++)
	{
		for(j=1;j<=net->n_uavs;j++)
		{
			ia[counter] = i;
			ja[counter] = j;
			ar[counter] = ( euclDistance(net->uavs[j],grnds[i]) <= range ? 1.0 : 0.0 );
		}
		counter++;
	}

	glp_load_matrix(lp, 9, ia, ja, ar);
	int res_solver=glp_simplex(lp, NULL);
	z = glp_get_obj_val(lp);
	
	/* Gather results */
	soln=(double*)calloc(net->n_uavs+1,sizeof(double));
	for(j=1;j<=net->n_uavs;j++)
		soln[j]= glp_get_col_prim(lp, j);
	
	glp_delete_prob(lp);

};
