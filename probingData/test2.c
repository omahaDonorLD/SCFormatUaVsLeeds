
#include <stdio.h>
#include <stdlib.h>
#include <glpk.h>

int main(void)
{

	//--=--------TEST START

		/* Build it first */
	glp_prob *lp1;
	int ia1[7],ja1[7];
	double ar1[25];
	
	lp1 = glp_create_prob();
	glp_set_prob_name(lp1, "test smallest set cover");
	glp_set_obj_dir(lp1, GLP_MIN);

	glp_add_rows(lp1, 2);

	glp_set_row_name(lp1, 1, "r1");
	glp_set_row_bnds(lp1, 1, GLP_LO, 1.0, 0.0);
	glp_set_row_name(lp1, 2, "r2");
	glp_set_row_bnds(lp1, 2, GLP_LO, 1.0, 0.0);

	glp_add_cols(lp1, 3);

	glp_set_col_name(lp1, 1, "x1");
	glp_set_col_bnds(lp1, 1, GLP_DB, 0.0, 1.0);
	glp_set_obj_coef(lp1, 1, 1.0);
	//glp_set_col_kind(lp, j, GLP_BV);
	//glp_set_col_kind(lp, j, GLP_BV);

	glp_set_col_name(lp1, 2, "x2");
	glp_set_col_bnds(lp1, 2, GLP_DB, 0.0, 1.0);
	glp_set_obj_coef(lp1, 2, 1.0);
	//glp_set_col_kind(lp, j, GLP_BV);
	//glp_set_col_kind(lp, j, GLP_BV);

	glp_set_col_name(lp1, 3, "x3");
	glp_set_col_bnds(lp1, 3, GLP_DB, 0.0, 1.0);
	glp_set_obj_coef(lp1, 3, 1.0);
	//glp_set_col_kind(lp, j, GLP_BV);
	//glp_set_col_kind(lp, j, GLP_BV);

	ia1[1] = 1; ja1[1] = 1; ar1[1] = 1;
	ia1[2] = 1; ja1[2] = 2; ar1[2] = 0;
	ia1[3] = 1; ja1[3] = 3; ar1[3] = 1;
	ia1[4] = 2; ja1[4] = 1; ar1[4] = 0;
	ia1[5] = 2; ja1[5] = 2; ar1[5] = 1;
	ia1[6] = 2; ja1[6] = 3; ar1[6] = 1;

	glp_load_matrix(lp1, 6, ia1, ja1, ar1);
	
	int ind1[7];
	double val1[7];
int j=1;
for(j=1;j<=3;j++)
{
ind1[j]=-1;
val1[j]=-1.;
}

	int len1=glp_get_mat_row(lp1, 2, ind1, val1);

printf(" TEST SET 2 ");
for(j=1;j<=3;j++)
if(ind1[j]>=0 || val1[j]>=0)
printf(" (%d,%f) ",ind1[j],val1[j]);
printf("\n");
	
	int res_solver1=glp_simplex(lp1, NULL);
	double z1 = glp_get_obj_val(lp1);
	
	/* Gather results */
	double* soln1=(double*)calloc(5,sizeof(double));
	for(j=1;j<=3;j++)
		soln1[j]= glp_get_col_prim(lp1, j);
	
printf("RES : ");
	for(j=1;j<=3;j++)
		if(soln1[j]>0)
			printf(" %d %f || ",j,soln1[j]);

printf("\n z== %f\n",z1);

	glp_delete_prob(lp1);
	
	//--------TEST END

return 0;

}
