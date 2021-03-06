
#include <stdio.h>
#include <stdlib.h>
#include <glpk.h>

int main(void)
{

	//--=--------TEST START

		/* Build it first */
	glp_prob *lp1;
	int ia1[25],ja1[25];
	double ar1[25];
	
	lp1 = glp_create_prob();
	glp_set_prob_name(lp1, "test set cover");
	glp_set_obj_dir(lp1, GLP_MIN);

	glp_add_rows(lp1, 6);

	glp_set_row_name(lp1, 1, "r1");
	glp_set_row_bnds(lp1, 1, GLP_LO, 2.0, 0.0);
	glp_set_row_name(lp1, 2, "r2");
	glp_set_row_bnds(lp1, 2, GLP_LO, 2.0, 0.0);
	glp_set_row_name(lp1, 3, "r3");
	glp_set_row_bnds(lp1, 3, GLP_LO, 2.0, 0.0);
	glp_set_row_name(lp1, 4, "r4");
	glp_set_row_bnds(lp1, 4, GLP_LO, 2.0, 0.0);
	glp_set_row_name(lp1, 5, "r5");
	glp_set_row_bnds(lp1, 5, GLP_LO, 2.0, 0.0);
	glp_set_row_name(lp1, 6, "r6");
	glp_set_row_bnds(lp1, 6, GLP_LO, 2.0, 0.0);

	glp_add_cols(lp1, 4);

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

	glp_set_col_name(lp1, 4, "x4");
	glp_set_col_bnds(lp1, 4, GLP_DB, 0.0, 1.0);
	glp_set_obj_coef(lp1, 4, 1.0);
	//glp_set_col_kind(lp, j, GLP_BV);
	//glp_set_col_kind(lp, j, GLP_BV);

	ia1[1] = 1; ja1[1] = 1; ar1[1] = 1;
	ia1[2] = 1; ja1[2] = 2; ar1[2] = 1;
	ia1[3] = 1; ja1[3] = 3; ar1[3] = 1;
	ia1[4] = 1; ja1[4] = 4; ar1[4] = 0;
	ia1[5] = 2; ja1[5] = 1; ar1[5] = 1;
	ia1[6] = 2; ja1[6] = 2; ar1[6] = 1;
	ia1[7] = 2; ja1[7] = 3; ar1[7] = 1;
	ia1[8] = 2; ja1[8] = 4; ar1[8] = 0;
	ia1[9] = 3; ja1[9] = 1; ar1[9] = 1;
	ia1[10] = 3; ja1[10] = 2; ar1[10] = 1;
	ia1[11] = 3; ja1[11] = 3; ar1[11] = 1;
	ia1[12] = 3; ja1[12] = 4; ar1[12] = 0;
	ia1[13] = 4; ja1[13] = 1; ar1[13] = 1;
	ia1[14] = 4; ja1[14] = 2; ar1[14] = 1;
	ia1[15] = 4; ja1[15] = 3; ar1[15] = 0;
	ia1[16] = 4; ja1[16] = 4; ar1[16] = 1;
	ia1[17] = 5; ja1[17] = 1; ar1[17] = 1;
	ia1[18] = 5; ja1[18] = 2; ar1[18] = 1;
	ia1[19] = 5; ja1[19] = 3; ar1[19] = 0;
	ia1[20] = 5; ja1[20] = 4; ar1[20] = 1;
	ia1[21] = 6; ja1[21] = 1; ar1[21] = 1;
	ia1[22] = 6; ja1[22] = 2; ar1[22] = 1;
	ia1[23] = 6; ja1[23] = 3; ar1[23] = 0;
	ia1[24] = 6; ja1[24] = 4; ar1[24] = 1;

	glp_load_matrix(lp1, 24, ia1, ja1, ar1);
	
	int ind1[5];
	double val1[5];
int j=1;
for(j=1;j<=4;j++)
{
ind1[j]=-1;
val1[j]=-1.;
}

	int len1=glp_get_mat_row(lp1, 3, ind1, val1);

printf(" TEST SET 2 ");
for(j=1;j<=4;j++)
if(ind1[j]>=0 || val1[j]>=0)
printf(" (%d,%f) ",ind1[j],val1[j]);
printf("\n");
	
	int res_solver1=glp_simplex(lp1, NULL);
	double z1 = glp_get_obj_val(lp1);
	
	/* Gather results */
	double* soln1=(double*)calloc(5,sizeof(double));
	for(j=1;j<=4;j++)
		soln1[j]= glp_get_col_prim(lp1, j);
	
printf("RES : ");
	for(j=1;j<=4;j++)
		if(soln1[j]>0)
			printf(" %d %f || ",j,soln1[j]);

printf("\n z== %f\n",z1);

	glp_delete_prob(lp1);
	
	//--------TEST END

return 0;

}
