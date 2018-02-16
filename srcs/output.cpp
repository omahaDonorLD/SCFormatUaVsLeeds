#include "../heads/output.hpp"

using namespace std;



/*!
 *	\brief Fonction initialisant les futurs "path" d'output.
 *
 *	Fonction permettant d'initialiser les différents "path" d'output, de les mettre à jour et de supprimer tout fichier existant dans les dossier de même type dans les dossier.
 *
 *	\param	in_generation_path 			chemin des générations à écrire
 *	\param	in_pareto_path 				chemin des rang1 à écrire
 *	\param	in_plot_path 				chemin des images à créer
 *	\param	in_optimal_solution_path 	chemin des solutions optimales si existantes
 *	\param 	nbrObj 						nombre d'objectif des solutions à écrire
 *
 */
Output::Output(string in_generation_path, string in_pareto_path, string in_plot_path, string in_optimal_solution_path, int nbrObj)
{
	generation_path = in_generation_path;
	pareto_path = in_pareto_path;
	plot_path = in_plot_path;
	optimal_solution_path = in_optimal_solution_path;
	nbr_obj = nbrObj ;


	string to_remove;
	int status=0;
	int i=0;

	string s_i;

	/* supression des anciens fichier de donnée de génération */
	while(status != -1)
	{
		ostringstream convert ;
		convert << i;
		s_i = convert.str();
		to_remove = generation_path+s_i+".dat";
		status=remove(to_remove.c_str());
		i++;
	}
	status=0;
	i=0;
	/* suprressions des anciens fichier de donnée de rang 1 */
	while(status != -1)
	{
		ostringstream convert ;
		convert << i;
		s_i = convert.str();
		to_remove = pareto_path+s_i+".dat";
		status=remove(to_remove.c_str());
		i++;
	}
	status=0;
	i=0;
	/* suppressions des anciennes images */
	while(status != -1)
	{
		ostringstream convert ;
		convert << i;
		s_i = convert.str();
		to_remove = plot_path+s_i+".png";
		status=remove(to_remove.c_str());
		i++;
	}
	/* suppressions de données supplémentaires */
}



Output::~Output()
{

}



/*!
 *	\brief Permet d'écrire les valeurs d'une population à ploter plus tard.
 *
 *	Permet à partir d'une population d'individu d'écrire un fichier .txt/.dat ( ou autre ) qui contiendra les valeurs d'objectifs de chaque individu permettant un plotage ultérieur.
 *
 *	\param population	Population d'individu
 *	\param output_path	Chemin d'accès où écrire le fichier
 *
 */
void Output::writeData(vector<Individu*> population , string output_path)
{
	ofstream out_data(output_path.c_str(), ios::out | ios::trunc);

	if(out_data)
	{
		/* écriture des valeurs des objectifs de chaque individu de la population */
		for(unsigned int sol = 0; sol < population.size(); ++sol)
		{
			/* écriture des valeurs des objectifs de l'individu en cours */
			for( int obj = 0; obj < nbr_obj; ++obj)
			out_data<<population[sol]->getObjectiveValue(obj)<<"\t";

			out_data<<endl;
		}
	}
	else
	{
		cerr<<"Erreur ouverture fichier"<<endl;
	}
}





/*!
 *	\brief Permet l'écriture de la population finale ( après épuration ) dans un fichier.
 *
 *	Fonction permettant l'écriture des individus de la solution finale dans un fichier pour un plotage éventuelle.
 *
 *	\param mySol	Population d'individu finale après épuration de celle-ci.
 *
 */
void Output::writeSol(vector<Individu*> mySol )
{
	string output_path = "results/mySol.dat" ;
	ofstream out_data(output_path.c_str(), ios::out | ios::trunc);

	if(out_data)
	{
		/* écriture des valeurs d'objectifs de chaque individu */
		for(unsigned int sol = 0; sol < mySol.size(); ++sol)
		{
			/* écriture des valeurs d'objectifs de l'individu en cours */
			for( int obj = 0 ; obj < nbr_obj ; ++obj )
			out_data<<mySol[sol]->getObjectiveValue(obj)<<"\t";

			out_data<<endl;
		}
	}
	else
	{
		cerr<<"Erreur ouverture fichier"<<endl;
	}
}



/*!
 *	\brief Permet l'écriture dans un fichier généralisé d'une génération d'individu.
 *
 *	Fonction permettant d'écrire de façon générique chaque génération d'individu dans un fichier dont le nom comportera l'indice de la génération.
 *
 *	\param	population		Population d'individu dont l'on écrit les valeurs d'objectif
 *	\param generation_num	indice de la génération en cours.
 *
 */
void Output::writeGeneration(vector<Individu*> population, int generation_num)
{
	cout<<"smthihn"<<generation_num<<" and "<< generation_path<<endl;
	/* converstion de la génération actuelle en string */
	ostringstream convert ;
    convert << generation_num ;
    string num = convert.str() ;
	string path = generation_path+num+".dat" ;

	/* écriture dans le fichier correspondant */
	writeData(population, path);
}



/*!
 *	\brief Permet l'écriture dans un fichier généralisé du rang1 d'une génération d'individu.
 *
 *	Fonction permettant d'écrire de façon générique chaque génération d'individu de rang1 dans un fichier dont le nom comportera l'indice de la génération.
 *
 *	\param rank1			Population de rang1 dont l'on écrit les valeurs d'objectif
 *	\param generation_num	indice de la génération en cours
 *
 */
void Output::writeRank1(vector<Individu*> rank1, int generation_num)
{
	/* convertion de la génération actuelle en string */
	ostringstream convert ;
    convert << generation_num ;
    string num = convert.str() ;
	string path = pareto_path+num+".dat";

	/* écriture dans le fichier correspondant */
	writeData(rank1, path);
}



/*!
 *	\brief Permet la création d'image correspondant au fichier de données.
 *
 *	Fonction permettant de réguler la création d'image correspondant aux données de générations, de rang1, d'optimum, de résultat épuré et de rang moyen à partir des fichiers correspondants.
 *
 *	\param nbr_gen	Nombre de génération ayant servi jusque là et ainsi nombre de fichier à ploter.
 *
 */
void Output::plotData(int nbr_gen)
{
	int objectif1, objectif2;
	int borders[4];
	ifstream read_pareto(optimal_solution_path.c_str(), ios::in);

	if(read_pareto)
	{
		read_pareto>>objectif1>>objectif2;
		borders[0]=objectif1;
		borders[1]=objectif2;
		borders[2]=objectif1;
		borders[3]=objectif2;

		/* Recherche des bords extrèmes des solutions exactes */
		while( !read_pareto.eof() )
		{
			read_pareto>>objectif1>>objectif2;

			if(objectif1<borders[0])
				borders[0]=objectif1;
			else if(objectif1>borders[2])
				borders[2]=objectif1;

			if(objectif2<borders[1])
				borders[1]=objectif2;
			else if(objectif2>borders[3])
				borders[3]=objectif2;
		}
		read_pareto.close();
	}
	else
	{
		cout<<"Pareto exact non trouvé, le rang 1 de la dernière génération sera utilisé comme objectif final."<<endl;
		ostringstream convert ;
		convert << (nbr_gen-1) ;
		string j = convert.str();
		optimal_solution_path = pareto_path + j +".dat";

		ifstream read_rank1(optimal_solution_path.c_str(), ios::in);

		read_rank1>>objectif1>>objectif2;
		borders[0]=objectif1;
		borders[1]=objectif2;
		borders[2]=objectif1;
		borders[3]=objectif2;

		/* Recherches des bords extrêmes des solutions de rang1 si exacte non trouvée */
		while( !read_rank1.eof() )
		{
			read_rank1>>objectif1>>objectif2;

			if(objectif1<borders[0])
				borders[0]=objectif1;
			else if(objectif1>borders[2])
				borders[2]=objectif1;

			if(objectif2<borders[1])
				borders[1]=objectif2;
			else if(objectif2>borders[3])
				borders[3]=objectif2;
		}
		read_rank1.close();
	}

	FILE *gnuplotpipe;
	gnuplotpipe=popen("gnuplot -persist","w");



	if (!gnuplotpipe)
	{
		cerr<< "Gnuplot not found !"<<endl;
	}
	else
	{
		for (int gen = 0; gen < nbr_gen ; ++gen)
		{
			ostringstream convert ;
		    convert << gen ;
		    string j = convert.str() ;

		    plotSolutions(gnuplotpipe, generation_path+j+".dat", pareto_path+j+".dat", optimal_solution_path, plot_path+j+".png", borders);

		}

		fprintf(gnuplotpipe,"exit\n");
		pclose(gnuplotpipe);
	}
}



/*!
 *	\brief Crée une image d'une population d'individu
 *
 *	Fonction permettant de créer une image comportant les optimums si existant ainsi que la population d'individu tout en mettant en avant les rangs 1.
 *
 *	\param gnuplotpipe		pipe gnuplot permettant la création de l'image
 *	\param datafilename		nom du fichier servant à la création des points de la population ( en rouge )
 *	\param currparetofile	nom du fichier servant à la création des points de rang 1 de la population ( en bleu )
 *	\param optparetofile	nom du fichier servant à la création des solutions exactes ou rang 1 final si non existante ( en vert )
 *	\param imagefilename	nom du fichier image que l'on va écrire
 *	\param borders			bord géométique de l'image
 *
 */
void Output::plotSolutions(FILE *gnuplotpipe, string datafilename, string currparetofile, string optparetofile, string imagefilename, int borders[])
{
	fprintf(gnuplotpipe, "set terminal png size 800,600 enhanced font \"Helvetica,10\"\n");
	fprintf(gnuplotpipe, "set output '%s'\n", imagefilename.c_str());


	fprintf(gnuplotpipe, " set xrange [%d : %d] \n", borders[0]-500 ,  borders[2]+500);
  	fprintf(gnuplotpipe, " set yrange [%d : %d] \n",  borders[1]-500,  borders[3]+500);

	fprintf(gnuplotpipe,
			"plot \"%s\" using 1:2 notitle w points ps 1, \"%s\" using 1:2 notitle w points ps 1, \"%s\" using 1:2 notitle w points ps 1 \n" ,
	 		datafilename.c_str(),
	 		optparetofile.c_str(),
	 		currparetofile.c_str());

	fflush(gnuplotpipe);
}



/*!
 *	\brief Crée une image représentant l'évolution des rangs moyens au fil des générations.
 *
 *	\param gnuplotpipe 		pipe gnuplot permettant la création de l'image
 *	\param rankmax			nom du fichier servant à la création des points de rang max taille N ( rouge )
 *	\param totalRank		nom du fichier servant à la création des points de rang max taille 2N ( vert )
 *	\param rankAve			nom du fichier servant à la création des points de rang moyen N et 2N ( bleu )
 *	\param					nom du fichier image que l'on va écrire
 *	\param					bord géométrique de l'image
 *
 */
void Output::plotRank(FILE *gnuplotpipe, string rankmax, string totalRank, string rankAve, string imagefilename, int borders[])
{
	fprintf(gnuplotpipe, "set terminal png size 800,600 enhanced font \"Helvetica,10\"\n");
	fprintf(gnuplotpipe, "set output '%s'\n", imagefilename.c_str());


	fprintf(gnuplotpipe, " set xrange [%d : %d] \n", borders[1],  borders[3]);
  	fprintf(gnuplotpipe, " set yrange [%d : %d] \n",  borders[0],  borders[2]);

	fprintf(gnuplotpipe,
			"plot \"%s\" using 1:2 notitle w points ps 1, \"%s\" using 1:2 notitle w points ps 1, \"%s\" using 1:2 notitle w points ps 1 \n" ,
	 		rankmax.c_str(),
	 		totalRank.c_str(),
	 		rankAve.c_str() );

	fflush(gnuplotpipe);
}



/*!
 *	\brief Crée une image représentant la solution finale après épuration
 *
 *	\param gnuplotpipe		pipe gnuplot permettant la création de l'image
 *	\param mysol			nom du fichier servant à la création des points de la solution épurée ( rouge )
 *	\param currparetofile	nom du fichier servant à la création des points de rang 1 afin de les mettres en avant ( vert )
 *	\param imagefilename	nom du fichier image que l'on va écrire
 *	\param borders			bord géométrique de l'image
 *
 *	\return
 *
 */
void Output::plotSol(FILE *gnuplotpipe, string mySol, string currparetofile, string imagefilename, int borders[])
{
	fprintf(gnuplotpipe, "set terminal png size 800,600 enhanced font \"Helvetica,10\"\n");
	fprintf(gnuplotpipe, "set output '%s'\n", imagefilename.c_str());


	fprintf(gnuplotpipe, " set xrange [%d : %d] \n", borders[0]-500 ,  borders[2]+500);
  	fprintf(gnuplotpipe, " set yrange [%d : %d] \n",  borders[1]-500,  borders[3]+500);

	fprintf(gnuplotpipe,
			"plot \"%s\" using 1:2 notitle w points ps 1, \"%s\" using 1:2 notitle w points ps 1 \n" ,
	 		mySol.c_str(),
	 		currparetofile.c_str() );

	fflush(gnuplotpipe);
}
