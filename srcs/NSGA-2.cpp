#include <vector>
#include <climits>
#include <iostream>
#include <stdlib.h>
#include <math.h>

#include <string>
#include <ctime>
#include <sstream>
#include <fstream>

#include <stdio.h>


#include "../heads/individuFactory.hpp"
#include "../heads/individu.hpp"
#include "../heads/output.hpp"

#define ERROR(x)  fprintf(stderr, x), fprintf(stderr, "\n"), exit(1)

using namespace std;



int Gener = 0 ;	 			// Number of Iteration
int n_Clone = 0 ;			// Number of clone removed in one generation.
int n_Objectives = 0 ;		// Number of objectives
string param_path = "paramsFiles/param.txt" ; 						// chemin d'accès au fichier paramètre par défaut



/* Can be read in a param file */
bool TESTMODE = 0 ;			// 0 If display infos, 1 otherwise ( still testing ).
int PopSize = 100 ; 	 	// Population Size not used here because it's time dependent
double split = 0.5 ; 		// Correspond au pourcentage de population graps => PopSize*k ici 50% si supérieur ou égale à 1 valeur brute d'individu par direction
int n_Delta = 3 ; 			// Correspond au nombre de Delta voulu par paire d'objectif
double alpha = 0.7 ; 		// Correspond à l'alpha nécessaire à la RCList du GRASP.
double mutaprob = 0.01 ; 	// Correspond à la chance de muter
double rank1maxsize = 0.5 ;	// Correspond au pourcentage de la population que peut représenter le rang 1, augmentation de la taille de la population si supérieur.
double crowdSave = 0.1 ;	// Correspond au pourcentage de la population sauver par la seconde valeur de crowding, si supérieur à 1 valeur brute d'individu sauvé
double crowdTotal = 0.5 ;	// Correspond au pourcentage de rang que l'on considère pour le calcul du crowding2, ici 50%.
int HyperVolCD = 10 ;
double PercentEvolution = 0.1 ;
double GRASPercent = 0.5 ;
int MaxGener = 1000 ;
int MaxTime = 1000 ;
string inst_path = "../instances/instanceTest250.txt" ; 	// chemin d'accès au fichier problème par défaut
string opt_path = "../optimum/optimum250.txt" ; 			// chemin d'accès à la réponse du problème par défaut
bool removeclone = 1 ;		// 0 si on garde les clones, 1 sinon.
bool parentarena = 1 ;		// 0 si on sélectionne les parents totalement aléatoirement, 1 sinon.
bool gnuplot = 1 ;			// 0 si on supprime l'affiche, 1 sinon.
bool results = 1 ;			// 0 si l'on ne veut pas les résultats, 1 sinon.


/* Variables utilisées pour les calcules d'hypervolume */
int size_Front2 ;
int nonGRASP_R1 = 0 ;

bool stillgoing = true ;
bool startHyperVol = false ;
vector< Individu* > lastRank1 ;



/*!
 *	\brief
 *
 *	Commentaire
 *
 *	\param
 *
 *	\return
 *
 */



/*

Déclaration de toutes les fonctions

*/
void generate(IndividuFactory usine, vector<Individu*> &myResult);
void addSol(int addsol, vector<Individu*> &myPop, vector< vector<Individu*> > &myRank );

void quicksortcrowd(vector<Individu*> &myCrowd, int begin, int end);
void quicksortn_Objectives(int n, vector<Individu*> &myRank, int begin, int end);
void crowding(vector<Individu*> &myRank, double inMaxMin[][2] );
void quicksortcrowdTotal(vector<Individu*> &myCrowd, int begin, int end);
void crowdingTotal(vector<Individu*> &myCrowd, double inMaxMin[][2] );
void saveTotalCrowd(	vector<Individu*> &myPop,
						vector< vector<Individu*> >&toReturn,
						double fMaxMin[][2]  );

void reCreate( 	vector< vector<Individu*> > &toReturn,
				vector<Individu*> &myPop,
				double fMaxMin[][2] );

int removeClone( vector< vector<Individu*> > &toReturn );

vector< vector<Individu*> > sortRank( vector<Individu*> &myPop );

Individu* selectParentRND( vector<Individu*> myParent );
bool readparam(string param_path) ;
void quality(string opt_path, vector< vector<Individu*> > &ranKing) ;

void extractSol( vector<Individu*> &myPop, vector<Individu*> &ranKing, Output my_ouput );

bool hyperVolumeDist(vector< Individu* > Rank1, int nbr_objectives);
double spread( vector< Individu* > result );



//methodes de calcul d'hypervolume de Zitzler
int  Dominates(double  point1[], double  point2[], int  noObjectives);
void  Swap(double  *front[], int  i, int  j);
int  FilterNondominatedSet(double  *front[], int  noPoints, int  noObjectives);
double  SurfaceUnchangedTo(double  *front[], int  noPoints, int  objective);
int  ReduceNondominatedSet(double  *front[], int  noPoints, int  objective, double  threshold);
double  CalculateHypervolume(double  *front[], int  noPoints,int  noObjectives);
int  ReadFront(double  **frontPtr[], FILE  *file, int  noObjectives);
int  MergeFronts(double  **frontPtr[], double  *front1[], int  sizeFront1, double*  front2[], int  sizeFront2, int  noObjectives);
void  DeallocateFront(double**  front, int  noPoints);
int  hypervolume(string pareto1, string pareto2, int nbr_objectives);

/*

Fin de la déclaration des fonctions

*/




/*!
 *	\brief Moteur de Recherche NSGA II
 *
 *	Permet de créer une population de taille interactive ( suivant l'utilisateur ) et de trouver le front pareto optimal après le nombre d'itération ou le temps choisi.
 *
 *	\param argc nombre de paramètre ajouté au lancement de NSGA2
 *	\param argv paramètre ajouté à NSGA2
 *
 *	\return 0 si tout s'est bien déroulé
 *
 */
int main(int argc,char** argv)
{
	time_t tbegin ; // Variable de début de temps, au lancement de l'algorithme, après que l'utilisateur ait entré les données nécessaires.
	time_t tend ; // Variable de fin de temps, mise à jour après chaque génération permettant l'arrêt de celles-ci

	bool test = false ;


	/* Reading Param File */
	if( argc > 1 )
	{
		param_path = argv[1] ;
		test = readparam( param_path ) ;

		if( !test )
		{
			return 0 ;
		}
	}
	else
	{
		test = readparam( "paramsFiles/param.txt" ) ;

		if( !test )
		{
			return 0 ;
		}
	}



	/* Initialization */
	srand(time(NULL)) ; // mise à jour de la seed random
	IndividuFactory usine(inst_path) ; // initialisation du problème
	n_Objectives = Individu::getNbrObjectives() ;

	vector<Individu*> myPopulation ;

	tbegin = time(NULL); // Commencement du temps de NSGA2

	Output my_ouput("data/solution",
					"pareto/pareto",
					"image/image",
					opt_path, n_Objectives); // écriture des données initiales par gnuplot



	/* Begin Generation of individuals */
	generate(usine, myPopulation ) ; // générations d'une population initiale

	cout << " End Generate " << endl ;

	tend = time(NULL) ;

	cout << " Population created in " << tend-tbegin << " seconds "  << endl ;
	cout << myPopulation.size() << endl ;


	/* First Ranking */
	vector< vector<Individu*> > ranKing ;

	ranKing = sortRank( myPopulation ) ; // ranking et crowding de la population initiale.

	/* Initialisation pour l'hypervol */
	for(unsigned int i = 0 ; i < ranKing[0].size() ; ++i )
		lastRank1.push_back( ranKing[0][i] );

	size_Front2 = lastRank1.size() ;



	/* Writing DATA */
	my_ouput.writeGeneration(myPopulation, 0); // écriture des données correspondantes par gnuplot sur la population
	my_ouput.writeRank1(ranKing[0], 0); // écriture des données correspondantes par gnuplot sur le front pareto optimal de notre population



	cout << " Begin Generation " << endl ;

	Gener = 1 ;	// Permet de savoir à quelle génération nous sommes rendus.
	tend = time(NULL) ; // initialisation du temps de fin, au cas où la génération, le ranking et crowding de la population initiale dépasse déjà le temps souhaité



	/* EVOLUTION PART */
	while( stillgoing ) // Tant que le temps passe est inférieur au temps voulu, ou que le nombre de generations est inférieur aux generations voulues on fait des générations
	{
		if( !TESTMODE )
			cout << "Generation : " << Gener << endl ;



		/* Création des enfants */
		for( int j = 0; j < PopSize; ++j) // Création de PopSize enfants.
		{

			Individu *dad ;
			Individu *mum ;

			if( parentarena )
			{
				dad = selectParentRND(myPopulation) ; // Sélection Random + Championnat de deux parents. Le meilleur en ressort vainqueur.
				mum = selectParentRND(myPopulation) ; // idem pour la sélection du second parent.
			}
			else
			{
				int rnd1 = (int) floor(rand() % PopSize) ;
				int rnd2 = (int) floor(rand() % PopSize) ;

				dad = myPopulation[rnd1] ;
				mum = myPopulation[rnd2] ;
			}

			Individu *child ;
			child = Individu::getChildFrom(dad,mum) ;

			double rnd3 = (int) rand() % 1000000 ;
			rnd3 = rnd3 / 1000000.0 ;

			if( rnd3 < mutaprob )
			{
				child->mutate() ;
			}

			myPopulation.push_back( child ); // Création de deux enfants à partir de ces deux parents, le meilleur en ressort vainqueur

		}



		/* RANKING NEW POPULATION */
		//cout << " myPopulation size = " << myPopulation.size() << endl ;

		ranKing.clear() ;

		//cout << " BEGIN RANKING " << endl ;

		ranKing = sortRank(myPopulation) ; // Mise à jour du ranking et crowding de la nouvelle population et reconstruction de celle-ci avec les PopSize meilleurs éléments.

		//cout << " END RANKING " << endl ;


		if( !TESTMODE )
		{
        	cout << " Number of solution in pareto optimal set = " << ranKing[0].size() << endl ;
        	cout << " Population Size : " << PopSize << endl ;
    	}



		/* WRITING DATA */
		if( gnuplot )
		{
      my_ouput.writeGeneration(myPopulation, Gener); // écriture par gnuplot des données de la génération i
			my_ouput.writeRank1(ranKing[0], Gener); // écriture par gnuplot des données du front pareto i
		}


		startHyperVol = startHyperVol || ( nonGRASP_R1 > ranKing[0].size()*(1-GRASPercent) ) ;


		/* Vérification de l'évolution de l'hypervolume */
		if( startHyperVol && Gener%HyperVolCD == 0 )
			stillgoing = hyperVolumeDist(ranKing[0], n_Objectives) ;



		// cout << Gener << " " << GenerCrowdTotal << " " << crowdTotal << endl ;

		tend = time(NULL) ; // mise à jour du temps de fin

		stillgoing = stillgoing && Gener < MaxGener && tend-tbegin < MaxTime ;

		++Gener ; // incrémentation du compteur de génération


	}

	// cout << " End " << endl ;

	cout << " Around "<< tend-tbegin <<" seconds have passed " << endl ;
	cout << " Number of solution in pareto optimal set = " << ranKing[0].size() << endl ;
	cout << " Number of gener : " << Gener << endl ;

    if( !gnuplot )
    {
    	Gener = 1 ;
    	my_ouput.writeGeneration(myPopulation, Gener); // écriture de la dernière population si gnuplot inactif
		my_ouput.writeRank1(ranKing[0], Gener); // écriture par gnuplot des données du front pareto de la dernière population si gnuplot inactif
		Gener++ ;
    }

    if( n_Objectives < 3 )
    {
    	my_ouput.plotData(Gener); // Création des images correspondant aux populations de chaque générations
   	}


   	// Permet d'écrire les résultats obtenu dans un fichier
   	if( results )
   	{
   		time_t t;
   		time(&t);

   		string stringtime = ctime(&t) ;
   		string resultpath = "results/"+stringtime+".dat" ;

   		ofstream out_data(resultpath.c_str(), ios::out | ios::trunc);

			if(out_data)
			{
				for(unsigned int sol = 0; sol < myPopulation.size(); ++sol)
				{
					out_data<< myPopulation[sol]->toString() << endl ;
				}
			}
			else
			{
				cerr<<"Erreur ouverture fichier function "<< __FUNCTION__ <<endl;
			}
   	}


   	quality(opt_path, ranKing) ; // Calcul de la qualité "euclidienne" de notre solution

   	double SpreadValue = spread(myPopulation);

   	cout << " Spread Quality " << SpreadValue << endl ;

    cout << " Good Bye World " << endl ;

	return 0;
}




/*!
 *  \brief Générateur d'individu.
 *
 *	Méthode servant à générer des invidus de façon aléatoire avec éventuellement une portion d'individu élite.
 *
 *	\return La liste d'individu ainsi créée.
 */
void generate(IndividuFactory usine, vector<Individu*> &myResult )
{

	bool clone ;
	int sum = 0 ;
	int n_Grasp = 0 ;

	vector<int*> Coef ;

	Individu* temp ;

	/*

	Création des directions GRASP

	*/
	if( n_Delta > 2 )
	{

		int delta[n_Delta-2] ; // Coefficient que l'on trouvera entre chaque pair d'objectif

		for( int i = 0 ; i < n_Delta-2 ; ++i ) // Calcul des Coefficients que l'on trouvera entre chaque pair d'objectif de 0 à 100 par saut de 100/n_Delta
		{
			delta[i] = (i+1)*100/(n_Delta-1)  ;
		}

		/*
		Ci après, Création des coefficients entre chaque objectifs
		*/

		sum = n_Objectives ;
		int verif ;

		for( int i = 1 ; i < n_Objectives ; ++i )
			sum += (n_Delta-2)*i ;

		if( n_Objectives > 2 )
		{
			verif = sum ;
			sum++ ;
		}
		else
		{
			verif = sum ;
		}

		int counti = 0 ;
		int countj = 1 ;
		int countk = 0 ;

		for( int i = 0 ; i < n_Objectives ; ++i )
		{
			Coef.push_back( new int[n_Objectives] );

			for( int j = 0 ; j < n_Objectives ;++j )
				Coef[i][j] = 0 ;

			Coef[i][i] = 100 ;
		}

		for( int i = n_Objectives ; i < verif ; ++i )
		{
			Coef.push_back( new int[n_Objectives] );

			for( int j = 0 ; j < n_Objectives ;++j )
				Coef[i][j] = 0 ;

			Coef[i][counti] = delta[countk] ;
			Coef[i][countj] = 100 - delta[countk];

			countk = (countk + 1) % (n_Delta-2) ;

			if( countk == 0 )
			{
				countj = (countj + 1) % n_Objectives ;
			}
			if( countj == 0 )
			{
				counti ++ ;
				countj = counti+1 ;
			}

		}

		if( n_Objectives > 2 ) // Création d'une direction "centrale" si le nombre d'objectif est supérieur à 2
		{
			Coef.push_back( new int[n_Objectives] ) ;

			for(int i = 0 ; i < n_Objectives ; ++i )
				Coef[sum-1][i] = 100/n_Objectives ;
		}

		cout << " Coef Size : " << Coef.size() << endl ;

	}
	else if( n_Delta == 2 )
	{
		sum = n_Objectives ;

		if( n_Objectives > 2 )
			sum++ ;

		for( int i = 0 ; i < n_Objectives ; ++i )
		{
			Coef.push_back( new int[n_Objectives] );

			for( int j = 0 ; j < n_Objectives ;++j )
				Coef[i][j] = 0 ;

			Coef[i][i] = 100 ;

		}

		if( n_Objectives > 2 ) // Création d'une direction "centrale" si le nombre d'objectif est supérieur à 2
		{
			Coef.push_back( new int[n_Objectives] ) ;

			for(int i = 0 ; i < n_Objectives ; ++i )
				Coef[sum-1][i] = 100/n_Objectives ;
		}
	}
	else if( n_Delta < 2 )
	{
		int verif ;

		for( int i = 0 ; i < n_Objectives ; ++i )
			sum += i ;

		if( n_Objectives > 2 )
		{
			verif = sum ;
			sum++ ;
		}
		else
		{
			verif = sum ;
		}

		int counti = 0 ;
		int countj = 1 ;

		for( int i = 0 ; i < verif ; ++i )
		{
			Coef.push_back( new int[n_Objectives] );

			for( int j = 0 ; j < n_Objectives ;++j )
				Coef[i][j] = 0 ;

			Coef[i][counti] = 50 ;
			Coef[i][countj] = 50 ;

			countj = (countj + 1) % n_Objectives ;

			if( countj == 0 )
			{
				counti ++ ;
				countj = counti+1 ;
			}

		}

		if( n_Objectives > 2 ) // Création d'une direction "centrale" si le nombre d'objectif est supérieur à 2
		{
			Coef.push_back( new int[n_Objectives] ) ;

			for(int i = 0 ; i < n_Objectives ; ++i )
				Coef[sum-1][i] = 100/n_Objectives ;
		}

	}

	if( split >= 1 )
	{
		n_Grasp = split ;
	}
	else
	{
		n_Grasp = (int) PopSize*split / Coef.size() ;
	}



	///////////////////////////////////////////////////////////////////////////////////////////////////////



	/*

	Création des individus GRASP

	*/
	for( int p = 0 ; p < sum ; ++p )
	{

		Individu::createGRASPFIT( Coef[p], alpha ) ;

		for( int i = 0 ; i < n_Grasp ; ++i )
		{
			clone = true ;
			while( clone ) // Vérification de la création d'un clône. Si clone, on re-essaye. Les clônes sont actuellement rare sur les problèmes de sac à dos multi objectifs
			{

				temp = usine.createGraspSolution(alpha) ; // création d'une solution GRASP

				clone = false ;

				for(unsigned int j = 0; j < myResult.size() && !clone ; ++j ) // Vérification de si notre solution est un clône ou non.
				{
					clone = myResult[j]->isCloneOf(temp) ;
				}

				if( !clone ) // Si ce n'est pas un clône, rajout dans la population.
				{
					myResult.push_back( temp ) ;
				}
			}
		}

	}

	cout << "Size after Grasp : " << myResult.size() << endl ;


	///////////////////////////////////////////////////////////////////////////////////////////////////////


	// Create RND
	while( myResult.size() < PopSize ) // Création du reste de la population de façon random
	{

		Individu *temp = usine.createRandomSolution() ; // Création random d'un individu

		clone = false ;

		for(unsigned int i = 0; i < myResult.size() && !clone ; ++i ) // Vérification d'un éventuel clône
		{
			clone = myResult[i]->isCloneOf(temp) ;
		}

		if( !clone ) // Si non clône, rajout de l'individu dans la population.
		{
			myResult.push_back( temp ) ;
		}
	}
}



/*!
 *  \brief Générateur d'individu muté à partir d'un individu actuel.
 *
 *	Méthode servant à générer des invidus muté de façon aléatoire lorsque l'on en a besoin de façon ponctuelle afin de les rajouter à la population actuelle.
 *
 */
void addSol(int addsol, vector<Individu*> &myPop, vector< vector<Individu*> > &myRank )
{
	int rnd1;
	int rnd2;
	bool clone;

	for(int i = 0; i < addsol ; ++i )
	{
		clone = true ;

		while( clone )
		{
			rnd1 = rand() % myRank.size();
			rnd2 = rand() % myRank[rnd1].size();

			Individu *temp = new Individu(myRank[rnd1][rnd2]) ;

			temp->mutate() ;

			if( !temp->isCloneOf(myRank[rnd1][rnd2]) )
			{
				myPop.push_back(temp);
				clone = false ;
			}
		}
	}
}



/*!
 *	\brief Tri du crowd en ~nlog(n)
 *
 *	Tri d'un rang par ordre croissant de crowding en ~nlog(n)
 *
 *	\param myCrowd  population à trier
 *  \param begin    iterateur de commencement de tri
 *  \param end      iterateur de fin de tri
 *
 */
void quicksortcrowd(vector<Individu*> &myCrowd, int begin, int end)
{
	if( end - begin > 0 )
	{
		int pointeur = begin ;
		for(int i = begin; i < end ; ++i)
		{
			if( myCrowd[i]->crowding_value < myCrowd[end]->crowding_value )
			{
				Individu *temp = myCrowd[pointeur] ;
				myCrowd[pointeur] = myCrowd[i] ;
				myCrowd[i] = temp ;
				pointeur ++ ;
			}
		}

		Individu *temp = myCrowd[pointeur] ;
		myCrowd[pointeur] = myCrowd[end] ;
		myCrowd[end] = temp ;

		quicksortcrowd( myCrowd, begin, pointeur-1) ;
		quicksortcrowd( myCrowd, pointeur+1, end) ;
	}
}



/*!
 *	\brief Tri de la Fitness de l'objectif n en ~nlog(n)
 *
 *	Permet de trier une population d'individu selon un objectif précis, et la fitness associé.
 *
 *  \param n_Objectives      Objectif selon lequel trier
 *	\param myRank   Population à trier
 *  \param begin    iterateur de commencement de tri
 *  \param end      iterateur de fin de tri
 *
 */
void quicksortn_Objectives(int n, vector<Individu*> &myRank, int begin, int end)
{
	if( end - begin > 0 )
	{
		int pointeur = begin ;

		for(int i = begin; i < end ; ++i)
		{
			if( myRank[i]->getObjectiveValue(n) < myRank[end]->getObjectiveValue(n) )
			{
				Individu *temp = myRank[pointeur] ;
				myRank[pointeur] = myRank[i] ;
				myRank[i] = temp ;
				pointeur ++ ;
			}
		}

		Individu *temp = myRank[pointeur] ;
		myRank[pointeur] = myRank[end] ;
		myRank[end] = temp ;

		quicksortn_Objectives(n, myRank, begin, pointeur-1) ;
		quicksortn_Objectives(n, myRank, pointeur+1, end) ;
	}
}



/*!
 *	\brief Calcul le Crowding des individus d'une population
 *
 *	Fonction permettant de calculer le crowding de chaque individu d'une population d'un même rang et de le mettre à jour
 *
 *	\param myRank   Population d'individu d'un même rang à trier
 *  \param inMaxMin Tableau d'entier comportant les max et min de chaque objectif
 *
 */
void crowding(vector<Individu*> &myRank, double inMaxMin[][2] )
{

    for(unsigned int i = 0 ; i < myRank.size() ; ++i )
    {
        myRank[i]->crowding_value = 0 ; // Remise à zéro de tous les crowdings
    }

	for( int i = 0 ; i < n_Objectives ; ++i )
	{
		quicksortn_Objectives(i, myRank, 0, myRank.size()-1 ) ; // Tri de chaque individu en fonction de sa fitness dans l'objectif i

		myRank[0]->crowding_value = INT_MAX ; // Permet de garder les extrémités.
		myRank[myRank.size()-1]->crowding_value = INT_MAX ; // Permet de garder les extrémités.

		for(unsigned int j = 1; j < myRank.size()-1 ; ++j ) // Ajout de la valeur actuelle du crowding pour chaque individu
		{
			myRank[j]->crowding_value += ( myRank[j+1]->getObjectiveValue(i)-myRank[j-1]->getObjectiveValue(i) ) / ( inMaxMin[i][0] - inMaxMin[i][1] ) ;
		}

	}
}



/*!
 *	\brief Tri du crowd2 en ~nlog(n)
 *
 *	Tri d'un rang par ordre croissant de crowding2 en ~nlog(n)
 *
 *	\param myCrowd  population à trier
 *  \param begin    iterateur de commencement de tri
 *  \param end      iterateur de fin de tri
 *
 */
void quicksortcrowdTotal(vector<Individu*> &myCrowd, int begin, int end)
{
	if( end - begin > 0 )
	{
		int pointeur = begin ;
		for(int i = begin; i < end ; ++i)
		{
			if( myCrowd[i]->crowding_Total < myCrowd[end]->crowding_Total )
			{
				Individu *temp = myCrowd[pointeur] ;
				myCrowd[pointeur] = myCrowd[i] ;
				myCrowd[i] = temp ;
				pointeur ++ ;
			}
		}

		Individu *temp = myCrowd[pointeur] ;
		myCrowd[pointeur] = myCrowd[end] ;
		myCrowd[end] = temp ;

		quicksortcrowdTotal( myCrowd, begin, pointeur-1) ;
		quicksortcrowdTotal( myCrowd, pointeur+1, end) ;
	}
}




/*!
 *	\brief Calcul le Crowding2 des individus d'une population
 *
 *	Fonction permettant de calculer le crowding2 de chaque individu d'une population et de le mettre à jour
 *
 *	\param myCrowd   Population d'individu à trier
 *  \param inMaxMin Tableau d'entier comportant les max et min de chaque objectif
 *
 */
void crowdingTotal(vector<Individu*> &myCrowd, double inMaxMin[][2] )
{
	for(unsigned int i = 0 ; i < myCrowd.size() ; ++i )
    {
        myCrowd[i]->crowding_Total = 0 ; // Remise à zéro de tous les crowdings
    }

	for( int i = 0 ; i < n_Objectives ; ++i )
	{
		quicksortn_Objectives(i, myCrowd, 0, myCrowd.size()-1 ) ; // Tri de chaque individu en fonction de sa fitness dans l'objectif i

		myCrowd[0]->crowding_Total = INT_MAX ; // Permet de garder les extrémités.
		myCrowd[myCrowd.size()-1]->crowding_Total = INT_MAX ; // Permet de garder les extrémités.

		for(unsigned int j = 1; j < myCrowd.size()-1 ; ++j ) // Ajout de la valeur actuelle du crowding pour chaque individu
		{
			myCrowd[j]->crowding_Total += ( myCrowd[j+1]->getObjectiveValue(i)-myCrowd[j-1]->getObjectiveValue(i) ) / ( inMaxMin[i][0] - inMaxMin[i][1] ) ;
		}

	}

}



/*!
 *	\brief Fonction calculant et sauvegardant les individus avec un crowding2 intéressant
 *
 *	La fonction permet à l'aide de la population triée en rang de retrouver les individus intéressants afin de permettre l'uniformité du front de solution.
 *
 *	\param myPop population à laquelle sera rajoutée les individus désignés comme intéressant par le crowding2
 *	\param toReturn ranking actuelle de la population de taille 2N dans lequel on ira chercher les individus pour le calcul du crowding2
 *	\param fMaxMin Tableau d'entier comportant les max et min de chaque objectif
 *
 *	\return
 *
 */
void saveTotalCrowd(	vector<Individu*> &myPop,
						vector< vector<Individu*> >&toReturn,
						double fMaxMin[][2]  )
{
	if( crowdSave > 0 )
	{

		// Calcul du nombre de rang à traiter

		int rankave =  (int) toReturn.size()*crowdTotal ;

		if( rankave >= toReturn.size() )
			rankave = toReturn.size()-1 ;

		// Copie des rangs à traiter dans la population vide
		for( int i = 0 ; i <= rankave ; i++ )
		{
			for(unsigned int j = 0 ; j < toReturn[i].size() ; j++ )
			{
				myPop.push_back( toReturn[i][j] );
			}
		}

		int inRank ;


		// Calcul du Crowding2 et sauvegarde

		if( myPop.size() > 0 )
		{
			/* Calcul du crowding2 et trie en fonction de celui-ci */
			crowdingTotal( myPop, fMaxMin );


			/* Suppression des rangs 1 pour qu'il ne fausse pas la sauvegarde des autres rang1 */
			for(unsigned int i = 0 ; i < myPop.size() ; ++i )
			{
				if( myPop[i]->getRank() == 1 )
				{
					myPop.erase(myPop.begin()+i) ;
					i = -1 ;
				}

			}

			quicksortcrowdTotal( myPop, 0 , myPop.size()-1 );



			/* Si crowdSave est un pourcentage épuration de la population jusqu'à obtenir le pourcentage voulu */
			if( crowdSave <= 1 )
			{
				/* épuration de la population en commençant par les individus au crowding le plus faible */
				while( myPop.size() > crowdSave * PopSize )
				{
					myPop.erase(myPop.begin()) ;
				}

				/* Retrait de chaque individu des rangs correspondant pour éviter les clones */
				for(unsigned int i = 0 ; i < myPop.size() ; ++i )
				{
					inRank = myPop[i]->getRank()-1 ;

					/* vérification et suppression du clone correspondant */
					for(unsigned int j = 0 ; j < toReturn[inRank].size() ; ++j )
					{
						if( toReturn[inRank][j]->isCloneOf( myPop[i] ) )
						{
							toReturn[inRank].erase(toReturn[inRank].begin()+j ) ;
	 						j = INT_MAX ;
						}
					}
				}


			}
			else
			{
				/* épuration de la population en commençant par les individus au crowding le plus faible jusqu'à obtenir crowdSave individu */
				while( myPop.size() > crowdSave )
				{
					myPop.erase(myPop.begin()) ;
				}

				/* Retrait de chaque individu des rangs correspondant pour éviter les clones */
				for(unsigned int i = 0 ; i < myPop.size() ; ++i )
				{
					inRank = myPop[i]->getRank()-1 ;

					/* vérification et suppression du clone correspondant */
					for(unsigned int j = 0 ; j < toReturn[inRank].size() ; ++j )
					{
						if( toReturn[inRank][j]->isCloneOf( myPop[i] ) )
						{
							toReturn[inRank].erase(toReturn[inRank].begin()+j ) ;
							j = INT_MAX ;
						}
					}
				}
			}
		}
	}
}



/*!
 *	\brief Permet de Recréer une population de taille PopSize
 *
 *	Fonction permettant de recréer la population à partir des rangs précalculer. On ajoute automatiquement les rangs tout en calculant le crowding jusqu'à arriver au rang intermédiaire. Une fois arrivée à celui-ci. On calcul le crowding, on trie les individus suivant le crowding et on garde les plus intéressant
 *
 *  \param toReturn     ensemble des rangs comportant les individus
 *  \param myPop        population d'individu actuellement traitée
 *  \param fMaxMin      permet de calculer le crowding
 *
 */
void reCreate( 	vector< vector<Individu*> > &toReturn,
				vector<Individu*> &myPop,
				double fMaxMin[][2] )
{
	// ici on s'arrête avant le rang intermédiaire ou le nombre de "population" restante après la suppression des clones

    int inlength = toReturn[0].size() + myPop.size() ; // initialisation de la taille actuelle de la population, celle-ci pouvant être déjà composé d'élément créé d'urgence lorsque la taille double, augmente. Mais alors, chaque individu du ranking se trouvera forcément dedans.

    int ranking = 0 ;

    while( inlength < PopSize ) // && inlength < ((1+rank1maxsize)*PopSize-n_Clone) ) // Tant qu'on ne dépasse pas PopSize ou qu'on ne dépasse pas 2PopSize moins le nombre de clône ( si celui-ci est supérieur à PopSize ) on rajoute nos rangs. La population a déjà été nettoyé des clônes.
    {
        for(unsigned int i = 0; i < toReturn[ranking].size() ; ++i ) // Rajout du rang actuellement traité
        {
            myPop.push_back(toReturn[ranking][i]) ;
        }

		if(parentarena && toReturn[ranking].size() > 0 )
        	crowding( toReturn[ranking], fMaxMin ); // Calcul du crowding nécéssaire à la sélection de parent pour chaque rang.

        ranking ++ ;

        inlength += toReturn[ranking].size() ; // rajout de la taille du rang prochain, pour vérifier s'il dépasse ou non PopSize avant de le rajouter.
    }

	//cout << " ranking " << ranking << " / " << toReturn.size() << endl ;
	//cout << " POPSIZE " << myPop.size() << " / " << PopSize << endl ;

	if( myPop.size() < PopSize )
	{
		crowding( toReturn[ranking], fMaxMin ); // Calcul du crowding du rang intermédiaire

		int i = toReturn[ranking].size()-1 ;

		quicksortcrowd( toReturn[ranking], 0, i ); // Tri des individus selon leur crownding

		while( myPop.size() < PopSize )
		{
			myPop.push_back( toReturn[ranking][i] ); // rajout des individus avec le plus grand crowding à la population
			--i ;
		}
	}

}



/*!
 *	\brief Permet de retirer les clones d'un ranking
 *
 *	Fonction permettant de retirer les clones appartenant au ranking actuel. La fonction vérifier rang par rang si des clones existent dans le rang. En effet, vérifier la population entière n'est pas nécessaire, des clones étant nécessairement dans le même rang. Cette fonction ne touche pas à la population qui sera par le futur recréer à partir des rangs.
 *
 *  \param myResult         ensemble des rangs dont chaque int represente l'iterateur dans la population permet la corrélation entre toReturn et les sets de dominance
 *  \param toReturn         ensemble des rangs comportant les individus
 *  \param myPop            population d'individu actuellement traitée
 *
 *	\return le nombre de clône supprimer
 *
 */
int removeClone( vector< vector<Individu*> > &toReturn )
{
    int inlength = 0 ;

    for(unsigned int ranking = 0 ; ranking < toReturn.size() ; ++ranking ) // Si il existe des clones, ils se trouveront dans le même rang, on traite donc la vérification rang par rang
    {


        for(unsigned int i=0; i < toReturn[ranking].size(); ++i )
        {
            for(unsigned int j = i+1 ; j < toReturn[ranking].size() ; ++j )
            {
                if( toReturn[ranking][i]->isCloneOf(toReturn[ranking][j]) ) // si il y a un clone, on le supprime du rang.
                {
                    toReturn[ranking].erase(toReturn[ranking].begin() + j ) ;
                    inlength ++ ;
                    j = i ;
                }
            }
        }
    }

    // cout << " NB of non clone : " << inlength << endl ;

    return inlength ;

}



/*!
 *	\brief Range une population d'individu en rang et recréer à partir du ranking une population
 *
 *	Fonction permettant le rangement d'une population en rang, de supprimer les clones et de recréer à partir du rangement une population de taille PopSize en favorisant l'élitisme
 *
 *	\param myPop    Population à ranger.
 *
 *	\return le rangement qui a permis de recréer la population, indépendemment de la nouvelle population.
 *
 */
vector< vector<Individu*> > sortRank(vector<Individu*> &myPop)
{

	vector< vector<int> > dominance ; // Set de dominance de chaque individu
	vector< vector<int> > myResult ; // permet la corrélation entre la population actuel et le reste
	vector< vector<Individu*> > toReturn ; // ranking à retourner.


	myResult.push_back( vector<int>() ) ;
	toReturn.push_back( vector<Individu*>() ) ;



	/*
	int domineCounter[ myPop.size() ] ; // pour chaque individu, compte combien de fois il est dominé

	for( int i = 0; i < myPop.size() ; ++i ) // initialisation des compteurs et des sets de dominances
	{
		dominance.push_back( vector<int>() ) ;
		domineCounter[i] = 0 ;
	}
	*/



	// Recherche des min et max selon chaque objectif
	double fMaxMin[n_Objectives][2] ;

	for( int i = 0 ; i < n_Objectives ; ++i )
	{
		fMaxMin[i][0] = 0 ;
		fMaxMin[i][1] = INT_MAX ;
	}


	for(unsigned int i = 0 ; i < myPop.size() ; ++i )
	{
		for(int j = 0 ; j < n_Objectives ; ++j )
		{
			if( myPop[i]->getObjectiveValue(j) > fMaxMin[j][0] )
			{
				fMaxMin[j][0] = myPop[i]->getObjectiveValue(j) ;
			}
			if( myPop[i]->getObjectiveValue(j) < fMaxMin[j][1] )
			{
				fMaxMin[j][1] = myPop[i]->getObjectiveValue(j) ;
			}
		}
	}
	// Fin de la recherche

	vector<bool> currentRank ;
	int rank = 0 ;

	while( myPop.size() > 0 )
	{
		toReturn.push_back( vector<Individu*>() );
		currentRank = vector<bool>(myPop.size(),1);

		for(unsigned int i = 0 ; i < myPop.size() ; ++i )
		{
			if( currentRank[i] )
			{
				for(unsigned int j = i+1 ; j < myPop.size() ; ++j )
				{
					if( currentRank[j] )
					{
						if( myPop[i]->domine(myPop[j]) )
						{
							currentRank[j] = 0 ;
						}
						else if( myPop[j]->domine(myPop[i]) )
						{
							currentRank[i] = 0 ;
						}
					}
				}
			}

			if( currentRank[i] )
			{
				toReturn[rank].push_back(myPop[i]);
				myPop[i]->setRank(rank+1);
			}
			else
			{
				toReturn[rank+1].push_back(myPop[i]);
			}
		}

		rank++ ;
		myPop.clear() ;
		myPop = toReturn[rank] ;
		toReturn[rank].clear() ;
		currentRank.clear() ;
	}

	toReturn.erase( toReturn.end() ) ;



	if( removeclone )
   		n_Clone = removeClone( toReturn ) ; // Suppression des clônes



	// Début de la recréation de la population à partir du nouveau ranking

    myPop = vector<Individu*>()  ; // Réinitialisation de la population

    // Début de la sauvegarde "crowding Gandibleux"

	saveTotalCrowd( myPop, toReturn, fMaxMin ) ;

	// Fin de la sauvegarde "crowding Gandibleux"


	if( n_Clone > PopSize ) // réparation de la population si jamais le nombre de clône était trop important.
		{
			addSol(n_Clone-PopSize ,myPop,toReturn) ;
		}

	int toReturnSize = toReturn[0].size() ;

	if( toReturn[0].size() > PopSize*rank1maxsize ) // Si la taille du rang 1 est supérieur à un pourcentage de la population, on augmente la taille de celle-ci.
    {
			if( toReturnSize*(1+rank1maxsize) > 2*PopSize-n_Clone )
	        {
	        	addSol( toReturnSize*(1/rank1maxsize+rank1maxsize) - 2*PopSize + n_Clone ,myPop,toReturn) ; // on rajoute à la population un nombre de solution correspondant au trou éventuelle que connaitrait celle-ci sans ce rajout.
			    PopSize =  toReturnSize*(1/rank1maxsize+rank1maxsize) ;
			}
	       	else
	       	{
	       		if( 2*PopSize-n_Clone < PopSize*(1+rank1maxsize) )
	       			addSol( PopSize*(1+rank1maxsize) - ( 2*PopSize - n_Clone ) ,myPop,toReturn) ;
	       		else if( Gener == 0 && split*n_Delta < PopSize )
	       			addSol( PopSize*rank1maxsize ,myPop,toReturn) ;

	       		PopSize = PopSize*(1+rank1maxsize);
			}

	}



	//cout << " BEGIN RECREATE " << endl ;

    reCreate( toReturn, myPop, fMaxMin ) ; // récréation de la population à partir du ranking

    //cout << " END RECREATE " << endl ;

	return toReturn ;
}



/*!
 *	\brief Permet de Sélectionner un parent
 *
 *	Fonction servant à sélectionner deux parents aléatoirement, les compare et garde le meilleur.
 *
 *	\param myParent Population dans laquelle on selectionne les parents aléatoirement
 *
 *	\return un parent
 *
 */
Individu* selectParentRND( vector<Individu*> myParent )
{
	int rnd1 = (int) floor(rand() % PopSize) ;
	int rnd2 = (int) floor(rand() % PopSize) ;

	/* Moins aggressif au centre et les bords se perdent moins avec seulement crowding sur les parents */

	if( Gener % 2 == 0 ){

	/* Choix en fonction du rang puis du crowding1 pour un même rang */
		if( myParent[rnd1]->getRank() > myParent[rnd2]->getRank() )
		{
			return myParent[rnd1] ;
		}
		else if( myParent[rnd2]->getRank() > myParent[rnd1]->getRank() )
		{
			return myParent[rnd2] ;
		}
		else if( myParent[rnd1]->crowding_value > myParent[rnd2]->crowding_value )
		{
			return myParent[rnd1] ;
		}
		else
		{
			return myParent[rnd2] ;
		}
	}
	else
	{
	/* Choix en fonction du crowding2 */
		if( myParent[rnd1]->crowding_Total > myParent[rnd2]->crowding_Total )
		{
			return myParent[rnd1] ;
		}
		else
		{
			return myParent[rnd2] ;
		}
	}
}



/*!
 *	\brief Permet de lire le fichier de paremètre
 *
 *	Fonction servant à lire un fichier de paramètre conçu d'une manière générique préalablement décidé.
 *
 *	\return un booléen, TRUE si le fichier a été lu, FALSE sinon.
 *
 */
bool readparam(string param_path)
{

	ifstream read_param(param_path.c_str(), ios::in);
	string temp ;

	if( read_param ) // si le fichier existe, on continue notre lecture.
	{
	 	/* PROBLEM_CRITERIA */

	 	read_param >> temp ;
	 	read_param >> temp ;
	 	read_param >> inst_path ;

	 	read_param >> temp ;
	 	read_param >> opt_path ;

	 	read_param >> temp ;
	 	read_param >> PopSize ;

	 	read_param >> temp ;
	 	read_param >> mutaprob ;

	 	read_param >> temp ;
	 	read_param >> rank1maxsize ;

	 	read_param >> temp ;
	 	read_param >> crowdSave ;

	 	read_param >> temp ;
	 	read_param >> crowdTotal ;

	 	/* STOP_CRITERIA */

	 	read_param >> temp ;
		read_param >> temp ;
	 	read_param >> MaxTime ;

	 	read_param >> temp ;
	 	read_param >> MaxGener ;

	 	read_param >> temp ;
	 	read_param >> HyperVolCD ;

	 	read_param >> temp ;
	 	read_param >> PercentEvolution ;

	 	read_param >> temp ;
	 	read_param >> GRASPercent ;

	 	/* GRASP_CRITERIA */

	 	read_param >> temp ;
	 	read_param >> temp ;
	 	read_param >> split ;

	 	read_param >> temp ;
	 	read_param >> n_Delta ;

	 	read_param >> temp ;
	 	read_param >> alpha ;

	 	/* NSGA_CRITERIA */

	 	read_param >> temp ;
	 	read_param >> temp ;
		read_param >> TESTMODE ;

	 	read_param >> temp ;
	 	read_param >> removeclone ;

	 	read_param >> temp ;
	 	read_param >> parentarena ;

	 	read_param >> temp ;
	 	read_param >> gnuplot ;

	 	read_param >> temp ;
	 	read_param >> results ;

	 	return 1 ;
	}
	else
	{
		cout << " Param File Not Found, Program Will Terminate. " << endl;
		return 0 ;
	}
}



/*!
 *	\brief Fonction Calculant la qualité
 *
 *	Fonction calculant la distance euclidienne entre la solution pareto optimal trouvé par l'algorithme et celle fournit à l'utilisateur ci-celle-ci existe.
 *
 *	\param opt_path Chemin d'accès vers le fichier contenant les résultats optimaux
 * 	\param &ranking Notre ranking final permettant de calculer la distance au rang 1
 *
 */
void quality(string opt_path, vector< vector<Individu*> > &ranKing)
{
	ifstream read_pareto(opt_path.c_str(), ios::in);

	if( read_pareto ) // si le fichier existe, on continue notre lecture.
	{
		vector<int> optimum ; // contiendra les points des optimum

		int temp ;
		int n_zero = 0 ;
		int dist ;
		int	distemp ;
		int pow2 ;

		double totaldist = 0; // distance total calculé, servira à calculer la moyenne

		double best = INT_MAX; // meilleur distance trouvée

		read_pareto >> temp ;

	 	/* Lecture des valeurs des solutions exactes */
		while( !read_pareto.eof() ) // lecture et mise en mémoire de chacun des points de l'optimum
		{
			optimum.push_back(temp) ;

			for( int i = 1 ; i < n_Objectives ; ++i )
			{
				 read_pareto >> temp ;
				 optimum.push_back(temp) ;
			}

			read_pareto >> temp ;

		}
		read_pareto.close();

		// Ci après, calcul de la distance moyenne entre le front pareto optimal et l'optimum
		for(unsigned int i = 0 ; i < ranKing[0].size() ; ++i )
		{
			dist = INT_MAX ;


			/* Recherche de la solution exacte la plus proche de l'individu */
			for(unsigned int j = 0 ; j < optimum.size()/n_Objectives ; ++j )
			{

				distemp = 0 ;

				for( int k = 0 ; k < n_Objectives ; ++k ) // Recherche du point optimal le plus proche de notre point actuellement traité
				{
					pow2 = (ranKing[0][i]->getObjectiveValue(k) - optimum[j*n_Objectives+k]) ;
					distemp += pow2*pow2 ;
				}

				if( dist > distemp ) // mise à jour de la distance
				{
					dist = distemp ;
				}

				if( best > distemp ) // mise à jour de la meilleur distance
				best = distemp ;

				if( distemp == 0 ) // mise à jour du nombre de solution optimale trouvé
				n_zero++ ;
			}
			//cout << " Distance object " << i << " equal " << sqrt(dist) << endl ;
			totaldist += sqrt(dist) ;
		}

		// cout << totaldist << endl ;
		best = sqrt(best) ;
		totaldist = totaldist / ranKing[0].size() ;

		cout << " Average Distance From Pareto : " << totaldist << endl ;
		cout << " Best Distance Found : " << best << endl ;
		cout << " Percentage of Optimum found : " << n_zero << "/" << optimum.size()/2 << " = " << (double) n_zero*100 / (optimum.size()/2) << "%" << endl ;

		ostringstream convert ;
		convert << (Gener-1);
		string s_Gener = convert.str();

		s_Gener = "pareto/pareto" + s_Gener + ".dat";



		cout<<"Hyper Volume Quality :"<<endl;

		/* Calcul d'hypervolume de Zitzler */
		hypervolume(opt_path, s_Gener, n_Objectives);
	}
	else
		cout << " Optimum Files not found/readable " << endl;
}



/*!
 *	\brief Calcul la différence d'hypervolume entre deux Rang1 de générations différentes.
 *
 *	Permet le calcul de l'hypervolume du Rang1 de la génération actuelle et de le comparer avec le dernier HyperVolume calculé afin d'en obtenir une qualité d'avancée de la population.
 *
 *	\param Rank1 Rang1 dont l'on va calculerl l'hypervolume.
 *	\param nbr_objectives nombre d'objectif sur lesquels il faut calculer l'hypervolume
 *
 *	\return un booléen nous apprenant si la population actuelle a ou non avancé d'un pourcentage significatif ( donné en paramètre ) par rapport au dernier essai.
 *
 */
bool hyperVolumeDist(vector< Individu* > Rank1, int nbr_objectives)
{
		double    **front1, **front2;
		int       size_Front1;
		double    volFront1, volFront2;
		vector< int > Nadir ;

		/* Initialisation des points Nadir */
		for( int i = 0 ; i < nbr_objectives ; ++i )
			Nadir.push_back( INT_MAX ) ;

		size_Front1 = Rank1.size() ;

		/* Initialisation des tableaux d'objectif du rank1 et lastrank1 */
		front1 = (double**)malloc(size_Front1 * sizeof(double *));
	  	for (int i = 0; i < size_Front1; i++)
			(front1)[i] = (double*)malloc(nbr_objectives * sizeof(double));

		front2 = (double**)malloc(size_Front2 * sizeof(double *));
	  	for (int i = 0; i < size_Front2; i++)
			(front2)[i] = (double*)malloc(nbr_objectives * sizeof(double));



		/* Récupération des valeurs d'objectif du rang 1 actuel */
	   	for( int i = 0 ; i < size_Front1 ; ++i )
	  	{
	  		for( int j = 0 ; j < nbr_objectives ; ++j )
	  			front1[i][j] = Rank1[i]->getObjectiveValue(j) ;
	  	}
	  	/* Récupération des valeurs d'objectif de l'ancien rang 1 */
	  	for( int i = 0 ; i < size_Front2 ; ++i )
	  	{
	  		for( int j = 0 ; j < nbr_objectives ; ++j )
	  			front2[i][j] = lastRank1[i]->getObjectiveValue(j) ;
	  	}


		/* Finding Nadir Point */
		for( int j = 0 ; j < nbr_objectives ; ++j )
		{
			/* Finding Nadir Point in Rank1 */
			for( int i = 0 ; i < size_Front1 ; ++i )
			{
				if( Nadir[j] > front1[i][j] )
					Nadir[j] = front1[i][j] ;
			}
			/* Finding Nadir Point in lastRank1 */
			for( int i = 0 ; i < size_Front2 ; ++i )
			{
				if( Nadir[j] > front2[i][j] )
					Nadir[j] = front2[i][j] ;
			}

		}



		/* Updating point with nadir */
		for( int i = 0 ; i < size_Front1 ; ++i )
		{
			for( int j = 0 ; j < nbr_objectives ; ++j )
				front1[i][j] += -1*Nadir[j] ;
		}

		for( int i = 0 ; i < size_Front2 ; ++i )
		{
			for( int j = 0 ; j < nbr_objectives ; ++j )
				front2[i][j] += -1*Nadir[j] ;
		}



		/* calculate dominated hypervolume */
		volFront1 = CalculateHypervolume(front1, size_Front1, nbr_objectives);
		volFront2 = CalculateHypervolume(front2, size_Front2, nbr_objectives);

		/* Preparing next time */
		lastRank1.clear() ;

		for(unsigned int i = 0 ; i < Rank1.size() ; ++i )
			lastRank1.push_back( Rank1[i] ) ;

		size_Front2 = size_Front1 ;

		/* retour de l'avancé du rang 1 actuel en fonction de PercentEvolution */
		if( (volFront1 - volFront2) / volFront2 > PercentEvolution || (volFront2 - volFront1) / volFront1 > PercentEvolution )
			return true ;
		else return false ;
}




/*!
 *	\brief Calcul la mesure de qualité du spreading
 *
 *	Permet d'obtenir une mesure de qualité vis à vis du spreading de la population.
 *
 *	\param result population d'individu que l'on aurait obtenu suite à un filtrage de la population totale.
 *
 *	\return un double correspondant à la qualité du spread. Si proche de 0, la qualité est meilleure ( selon Deb ).
 *
 */
double spread( vector< Individu* > result )
{
	double fMaxMin[n_Objectives][2] ;
	double AveCrowd = 0.0 ;
	double spread = 0.0 ;

	for( int i = 0 ; i < n_Objectives ; ++i )
	{
		fMaxMin[i][0] = 0 ;
		fMaxMin[i][1] = INT_MAX ;
	}

	for(unsigned	 int i = 0 ; i < result.size() ; ++i )
	{
		for(int j = 0 ; j < n_Objectives ; ++j )
		{
			if( result[i]->getObjectiveValue(j) > fMaxMin[j][0] )
			{
				fMaxMin[j][0] = result[i]->getObjectiveValue(j) ;
			}
			if( result[i]->getObjectiveValue(j) < fMaxMin[j][1] )
			{
				fMaxMin[j][1] = result[i]->getObjectiveValue(j) ;
			}
		}
	}

	crowdingTotal(result, fMaxMin );
	quicksortcrowdTotal(result, 0, result.size()-1 );

	for(unsigned int i = 0 ; i < result.size()-(n_Objectives*2) ; ++i )
	{
		AveCrowd += result[i]->crowding_Total ;
	}

	AveCrowd = AveCrowd / ( result.size()-(n_Objectives*2) ) ;



	for(unsigned int i = 0 ; i < result.size()-(n_Objectives*2) ; ++i )
	{
		if( result[i]->crowding_Total - AveCrowd > 0 )
			spread += result[i]->crowding_Total - AveCrowd ;
		else
			spread += AveCrowd - result[i]->crowding_Total ;
	}

	spread = spread / ( ( result.size() - (n_Objectives*2) ) * AveCrowd ) ;


	return spread ;
}



/* Les fonctions qui suivent ont été copiées puis adaptées à notre code à partir des méthodes de différence d'hypervolume de M. Eckart Zitzler */


/*---------------------------------------------------------------------------*/
/* Hypervolume Metric Calculation                                            */
/*---------------------------------------------------------------------------*/
/* This program calculates for a given set of objective vectors the volume   */
/* of the dominated space, enclosed by the nondominated points and the       */
/* origin. Here, a maximization problem is assumed, for minimization         */
/* or mixed optimization problem the objective vectors have to be trans-     */
/* formed accordingly. The hypervolume metric has been proposed in:          */
/*                                                                           */
/* 1. E. Zitzler and L. Thiele. Multiobjective Optimization Using            */
/*    Evolutionary Algorithms - A Comparative Case Study. Parallel Problem   */
/*    Solving from Nature - PPSN-V, September 1998, pages 292-301.           */
/*                                                                           */
/* A more detailed description and extensions can be found in:               */
/*                                                                           */
/* 2. E. Zitzler. Evolutionary Algorithms for Multiobjective Optimization:   */
/*    Methods and Applications. Swiss Federal Institute of Technology (ETH)  */
/*    Zurich. Shaker Verlag, Germany, ISBN 3-8265-6831-1, December 1999.     */
/*                                                                           */
/*---------------------------------------------------------------------------*/
/* Usage (command line parameters):                                          */
/*                                                                           */
/* 1) number of objectives                                                   */
/* 2) name of the file containing the first set of objective vectors         */
/* 3) name of the file containing the second set of objective vectors (this  */
/*    parameter is optional)                                                 */
/*                                                                           */
/* The file format is as follows: Each line in the file describes one point  */
/* of the trade-off front and contains a sequence of numbers, separated by   */
/* blanks. Per line, the first number corresponds to the first objective,    */
/* the second number to the second objective, and so forth.                  */
/*                                                                           */
/* Output:                                                                   */
/*                                                                           */
/* 1) volume of the space dominated by the first set of objective vectors    */
/* 2) volume of the space dominated by the second set of objective vectors   */
/*    (if the third command line parameter has been specified)               */
/* 3) volume of the space which is dominated by the first set of objective   */
/*    but not by the second set of objective vectors (if the third command   */
/*    line parameter has been specified)                                     */
/* 4) volume of the space which is dominated by the second set of objective  */
/*    but not by the first set of objective vectors (if the third command    */
/*    line parameter has been specified)                                     */
/*                                                                           */
/* Outputs 1+2) refer to the S metric and outputs 3+4) to the D metric as    */
/* described in reference 2 (see above).                                     */
/*---------------------------------------------------------------------------*/
/* Eckart Zitzler                                                            */
/* Computer Engineering and Networks Laboratory (TIK)                        */
/* Swiss Federal Institute of Technology (ETH) Zurich, Switzerland           */
/* (c)2001                                                                   */
/*---------------------------------------------------------------------------*/

int  Dominates(double  point1[], double  point2[], int  noObjectives)
     /* returns true if 'point1' dominates 'points2' with respect to the
	to the first 'noObjectives' objectives */
{
  int  i;
  int  betterInAnyObjective;

  betterInAnyObjective = 0;
  for (i = 0; i < noObjectives && point1[i] >= point2[i]; i++)
      if (point1[i] > point2[i])
	betterInAnyObjective = 1;
  return (i >= noObjectives && betterInAnyObjective);
} /* Dominates */

void  Swap(double  *front[], int  i, int  j)
{
  double  *temp;

  temp = front[i];
  front[i] = front[j];
  front[j] = temp;
} /* Swap */

int  FilterNondominatedSet(double  *front[], int  noPoints, int  noObjectives)
     /* all nondominated points regarding the first 'noObjectives' dimensions
	are collected; the points referenced by 'front[0..noPoints-1]' are
	considered; 'front' is resorted, such that 'front[0..n-1]' contains
	the nondominated points; n is returned */
{
  int  i, j;
  int  n;

  n = noPoints;
  i = 0;
  while (i < n) {
    j = i + 1;
    while (j < n) {
      if (Dominates(front[i], front[j], noObjectives)) {
	/* remove point 'j' */
	n--;
	Swap(front, j, n);
      }
      else if (Dominates(front[j], front[i], noObjectives)) {
	/* remove point 'i'; ensure that the point copied to index 'i'
	   is considered in the next outer loop (thus, decrement i) */
	n--;
	Swap(front, i, n);
	i--;
	break;
      }
      else
	j++;
    }
    i++;
  }
  return n;
} /* FilterNondominatedSet */


double  SurfaceUnchangedTo(double  *front[], int  noPoints, int  objective)
     /* calculate next value regarding dimension 'objective'; consider
	points referenced in 'front[0..noPoints-1]' */
{
  int     i;
  double  minValue, value;

  if (noPoints < 1)  ERROR("run-time error");
  minValue = front[0][objective];
  for (i = 1; i < noPoints; i++) {
    value = front[i][objective];
    if (value < minValue)  minValue = value;
  }
  return minValue;
} /* SurfaceUnchangedTo */

int  ReduceNondominatedSet(double  *front[], int  noPoints, int  objective,
			   double  threshold)
     /* remove all points which have a value <= 'threshold' regarding the
	dimension 'objective'; the points referenced by
	'front[0..noPoints-1]' are considered; 'front' is resorted, such that
	'front[0..n-1]' contains the remaining points; 'n' is returned */
{
  int  n;
  int  i;

  n = noPoints;
  for (i = 0; i < n; i++)
    if (front[i][objective] <= threshold) {
      n--;
      Swap(front, i, n);
    }
  return n;
} /* ReduceNondominatedSet */

double  CalculateHypervolume(double  *front[], int  noPoints,
			     int  noObjectives)
{
  int     n;
  double  volume, distance;

  volume = 0;
  distance = 0;
  n = noPoints;
  while (n > 0) {
    int     noNondominatedPoints;
    double  tempVolume, tempDistance;

    noNondominatedPoints = FilterNondominatedSet(front, n, noObjectives - 1);
    tempVolume = 0;
    if (noObjectives < 3) {
      if (noNondominatedPoints < 1)  ERROR("run-time error");
      tempVolume = front[0][0];
    }
    else
      tempVolume = CalculateHypervolume(front, noNondominatedPoints,
					noObjectives - 1);
    tempDistance = SurfaceUnchangedTo(front, n, noObjectives - 1);
    volume += tempVolume * (tempDistance - distance);
    distance = tempDistance;
    n = ReduceNondominatedSet(front, n, noObjectives - 1, distance);
  }
  return volume;
} /* CalculateHypervolume */

int  ReadFront(double  **frontPtr[], FILE  *file, int  noObjectives)
{
  int     noPoints;
  int     i;
  double  value;

  /* check file and count points */
  noPoints = 0;
  while (!feof(file)) {
    for (i = 0; i < noObjectives && fscanf(file, "%lf", &value) != EOF; i++);
    if (i > 0 && i < noObjectives)  ERROR("data in file incomplete");
    noPoints++;
  }
  /* allocate memory */
  *frontPtr =(double**) malloc(noPoints * sizeof(double *));
  if (*frontPtr == NULL)  ERROR("memory allocation failed");
  for (i = 0; i < noPoints; i++) {
    (*frontPtr)[i] = (double*)malloc(noObjectives * sizeof(double));
    if ((*frontPtr)[i] == NULL)  ERROR("memory allocation failed");
  }
  /* read data */
  rewind(file);
  noPoints = 0;
  while (!feof(file)) {
    for (i = 0; i < noObjectives; i++) {
      if (fscanf(file, "%lf", &value) != EOF)
	(*frontPtr)[noPoints][i] = value;
      else
	break;
    }
    if (i > 0 && i < noObjectives)  ERROR("data in file incomplete");
    noPoints++;
  }
  if (noPoints < 1)  ERROR("file contains no data");
  return noPoints;
} /* ReadFront */

int  MergeFronts(double  **frontPtr[], double  *front1[], int  sizeFront1,
		 double*  front2[], int  sizeFront2, int  noObjectives)
{
  int  i, j;
  int  noPoints;

   /* allocate memory */
  noPoints = sizeFront1 + sizeFront2;
  *frontPtr = (double**)malloc(noPoints * sizeof(double *));
  if (*frontPtr == NULL)  ERROR("memory allocation failed");
  for (i = 0; i < noPoints; i++) {
    (*frontPtr)[i] = (double*)malloc(noObjectives * sizeof(double));
    if ((*frontPtr)[i] == NULL)  ERROR("memory allocation failed");
  }
  /* copy points */
  noPoints = 0;
  for (i = 0; i < sizeFront1; i++) {
    for (j = 0; j < noObjectives; j++)
      (*frontPtr)[noPoints][j] = front1[i][j];
    noPoints++;
  }
  for (i = 0; i < sizeFront2; i++) {
    for (j = 0; j < noObjectives; j++)
      (*frontPtr)[noPoints][j] = front2[i][j];
    noPoints++;
  }

  return noPoints;
} /* MergeFronts */

void  DeallocateFront(double**  front, int  noPoints)
{
  int  i;

  if (front != NULL) {
    for (i = 0; i < noPoints; i++)
      if (front[i] != NULL)
	free(front[i]);
    free(front);
  }
} /* DeallocateFront */

int hypervolume(string pareto1, string pareto2, int nbr_objectives)
{
  FILE      *file1, *file2;
  double    **front1, **front2, **front3;
  int       sizeFront1, sizeFront2, sizeFront3;
  int       redSizeFront1, redSizeFront2, redSizeFront3;
  double    volFront1, volFront2, volFront3;
  int       noObjectives;

  /* check parameters */

  noObjectives = nbr_objectives;
  file1 = fopen(pareto1.c_str(), "r");
  file2 = fopen(pareto2.c_str(), "r");

  if(!file2)
  	cout<<"Erreur: fichier "<<pareto2<<endl;



  /* read in data */
  sizeFront1 = ReadFront(&front1, file1, noObjectives);
  fclose(file1);

  sizeFront2 = ReadFront(&front2, file2, noObjectives);
  sizeFront3 = MergeFronts(&front3, front1, sizeFront1, front2, sizeFront2,
				 noObjectives);
  fclose(file2);
  /* calculate dominated hypervolume */
  redSizeFront1 = FilterNondominatedSet(front1, sizeFront1, noObjectives);
  volFront1 = CalculateHypervolume(front1, redSizeFront1, noObjectives);
  printf("%.10f ", volFront1);


  redSizeFront2 = FilterNondominatedSet(front2, sizeFront2, noObjectives);
  volFront2 = CalculateHypervolume(front2, redSizeFront2, noObjectives);
  printf("%.10f ", volFront2);
  redSizeFront3 = FilterNondominatedSet(front3, sizeFront3, noObjectives);
  volFront3 = CalculateHypervolume(front3, redSizeFront3, noObjectives);
  printf("%.10f %.10f", volFront3 - volFront2,
	   volFront3 - volFront1);

  DeallocateFront(front1, sizeFront1);
  DeallocateFront(front2, sizeFront2);
  DeallocateFront(front3, sizeFront3);
  printf("\n");
}
