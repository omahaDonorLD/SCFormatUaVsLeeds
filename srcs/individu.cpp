
#include "../heads/individu.hpp"
#include "kmeans.c"

using namespace std;

// Creates randomly a solution
sln* randomize()
{
	sln rand_sln;
	return rand_sln;
}


// creates individuals using the clustering method kmeans to partition the ground nodes
sln* clustering()
{
	sln sln_kmeans;
}


// creates individuals that are equally spread on the search space : Width/(2*radius-constant). With "cst" so that pairs of UaVs share a space equal to that constant (kind of ;) )
sln* uniform()
{
	sln uniform_sln;
	return uniform_sln;
}


/**
 * \brief       Méthode de classe qui initialise le problème.
 * \details    Cette méthode lit un fichier de données formaté pour initialiser les paramètres du problème, comme le nombre d'items, les contraintes, les objectifs, etc.
 * \param   	instance_path		Le chemin du fichier de données.
 *
 */
void Individu::initProblem(string instance_path) //parse le fichier d'instance de Zitzler
{
	ifstream read(instance_path.c_str(), ios::in);

	nbr_cnstrts = 2;  //attention ici seulement pour notre type de sac à dos

	if(read)
	{
		string value;

		for(int i = 0; i<4; ++i)
		{//lit les premiers mots
			read>>value;
		}

		read>>nbr_objs; //lit le nombre d'objectifs

		cout<<"nbr Objectif : "<<nbr_objs<<endl;

		read>>value;

		read>>nbr_uavs; //lit le nombre d'item
		cout<<"nbr items : "<<nbr_uavs<<endl;

		read>>value;
		read>>value;

		for(int curr_obj = 0; curr_obj < nbr_objs; ++curr_obj)
		{
			//cout<<"Objectif : "<<curr_obj<<endl;
			read>>value; //lit le =
			read>>value; read>>value;//lit le "knapsack n "
			read>>value;// lit le "capacity:"

			int nbr;

			read>>nbr; //lit la capacité
			cout<<"capacity : "<<nbr<<endl;

			if (curr_obj < 2) //on se limite ici à deux contraintes
				capacities.push_back(nbr);

			for(int curr_item = 0; curr_item < nbr_uavs; ++curr_item)
			{
				//cout<<"Item : "<<curr_item<<endl;
				read>>value; read>>value;//lit le "item n "

				read>>value; //lit le "weight:"
				read>>nbr; //lit la valeur du poids

				if(curr_obj < 2)//on se limite ici à deux contraintes
					weights.push_back(nbr);

				read>>value; //lit le "profits:"
				read>>nbr; //lit la valeur du profit

				profits.push_back(nbr);
			}
		}

	}
	else
	{
		cout<<"Erreur ouverture fichier function "<< __FUNCTION__ <<endl;
	}

	read.close();

	cout<<"End Init"<<endl;
}



/**
 * \brief       Retourne la valeur d'un objectif de l'individu.
 * \details    Calcule la valeur d'un objectif de l'individu à partir des paramètres du problèmes et des items sélectionnés.
 * \param   	num_objective	 le numéro de l'objectif calculé. (Commence à 0)
 * \return    Un \e int représentant l'objectif calculée.
 */
int Individu::getObjectiveValue(int num_objective)
{
	return current_objective_value[num_objective];
}



/**
 * \brief       Retourne la valeur d'une contrainte de l'individu.
 * \details    Calcule la valeur d'une contrainte de l'individu à partir des paramètres du problèmes et des items sélectionnés.
 * \param   	num_contraint	 le numéro de la contrainte calculé. (Commence à 0)
 * \return    Un \e int représentant l'objectif calculée.
 */
int Individu::getContraintValue(int num_contraint)
{
	return current_constraint_value[num_contraint];
}



/*!
 *	\brief Calcule la valeur de chaque Objectif et la met à jour.
 *
 *	Fonction permettant de mettre à jour après différent changement éventuellement aléatoire les valeurs de chaque objectifs.
 *
 */
void Individu::computeObjectives()
{
	/* Calcule la valeur de chaque objectif */
	for(int curr_obj = 0; curr_obj<nbr_objs; ++curr_obj)
	{
		int curr_total = 0;

		/* Calcule la valeur de l'objectif en cours */
		for(int curr_item = 0; curr_item < nbr_uavs; ++curr_item)
		{
			curr_total+= picked_objects[curr_item] * getProfit(curr_item, curr_obj);
		}

		current_objective_value[curr_obj] = curr_total;
	}
}



 /*!
 *	\brief Calcule la valeur de chaque contraintes et la met à jour.
 *
 *	Fonction permettant de mettre à jour après différent changement éventuellement aléatoire les valeurs de chaque contraintes.
 *
 */
void Individu::computeConstraints()
{
	/* Calcule la valeur de chaque contrainte */
	for(int curr_obj = 0; curr_obj<nbr_cnstrts; ++curr_obj)
	{
		int curr_total = 0;

		/* Calcule la valeur de la contrainte en cours */
		for(int curr_item = 0; curr_item < nbr_uavs; ++curr_item)
		{
			curr_total+= picked_objects[curr_item] * getWeight(curr_item, curr_obj);
		}

		current_constraint_value[curr_obj] = curr_total;
	}
}



/**
 * \brief       Calcule la dominance entre deux solutions
 * \details    Les valeurs des objectifs sont comparées pour obtenir le résultat de la domainance stricte.
 * \param    to_compare    Solution avec laquelle la solution courante est comparée.
 * \return    Un \e bool indiquant si la solution courante domine ou pas celle fournie en paramètre.
 */
bool Individu::domine(Individu* to_compare)
{
	bool strictly = false;

	/* Vérification de la dominance suivant chaque objectif */
	for(int obj = 0; obj<nbr_objs; ++obj)
	{
		if( this->getObjectiveValue(obj) < to_compare->getObjectiveValue(obj) )
			return false; 		// Si il est dominé selon un objectif, l'on renvoie faux
		else if( this->getObjectiveValue(obj) > to_compare->getObjectiveValue(obj) )
			strictly = true;	// Si il domine strictement suivant un objectif, on se prépare à renvoyer vrai.
	}

	return strictly;
}



/**
 * \brief       Indique si deux solution sont identiques.
 * \details   	Les items sélectionnés sont comparés.
 * \param    inIndividu		Solution avec laquelle la solution courante est comparée.
 * \return    Un \e bool indiquant si les deux solutions sont identiques.
 */
bool Individu::isCloneOf(Individu* inIndividu)
{
	/* Vérification de la valeur des objectifs */
	for( int i = 0; i < nbr_objs ; ++i )
		if( this->getObjectiveValue(i) != inIndividu->getObjectiveValue(i) )
			return false ;
	/* Si la valeur des objectifs n'a pu déterminer qu'il n'était pas un clône, on vérifie chaque valeur d'item ( comparaison des codes génétiques ) */
	for(int curr_item = 0 ; curr_item < nbr_uavs; ++curr_item)
		if(this->picked_objects[curr_item] != inIndividu->picked_objects[curr_item])
			return false;

	return true;
}



/**
 * \brief       Fournit un enfant issu de deux solutions parents
 * \details    Deux enfants sont produits par un crossover à un point. Si l'un des deux enfant domine l'autre, celui-ci est sélectionné. Sinon, le second est pruis paer défaut.
 				A l'issue du cross-over, les deux enfants sont rendu faisables.
 * \param    parent1         Le premier parent.
 * \param    parent2         Le second parent.
 * \return    L'enfant sélectionné.
 */
Individu* Individu::getChildFrom(Individu* parent1, Individu* parent2)//crossOver à un point
{
	int cross_point = (rand()% nbr_uavs);

	Individu* child1 = new Individu(parent1);

	Individu* child2 = new Individu(parent2);

	/* Copie des items du parents opposé pour le cross-over */
	for(int curr_item = 0; curr_item < cross_point; ++curr_item)
	{
		child1->picked_objects[curr_item] = parent2->picked_objects[curr_item];
		child2->picked_objects[curr_item] = parent1->picked_objects[curr_item];
	}


	/* On rend les enfants viables */
	child1->makeFeasible();
	child2->makeFeasible();


	/* Comparaison des enfants */
	if(child1->domine(child2))
	{
		return child1;
	}
	else
	{
		return child2;
	}

}



/**
 * \brief       Mute une solution
 * \details    Cette méthode est  pour l'instant vide.
 */
void Individu::mutate()
{
	this->restart2() ;

	/* Problem with restart3() */
}



/*!
 *	\brief Mutation 1
 *
 *	Mutation consistant à prendre un 1 aléatoirement et à le placer à 0.
 *
 */
void Individu::restart()
{
	vector<int> Picked1 ;

	/* Recherche des items à 1 que l'on peut passer à 0 */
	for( int curr_item = 0 ; curr_item < nbr_uavs ; curr_item++ )
	{
		if( picked_objects[curr_item] )
			Picked1.push_back(curr_item);
	}

	/* Sélection et suppression d'un item aléatoire */
	int rnd = rand() % Picked1.size() ;

	this->removeObject( Picked1[rnd] ) ;
}



/*!
 *	\brief Mutation 2
 *
 *	Mutation consistant à rajouter tous les items disponibles jusqu'à ce le sac à dos soit plein.
 *
 */
void Individu::restart2()
{
	vector<int> Picked0 = vector<int>() ;
	vector<int> topick = vector<int>() ;

	/* Sélection des items à 0 à rajouter */
	for( int curr_item = 0 ; curr_item < nbr_uavs ; curr_item++ )
	{

		/* Vérification que l'item en cours est à 0 */
		if( !(this->picked_objects[curr_item]) )
		{
			bool pickable = true ;

			/* Vérification qu'il ne dépasse aucune contrainte si rajouté */
			for( int curr_dim = 0 ; curr_dim < nbr_cnstrts && pickable ; curr_dim++)
				pickable = ( this->current_constraint_value[curr_dim] + getWeight(curr_item, curr_dim) ) < capacities[curr_dim]  ;

			/* Si rajoutable on le rajoute dans les items à rajouter  */
			if( pickable )
			{
				Picked0.push_back(curr_item);
				topick.push_back(curr_item);
			}
		}
	}

	/* Si il existe des items à rajouter, sélection aléatoire d'un et rajout de celui-ci */
	if( Picked0.size() > 0 )
	{
		int rnd = rand() % Picked0.size() ;

		this->addObject( Picked0[rnd] ) ;
	}

	/* à partir de la liste d'item précédente et tant qu'elle ne sera pas vide, recherche d'item à rajouter */
	while( !(Picked0.size() == 0) )
	{
		Picked0 = vector<int>() ;

		for(unsigned int i = 0 ; i < topick.size() ; ++i )
		{

			if( !(this->picked_objects[topick[i]]) )
			{
				bool pickable = true ;

				for( int curr_dim = 0 ; curr_dim < nbr_cnstrts && pickable; curr_dim++)
					pickable = ( this->current_constraint_value[curr_dim] + getWeight(topick[i], curr_dim) ) < capacities[curr_dim] ;

				if( pickable )
				{
					Picked0.push_back(topick[i]);
				}
			}
		}

		if( Picked0.size() > 0 )
		{
			int rnd = rand() % Picked0.size() ;

			this->addObject( Picked0[rnd] ) ;
		}
	}
}



/*!
 *	\brief Mutation 3
 *
 *	Mutation consistant à prendre aléatoirement un item à 1 et un item à 0 et échanger leur valeur ( avec vérification préalable de la possibilité de cette échange ).
 *
 */
void Individu::restart3()
{
	vector<int> Picked1 ;

	/* Recherche des items à 1 pour l'échange 1-1 */
	for( int curr_item = 0 ; curr_item < nbr_uavs ; curr_item++ )
	{
		if( picked_objects[curr_item] )
			Picked1.push_back(curr_item);
	}

	vector< vector<int> > Couple;

	/* Recherche des items à 0 pouvant s'échanger avec les items à 1 */
	for(unsigned int i = 0 ; i < Picked1.size() ; ++i )
	{
		Couple.push_back( vector<int>() ); // liste de possibilité des items à 1

		/* Recherche des items à 0 pouvant s'échanger avec UN item à 1 */
		for( int curr_item = 0 ; curr_item < nbr_uavs ; curr_item++ )
		{
			if( !picked_objects[curr_item] )
			{
				bool pickable = true ;

				/* Si l'échange ne brise pas les contraintes rajout de l'item à 0 dans la liste de possibilité de l'item à 1 */
				for( int curr_dim = 0; curr_dim < nbr_cnstrts ; curr_dim++)
					pickable = pickable && ( ( current_constraint_value[curr_dim] + getWeight(curr_item, curr_dim) - getWeight(Picked1[i], curr_dim) ) < capacities[curr_dim] ) ;

				if( pickable )
					Couple[i].push_back(curr_item);
			}
		}
	}

	vector<int> test ;

	/* vérification des items à 1 dont l'on a trouvé un échange possible */
	for(unsigned int i = 0 ; i < Picked1.size() ; ++i )
	{
		if( Couple[i].size() > 0 )
		{
			test.push_back(i) ;
		}
	}

	/* si des items à 1 sont échangeable, sélection aléatoire d'un item à 1 et d'un 0 correspondant */
	if( test.size() > 0 )
	{
		int rnd1 = rand() % test.size() ;
		int rnd2 = rand() % Couple[Picked1[test[rnd1]]].size();

		this->removeObject(Picked1[rnd1]) ;
		this->addObject(Couple[rnd1][rnd2]) ;
	}
}



/**
 * \brief      Récupère le poids d'un item.
 * \details    Calcule le poids d'un item donné dans une contrainte donnée. Cette valeur est une donnée du problème, pas d'une solution particulière.
 * \param    num_item       Le numéro de l'item dont on veut connaître le poids.
 * \param    dimension      Le numéro de la contrainte pour laquelle on veut connaître le poids.
 * \return    Un \e int donnant ce poids.
 */
int Individu::getWeight(int num_item, int dimension)
{
	return weights[dimension * nbr_uavs  + num_item];
}



/**
 * \brief       Récupère le profit (valeur) d'un item.
 * \details    Calcule le profit d'un item donné dans un objectif donnée. Cette valeur est une donnée du problème, pas d'une solution particulière.
 * \param    num_item       Le numéro de l'item dont on veut connaître le profit.
 * \param    num_objective      Le numéro de l'objectif pour laquelle on veut connaître le profit.
 * \return     Un \e int donnant ce profit.
 */
int Individu::getProfit(int num_item, int num_objective)
{
	return profits[num_objective * nbr_uavs  + num_item];
}



/**
 * \brief       Calcule le poids total d'une solution
 * \details    La somme des poids dans une contrainte donnée est calculée et renvoyée. Cette valeur dépend de la solution, seuls les pouids des items sélectionnés sont sommés.
 * \param    dimension      Le numéro de la contrainte pour laquelle on veut connaître le poids.
 * \return    Un \e int représentant la somme des poids calculée.
 */
int Individu::getTotalWeight(int dimension)
{
	/*int total_weight = 0;

	for(int curr_item = 0; curr_item < nbr_uavs; ++curr_item)
	{
		total_weight+= picked_objects[curr_item] * getWeight(curr_item, dimension);
	}

	return total_weight;*/

	return getContraintValue(dimension) ;
}



/**
 * \brief      Vérifie si une contrainte d'une solution proposée est respectée.
 * \param    dimension         le numéro de la contrainte dont on veut vérifier si elle est respectée.
 * \return    Un \e bool représentant si la contrainte est respectée ou pas.
 */
bool Individu::constraintRespected(int dimension)
{
	return (this->current_constraint_value[dimension] < capacities[dimension]);
}



/**
 * \brief       Calcule si la solution proposée est faisable ou pas.
 * \details    Cette méthode vérifie que chaque contrainte est bien respectée.
 * \return    Un \e bool représentant si la solution est faisable ou pas.
 */
bool Individu::isFeasible()
{
	/* Vérification que l'individu est valide ( ne dépasse pas les contraintes ) */
	for(int curr_constraint = 0; curr_constraint < nbr_cnstrts; ++curr_constraint)
		if( !(this->constraintRespected(curr_constraint)) )
			return false;

	return true;
}



/**
 * \brief       Rend une solution faisable.
 * \details    Cette méthode procéde contrainte par contrainte. Elle choisit l'objet dont le poids est le plus proche du débordement de la contrainte et les retire. Si ce n'est pas suffisant, elle en retire un autre, etc, et ce pour chaque contrainte une par une tant que la solution n'est pas faisable.
  */
void Individu::makeFeasible() //rend la solution faisable en retirant les éléments les plus lourds ou les plus proches du débordement
{



	/* Calcule des objectifs et contraintes de l'individu */

	this->computeConstraints() ;
	this->computeObjectives();


	int overflow ;
	int overflow2 ;
	int selObject ;

	while(!this->isFeasible()) //rend la solution faisable si elle ne l'est pas.
	{
		selObject = 0 ;
		overflow = 0 ;

		for( int curr_item = 1 ; curr_item < nbr_uavs && !picked_objects[curr_item-1] ; curr_item++ )
		{
			selObject = curr_item ;
		}

		for( int curr_constraint = 0 ; curr_constraint < nbr_cnstrts ; curr_constraint++ )
		{
			if( !constraintRespected(curr_constraint) )
			overflow += (current_constraint_value[curr_constraint]-capacities[curr_constraint])-getWeight(selObject, curr_constraint ) ;
		}

		for( int curr_item = selObject ; curr_item < nbr_uavs ; curr_item++ )
		{
			if( picked_objects[curr_item] )
			{
				overflow2 = 0 ;

				for( int curr_constraint = 0 ; curr_constraint < nbr_cnstrts ; curr_constraint++ )
				{
					if( !constraintRespected(curr_constraint) )
					overflow2 += (current_constraint_value[curr_constraint]-capacities[curr_constraint])-getWeight(curr_item, curr_constraint ) ;
				}

				if( overflow > 0 && overflow2 < overflow )
				{
					selObject = curr_item ;
					overflow = overflow2 ;
				}
				else if( overflow < 0 && overflow2 < 0 && overflow2 > overflow )
				{
					selObject = curr_item ;
					overflow = overflow2 ;
				}
			}
		}

		/*
		for( int curr_item = 0 ; curr_item < nbr_uavs ; curr_item++ )
		{
			overflow[curr_item] = 0 ;
		}

		for( int curr_item = 0 ; curr_item < nbr_uavs ; curr_item++ )
		{
			if( picked_objects[curr_item] )
			{
				for( int curr_constraint = 0 ; curr_constraint < nbr_cnstrts ; curr_constraint++ )
				{
					if( !constraintRespected(curr_constraint) )
					overflow[curr_item] += (current_constraint_value[curr_constraint]-capacities[curr_constraint])-getWeight(curr_item, curr_constraint ) ;
				}
			}
			else
			overflow[curr_item] = INT_MAX ;

		}

		for( int curr_item = 1 ; curr_item < nbr_uavs ; curr_item++ )
		{
			if( picked_objects[curr_item] )
			{
				if( overflow[selObject] > 0 && overflow[curr_item] < overflow[selObject] )
					selObject = curr_item ;
				else if( overflow[selObject] < 0 && overflow[curr_item] < 0 && overflow[curr_item] > overflow[selObject] )
					selObject = curr_item ;
			}
		}*/

		removeObject(selObject) ;
	}


	/*
	if(!this->isFeasible()) //rend la solution faisable si elle ne l'est pas.
	{
		int min_difference;
		int heaviest_weight;
		int overflow;
		int closest_item;
		int heaviest_item;
		int curr_constraint = 0;
		do
		{
			while(!this->constraintRespected(curr_constraint)) //tant que la contrinte n'est pas respectée.
			{
				closest_item = -1; //L'item le plus proche du débordement
				heaviest_item = -1; //l'item le plus lourd.
				heaviest_weight = 0;// le poids le plus lourd
				min_difference = 32767; //la différence minimale entre un poids et le débordement.

				overflow = this->getContraintValue(curr_constraint) - capacities[curr_constraint]; // Le débordement.

				for(int curr_item = 0; curr_item < nbr_uavs; ++curr_item) // recherche l'objet le plus lourd et le plus proche.
				{
					if( picked_objects[curr_item])
					{
						int curr_weight = getWeight(curr_item, curr_constraint);
						if(curr_weight >= overflow && (curr_weight-overflow)< min_difference)
						{
							min_difference = curr_weight-overflow;
							closest_item = curr_item;
						}
						if(curr_weight>heaviest_weight)
						{
							heaviest_weight = curr_weight;
							heaviest_item = curr_item;
						}
					}
				}

				if(closest_item == -1) //Si aucun objet ne suffit à retirer le débordement
				{
					this->removeObject(heaviest_item) ;//retire l'objet le plus lourd
				}
				else
				{
					this->removeObject(closest_item) ; //retire l'objet le plus proche.
				}
			}
			curr_constraint++; //passe à la contrainte suivante.
		}while(curr_constraint<nbr_cnstrts && !this->isFeasible()); //tant que la solution n'est pas faisable.
	}*/



}



/*!
 *	\brief Ajoute un objet à un individu
 *
 *	Permet la gestion de l'ajout d'un objet à un individu et par cela la mise à jour des objectifs et contraintes.
 *
 *	\param num_item numéro de l'item à rajouter
 *
 */
void Individu::addObject(int num_item)
{
	/* Si l'object n'est pris on le rajoute avec mise à jour des objectifs et contraintes */
	if( !picked_objects[num_item] )
	{
		picked_objects[num_item] = 1 ;

		/* Mise à jour des objectifs */
		for( int curr_obj = 0 ; curr_obj < nbr_objs ; curr_obj++ )
		{
			current_objective_value[curr_obj] += getProfit(num_item, curr_obj) ;
		}
		/* Mise à jour des contraintes */
		for( int curr_dim = 0 ; curr_dim < nbr_cnstrts ; curr_dim++ )
		{
			current_constraint_value[curr_dim] += getWeight(num_item, curr_dim) ;
		}
	}
}



/*!
 *	\brief Retire un objet à un individu
 *
 *	Permet la gestion de  la suppression d'un objet à un individu et par cela la mise à jour des objectifs et contraintes.
 *
 *	\param num_item numéro de l'item à supprimer
 *
 */
void Individu::removeObject(int num_item)
{
	/* Si l'object est pris on le supprime et on met à jour objectif et contrainte */
	if( picked_objects[num_item] )
	{
		picked_objects[num_item] = 0 ;
		/* Mise à jour des objectifs */
		for( int curr_obj = 0 ; curr_obj < nbr_objs ; curr_obj++ )
		{
			current_objective_value[curr_obj] += -1*getProfit(num_item, curr_obj) ;
		}
		/* Mise à jour des contraintes */
		for( int curr_dim = 0 ; curr_dim < nbr_cnstrts ; curr_dim++ )
		{
			current_constraint_value[curr_dim] += -1*getWeight(num_item, curr_dim) ;
		}
	}
}



/**
 * \brief       Choisi aléatoirement une solution faisable.
 * \details    Les items sont sélmectionnés aléatoirement, et la solution est ensuite rendu faisable par
 */
void Individu::randomize()
{
	isGRASP = 0 ;
	/* Pour chaque item, on tire aléatoirement entre 0 et 1 */
	for(int curr_item = 0; curr_item < nbr_uavs; ++curr_item)
	{
		picked_objects[curr_item] = (rand()%2);
	}

	/* On rend valide s'il ne l'est pas l'individu */
	makeFeasible();
}



/*!
 *	\brief Permet de créer les fitness des items suivant une famille de direction
 *
 *	Cette fonction permet la création des fitness de tous les items suivant les directions données et de créer la première RCList correspondante à partir du coefficient glouton/aléatoire alpha.
 *
 *	\param Coef coefficient des directions leur permettant d'avoir un "poids"
 *	\param alpha coefficient glouton/aléatoire permettant de créer la première RCList.
 *
 */
void Individu::createGRASPFIT(int* Coef, double alpha)
{

	//cout << " Start GRASPFIT " << endl ;
	double fitness[nbr_uavs] ;
	int fitweight[nbr_uavs];
	RCList.clear() ;
	GRASPFIT = vector<double>() ;

	/* Calcul des valeurs GRASP de chaque item en fonction des coefficients de direction */
	for(int curr_item = 0; curr_item < nbr_uavs ; ++curr_item)
	{
		fitness[curr_item] = 0 ;
		fitweight[curr_item] = 0 ;

		/* Calcul de la valeur de profit de l'item selon les coefficients de direction */
		for( int curr_coef = 0; curr_coef < nbr_objs ; ++curr_coef)
		{
			fitness[curr_item] += getProfit(curr_item,curr_coef) * Coef[curr_coef];
			// somme coef = 100
		}
		/* Calcul du poid de l'item selon les coefficients de direction */
		for( int curr_weight = 0; curr_weight < nbr_cnstrts ; ++ curr_weight )
		{
			fitweight[curr_item] += getWeight(curr_item, curr_weight);
		}
		/* Calcul du fitness final de l'item */
		fitness[curr_item] = fitness[curr_item] / fitweight[curr_item];
	}

	/* Report des valeurs de fitness de chaque item dans GRASPFIT */
	for( int i = 0 ; i < nbr_uavs ; ++i )
	{
		GRASPFIT.push_back( fitness[i] );
	}

	int cmin = INT_MAX ;
	int cmax = 0 ;

	/* Recherche des fitness min et max */
	for( int i = 0; i < nbr_uavs ; ++i )
	{
			if( GRASPFIT[i] > cmax )
			{
				cmax = GRASPFIT[i] ;
			}
			else if( GRASPFIT[i] < cmin )
			{
				cmin = GRASPFIT[i] ;
			}
	}

	/* Création de la première RCList GRASP en fonction du coefficient glouton aléatoire alpha */
	for( int i = 0; i < nbr_uavs ; ++i )
	{
		if( GRASPFIT[i] >= ( cmin + alpha*(cmax-cmin) ) )
		{
			RCList.push_back(i) ;
		}
	}

//cout << " End GRASPFIT " << endl ;
}



/**
 * \brief      Algorithme Glouton
 * \details    Algorithme Glouton prenant en priorité les éléments dont la fonction d'utilité est la plus forte
 *
 * \param    *Coef         	Tableau des Coefficients à utiliser pour le calcul de la fonction d'utilité.
 *
 */
void Individu::GRASP( double alpha )
{
	isGRASP = 1 ;

    vector<int> RCList2 ;

    int rnd ;

    /* Rajout d'un item de la RCList */
    rnd = rand() % RCList.size() ;
	this->addObject( RCList[rnd] ) ;

    /* Tant que la solution est faisable on rajoute des items des RCLists successives */
	while(this->isFeasible())
    {

    	RCList2.clear() ;

    	int cmin = INT_MAX ;
    	int cmax = 0 ;

    	/* Recherche des fitness min et max en fonction des items non pris */
    	for( int i = 0; i < nbr_uavs ; ++i )
    	{
    		if( picked_objects[i] == 0)
    		{
    			if( GRASPFIT[i] > cmax )
    			{
    				cmax = GRASPFIT[i] ;
    			}
    			else if( GRASPFIT[i] < cmin )
    			{
    				cmin = GRASPFIT[i] ;
    			}
    		}
    	}

    	/* Calcul des nouvelles RCList en fonction des items non pris */
    	for( int i = 0; i < nbr_uavs ; ++i )
    	{
    		if( GRASPFIT[i] >= ( cmin + alpha*(cmax-cmin) ) && picked_objects[i]==0 )
    		{
    			RCList2.push_back(i) ;
    		}
    	}

        /* rajout aléatoire d'un objet */
        rnd = rand() % RCList2.size() ;

        this->addObject( RCList2[rnd] ) ;
    }

	/* Suppression du dernier item rendant l'individu non valide */
    this->removeObject( RCList2[rnd] );
	/* Intensification de l'individu via restart2 */
	this->restart2();
}



/**
 * \brief       Affiche les caractéristiques de la solution dans le terminal
 */
void Individu::printToScreen()
{
	cout<<"Items: "<<endl;
	for(int curr_item = 0; curr_item < nbr_uavs; ++curr_item)
	{
		cout<<"\tItem nb "<<curr_item<<" : "<<picked_objects[curr_item]<<endl;
	}

	cout<<"Constraints: "<<endl;
	for(int curr_constraint = 0; curr_constraint < nbr_cnstrts; ++curr_constraint)
	{
		cout<<"\tConstraint nb "<<curr_constraint<<" : "<<getContraintValue(curr_constraint)<<" out of "<<capacities[curr_constraint]<<endl;
	}
	cout<<"Objectives: "<<endl;
	for(int curr_fitness = 0; curr_fitness < nbr_objs; ++curr_fitness)
	{
		cout<<"\tObjective nb "<<curr_fitness<<" : "<<getObjectiveValue(curr_fitness)<<endl;
	}
}



/*!
 *	\brief Retourne la génétique de l'individu
 *
 *	Permet la création d'une chaîne de caractère correspondant à la génétique de l'individu.
 *
 *	\return un string correspondant à la génétique de l'individu.
 *
 */
string Individu::toString()
{

	string toReturn = "" ;

	/* Lecture du code génétique et écriture de celui-ci */
	for(int curr_item = 0; curr_item < nbr_uavs; ++curr_item)
	{
		if( picked_objects[curr_item] )
			toReturn += "1" ;
		else
			toReturn += "0" ;
	}

	return toReturn ;

}



/**
 * \brief       Fonction de test
 * \details    Teste l'initialisation du problème.
 */
void Individu::testData()
{
	cout<<"Capacité 0 = "<< capacities[0]<<endl;
	cout<<"Profit 0 0 = "<< getProfit(0,0)<<endl;
	cout<<"weight 0 0 = "<< getWeight(0,0)<<endl;
	cout<<"Profit 0 249 = "<< getProfit(249,0)<<endl;
	cout<<"weight 0 249 = "<< getWeight(249,0)<<endl;


	cout<<"Capacité 0 = "<< capacities[1]<<endl;
	cout<<"Profit 0 0 = "<< getProfit(0,1)<<endl;
	cout<<"weight 0 0 = "<< getWeight(0,1)<<endl;
	cout<<"Profit 0 249 = "<< getProfit(249,1)<<endl;
	cout<<"weight 0 249 = "<< getWeight(249,1)<<endl;
}
