/*
 *  simulation.h
 *  SingleSiteSorter
 *
 *  Created by David Bryant on 29/08/08.
 *  Copyright 2008 David Bryant
 *
 */

#ifndef SORTING_SIMULATION_H
#define SORTING_SIMULATION_H

#include "phylib.h"

using namespace Phylib;

extern string g_simtree;

/**
 Node type for gene trees 
 **/
class geneTreeNode : public basic_newick {
public: 
        int state;
	int species;
	double height;
	double P[2][2]; //Transition probabilities.
};

/**
 Node type for species trees
 **/
class simNodeData : public basic_newick {
public:
        double theta;
	//MOI
	double theta_right;
	list< phylo<geneTreeNode> > lineages_right;
	double gamma_right;
	int numberCoalescences_right;
	//END MOI
	double gamma;	
	list< phylo<geneTreeNode> > lineages;
	int numberCoalescences;
	int numberCoalescencesTotal;
	bool hasBeenTreated;

	
 simNodeData() : basic_newick(), gamma(1.0) , gamma_right(1.0), numberCoalescences(0), numberCoalescences_right(0), hasBeenTreated(0) {}
	//FIX IT to intinialize gamma_right 
                simNodeData(int id, double len, double gamma_val): basic_newick(id, len), gamma(gamma_val), numberCoalescences(0), numberCoalescences_right(0) {}
	        simNodeData(const basic_newick& data) {

		basic_newick::copy(data);
		
		
		size_t equalPos = data.meta_data.find_first_of("=");
		if (equalPos!=string::npos) {
			if (data.meta_data.compare(0,5,"theta") == 0)
				theta = atof((data.meta_data.substr(equalPos+1)).c_str());
			else if (data.meta_data.compare(0,15,"coalescenceRate")==0)
				theta = 2.0 / atof((data.meta_data.substr(equalPos+1)).c_str());
			else	{		
				theta = atof((data.meta_data.substr(equalPos+1)).c_str());
				cerr<<"Error reading in tree"<<endl;
			}
		}
		else
			theta = atof(data.meta_data.c_str());
		
		//TODO: Make this more robust.
		
		gamma = 0.0;
		numberCoalescences = 0;
	}
		
	void copy(const simNodeData& data) {
		basic_newick::copy(data);
		gamma = data.gamma;
		theta = data.theta;
		gamma_right = data.gamma_right;
		theta_right = data.theta_right;
		numberCoalescences = data.numberCoalescences;
		
		lineages.clear();
		lineages.insert(lineages.begin(),data.lineages.begin(),data.lineages.end());
		lineages_right.clear();
		lineages_right.insert(lineages_right.begin(),data.lineages_right.begin(),data.lineages_right.end());
		
	}
	
	simNodeData& operator=(const simNodeData& data) {
		if (this!=&data) 
			copy(data);
		return *this;
	}
	
};

/**
 Copies the tree structure and initialises the fields for gene tree simulation.
 **/
phylo<simNodeData> initialiseSimTree(phylo<basic_newick>& tree);


/*
 Copies the tree structure and initialises the fields for gene tree simulation in presence of network */
phylo<simNodeData> initialiseSimTreeNetwork(phylo<basic_newick>& tree, double rate, const vector<int>& retic_id, const vector<double>& retic_proba, list<phylo<simNodeData>::iterator > &myListNodesNetwork);

/**
 Simulates a gene tree. 
 **/
phylo<geneTreeNode> simulateGeneTree(phylo<simNodeData>& speciesTree, const vector<uint>& sampleSizes);


/**
Simulate a gene tree from a network
**/
void simulateSubGeneTreeAboveReticulateNode(phylo<simNodeData>::iterator& s);
void simulateSubGeneTreeWithinSingleBranch(phylo<simNodeData>::iterator& s);
phylo<geneTreeNode> simulateGeneTreeFromNetwork(phylo<simNodeData>& speciesTree, const vector<uint>& sampleSizes, list<phylo<simNodeData>::iterator > &myListNodesNetwork);

/**
 Simulate a single SNP. if rejectConstant = true, a non-constant site (polymorphism) will be generated.
 **/



/**
 Simulate multiple unlinked sites.
 
 If rejectConstant is true, only non Constant sites are returned. We use a simple rejection algorithm (any better ideas?)
 **/

void simulateMultipleSites(phylo<basic_newick>& tree, double u, double v, const vector<uint>& sampleSizes, int nSites, bool rejectConstant, bool onlyRootMutation, bool hasDominantMarkers, vector<vector<uint> >& redCounts, bool outputTrees = false,const vector<int>& retic_id=vector<int>(),const vector<double>& retic_proba=vector<double>());
void simulateMultipleSites(phylo<basic_newick>& tree, double u, double v, const vector<uint>& sampleSizes, int nSites, bool rejectConstant, bool onlyRootMutation, bool hasDominantMarkers, vector<vector<uint> >& redCounts, double& proportionConstant, bool outputTrees = false,const vector<int>& retic_id=vector<int>(), const vector<double>& retic_proba=vector<double>());

#endif
