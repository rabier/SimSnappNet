/*
 *  simulation.cpp
 *  SingleSiteSorter
 *
 *  Created by David Bryant on 29/08/08.
 *  Copyright 2008 David Bryant
 *  Modified by Charles-Elie Rabier to handle networks 
 */

#include "sortingSimulation.h"
//#include "UsefulForReticulations.h"
string g_simtree;

using namespace Phylib;
/**
 Copies the tree structure and initialises the fields for gene tree simulation.
 
 rate is the mutation rate, needed to convert theta values into gamma values
 **/
phylo<simNodeData> initialiseSimTree(phylo<basic_newick>& tree, double rate) {
	typedef phylo<simNodeData>::iterator S_ITER;	
	
	phylo<simNodeData> simTree;
	Phylib::copy(tree,simTree);
	
	//Convert the theta values into gamma values
	// The expected divergence between two individuals is 
	// theta = 2*rate/gamma.
	for(S_ITER s = simTree.leftmost_leaf();!s.null(); s = s.next_post()) {
		s->gamma = 2.0*rate/(s->theta);
		s->numberCoalescences = 0;
		s->numberCoalescencesTotal = 0;
	}
	
	
	return simTree;
}


// Function initialiseSimTreeNetwork is the same as initialiseSimTree for Network now !!!
// Transform tree into network structure and initialises the fields for gene tree simulation.

phylo<simNodeData> initialiseSimTreeNetwork(phylo<basic_newick>& tree, double rate, const vector<int>& retic_id, const vector<double>& retic_proba,list<phylo<simNodeData>::iterator > &myListNodesNetwork) {
	typedef phylo<simNodeData>::iterator S_ITER;
	
	phylo<simNodeData> simTree;
	Phylib::copy(tree,simTree);

	//Look for reticulation nodes ... the ghost will be the leaf with the #, the true will be the internal node with the #
	phylo<simNodeData>::iterator  GhostReticulate;
	phylo<simNodeData>::iterator  TrueReticulate;

	int MyBoolGhostRetic=0; 

	for(phylo<simNodeData>::iterator  MyPointer= simTree.leftmost_leaf();!MyPointer.null();MyPointer=MyPointer.next_post()) {	  
	    
	  MyBoolGhostRetic=0;
	  for(int i=0; i<retic_id.size(); i++){
	    if ((MyPointer->id==retic_id[i]) && (MyPointer.leaf())){
	       MyBoolGhostRetic=1;
	      break;}
	  }

	  if (MyBoolGhostRetic==0) {
	    myListNodesNetwork.push_back(MyPointer);
	  }
	}
	// we now have the list of nodes of interest

	std::list<phylo<simNodeData>::iterator >::iterator it=myListNodesNetwork.begin();
	//create a list of reticulation nodes (all the ghost and true retic nodes)
	list<phylo<simNodeData>::iterator >  myListNodesRetic;

	for(int i=0; i<retic_id.size(); i++){
	  // deal with ith reticulation node  
	  // look for true retic and ghost retic
	  for(phylo<simNodeData>::iterator  MyPointer= simTree.leftmost_leaf();!MyPointer.null();MyPointer=MyPointer.next_post()) {
	    if (MyPointer->id==retic_id[i]){
		myListNodesRetic.push_back(MyPointer);
	      }
	  }
	}
	  
	// the list size is myListNodesRetic.size()

	int i=0;// useful for initializing reticulation probabilities
	//regraft each true reticulation node to his second parent and handle all the pointer changes.

	//we use the fact that each true retic node is immediately followed by his ghost analogue in myListNodesRetic
	while(myListNodesRetic.size()>0) {	  
	  //take the first two elements
	  it=myListNodesRetic.begin();
	  TrueReticulate=*it;
	  //initialize the probability left	  
	  TrueReticulate.initialise_proba_left(retic_proba[i]);
	  advance(it,1);
	  GhostReticulate=*it;	  	  	 
	  myListNodesRetic.pop_front();
	  myListNodesRetic.pop_front();
	  simTree.graft_reticulate_node(GhostReticulate.par(), TrueReticulate);
	  i++;
	}

	GhostReticulate.set_null();
	TrueReticulate.set_null();
	//Be careful I have not handled numberCoalescences and numberCoalescencesTotal in the code
	// FIX IT Later, but we do not care at this time

	phylo<simNodeData>::iterator  s;
	it=myListNodesNetwork.begin();
	for(int i=0; i<myListNodesNetwork.size(); i++){
	  s=*it;	  
	  s->gamma = 2.0*rate/(s->theta);
	  s->numberCoalescences = 0;
	  s->numberCoalescencesTotal = 0;

	  if (s->theta_right!=0){
	    // handle reticulation node
	    // intialize gamma for branch on the right side of retic node
	    s->gamma_right = 2.0*rate/(s->theta_right);	    	    
	  }
	  advance(it,1);
	}

	return simTree;
}



/**
 Simulates a gene tree. 
 We use GIllespie's algorithm up each branch of the species tree to simulate the coalescent.
 **/
phylo<geneTreeNode> simulateGeneTree(phylo<simNodeData>& speciesTree, const vector<uint>& sampleSizes) {
	typedef phylo<simNodeData>::iterator S_ITER;
	typedef phylo<geneTreeNode>::iterator G_ITER;
	typedef list<phylo<geneTreeNode> >::iterator L_ITER;
	
	uint nleaves = 0;
	
	for(S_ITER s = speciesTree.leftmost_leaf(); !s.null(); s = s.next_post()) {
		if (s.leaf()) {
			//Create new lineages at the base of the branch, one for each individual in the sample.
			int count=sampleSizes[s->id];
			s->lineages.clear();
			for(int i=0;i<count;i++) {
				phylo<geneTreeNode> leafNode;
				G_ITER newLeaf = leafNode.insert_child(leafNode.header());
				newLeaf->id = nleaves;
				newLeaf->length = 0.0;
				newLeaf->species = s->id;
				nleaves++;
				s->lineages.push_back(leafNode);
			}
		

		} else {
			//Splice the lineages from the children to the bottom of this branch.
			s->lineages.clear();
			for(S_ITER c = s.left(); !c.null(); c=c.right()) 
				s->lineages.splice(s->lineages.end(),c->lineages);
		}
		
		//cerr<<"This gamma = "<<s->gamma<<endl;
		
		//Evolve the lineages up the branch.
		double height_in_branch = 0.0;
		int k = s->lineages.size();
		double branch_length = s->length;
		
		for( ; ; ) {
			
			
			if (k==1) { //All lineages coalesced
				if(s.root()) {//If we're at the root of the species tree, and coalesced to one lineage, we're done.
					phylo<geneTreeNode> geneTree;
					geneTree.swap(s->lineages.front());
					s->lineages.clear();
					return geneTree;
				}
				else {
					//Skip straight to the top of the branch
					for(L_ITER node = s->lineages.begin();node!=s->lineages.end();node++) 
						(*node).root()->length+=branch_length - height_in_branch;
					break;
				}	
			}
			
			//Waiting time until next coalescent event.
			
			double wait = random_exp(2.0 / ( (double)k*(k-1.0) * s->gamma));
			
			if (!s.root() && height_in_branch+wait>=branch_length) {
				//Reached the top of the branch
				for(L_ITER node = s->lineages.begin();node!=s->lineages.end();node++) 
					(*node).root()->length+=branch_length - height_in_branch;
				break;
			}
			
			s->numberCoalescences++; 
			for(L_ITER node = s->lineages.begin();node!=s->lineages.end();node++) {
				(*node).root()->length+=wait;
			}
			
			height_in_branch+=wait;
			
			
			
			
			//Choose a random pair (a,b),  a<b, to coalesce.
			int a = random_num(k);
			int b = random_num(k-1);
			if (b>=a) 
				b++;
			else { //Swap a and b so that a<b.
				int c=a; a=b; b = c;		
			}	
			
			//Update the gene tree and list of lineages.
			phylo<geneTreeNode> newNode;
			newNode.insert_child(newNode.header());
			newNode.root()->length = 0.0;
			
			L_ITER ptr = s->lineages.begin();
			std::advance(ptr,a);
			phylo<geneTreeNode>& child_a = *ptr;
			std::advance(ptr,b-a);
			phylo<geneTreeNode>& child_b = *ptr;
			
			newNode.graft_child(newNode.root(),child_a);
			newNode.graft_child(newNode.root(),child_b);
			ptr = s->lineages.begin();
			advance(ptr,b);
			s->lineages.erase(ptr);
			ptr = s->lineages.begin();
			advance(ptr,a);
			s->lineages.erase(ptr);
			s->lineages.push_back(newNode);
			k--;
			
			//TODO: Clean up this list pointer rubbish.
		}
	}
	throw PhylibException("Error in gene tree simulation");
	phylo<geneTreeNode> nullTree;
	return nullTree;
}



void simulateSubGeneTreeWithinSingleBranch(phylo<simNodeData>::iterator& s) {

  typedef list<phylo<geneTreeNode> >::iterator L_ITER;
  //Evolve the lineages up the branch.
  double height_in_branch = 0.0;
  int k = s->lineages.size();
  double branch_length = s->length;

  for( ; ; ) {
			
    if (k==0) { //no lineages quit
      break;}
		      	
    if (k==1) { //All lineages coalesced
			
      //Skip straight to the top of the branch
      for(L_ITER node = s->lineages.begin();node!=s->lineages.end();node++) 
	(*node).root()->length+=branch_length - height_in_branch;
      break;	
    }
			
    //Waiting time until next coalescent event.			
    double wait = random_exp(2.0 / ( (double)k*(k-1.0) * s->gamma));

    if ((height_in_branch+wait)>=branch_length) {
      //Reached the top of the branch
      for(L_ITER node = s->lineages.begin();node!=s->lineages.end();node++) {
	(*node).root()->length+=branch_length - height_in_branch;
      }
      break;
    }
			
    s->numberCoalescences++; 
    for(L_ITER node = s->lineages.begin();node!=s->lineages.end();node++) {
      (*node).root()->length+=wait;
    }
			
    height_in_branch+=wait;    
			
    //Choose a random pair (a,b),  a<b, to coalesce.
    int a = random_num(k);
    int b = random_num(k-1);
    if (b>=a) 
      b++;
    else { //Swap a and b so that a<b.
      int c=a; a=b; b = c;		
    }	
			
    //Update the gene tree and list of lineages.
    phylo<geneTreeNode> newNode;
    newNode.insert_child(newNode.header());
    newNode.root()->length = 0.0;
    
    L_ITER ptr = s->lineages.begin();
    std::advance(ptr,a);
    phylo<geneTreeNode>& child_a = *ptr;
    std::advance(ptr,b-a);
    phylo<geneTreeNode>& child_b = *ptr;
    
    newNode.graft_child(newNode.root(),child_a);
    newNode.graft_child(newNode.root(),child_b);
    ptr = s->lineages.begin();
    advance(ptr,b);
    s->lineages.erase(ptr);
    ptr = s->lineages.begin(); 
    advance(ptr,a);
    s->lineages.erase(ptr);
    s->lineages.push_back(newNode);
    k--;
    // there are  s->lineages.size() lineages going above this branch !!!
    //TODO: Clean up this list pointer rubbish.
 			
		}
}




//handle reticulation node
void simulateSubGeneTreeAboveReticulateNode(phylo<simNodeData>::iterator& s) {
  typedef list<phylo<geneTreeNode> >::iterator L_ITER;
	
  //cerr<<"This gamma = "<<s->gamma<<endl;
		
  // Handle the left branch

  //Evolve the lineages up the branch.
  double height_in_branch = 0.0;
  int k = s->lineages.size();
  double branch_length = s->length;		

  //have to handle k lineages on the branch of size branch_length
  for( ; ; ) {
 					   
    if (k==0) {
      // no lineage on the left branch above reticulation node    
      break;}
    	
    if (k==1) {
      // only one lineage remaining on the left branch above reticulation node
      //Skip straight to the top of the branch
      for(L_ITER node = s->lineages.begin();node!=s->lineages.end();node++) 
	(*node).root()->length+=branch_length - height_in_branch;
      break;
    }	
 
    // more than 1 lineage remaining on the left branch   
    //Waiting time until next coalescent event.

    double wait = random_exp(2.0 / ( (double)k*(k-1.0) * s->gamma));
    // recall that s->gamma is the gamma on the left side of reticulation node

    if ( (height_in_branch+wait)>=branch_length) {
      //Reached the top of the branch
      for(L_ITER node = s->lineages.begin();node!=s->lineages.end();node++) 
	(*node).root()->length+=branch_length - height_in_branch;
      break;
    }
			
    s->numberCoalescences++; 
    for(L_ITER node = s->lineages.begin();node!=s->lineages.end();node++) {
      (*node).root()->length+=wait;
    }

    height_in_branch+=wait;
    // there is one coalescent event
    //Choose a random pair (a,b),  a<b, to coalesce.
    int a = random_num(k);
    int b = random_num(k-1);
    if (b>=a) 
      b++;
    else { //Swap a and b so that a<b.
      int c=a; a=b; b = c;		
    }	
			
    //Update the gene tree and list of lineages.
    phylo<geneTreeNode> newNode;
    newNode.insert_child(newNode.header());
    newNode.root()->length = 0.0;

    L_ITER ptr = s->lineages.begin();
    std::advance(ptr,a);
    phylo<geneTreeNode>& child_a = *ptr;
    std::advance(ptr,b-a);
    phylo<geneTreeNode>& child_b = *ptr;

    newNode.graft_child(newNode.root(),child_a);
    newNode.graft_child(newNode.root(),child_b);
    ptr = s->lineages.begin();
    advance(ptr,b);
    s->lineages.erase(ptr);
    ptr = s->lineages.begin(); 
    advance(ptr,a);
    s->lineages.erase(ptr);
    s->lineages.push_back(newNode);
    k--;}
			
  //TODO: Clean up this list pointer rubbish.
  // there are s->lineages.size() that will go above the left branch of reticulation nod

  // Let us now handle the right branch of reticulation node
  //Evolve the lineages up the branch.
  double height_in_branch_right = 0.0;
  int k_right = s->lineages_right.size();
  double branch_length_right = s->length_right;
  // there are k_right lineages to handle on the right on a branch of size branch_length_right

  for( ; ; ) {
     
    if (k_right==0) {
      // there are no lineages going on the rightside
      // case where no lineage within a species
      break;}
    if (k_right==1) {
      // only one lineage remaining
      //Skip straight to the top of the branch
      for(L_ITER node = s->lineages_right.begin();node!=s->lineages_right.end();node++) 
	(*node).root()->length+=branch_length_right - height_in_branch_right;
      break;
    }	
 
    // more than 1 lineage remaining 
    //Waiting time until next coalescent event.
    double wait_right = random_exp(2.0 / ( (double)k_right*(k_right-1.0) * s->gamma_right));
 
    if ((height_in_branch_right+wait_right) >= branch_length_right) {
      //Reached the top of the branch
      for(L_ITER node = s->lineages_right.begin();node!=s->lineages_right.end();node++) 
	(*node).root()->length+=branch_length_right - height_in_branch_right;
      break;
    }
    
    s->numberCoalescences_right++; 
    for(L_ITER node = s->lineages_right.begin();node!=s->lineages_right.end();node++) {
      (*node).root()->length+=wait_right;
    }
    
    height_in_branch_right+=wait_right;
    // one coalescent event
    //Choose a random pair (a,b),  a<b, to coalesce.
    int a = random_num(k_right);
    int b = random_num(k_right-1);
    if (b>=a) 
      b++;
    else { //Swap a and b so that a<b.
      int c=a; a=b; b = c;		
    }	
			
    //Update the gene tree and list of lineages.
    phylo<geneTreeNode> newNode;
    newNode.insert_child(newNode.header());
    newNode.root()->length = 0.0;
    
    L_ITER ptr = s->lineages_right.begin();
    std::advance(ptr,a);
    phylo<geneTreeNode>& child_a = *ptr;
    std::advance(ptr,b-a);
    phylo<geneTreeNode>& child_b = *ptr;
    
    newNode.graft_child(newNode.root(),child_a);
    newNode.graft_child(newNode.root(),child_b);
    ptr = s->lineages_right.begin();
    advance(ptr,b);
    s->lineages_right.erase(ptr);
    ptr = s->lineages_right.begin(); 
    advance(ptr,a);
    s->lineages_right.erase(ptr);
    s->lineages_right.push_back(newNode);
    k_right--;
  }
  // there are s->lineages_right.size() that will go above the right branch of reticulation node
  // Recall that there are s->lineages.size() that will go above the left branch of reticulation node

}



/**
 Simulates a gene tree coming from a network  using the list of nodes of the network
 **/
phylo<geneTreeNode> simulateGeneTreeFromNetwork(phylo<simNodeData>& speciesTree, const vector<uint>& sampleSizes, list<phylo<simNodeData>::iterator > &myListNodesNetwork) {
  typedef phylo<simNodeData>::iterator S_ITER;
  typedef phylo<geneTreeNode>::iterator G_ITER;
  typedef list<phylo<geneTreeNode> >::iterator L_ITER;
	
  phylo<simNodeData>::iterator speciesTreeBis;
  uint nleaves = 0;    

  S_ITER s = speciesTree.leftmost_leaf();
  std::list<phylo<simNodeData>::iterator >::iterator myFirstSNode;
  bool dealWithMe;

  while(myListNodesNetwork.size()>0) {
    dealWithMe=0;    

    while(dealWithMe==0) {
      myFirstSNode=myListNodesNetwork.begin();    
      s=*myFirstSNode;
      //is the first Node of the List
    
      // check if the children have been treated 
      if (s.leaf()){
	dealWithMe=1;
      } else if (s.reticulate_node()){
	//s is a reticulation
	//it has one child
	//be careful s.numChildren() can be greater than one because of half sibs
	// that is why I am not using  s.numChildren() here	
	dealWithMe=s.left()->hasBeenTreated;
      }  else if (s.numChildren()>1){ 
	// s is not a reticulation and has at least 2 children
	// note that s.numChildren() can be equal to 3 but we take just the first two, because of half sibs
	dealWithMe=(s.left()->hasBeenTreated && s.left().right()->hasBeenTreated);
      }

      //take next node of the list when dealWithMe=0
      if (dealWithMe==0){
	// put the node at the end of the list
	myListNodesNetwork.push_back(s);	
	myListNodesNetwork.pop_front();
      }
    }


    if (s.leaf()) {
      //Create new lineages at the base of the branch, one for each individual in the sample.
      int count=sampleSizes[s->id];
      s->lineages.clear();
      for(int i=0;i<count;i++) {
	phylo<geneTreeNode> leafNode;
	G_ITER newLeaf = leafNode.insert_child(leafNode.header());
	newLeaf->id = nleaves;
	newLeaf->length = 0.0;
	newLeaf->species = s->id;
	nleaves++;
	s->lineages.push_back(leafNode);
      }

      // the sampleSizes[s->id] lineages for species s->id have been created
      
      if (count>0){
	// simulate coalescent process within the branch
	speciesTreeBis=s;
       	simulateSubGeneTreeWithinSingleBranch(speciesTreeBis);}
      // after the coalescent process, s->lineages.size() will be remaining      

    } else if (s.reticulate_node())  {
      // handle a reticulation node
      s->lineages.clear();
      s->lineages_right.clear();
      // have to handle s.left()->lineages.size() lineages 
      // species associated to this reticulation node is  s->id
      // reticul probability is s.retic_proba_left()
      L_ITER myPtr = s.left()->lineages.begin();
      // have to handle s.left()->lineages.size() lineages 
      // species associated to this reticulation node is  s->id
      // reticul probability is s.retic_proba_left()

      for(int i=0;i<s.left()->lineages.size();i++) {		   
	if (randu()<s.retic_proba_left()){// go left
	  //the lineage goes left
	  s->lineages.push_back(*myPtr);}
	else{ //the lineage goes  right
	  s->lineages_right.push_back(*myPtr);}	  
	advance(myPtr,1);
      }

      //Have now on the left s->lineages.size()
      //Have on the right s->lineages_right.size()
      speciesTreeBis=s;
      simulateSubGeneTreeAboveReticulateNode(speciesTreeBis);
      //it will simulate coalescent processes above the reticulation node (left and right branches)
      
    }  else if (s.left().reticulate_node()){ 
      // s is a parent of reticulation node, whose left child is a reticulation node    
      // We are on the internal branch going up from s->id
      // First handle right node of reticulation node, ie s
      //Splice the lineages from the children to the bottom of this branch.
       s->lineages.clear();
       S_ITER c = s.left();       
       s->lineages.splice(s->lineages.end(),c->lineages_right);
       c = c.right();      
       s->lineages.splice(s->lineages.end(),c->lineages); 
       speciesTreeBis=s;
       //simulate coalescent process in internal branch going up from s->id
       simulateSubGeneTreeWithinSingleBranch(speciesTreeBis);
     
    }else if (!s.root())  { 
      // we are on an internal branch going up from  s->id
      // s is an internal node, and not a reticulation node, and not a node whose left child is a reticulation node 
      //Splice the lineages from the children to the bottom of this branch.
      s->lineages.clear();
      // s.left()->lineages.size() lineages are coming from the left
      // s.left().right()->lineages.size() are coming from the right
      s->lineages.splice(s->lineages.end(),s.left()->lineages);
      s->lineages.splice(s->lineages.end(),s.left().right()->lineages);
      speciesTreeBis=s;
      //simulate coalescence process within the branch matching this species
      simulateSubGeneTreeWithinSingleBranch(speciesTreeBis);

    }else{
      // s is the root
      // we deal with the branch above the root			  
    	//Splice the lineages from the children to the bottom of this branch.
	s->lineages.clear();
	for(S_ITER c = s.left(); !c.null(); c=c.right()){	  
	  s->lineages.splice(s->lineages.end(),c->lineages);
	}

	//have to handle s->lineages.size() lineages
	//cerr<<"This gamma = "<<s->gamma<<endl;
		
	//Evolve the lineages up the branch.
	double height_in_branch = 0.0;
	int k = s->lineages.size();
	double branch_length = s->length;
	
		for( ; ; ) {
			
		      	
		       if (k==1) {
			 // we're at the root of the species tree, and/All lineages coalesced, we're done.
					phylo<geneTreeNode> geneTree;
					geneTree.swap(s->lineages.front());		
					speciesTreeBis=s;
					s->lineages.clear();				       
					return geneTree;
			       				
			}
			
			//Waiting time until next coalescent event.
			
			double wait = random_exp(2.0 / ( (double)k*(k-1.0) * s->gamma));					
			
			s->numberCoalescences++; 
			for(L_ITER node = s->lineages.begin();node!=s->lineages.end();node++) {
				(*node).root()->length+=wait;
			}
			
			height_in_branch+=wait;
			
			
			//Choose a random pair (a,b),  a<b, to coalesce.
			int a = random_num(k);
			int b = random_num(k-1);
			if (b>=a) 
				b++;
			else { //Swap a and b so that a<b.
				int c=a; a=b; b = c;		
			}	
			
			//Update the gene tree and list of lineages.
			phylo<geneTreeNode> newNode;
			newNode.insert_child(newNode.header());
			newNode.root()->length = 0.0;
			
			L_ITER ptr = s->lineages.begin();
			std::advance(ptr,a);
			phylo<geneTreeNode>& child_a = *ptr;
			std::advance(ptr,b-a);
			phylo<geneTreeNode>& child_b = *ptr;
			
			newNode.graft_child(newNode.root(),child_a);
			newNode.graft_child(newNode.root(),child_b);
			ptr = s->lineages.begin();
			advance(ptr,b);
			s->lineages.erase(ptr);
			ptr = s->lineages.begin(); 
			advance(ptr,a);
			s->lineages.erase(ptr);
			s->lineages.push_back(newNode);
			k--;
			
			//TODO: Clean up this list pointer rubbish.
			
		}
    }

    s->hasBeenTreated=1;
    myListNodesNetwork.pop_front();
  }
  //end while
  
  throw PhylibException("Error in gene tree simulation");
  phylo<geneTreeNode> nullTree;
  return nullTree;
	

}






typedef phylo<geneTreeNode>::iterator G_ITER;

/**
 Pre-compute the mutation probabilities using the 2-state model
 **/
static void computeTransitionProbs(phylo<geneTreeNode>& G, double u, double v) {
	double pi_0 = u/(u+v), pi_1 = 1.0-pi_0;  //0 = green, 1 = red.
	
	//	printf("je suis dans compute transition probs\n");
	for(G_ITER p = G.leftmost_leaf();p!=G.root();p=p.next_post()) {
		double t = p->length;
		double x = 1.0 - std::exp(-(u+v)*t);
		p->P[1][0] = pi_0*x;
		p->P[1][1] = 1.0 - pi_0*x;
		p->P[0][1] = pi_1*x;
		p->P[0][0] = 1.0 - pi_1*x;
	      
		/*
		printf("I am node %d\n", p->id);
		printf("Transition probability from green to green is %f\n", p->P[0][0]);
		printf("Transition probability from green to red is %f\n", p->P[0][1]);
		printf("Transition probability from red to red is %f\n", p->P[1][1]);
		printf("Transition probability from red to green is %f\n", p->P[1][0]);
		*/
	}
}

/**
 Simulate a two state character on the tree
 **/
static void simulateCharacter(phylo<geneTreeNode>& G, double pi_0) {
	G.root()->state = (randu()<pi_0)?0:1;
	
	for(G_ITER p = G.root().next_pre();!p.null();p=p.next_pre()) {
		uint i = p.par()->state;
		p->state = (randu()<p->P[i][0])?0:1;
		//	printf("voici mon etat %d\n",p->state);
	}
}

/**
 Check to see whether a character (given by allele counts) is constant or not.
 **/
static bool checkIfConstant(const vector<uint>& redCount, const vector<uint>& sampleSizes) {
	bool allGreen, allRed;
	allGreen = allRed = true;
	for(uint i=0;i<redCount.size();i++) {
		uint r = redCount[i];
		if (r>0)
			allGreen = false;
		if (r<sampleSizes[i])
			allRed = false;
		if (!allGreen && !allRed)
			break;
	}
	//cerr<<((allRed||allGreen)?"constant":"variable")<<endl;
	return (allRed || allGreen);
}


/**
 Simulate a gene tree and then simulate a single binary character on it. Returns allele counts for each species.
 Uses a rejection algorithm to simulate non-constant characters. The numberAttempts is the number of characters we had to simulate to a polymorphic one.
 **/
void simulateSingleSite(phylo<simNodeData>& speciesTree, double u, double v, const vector<uint>& sampleSizes, vector<uint>& redCounts, bool rejectConstant, bool onlyRootMutation, bool hasDominantMarkers, double& numberAttempts, bool outputTree, int site, list<phylo<simNodeData>::iterator >& myListNodesNetwork) {

        //CE : I assume that onlyRootMutation is set to false, ie default
        //since I did not handle the RC model
	
	//Now evolve the markers.
	double pi_0 = u/(u+v);
	uint nSpecies = sampleSizes.size();
	numberAttempts = 0;
	
	bool allConstant = true;
	
	vector<uint> effectiveSampleSizes(nSpecies);

	//need to copy myListNodesNetwork, otherwise it will get empty each time I generate a gene tree
	list<phylo<simNodeData>::iterator > myListNodesNetworkSaved;
	std::list<phylo<simNodeData>::iterator >::iterator it=myListNodesNetwork.begin();

		
	for(uint i=0;i<nSpecies;i++)
		if (hasDominantMarkers)
			effectiveSampleSizes[i] = 2*sampleSizes[i];
		else
			effectiveSampleSizes[i] = sampleSizes[i];
			
	do {
		numberAttempts++;		
		// we copy my  myListNodesNetwork into myListNodesNetworkSaved
		it=myListNodesNetwork.begin();
		for(int i=0; i<myListNodesNetwork.size(); i++){
		  myListNodesNetworkSaved.push_back(*it);
		  advance(it,1);
		}
		
		//First generate the gene tree
		phylo<geneTreeNode> G = simulateGeneTreeFromNetwork(speciesTree, effectiveSampleSizes, myListNodesNetworkSaved);
		
		computeTransitionProbs(G,u,v);		
		simulateCharacter(G, pi_0);
		redCounts.resize(nSpecies);
		std::fill(redCounts.begin(),redCounts.end(),0);
		vector<int> lastAllele(nSpecies);
		std::fill(lastAllele.begin(),lastAllele.end(),-1);		
		
	       		
		for(G_ITER p = G.root();!p.null();p=p.next_pre()) {
			if (p.leaf()) {
				if (hasDominantMarkers) {
					uint species = p->species;
					if (lastAllele[species]<0)
						lastAllele[species] = (int)p->state;
					else {
						if (lastAllele[species]==1 || p->state == 1)
							redCounts[p->species]++;
						lastAllele[species]=-1;
					}
				} else {
					if (p->state==1) {
						redCounts[p->species]++;
					}
				}
			}
		}
		
				
		
		if (rejectConstant) 
			allConstant = checkIfConstant(redCounts, sampleSizes);	
		
		
		G.clear();
		myListNodesNetworkSaved.clear();	
	} while(allConstant && rejectConstant);	


}

/**
 Simulate multiple sites.
 
 For each site, simulates a gene tree according to the multi-species network coalescent, then simulates a single binary character on that gene tree.
 If rejectConstant is set to true, the process will repeat at each site until a non-constant character is obtained. This could cause problems
 if the mutation rate is really low.
 **/

/**wrapper**/
void simulateMultipleSites(phylo<basic_newick>& tree, double u, double v, const vector<uint>& sampleSizes, int nSites, bool rejectConstant, bool onlyRootMutation, bool hasDominantMarkers, vector<vector<uint> >& redCounts, bool outputTrees,const vector<int>& retic_id,const vector<double>&  retic_proba) {
	double proportionConstant;
		
	simulateMultipleSites(tree,u,v,sampleSizes,nSites,rejectConstant, onlyRootMutation, hasDominantMarkers, redCounts,proportionConstant, outputTrees, retic_id, retic_proba);
	// recall that the size of retic_id is retic_id.size()
	// recall that the size of retic_proba is retic_proba.size()

}


void simulateMultipleSites(phylo<basic_newick>& tree, double u, double v, const vector<uint>& sampleSizes, int nSites, bool rejectConstant, bool onlyRootMutation, bool hasDominantMarkers, vector<vector<uint> >& redCounts, double& proportionConstant, bool outputTrees,const vector<int>&  retic_id, const vector<double>&  retic_proba) {
	typedef phylo<simNodeData>::iterator S_ITER;
	list<phylo<simNodeData>::iterator > myListNodesNetwork;
	
	redCounts.resize(nSites);
	double rate = 2.0*u*v/(u+v); 

	// Network initialization
	phylo<simNodeData> simTree = initialiseSimTreeNetwork(tree,rate,retic_id,retic_proba,myListNodesNetwork);	

	// for me : to see if it generates gene tree correctly
	//	phylo<geneTreeNode>  GUsingList=simulateGeneTreeFromNetwork(simTree, sampleSizes, myListNodesNetwork);
	//end for me
	// in fact genetrees are generated in function simulateSingleSite, see below
	//	uint numberAttempts;
	double numberAttempts;
	double numberAttemptsTotal;
	numberAttemptsTotal = 0;
	numberAttempts = 0;
	
	for(int i=0;i<nSites;i++) {
	 simulateSingleSite(simTree, u, v, sampleSizes, redCounts[i], rejectConstant, onlyRootMutation, hasDominantMarkers,numberAttempts, outputTrees,i+1,myListNodesNetwork);
	 numberAttemptsTotal+=numberAttempts;
	}

	// We are going to erase the network !!!
	phylo<simNodeData>::iterator TrueReticulate;

	// First we remove reticulation nodes
	std::list<phylo<simNodeData>::iterator >::iterator it=myListNodesNetwork.begin();

	for(int i=0; i< myListNodesNetwork.size(); i++){
	  	 
	  if ((*it).reticulate_node()){
	    // it is a reticulation node and not a leaf 
	    TrueReticulate=*it;
	    simTree.Remove_Reticulate_Node(TrueReticulate);	   
	  } 
	   advance(it,1);
	}
	// we have now removed reticulation nodes
	myListNodesNetwork.clear();

	// we can now erase the tree safetly since it is not a network anymore
       	simTree.clear();  

	if (rejectConstant)
		proportionConstant = 1.0-(double)nSites/(double)numberAttemptsTotal;
	else
		proportionConstant = -1.0;
	
	//cerr << "Number of constant sites " <<  (numberAttemptsTotal - nSites) << endl;
		
}

















