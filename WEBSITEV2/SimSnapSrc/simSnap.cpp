/*
 *  simSnap.cpp
 *  SimSnap
 *
 *  Created by David Bryant on 8/03/10.
 *  Copyright 2010 University of Auckland. All rights reserved.
 *
 *  Modified by Charles-Elie Rabier to handle networks
 */




#include <unistd.h>
#include "sortingSimulation.h"
#include "characterData.h"
//#include "posteriorCheck.h"

//#define DEBUG_2_SPECIES


void printUsage(ostream& os) {
	os<<"SimSnap\n\nSimulates SNPs on a species network.\n";
	os<<"Usage:\n\n\tSimSnap [-c]  [nsites] [seed] filename\n\n";
	os<<"\tFlags are:\n";
 
	os<<"\t-c \tInclude constant sites (default is to simulate only polymorphic sites)\n";
	//os<<"\t-r \tCondition on all mutation occurring at the root\n";

	//fixit later , ie parameters to fill xml file
	//os<<"\t-t \tOutput the gene trees used to generate each site\n";
	//os<<"\t-i \tStart chain at the values used for simulation\n";
	//os<<"\t\tLENGTH\tChain length\n";
	//os<<"\t\tPRIORBURN\tNumber of iterations with prior only\n";
	//os<<"\t\tALPHA BETA\tParameters of gamma prior\n";
	//os<<"\t\tLAMBDA\tParameter of Yule prior\n";
	//os<<"\t\tTREEWEIGHT\tWeight assigned to NodeSwapping operator\n";
	//os<<"\t\tFILEROOT\tName used for output files: .log, .trees\n\n";

	       
	os<<"INPUT FILE FORMAT:\n\n";
	os<<"<number of species>\n";
	os<<"<species1-name>\t<sample size species 1>\n";
	os<<"<species2-name>\t<sample size species 2>\n";
	os<<"...\n";
	os<<"<species_n-name>\t<sample size species n>\n";
	os<<"<mutation rate u> <mutation rate v>\n";
	os<<"<number of networks>\n";
	os<<"Networks in extended newick notation, with branch lengths. Theta values indicated in square brackets following node.\n";
	os<<"Be careful, the reticulation node is considered as a species, its name has to begin with #, followed by one letter, and the probability of going left. The sample size must be set to 0 for that node (see examples below!)\n";
	os<<"Optionall insert <scale> before the network to specify a branch length scaler.\n\n";
	os<<"EXAMPLE 5 SPECIES 1 RETICULATION\n\n";
        os<<"6\n";
        os<<"C 2\nR 2\nQ 2\nA 2\n#H0.3 0\nL 2\n";
        os<<"1.0 1.0\n1\n";
	os<<"<1.0>(C[0.005]:0.08,((R[0.005]:0.007,(Q[0.005]:0.004)#H0.3[0.005]:0.003)[0.005]:0.035,((A[0.005]:0.006,#H0.3[0.005]:0.002)[0.005]:0.016,L[0.005]:0.022)[0.005]:0.02)[0.005]:0.038)[0.005]:0.1;\n\n\n"<<endl;
	exit(1);
}


class ArgumentParser {
public: 
	bool outputXML;
	bool excludeConst;
	bool outputTrees;
	bool onlyRootMutation; 
	bool hasDominantMarkers;
	bool simulationOutput;
	bool initialiseAtTrue;
	bool snappDominant;
	bool snappNoMutation;
        bool useSNAPPMCMC;
	
	string inputfile;
	int nsites;
        int seedValue;
	ArgumentParser(int argc, char* argv[]) {


	        //CE : I did not handle hasDominantMarkers and snappDominant, and onlyRootMutation	  
                //CE : so leave those parameters to false
        
	        outputXML=true;
		excludeConst = true;
		outputTrees = false;
		simulationOutput = false;
		initialiseAtTrue = false;
		onlyRootMutation = false;
		hasDominantMarkers = false;
		snappDominant = false;
		snappNoMutation = false;
		useSNAPPMCMC = true;
		
		nsites = 0;
		seedValue=1;
		inputfile = "";
		
		if (argc<3)
			printUsage(cerr);
		
		int arg = 1;
		
		//First read in the flags.
		string flags = string(argv[arg]);
		if (flags[0]=='-') {
			if (flags.find_first_not_of("-nirdcstMDb")!=string::npos)
				printUsage(cerr);
			
			if (flags.find_first_of('n')!=string::npos)
				outputXML = false;
			if (flags.find_first_of('r')!=string::npos)
				onlyRootMutation = true;
			if (flags.find_first_of('c')!=string::npos)
				excludeConst = false;
			if (flags.find_first_of('t')!=string::npos)
				outputTrees = true;
			if (flags.find_first_of('s')!=string::npos)
				simulationOutput = true;
			if (flags.find_first_of('i')!=string::npos)
				initialiseAtTrue = true;
			if (flags.find_first_of('d')!=string::npos)
				hasDominantMarkers = true;
			if (flags.find_first_of('D')!=string::npos)
				snappDominant = true;
			if (flags.find_first_of('M')!=string::npos)
				snappNoMutation = true;
			if (flags.find_first_of('b')!=string::npos)
				useSNAPPMCMC = false;

			
			arg++;
			if (argc < 4)
				printUsage(cerr);
		}
		
		nsites = atoi(argv[arg]); 
		arg++;

		//CE have added a parameter for the seed
		seedValue = atoi(argv[arg]); 
		arg++;
		//END CE

		inputfile = string(argv[arg]);
	}
};





//CE  generates a xml for SnappNet
void output_xml_SnappNet(ostream& os, const vector<string>& taxa, phylo<basic_newick>& tree, const vector<uint>& sampleSizes, double u, double v, const vector<vector<uint> >&alleleCounts, bool simulationOutput, bool polymorphicOnly, bool initialiseFromTruth, bool snappDominant, bool snappNoMutation, bool useSNAPPMCMC, const string fileroot="test") {
	
       	os<<"<!-- Generated with SimSnapNet -->\n";
	os<<"<!-- -->\n";
	os<<"<!-- u = "<<u<<" v = "<<v<<"   -->\n";
	
	
	//TRanslate the meta-data in the tree. 
	phylo<basic_newick> ratetree;
	copy(tree,ratetree);
	for(phylo<basic_newick>::iterator p = ratetree.root();!p.null();p=p.next_pre()) {
		//cout<<(*p).meta_data<<endl;
		double rate;
		if ((*p).meta_data.compare(0,5,"theta") == 0) {
			size_t equalPos = (*p).meta_data.find_first_of("=");
			rate = 2.0/atof(((*p).meta_data.substr(equalPos+1)).c_str());
			stringstream s;
			s<<"rate="<<rate;;
			(*p).meta_data = s.str();
		}
	}
	
	
	//os<<"<?xml version='1.0' encoding='UTF-8'?>";
	os<<"\n";
	os<<"<beast namespace='beast.core:beast.core.util:beast.evolution.alignment:beast.evolution.tree.coalescent:beast.evolution.operators:beast.evolution.sitemodel:beast.evolution.substitutionmodel:beast.evolution.branchratemodel:beast.evolution.likelihood' version='2.4'>\n";


	os<<"<map name='Uniform'>beast.math.distributions.Uniform</map>\n";
	os<<"<map name='Exponential'>beast.math.distributions.Exponential</map>\n";
	os<<"<map name='LogNormal'>beast.math.distributions.LogNormalDistributionModel</map>\n";
	os<<"<map name='Normal'>beast.math.distributions.Normal</map>\n";
	os<<"<map name='Beta'>beast.math.distributions.Beta</map>\n";
	os<<"<map name='Gamma'>beast.math.distributions.Gamma</map>\n";
	os<<"<map name='LaplaceDistribution'>beast.math.distributions.LaplaceDistribution</map>\n";
	os<<"<map name='InverseGamma'>beast.math.distributions.InverseGamma</map>\n";
	os<<"<map name='OneOnX'>beast.math.distributions.OneOnX</map>\n";
	os<<"<map name='prior'>beast.math.distributions.Prior</map>\n";	


	//Compute largest sample size.
	int statecount = *max_element(sampleSizes.begin(),sampleSizes.end()); //if k is the largest sample size, states are 0...k, so k+1 states
	os<<"\t<data spec='snappNetForSimSnappNet.core.SnapData' id='snapalignment' name='alignment' dataType='integer'  statecount='" <<statecount + 1<<"'>\n";
	os<<"\t<rawdata  id='rawdata' spec='Alignment' dataType='integer'>\n";


	int taxaNumberWithoutRet =0;
	for(uint i=0;i<taxa.size();i++) {
	  //
	  char b = taxa[i].at(0);
	  //  printf("tt letter of taxon%d is  %c\n",i,b);
	  // do not consider species beginning with # since we use # for hybrids
	  if (b!='#'){
	  taxaNumberWithoutRet++;
	  os<<"\t\t<sequence taxon='"<<taxa[i]<<"1' totalcount='"<<sampleSizes[i] + 1<<"'>\n";
	  for(uint j=0;j<alleleCounts.size();j++)
	    os<<alleleCounts[j][i]<<",";
	  os<<"\n\t\t</sequence>\n";}
	}
	os<<"\t</rawdata>\n\n\n";
	

	for(uint i=0;i<taxa.size();i++) {
	  //
	  char b = taxa[i].at(0);
	  //  printf("tt letter of taxon%d is  %c\n",i,b);
	  // do not consider species beginning with # since we use # for hybrids
	  if (b!='#'){
	  taxaNumberWithoutRet++;
	  os<<"\t\t<taxonset id='"<<taxa[i]<<"' spec='TaxonSet'>\n"; 
	  os<<"\t\t\t<taxon id='"<<taxa[i]<<"1' spec='Taxon'/>\n";	  
	  os<<"\t\t</taxonset>\n";}
	}

	os<<"\t</data>\n\n\n";


	//////////////////////////////////////////////////////////////////////////////

	//CE init parameter
	os<<"\t<init spec='beast.util.TreeParser' id='newick:species' IsLabelledNewick='true' adjustTipHeights='false' newick='(((C:0.05,R:0.05):0.05,((A:0.05,L:0.05):0.025,Q:0.075):0.025):0)'/>\n\n\n";
	os<<"\t\t<run id='mcmc' spec='MCMC' chainLength='8000000' storeEvery='1000'>\n";

	os<<"\t\t\t<state id='state'>\n";
        os<<"\t\t\t<stateNode id='network:species' spec='snappNetForSimSnappNet.core.NetworkParser' tree='@newick:species'>\n";	 
	os<<"\t\t\t</stateNode>\n";
	os<<"\t\t\t<parameter id='originTime:species' lower='0.0' name='stateNode'>0.1</parameter>\n";
os<<"\t<!-- be careful, originTime:species is a scaling factor, if you do not want to scale anything it has to be set to the same value as origin (the height of) the network -->\n";

	os<<"\t\t\t<parameter id='netDivRate:species' lower='0.0' name='stateNode'>2.0</parameter>\n";
	os<<"\t\t\t<parameter id='turnOverRate:species' lower='0.0' upper='1.0' name='stateNode'>0.5</parameter>\n";
	//os<<"\t\t\t<parameter id='u' lower='0.0' upper='1.0' name='stateNode'>1.0</parameter>\n";	
	os<<"\t\t\t<parameter id='u' lower='0.0' upper='1.0' name='stateNode'>"
	  <<u<<"</parameter>\n";
	os<<"\t\t\t<parameter id='v' lower='0.0' upper='1.0' name='stateNode'>"
          <<v<<"</parameter>\n";
	//	os<<"\t\t\t<parameter id='v' lower='0.0' upper='1.0' name='stateNode'>1.0</parameter>\n";
	os<<"\t\t\t<parameter id='coalescenceRate' lower='0.0' upper='1.0' name='stateNode'>400</parameter>\n";
	os<<"\t\t\t</state>\n\n";
	
    	os<<"\t<init id='SNI' spec='snappNetForSimSnappNet.core.MySpeciesNetworkInitializerWithoutEmbedding' estimate='false' method='random' speciesNetwork='@network:species' origin='@originTime:species'>\n";
        os<<"\t</init>\n\n";

	os<<"\t<distribution id='posterior' spec='util.CompoundDistribution'>\n";
	os<<"\t\t<distribution id='prior' spec='util.CompoundDistribution'>\n";
	os<<"\t\t\t<distribution id='networkPrior' spec='snappNetForSimSnappNet.core.BirthHybridizationModel' network='@network:species' netDiversification='@netDivRate:species' turnOver='@turnOverRate:species' betaShape='1.0'/>\n";
	os<<"\t\t\t <prior id='networkOrigin' name='distribution' x='@originTime:species'>\n";
	os<<"\t\t\t\t <Exponential id='exponential.0' name='distr' mean='0.1'/>\n";
        os<<"\t\t\t </prior>\n";
	os<<"\t\t\t<prior id='netDivPrior' name='distribution' x='@netDivRate:species'>\n";
        os<<"\t\t\t\t<Exponential id='exponential.01' name='distr' mean='10.0'/>\n";
        os<<"\t\t\t</prior>\n";
	os<<"\t\t\t<prior id='turnOverPrior' name='distribution' x='@turnOverRate:species'>\n";
        os<<"\t\t\t\t<Beta id='betadistr.01' name='distr' alpha='1.0' beta='1.0'/>\n";
    os<<"\t\t\t</prior>\n";
 	os<<"\t\t\t<prior id='uPrior' name='distribution' x='@u'>\n";
        os<<"\t\t\t\t<OneOnX id='OneOnX.2' name='distr'/>\n";
    os<<"\t\t\t</prior>\n";
    os<<"\t\t\t<prior id='vPrior' name='distribution' x='@v'>\n";
        os<<"\t\t\t\t<OneOnX id='OneOnX.3' name='distr'/>\n";
    os<<"\t\t\t</prior>\n";
	os<<"\t\t\t<distribution spec='snappNetForSimSnappNet.core.SnappNetPrior' name='distribution' id='snapprior' coalescenceRate='@coalescenceRate'>\n";	
		os<<"\t\t\t\t<parameter id='alpha' estimate='false' lower='0.0' name='alpha'>1.0</parameter>\n";
        os<<"\t\t\t\t<parameter id='beta' estimate='false' lower='0.0' name='beta'>200.0</parameter>\n";
	os<<"\t\t\t</distribution>\n";
	os<<"\t\t</distribution>\n";
	os<<"\t\t<distribution id='likelihood' spec='snappNetForSimSnappNet.core.SnappNetNetworkLikelihood' data='@snapalignment' speciesNetwork='@network:species' mutationRateU='@u'  mutationRateV='@v' coalescenceRate='@coalescenceRate'>\n";
	os<<"\t\t</distribution>\n";
	os<<"\t</distribution>\n\n";

	os<<"\t <operator id='divrRateScale:species' spec='ScaleOperator' parameter='@netDivRate:species' scaleFactor='0.5' weight='10.0'/>\n";
	os<<"\t <operator id='turnOverScale:species' spec='ScaleOperator' parameter='@turnOverRate:species' scaleFactor='0.5' weight='10.0'/>\n";
 	os<<"\t <operator id='gammaProbUniform:species' spec='snappNetForSimSnappNet.operators.InheritanceProbUniform' speciesNetwork='@network:species' weight='10.0'/>\n";
	os<<"\t <operator id='gammaProbRndWalk:species' spec='snappNetForSimSnappNet.operators.InheritanceProbRndWalk' speciesNetwork='@network:species' weight='10.0'/>\n";
	os<<"\t <operator id='originMultiplier:species' spec='snappNetForSimSnappNet.operators.OriginMultiplier' speciesNetwork='@network:species' origin='@originTime:species' weight='0.0'/>\n";	
	os<<"\t <operator id='addReticulation:species' spec='snappNetForSimSnappNet.operators.AddReticulation' speciesNetwork='@network:species' weight='10.0' coalescenceRate='@coalescenceRate'/>\n";
	os<<"\t<operator id='deleteReticulation:species' spec='snappNetForSimSnappNet.operators.DeleteReticulation' speciesNetwork='@network:species' weight='10.0'  coalescenceRate='@coalescenceRate'/>\n";
	os<<"\t<operator id='networkMultiplier:species' spec='snappNetForSimSnappNet.operators.NetworkMultiplier' speciesNetwork='@network:species' origin='@originTime:species' weight='5.0'/>\n";
	os<<"\t <operator id='flipReticulation:species' spec='snappNetForSimSnappNet.operators.FlipReticulation' speciesNetwork='@network:species' weight='10.0'/>\n";
	os<<"\t <operator id='relocateBranch:species' spec='snappNetForSimSnappNet.operators.RelocateBranchNarrow' speciesNetwork='@network:species' weight='30.0'/>\n";
	os<<"\t <operator id='nodeSlider:species' spec='snappNetForSimSnappNet.operators.NodeSlider' speciesNetwork='@network:species' origin='@originTime:species' isNormal='true' sigma='0.005' weight='10.0'/>\n";
 	os<<"\t <operator id='NodeUniform:species' spec='snappNetForSimSnappNet.operators.NodeUniform' speciesNetwork='@network:species' weight='10.0'/>\n";
	os<<"\t <operator id='relocateBranchNarrow:species' spec='snappNetForSimSnappNet.operators.RelocateBranchNarrow' speciesNetwork='@network:species' weight='10.0'/>\n";
	os<<"\t <operator id='ChangeUAndV' spec='snappNetForSimSnappNet.operators.ChangeUAndV' u='@u' v='@v' window='0.1'  weight='10'/>\n";
	os<<"\t <operator id='ChangeGamma' spec='snappNetForSimSnappNet.operators.ChangeGamma'  scale='0.5' weight='10' coalescenceRate='@coalescenceRate'/>\n"; 
	os<<"\t <operator id='ChangeAllGamma' spec='snappNetForSimSnappNet.operators.ChangeAllGamma'  scale='0.5' weight='10' coalescenceRate='@coalescenceRate'/>\n\n";

	os<<"\t<logger id='screenlog' logEvery='1000' model='@posterior'>\n";
        os<<"\t\t <log idref='posterior'/>\n";
        os<<"\t\t <log id='ESS.0' spec='util.ESS' arg='@posterior'/>\n"; 
        os<<"\t\t <log idref='likelihood'/>\n";
        os<<"\t\t <log idref='prior'/>\n";         
        os<<"\t</logger>\n";

	os<<"\t<logger id='tracelog' fileName='test_3s.xml.trace.log' logEvery='1000' model='@posterior' sort='smart'>\n";
        os<<"\t\t <log idref='posterior'/>\n";
        os<<"\t\t <log idref='likelihood'/>\n";
        os<<"\t\t <log idref='prior'/>\n"; 
        os<<"\t\t <log idref='netDivRate:species'/>\n"; 
        os<<"\t\t <log idref='turnOverRate:species'/>\n"; 
        os<<"\t\t <log idref='originTime:species'/>\n";
        os<<"\t\t <log id='height:species' speciesNetwork='@network:species' spec='snappNetForSimSnappNet.core.NetworkStatLogger'/>\n";
	os<<"\t\t <log idref='u'/>\n";
	os<<"\t\t <log idref='v'/>\n";
        os<<"\t </logger>\n";

	os<<"\t<logger id='specieslog' fileName='test_3s.xml.species.trees' logEvery='1000' mode='tree'>\n";
        os<<"\t\t<log id='networkLogger:species' spec='snappNetForSimSnappNet.core.NetworkWithMetaDataLogger' speciesNetwork='@network:species' coalescenceRate='@coalescenceRate'/>\n";
	os<<"\t </logger>\n";

os<<"</run>\n";
os<<"\n";
os<<"\n";
os<<"</beast>\n";

}

//End CE generates xml for SnappNet


int main(int argc, char* argv[]) {
	
         //seed_random();
        //srand(2);	

	ArgumentParser ap(argc,argv);	
	//setting seed
	srand(ap.seedValue);

	//Read in the input file.
	//First line is number of taxa
	int nspecies;
	// nr of trees to generate sequences for
	int ntrees;
	ifstream is(ap.inputfile.c_str());
	if (!is) {
		cerr<<"Error reading in file "<<ap.inputfile<<"\n\n";
		printUsage(cerr);
	}
	is>>nspecies;
	if (nspecies<=0)
		printUsage(cerr);
	
	//Get the prefix from the filename.
	size_t dotpos = ap.inputfile.find_last_of('.');
	string fileroot;
	if (dotpos != string::npos) 
		fileroot = ap.inputfile.substr(0,dotpos);
	else
		fileroot = ap.inputfile;
	
	
	//Now read in the taxon names and sample sizes
	vector<uint> sampleSizes(nspecies);
	vector<string> species(nspecies);
	for(int i=0;i<nspecies;i++) {
		is>>species[i];
		is>>sampleSizes[i];
	}
	double u,v;
	is>>u;
	is>>v;
	
	is>>ntrees;
	if (ntrees<=0)
		printUsage(cerr);
	
	cerr<<"Simulating data from "<<ntrees<<" species trees\n";
	cerr<<"Mutation rates u = "<<u<<" v = "<<v<<"\n";
	cerr<<"Number of sites = "<<ap.nsites<<"\n";
	cerr<<"Number of trees = "<<ntrees<<"\n";
	cerr<<"Output to "<<((ap.outputXML)?"XML":"nexus")<<"\n";
	cerr<<((ap.excludeConst)?"Exclude":"Include")<<" constant characters\n\n";
	cerr<<((ap.outputTrees)?"Output":"Don't output")<<" trees\n";
	if (ap.onlyRootMutation)
		cerr<<"Mutation only at root\n";
	if (ap.hasDominantMarkers)
		cerr<<"Simulate dominant markers\n";
	//CE: until now, we leave  ap.hasDominantMarkers to default value false since we did not handle it
	//fixit later to handle dominant markers
	
	cerr<<"Species and sample sizes:\n";
	for(int i=0;i<nspecies;i++)
		cerr<<(i+1)<<" "<<species[i]<<"\t"<<sampleSizes[i]<<"\n";
	cerr<<endl;
	
	cerr<<"Fileroot = "<<fileroot<<endl;
	
	
	for (int iTree = 0; iTree < ntrees; iTree++) {
		//Now read in tree string
		string treeString;
		is>>treeString;
		//Check for scaling
		double scaling = 1.0;
		size_t pos = treeString.find_first_of("<");
		
		if (pos!=string::npos) {
			size_t pos2 = treeString.find_first_of(">");
			string scaleString = treeString.substr(pos+1,pos2-1);
			treeString = treeString.substr(pos2+1);
			//cerr<<"Scale String is "<<scaleString<<endl;
			scaling = strtod(scaleString.c_str(),NULL);
		}
		
		//cerr<<"Tree String is "<<treeString<<" scaling = "<<scaling<<endl;
		//TODO: read in semicolons, or at least check for them.
		
		phylo<basic_newick> tree;
		read_newick(treeString, tree, species, 0.0);
		
		//Scale tree here.
		
		for(phylo<basic_newick>::iterator p = tree.root();!p.null();p=p.next_pre()) {
			(*p).length *= scaling;
		}
		
		print_newick(cerr,tree,true,true);
		cerr<<endl;
		
		// MOIII

		// Write the newick with species names
		printf("\n");
		print_newick(cerr,tree,species,true,true);

		vector<string> taxa_names=species;
		//	printf("there is a total of  %lu species names\n",taxa_names.size());

		vector<int> retic_id;
		vector<double> retic_proba;
		int nbRetic=0;
				
		for(int i=0; i<taxa_names.size(); i++)
		  {
		    printf("The name of taxon%d is  %s\n",i,(taxa_names[i]).c_str());		    
		    char b = taxa_names[i].at(0);
		    // printf("First letter of taxon%d is  %c\n",i,b);

		    if (b=='#'){
		      printf(" this is a reticulation node\n");
		      // its id is i by definition
		      retic_id.push_back(i);
		      printf("the interesting part of the reticulation node is %s\n",taxa_names[i].substr(2,taxa_names[i].length()-2).c_str());

		      //We suppose that the reticulation name has only one letter e.g. #B, not #BB 
		      printf("the probability associated with this retic node is %f\n",atof(taxa_names[i].substr(2,taxa_names[i].length()-2).c_str()));
		      retic_proba.push_back(atof(taxa_names[i].substr(2,taxa_names[i].length()-2).c_str()));
		      nbRetic++;
		    }
		    
		  }


		printf("\n\n");
		printf("There are %d reticulation nodes\n",nbRetic);
		//printf("size of retic_id %ld is!!!\n",retic_id.size());
		//printf("size of retic_proba %ld is!!!\n",retic_proba.size());
		
		for(int i=0; i<nbRetic; i++){
		  printf("The id of the  %d th reticulation node is %d\n",i,retic_id[i]);
		  printf("the probability associated with this retic node is %f\n",(float) retic_proba[i]);
		}
		
		
		ostringstream s1;
		//s1 << fileroot<<"_tree_"<<(iTree+1);
		
		//CE changes to 
		s1 << fileroot<<"_seed_"<<ap.seedValue<<"_nsites_"<<ap.nsites;
		//END CE

		string shortFileName = s1.str();


		//fix it later, see if we keep the if
		if (ap.outputXML)
			s1<< ".xml";
		else
			s1 << ".nex";		
		
		string sFile = s1.str();
		
		cerr << "Writing " << sFile << endl;		
		ofstream * os =  new ofstream(sFile.c_str());
		vector<vector<uint> > alleleCounts;		
		
		printf("Let us simulate data !!!\n");

		simulateMultipleSites(tree, u, v, sampleSizes, ap.nsites, ap.excludeConst, ap.onlyRootMutation, ap.hasDominantMarkers, alleleCounts, ap.outputTrees,retic_id,retic_proba);
					
		if (ap.outputXML) {
        		   
       		output_xml_SnappNet(*os,species,tree,sampleSizes,u,v,alleleCounts,ap.simulationOutput,ap.excludeConst,ap.initialiseAtTrue,ap.snappDominant,ap.onlyRootMutation,ap.useSNAPPMCMC,shortFileName);

		} 

		(*os).close();
		
		//Remove scaling from tree
		for(phylo<basic_newick>::iterator p = tree.root();!p.null();p=p.next_pre()) {
			(*p).length /= scaling;
		}
	}
	
	
}

