/*
 * main.cpp
 *
 *  Created on: Dec 12, 2019
 *      Author: Quentin Marcou, Carlos Olivares
 *
 *  This source code is distributed as part of the IGoR software.
 *  IGoR (Inference and Generation of Repertoires) is a versatile software to analyze and model immune receptors
 *  generation, selection, mutation and all other processes.
 *   Copyright (C) 2017  Quentin Marcou
 *
 *   This program is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   This program is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.

 *   You should have received a copy of the GNU General Public License
 *   along with this program.  If not, see <https://www.gnu.org/licenses/>.
 *
 */

#include "../config.h"
#include "Deletion.h"
#include "Insertion.h"
#include "Genechoice.h"
#include "Model_Parms.h"
#include "Rec_Event.h"
#include "Singleerrorrate.h"
#include "Model_marginals.h"
#include <iostream>
#include "Aligner.h"
#include "GenModel.h"
#include "Dinuclmarkov.h"
#include "Counter.h"
#include "Coverageerrcounter.h"
#include "Bestscenarioscounter.h"
#include "Pgencounter.h"
#include "Errorscounter.h"
#include "Utils.h"
#include <chrono>
#include <set>


using namespace std;
// TODO: Given command line options we can construct a model for IGoR or OLGA maybe?
//

/*
class IOModel{
public:
  IOModel();
  virtual ~IOModel();
  vector Genechoice;


}
*/

int main(int argc , char* argv[]){

	//Command line argument iterator
	size_t carg_i = 1;
	//cout<<argv[argc-1]<<endl;

	// Test if some commands were supplied to IGoR
	//if(argc<2){
	//	return terminate_IGoR_with_error_message("The user did not supply IGoR any command.");
	//}

	//Task variables
	bool run_demo = false;
	bool align = false;
	bool infer = false;
	bool evaluate = false;
	bool generate = false;
	bool custom = false;

	//Common vars
	string batchname="";


	//Working directory vars
	bool wd = false;
	string cl_path;

	//Chains vars
	bool chain_provided = false;
	string chain_arg_str;
	string chain_path_str;
	bool has_D = false;

	//Species vars
	bool species_provided = false;
	string species_str = "";

	//Custom genomic templates loading variables
	bool custom_v = false;
	string custom_v_path;
	bool custom_d = false;
	string custom_d_path;
	bool custom_j = false;
	string custom_j_path;

	//Custom gene anchors loading variables
	bool custom_v_anchors = false;
	string custom_v_anchors_path;
	bool custom_j_anchors = false;
	string custom_j_anchors_path;

	//Genomic templates list and aligns parms
	// vector<pair<string,string>> v_genomic;
	// vector<pair<string,string>> d_genomic;
	// vector<pair<string,string>> j_genomic;
	unordered_map<string,size_t> v_CDR3_anchors;
	unordered_map<string,size_t> j_CDR3_anchors;

	//Model parms and marginals
	bool load_last_inferred_parms = false;
	bool custom_cl_parms = false;
	Model_Parms cl_model_parms;
	Model_marginals cl_model_marginals;

  
  /*
	if(chain_provided){
		if(chain_arg_str == "alpha"){
			has_D = false;
			chain_path_str = "tcr_alpha";
			try{
				v_genomic = read_genomic_fasta(string(IGOR_DATA_DIR) + "/models/"+species_str+"/"+chain_path_str+"/ref_genome/genomicVs.fasta");
			}
			catch(exception& e){
				return terminate_IGoR_with_error_message("Exception caught while reading TRA V genomic templates.",  e);
			}

			try{
				j_genomic = read_genomic_fasta(string(IGOR_DATA_DIR) + "/models/"+species_str+"/"+chain_path_str+"/ref_genome/genomicJs.fasta");
			}
			catch(exception& e){
				return terminate_IGoR_with_error_message("Exception caught while reading TRA J genomic templates.",  e);
			}
		}
		else if(chain_arg_str == "beta"){
			has_D = true;
			chain_path_str = "tcr_beta";
			try{
				v_genomic = read_genomic_fasta(string(IGOR_DATA_DIR) + "/models/"+species_str+"/"+chain_path_str+"/ref_genome/genomicVs.fasta");
			}
			catch(exception& e){
				return terminate_IGoR_with_error_message("Exception caught while reading TRB V genomic templates.",  e);
			}

			try{
				d_genomic = read_genomic_fasta(string(IGOR_DATA_DIR) + "/models/"+species_str+"/"+chain_path_str+"/ref_genome/genomicDs.fasta");
			}
			catch(exception& e){
				return terminate_IGoR_with_error_message("Exception caught while reading TRB D genomic templates.",  e);
			}

			try{
				j_genomic = read_genomic_fasta(string(IGOR_DATA_DIR) + "/models/"+species_str+"/"+chain_path_str+"/ref_genome/genomicJs.fasta");
			}
			catch(exception& e){
				return terminate_IGoR_with_error_message("Exception caught while reading TRB J genomic templates.",  e);
			}

		}
		else if(chain_arg_str == "light"){
			forward_list<string> error_messages;
			error_messages.emplace_front("If you wish to use IGoR on light chains please contact us so we can work on incorporating a light chain model to IGoR.");
			error_messages.emplace_front("Support for light chains in command line is not ready yet due to the lack of genomic templates and suitable model.");
			error_messages.emplace_front("Light chains support does not exist yet for command line!");
			return terminate_IGoR_with_error_message(error_messages);
		}
		else if( (chain_arg_str == "heavy_naive") or (chain_arg_str == "heavy_memory") ){
			has_D = true;
			chain_path_str = "bcr_heavy";
			try{
				v_genomic = read_genomic_fasta(string(IGOR_DATA_DIR) + "/models/"+species_str+"/"+chain_path_str+"/ref_genome/genomicVs.fasta");
			}
			catch(exception& e){
				return terminate_IGoR_with_error_message("Exception caught while reading IGH V genomic templates.",  e);
			}

			try{
				d_genomic = read_genomic_fasta(string(IGOR_DATA_DIR) + "/models/"+species_str+"/"+chain_path_str+"/ref_genome/genomicDs.fasta");
			}
			catch(exception& e){
				return terminate_IGoR_with_error_message("Exception caught while reading IGH D genomic templates.",  e);
			}

			try{
				j_genomic = read_genomic_fasta(string(IGOR_DATA_DIR) + "/models/"+species_str+"/"+chain_path_str+"/ref_genome/genomicJs.fasta");
			}
			catch(exception& e){
				return terminate_IGoR_with_error_message("Exception caught while reading IGH J genomic templates.",  e);
			}
			if( chain_arg_str == "heavy_naive" ){
				//Use a single error rate
			}
			else{
				//Memory

				//TODO infer only \mu for the hypermutation model
			}
		}
		//Read CDR3 anchors(cystein, tryptophan/phenylalanin indices)
		try{
			v_CDR3_anchors = read_gene_anchors_csv(string(IGOR_DATA_DIR) + "/models/"+species_str+"/"+chain_path_str+"/ref_genome/V_gene_CDR3_anchors.csv");
		}
		catch(exception& e){
			return terminate_IGoR_with_error_message("Exception caught while reading V CDR3 anchors.",  e);
		}

		try{
			j_CDR3_anchors = read_gene_anchors_csv(string(IGOR_DATA_DIR) + "/models/"+species_str+"/"+chain_path_str+"/ref_genome/J_gene_CDR3_anchors.csv");
		}
		catch(exception& e){
			return terminate_IGoR_with_error_message("Exception caught while reading J CDR3 anchors.",  e);
		}
	}

	//Read custom genomic templates if some custom ones were specified
	if(custom_v){
		try{
			v_genomic = read_genomic_fasta(custom_v_path);
		}
		catch(exception& e){
			return terminate_IGoR_with_error_message("Exception caught while reading user's custom V genomic templates.",  e);
		}
	}
	if(custom_d){
		has_D = true;
		try{
			d_genomic = read_genomic_fasta(custom_d_path);
		}
		catch(exception& e){
			return terminate_IGoR_with_error_message("Exception caught while reading user's custom D genomic templates.",  e);
		}
	}
	if(custom_j){
		try{
			j_genomic = read_genomic_fasta(custom_j_path);
		}
		catch(exception& e){
			return terminate_IGoR_with_error_message("Exception caught while reading user's custom J genomic templates.",  e);
		}
	}

	//Read custom CDR3
	if(custom_v_anchors){
		try{
			v_CDR3_anchors = read_gene_anchors_csv(custom_v_anchors_path);
		}
		catch(exception& e){
			return terminate_IGoR_with_error_message("Exception caught while reading user's custom V CDR3 anchors",  e);
		}
	}
	if(custom_j_anchors){
		try{
			j_CDR3_anchors = read_gene_anchors_csv(custom_j_anchors_path);
		}
		catch(exception& e){
			return terminate_IGoR_with_error_message("Exception caught while reading user's custom J CDR3 anchors",  e);
		}
	}

	//Make sure passed arguments are unambiguous
	if(custom_cl_parms and load_last_inferred_parms){
		return terminate_IGoR_with_error_message("Setting a custom model and loading the last inferred model in the same command is ambiguous!");
	}

	//Load last inferred model
	if(load_last_inferred_parms){
		clog<<"Loading last inferred model..."<<endl;
		try{
			cl_model_parms.read_model_parms(cl_path +  batchname + "inference/final_parms.txt");
			cl_model_marginals = Model_marginals(cl_model_parms);
			cl_model_marginals.txt2marginals(cl_path +  batchname + "inference/final_marginals.txt",cl_model_parms);

			//Check if the model contains a D gene event in order to load the alignments
			auto events_map = cl_model_parms.get_events_map();
			if(events_map.count(tuple<Event_type,Gene_class,Seq_side>(GeneChoice_t,D_gene,Undefined_side))>0){
				has_D = true;
			}
		}
		catch(exception& e){
			return terminate_IGoR_with_error_message("Exception caught while loading last inferred model, please check that the model exists",e);
		}
	}
  */



	/*
	 * If some custom genomic templates were supplied, two possible cases here:
	 * - if any supplied genomic template is absent from the model, or if its actual sequence is different the marginals will be re-initialized
	 * - if all the supplied templates were already contained in the model the missing one will be set to 0 probability,
	 * 	 the others will keep their probability ratio.
	 *
	 * This will be executed whether using a supplied model, a custom model or the last inferred one.
	 */
	if((infer or evaluate or generate)){
		bool any_custom_gene = false;
		unordered_map<tuple<Event_type,Gene_class,Seq_side>,shared_ptr<Rec_Event>> tmp_events_map = cl_model_parms.get_events_map();
	}




  /*Run this sample demo code
   *
   * Outline:
   *
   * Read TCRb genomic templates
   *
   * Align the sequences contained in the /demo/murugan_naive1_noncoding_demo_seqs.txt file to those templates
   *
   * Create a TCRb model, a simple error rate and the corresponding marginals
   *
   * Show reads and write functions for the model and marginals
   *
   * Infer a model from the sequences (perform 10 iterations of EM)
   *
   * Generate sequences from the obtained model
   *
   */

  clog<<"Reading genomic templates"<<endl;

  std::string fln_V_GenomicFasta = "";
  std::string fln_D_GenomicFasta = "";
  std::string fln_J_GenomicFasta = "";
  std::cout << "IGOR_DATA_DIR : " << string(IGOR_DATA_DIR) << std::endl;

	//v_genomic = read_genomic_fasta(string(IGOR_DATA_DIR) + "/models/"+species_str+"/"+chain_path_str+"/ref_genome/genomicVs.fasta");
  species_str = "human";
	chain_path_str = "tcr_beta";
	string str_path_ref_genome = string(IGOR_DATA_DIR) + "/models/"+species_str+"/"+chain_path_str+"/ref_genome/"; // genomicVs.fasta";
	fln_V_GenomicFasta = str_path_ref_genome + "genomicVs.fasta";
	fln_D_GenomicFasta = str_path_ref_genome + "genomicDs.fasta";
	fln_J_GenomicFasta = str_path_ref_genome + "genomicJs.fasta";
	
  vector<pair<string,string>> v_genomic = read_genomic_fasta( fln_V_GenomicFasta );
  vector<pair<string,string>> d_genomic = read_genomic_fasta( fln_D_GenomicFasta );
  vector<pair<string,string>> j_genomic = read_genomic_fasta( fln_J_GenomicFasta );


  clog<<"Construct the model"<<endl;
  //Construct a TCRb model
  Gene_choice v_choice(V_gene,v_genomic);
  v_choice.set_nickname("v_choice");
  v_choice.set_priority(7);
  Gene_choice d_choice(D_gene,d_genomic);
  d_choice.set_nickname("d_gene");
  d_choice.set_priority(6);
  Gene_choice j_choice(J_gene,j_genomic);
  j_choice.set_nickname("j_choice");
  j_choice.set_priority(7);

  
	Deletion v_3_del(V_gene,Three_prime,make_pair(-4,16));//16
  v_3_del.set_nickname("v_3_del");
  v_3_del.set_priority(5);
  Deletion d_5_del(D_gene,Five_prime,make_pair(-4,16));
  d_5_del.set_nickname("d_5_del");
  d_5_del.set_priority(5);
  Deletion d_3_del(D_gene,Three_prime,make_pair(-4,16));
  d_3_del.set_nickname("d_3_del");
  d_3_del.set_priority(5);
  Deletion j_5_del(J_gene,Five_prime,make_pair(-4,18));
  j_5_del.set_nickname("j_5_del");
  j_5_del.set_priority(5);

  Insertion vd_ins(VD_genes,make_pair(0,30));
  vd_ins.set_nickname("vd_ins");
  vd_ins.set_priority(4);
  Insertion dj_ins(DJ_genes,make_pair(0,30));
  dj_ins.set_nickname("dj_ins");
  dj_ins.set_priority(2);

  Dinucl_markov markov_model_vd(VD_genes);
  markov_model_vd.set_nickname("vd_dinucl");
  markov_model_vd.set_priority(3);

  Dinucl_markov markov_model_dj(DJ_genes);
  markov_model_dj.set_nickname("dj_dinucl");
  markov_model_dj.set_priority(1);


  Model_Parms parms;

  //Add nodes to the graph
  parms.add_event(&v_choice);
  parms.add_event(&d_choice);
  parms.add_event(&j_choice);

  parms.add_event(&v_3_del);
  parms.add_event(&d_3_del);
  parms.add_event(&d_5_del);
  parms.add_event(&j_5_del);

  parms.add_event(&vd_ins);
  parms.add_event(&dj_ins);

  parms.add_event(&markov_model_vd);
  parms.add_event(&markov_model_dj);


  //Add correlations
  parms.add_edge(&v_choice, &v_3_del);
  parms.add_edge(&j_choice, &j_5_del);
  parms.add_edge(&d_choice, &d_3_del);
  parms.add_edge(&d_choice, &d_5_del);
  parms.add_edge(&d_5_del, &d_3_del);
  parms.add_edge(&j_choice, &d_choice);


  //Create the corresponding marginals
  Model_marginals model_marginals(parms);
  model_marginals.uniform_initialize(parms); //Can also start with a random prior using random_initialize()

  //Instantiate an error rate
  Single_error_rate error_rate(0.001);

  parms.set_error_ratep(&error_rate);

  clog<<"Write and read back the model"<<endl;
  //Write the model_parms into a file
  //parms.write_model_parms(string(cl_path + "/demo_write_model_parms.txt"));
  string flnModel_parms = "tmp_model_parms.txt";
  parms.write_model_parms(flnModel_parms);

  //Write the marginals into a file
  string flnModel_marginals = "tmp_model_marginals.txt";
  model_marginals.write2txt(flnModel_marginals, parms);

  //Read a model and marginal pair
  Model_Parms read_model_parms;
  read_model_parms.read_model_parms(flnModel_parms);
  Model_marginals read_model_marginals(read_model_parms);
  read_model_marginals.txt2marginals(flnModel_marginals, read_model_parms);


	return EXIT_SUCCESS;

}



