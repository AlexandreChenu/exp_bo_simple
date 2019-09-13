#ifndef __TOOLS_HPP__
#define __TOOLS_HPP__

#include <csignal>

#include <iostream>
#include <fstream> 
#include <sstream>
#include <string>

#include <vector>
#include <map>
#include <assert.h>
#include <stdio.h>      /* printf */

#include "params.hpp"




std::map< std::vector<double>, Params_ad::archiveparams::elem_archive, Params_ad::archiveparams::classcomp> load_archive(std::string archive_name, int n_it)
{

    int dim_in = Params_ad::eval_params::bd_dim; //accessible depuis le fichier eval? 

    std::cout << "Loading text file..." << std::endl;

    std::ifstream dic_file;
    dic_file.open(archive_name);

    assert(dic_file.is_open());

    std::map< std::vector<double>, Params_ad::archiveparams::elem_archive, Params_ad::archiveparams::classcomp> archive;

    std::string line;
    int i = 0;

    std::cout << "Loading dictionary" << std::endl;

    while(getline(dic_file, line)){

    	Params_ad::archiveparams::elem_archive elem;

    	std::stringstream ss(line);

    	std::vector<double> zones(dim_in);

        ss >> elem.filename;
        ss >> elem.fit;
        ss >> zones[0];
        ss >> zones[1];
        ss >> zones[2];

        elem.zones = zones;
	//elem.set_chosen(false);
	//
	//Params_ad::archiveparams::classcomp classcomp;

        //archive.insert({zones, elem}); //on ajoute chaque élément de chaque ligne dans le dictionaire
	archive[zones] = elem;
        i++;
    }	

    if (archive.size() == 0){
      std::cerr << "ERROR: Could not load the archive." << std::endl;
      return archive;
    } 
   
   std::cout << "i: " << i << std::endl; 

    std::cout << archive.size() << " elements loaded" << std::endl;
    return archive;
}

#endif
