#ifndef EXHAUSTIVE_SEARCH_ARCHIVE_HPP_
#define EXHAUSTIVE_SEARCH_ARCHIVE_HPP_
namespace limbo {
    namespace opt {
        template <typename Params>
        struct ExhaustiveSearchArchive {

            ExhaustiveSearchArchive() {}
            template <typename F>
            Eigen::VectorXd operator()(const F& f, const Eigen::VectorXd& init, bool bounded) const
            {
                float best_acqui = -INFINITY;
                Eigen::VectorXd result;
		//std::vector<double> best_zones;
                typedef typename Params::archiveparams::archive_t::const_iterator archive_it_t;

                for (archive_it_t it = Params::archiveparams::archive.begin(); it != Params::archiveparams::archive.end(); ++it) {

                    Eigen::VectorXd temp(it->first.size());
		//    std::vector<double> temp_zones(it->second.zones.size());

                    for (size_t i = 0; i < it->first.size(); i++)
                        temp[i] = it->first[i]; //ici on fait que convertir un vecteur en eigen encore -> on récupère le behavior descriptor dans notre cas

                    float new_acqui = eval(f, temp); //test du temp
                    //float new_acqui = it->second.fit; //on se base sur la fitness déjà calculée pendant l'apprentissage? 
                    
		  //  for(int i=0; it->first.size(); i++)
		//	    temp_zones[i] = it->first[i];

                    if (best_acqui < new_acqui || it == Params::archiveparams::archive.begin()) {
                        best_acqui = new_acqui;
                        result = temp;
		//	for(int i = 0; i < temp_zones.size(); i++)
		//		best_zones[i] = temp_zones[i];
                    }
                }

		//Params::archiveparams::archive(best_zones);
		
                //return result;
		//Eigen::VectorXd result(Params::archiveparams::n_best[0].size());

		//for (int i= 0; i<Params::archiveparams::n_best[0].size(); i++){
		//	result[i] = Params::archiveparams::n_best[0][i];}
		
		//std::vector<std::vector<double>>::const_iterator it_t;

		//Params::archiveparams::n_best.erase(Params::archiveparams::n_best.begin());

		return result;

        }
    };
}

}
#endif
