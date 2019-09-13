#ifndef __PARAMS__ADAPTATION_HPP__
#define __PARAMS__ADAPTATION_HPP__

#include <Eigen/Core>
#include <limbo/limbo.hpp>

using namespace limbo;

struct Params_ad {


    struct bayes_opt_boptimizer : public defaults::bayes_opt_boptimizer {
    };

    struct bayes_opt_bobase : public defaults::bayes_opt_bobase {
        BO_PARAM(int, stats_enabled, true);
    };

    // no noise
    struct kernel : public defaults::kernel {
        BO_PARAM(double, noise, 0.5);
    };

    struct kernel_maternfivehalves : public defaults::kernel_maternfivehalves {
        BO_PARAM(double, l, 0.4);
    };

    struct stop_maxiterations {
        BO_DYN_PARAM(int, iterations);
    };

    struct init_randomsampling {
        BO_PARAM(int, samples, 10);
    };

    struct stop_maxpredictedvalue: public defaults::stop_maxpredictedvalue{
    };
    
    struct acqui_ucb : public defaults::acqui_ucb {
        BO_PARAM(double, alpha, 0.2);
    };

    struct archiveparams {
        struct elem_archive {
            float fit; //fitness - double
            std::vector<double> zones; //behavior descriptor - 3 zones
            std::string filename; //name of the corresponding model - string
            
	    };

        struct classcomp {
             bool operator()(const std::vector<double>& lhs, const std::vector<double>& rhs) const
             {
                 assert(lhs.size() == 3 && rhs.size() == 3);
                 int i = 0;
                 while (i < 2 && std::round(lhs[i] *10) == std::round(rhs[i] * 10)) //lhs[i]==rhs[i])
                     i++;
                 return std::round(lhs[i] * 10) < std::round(rhs[i] * 10); //lhs[i]<rhs[i];
             }
        };
		
         typedef std::map<std::vector<double>, elem_archive, classcomp> archive_t;
         static std::map<std::vector<double>, elem_archive, classcomp> archive; //archive

        //typedef std::map< std::vector<double>, elem_archive> archive_t;
        //static std::map< std::vector<double>, elem_archive> archive;

    };
    

    struct eval_params {

        static Eigen::Vector2d pos; //position de l'obstacle
        static std::string filename; //chemin jusqu'au fichiers
        static int bd_dim;
};

};

#endif
