#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <map>
#include <assert.h>
#include <stdio.h>      /* printf */
#include <math.h>  

#include <cmath>
#include <algorithm>
#include <cstdlib>
#include <Eigen/Core>

#include <limbo/limbo.hpp>

#include "exhaustive_search_archive.hpp"
#include "mean_archive.hpp"
#include "eval.hpp"

#include <boost/test/unit_test.hpp>
#include <boost/serialization/shared_ptr.hpp>
#include <boost/serialization/nvp.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/fusion/include/vector.hpp> 

#include <modules/nn2/mlp.hpp>
#include <modules/nn2/gen_dnn.hpp>
#include <modules/nn2/phen_dnn.hpp>
#include <modules/nn2/gen_dnn_ff.hpp>

#include "fit_behav.hpp"

#include "params.hpp"

#include "tools.hpp"



Params_ad::archiveparams::archive_t Params_ad::archiveparams::archive;
Eigen::Vector2d Params_ad::eval_params::pos;
std::string Params_ad::eval_params::filename;
int Params_ad::eval_params::bd_dim;
BO_DECLARE_DYN_PARAM(int, Params_ad::stop_maxiterations, iterations);

int main(int argc, char** argv)
{   
    using namespace sferes;
    using namespace nn;
    using namespace limbo;

    typedef nn_mlp<Params> fit_t; 
    typedef phen::Parameters<gen::EvoFloat<1, Params>, fit::FitDummy<>, Params> weight_t;
    //typedef phen::Parameters<gen::EvoFloat<1, Params>, fit::FitDummy<>, Params> bias_t;
    typedef PfWSum<weight_t> pf_t;
    typedef AfSigmoidNoBias<> af_t;
    typedef sferes::gen::DnnFF<Neuron<pf_t, af_t>,  Connection<weight_t>, Params> gen_t; // TODO : change by DnnFF in order to use only feed-forward neural networks
                                                                                       // TODO : change by hyper NN in order to test hyper NEAT 
    typedef phen::Dnn<gen_t, fit_t, Params> phen_t;


        // //limbo types definition
    typedef kernel::MaternFiveHalves<Params_ad> Kernel_t;
    typedef opt::ExhaustiveSearchArchive<Params_ad> InnerOpt_t;
    typedef boost::fusion::vector<stop::MaxIterations<Params_ad>> Stop_t;
    //typedef boost::fusion::vector<stop::MaxIterations<Params_ad>, stop::MaxPredictedValue<Params_ad> > Stop_t;
    typedef mean::MeanArchive<Params_ad> Mean_t;
    typedef boost::fusion::vector<limbo::stat::Samples<Params_ad>, limbo::stat::BestObservations<Params_ad>> Stat_t;

    typedef init::NoInit<Params_ad> Init_t;
    typedef model::GP<Params_ad, Kernel_t, Mean_t> GP_t;
    typedef acqui::UCB<Params_ad, GP_t> Acqui_t;

//   typedef init::RandomSampling<Params_ad> Init_t;

    bayes_opt::BOptimizer<Params_ad, modelfun<GP_t>, initfun<Init_t>, acquifun<Acqui_t>, acquiopt<InnerOpt_t>, statsfun<Stat_t>, stopcrit<Stop_t>> opt;


    //init params
    int n_it = 20;
    Params_ad::eval_params::pos[0] = 0.2;
    Params_ad::eval_params::pos[1] = 0.25;
    Params_ad::eval_params::filename = "/git/sferes2/exp/ex_data/2019-08-20_13_59_03_1913/";
    Params_ad::eval_params::bd_dim = 3;
    Params_ad::archiveparams::archive = load_archive(Params_ad::eval_params::filename + "dict_models.txt", n_it); //TODO hard_code path
    Params_ad::stop_maxiterations::set_iterations(n_it);  

    std::cout << "start optimization" << std::endl;

    opt.optimize(Eval<phen_t>());

    auto val = opt.best_observation();
    Eigen::VectorXd result = opt.best_sample().transpose();
    
    std::cout << val << " res  " << result.transpose() << std::endl;
    
    std::vector<double> final_bd(Params_ad::eval_params::bd_dim);

    for (int i = 0; i < Params_ad::eval_params::bd_dim; i++)
	    final_bd[i] = result[i];

    std::cout << "best model after bo: " << Params_ad::archiveparams::archive.at(final_bd).filename << std::endl;

    return 0;
}
