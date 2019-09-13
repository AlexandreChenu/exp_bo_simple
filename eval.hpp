#ifndef __EXPERIMENT__EVAL_HPP__
#define __EXPERIMENT__EVAL_HPP__

#ifndef EIGEN3_ENABLED 
#define EIGEN3_ENABLED

//#include "desc_hexa.hpp"

#include <modules/nn2/mlp.hpp>
#include <modules/nn2/gen_dnn.hpp>
#include <modules/nn2/phen_dnn.hpp>
#include <modules/nn2/gen_dnn_ff.hpp>


#include <boost/test/unit_test.hpp>
#include <boost/serialization/shared_ptr.hpp>
#include <boost/serialization/nvp.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>

#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <map>
#include <assert.h>
#include <stdio.h>      /* printf */
#include <math.h>  

#include <typeinfo>
#include <cmath>
#include <algorithm>
#include <cstdlib>
#include <Eigen/Core>

#include "params.hpp"

template<typename phen_t>
struct Eval {
    BO_PARAM(size_t, dim_in, 3); //behavior descriptor 
    BO_PARAM(size_t, dim_out, 1); //fitness

    // the function to be optimized
  Eigen::VectorXd operator()(const Eigen::VectorXd& x) const//x c'est le behavior descriptor 
  {
    
    std::cout << "operator eval" << std::endl; 
    std::vector<double> key(x.size(), 0);  
    Eigen::VectorXd::Map(key.data(), key.size()) = x; //a quoi sert cette ligne??

    std::cout << "behavior descriptor evaluated: " << key[0] << " / " << key[1] << std::endl;
    std::string model_name = Params_ad::archiveparams::archive.at(key).filename; //si la key est out of range c'est surement que le bd n'est pas reconnu

    std::cout << "model name: " << model_name << std::endl;
    Eigen::Vector2d obst = Params_ad::eval_params::pos; //

    float real_fitness = eval(obst, model_name);

    Eigen::VectorXd result(1);
    result[0] = real_fitness;
    
    std::cout<<"fit in archive: "<<Params_ad::archiveparams::archive.at(key).fit<<std::endl;
    std::cout<<"real fitness: "<< real_fitness << std::endl;
    return result;
     
  }

  float eval(Eigen::Vector2d& obst,  std::string model_name) const 
  { 
    std::cout << "function eval" << std::endl;
    //init model 
    typedef boost::archive::binary_iarchive ia_t;
    phen_t model; 

    const std::string modelname = Params_ad::eval_params::filename + model_name + ".bin";

    std::cout << modelname << std::endl;

    std::cout << "model...loading" << std::endl;
    {
    std::ifstream ifs(modelname , std::ios::binary);

    ia_t ia(ifs);
    ia >> model;
    }

    std::cout << "model...loaded" << std::endl;

    model.develop();
    std::cout << "model developed" << std::endl;

    //init variables
    double _vmax = 1;
    double _delta_t = 0.1;
    double _t_max = 10; //TMax guidé poto
    Eigen::Vector3d robot_angles;
    Eigen::Vector3d target;
    double dist = 0;
    std::ofstream logfile;
    Eigen::Vector3d prev_pos; //compute previous position
    Eigen::Vector3d output;
    Eigen::Vector3d pos_init;

    robot_angles = {0,M_PI,M_PI}; //init everytime at the same place
    target = {-0.211234, 0.59688,0.0}; //target fixe dans notre cas

    //open logfile
    std::cout << "open logfile" << std::endl;
    std::string filename = Params_ad::eval_params::filename + model_name + ".txt";
    //open logfile
    logfile.open(filename);

    std::cout << "logfile opened" << std::endl;

    //get gripper's position
    std::vector<double> rfm_pos = real_forward_model(robot_angles, obst);

    for (int i=0; i<3; i++){
      prev_pos[i] = rfm_pos[i];
      pos_init[i] = rfm_pos[i];}

    std::vector<float> inputs(5);

    int t = 0;
    bool terminal = false;
    bool collision = false;

    while(!terminal){
          
          //compute nn inputs
          inputs[0] = target[0] - prev_pos[0]; //get side distance to target
          inputs[1] = target[1] - prev_pos[1]; //get front distance to target
          inputs[2] = robot_angles[0];
          inputs[3] = robot_angles[1];
          inputs[4] = robot_angles[2];

          //save nn inputs
          logfile << inputs[0] << "    " << inputs[1] << "    ";
	 
          for (int j = 0; j < 10   + 1; ++j){
            model.nn().step(inputs);}

          for (int indx = 0; indx < 3; ++indx){
            output[indx] = 2*(model.nn().get_outf(indx) - 0.5)*_vmax; //Remap to a speed between -v_max and v_max (speed is saturated)
            robot_angles[indx] += output[indx]*_delta_t; //Compute new angles
          }
          logfile << output[0] << "    " << output[1] << "    ";

          //Eigen::Vector3d new_pos;
          std::vector<double> real_out(prev_pos.size()+1,0);
          real_out = real_forward_model(robot_angles, obst); //compute simulation

          for (int i = 0; i<prev_pos.size(); i++){
            prev_pos[i] = real_out[i]; //replace previous positions with the new ones
          }

	  //std::cout << "position :\n" << prev_pos << std::endl;

          //write data into logfile
          logfile << robot_angles[0] << "    " << robot_angles[1] << "    " << robot_angles[2] << "    ";
          logfile << prev_pos[0] << "    " << prev_pos[1] << "    ";
          logfile << target[0] << "    " << target[1] << "\n";
          logfile << obst[0] << "    " << obst[1] << "\n";

          //iterate through time
          t++;

  	  if (sqrt((target[0] - prev_pos[0])*(target[0] - prev_pos[0]) + (target[1] - prev_pos[1])*(target[1] - prev_pos[1])) < 0.02){
		dist -= sqrt((target[0] - prev_pos[0])*(target[0] - prev_pos[0]) + (target[1] - prev_pos[1])*(target[1] - prev_pos[1])); 
		}	
          else {
                dist -= log(1+t) + sqrt((target[0] - prev_pos[0])*(target[0] - prev_pos[0]) + (target[1] - prev_pos[1])*(target[1] - prev_pos[1]));
          } 

          //check if simulation ended
          if (t >= _t_max/_delta_t - 1 || real_out[3] == 1){
            terminal = true;

	    if (real_out[3] == 1)
		    collision = true;
          }

        }
    logfile.close();

	if (sqrt((target[0] - prev_pos[0])*(target[0] - prev_pos[0]) + (target[1] - prev_pos[1])*(target[1] - prev_pos[1])) < 0.02){
           return( 1.0 + dist/500); // -> 1
        }

        if (collision){
          return( -1 + dist/500);} // -> -1

	else 
		return(dist/500);
  }

  std::vector<double> real_forward_model(Eigen::VectorXd a, Eigen::Vector2d obst) const
  
  {
    
    Eigen::VectorXd _l_arm=Eigen::VectorXd::Ones(a.size()+1);
    _l_arm(0)=0;
    _l_arm = _l_arm/_l_arm.sum();
    
    Eigen::Matrix4d mat=Eigen::Matrix4d::Identity(4,4);
    Eigen::MatrixXd joints = Eigen::MatrixXd::Zero(a.size()+1,3);
    
    for(size_t i=0; i<a.size();i++){
        
        Eigen::Matrix4d submat;
        submat<<cos(a(i)), -sin(a(i)), 0, _l_arm(i), sin(a(i)), cos(a(i)), 0 , 0, 0, 0, 1, 0, 0, 0, 0, 1;
        mat=mat*submat;

        Eigen::VectorXd joint = mat*Eigen::Vector4d(0,0,0,1);
        
        joints(i, 0) = joint[0];
        joints(i, 1) = joint[1];
        joints(i, 2) = joint[2];
    }

    Eigen::Matrix4d submat;
    submat<<1, 0, 0, _l_arm(a.size()), 0, 1, 0 , 0, 0, 0, 1, 0, 0, 0, 0, 1;
    mat=mat*submat;
    Eigen::VectorXd v=mat*Eigen::Vector4d(0,0,0,1); //gripper's position
    
    joints(a.size(),0)= v[0];
    joints(a.size(),1)= v[1];
    joints(a.size(),2)= v[2];
    
    Eigen::Vector2d v_joint_obst; //vector joint to obstacle
    Eigen::Vector2d v_joint_njoint; //vector joint to next joint
    
    bool collision = false;

    for(size_t i=0; i<a.size();i++){
        
        v_joint_obst[0] = obst[0] - joints(i,0);
        v_joint_obst[1] = obst[1] - joints(i,1);
        
        v_joint_njoint[0] = joints(i+1,0) - joints(i,0);
        v_joint_njoint[1] = joints(i+1,1) - joints(i,1);
        
        double angle = v_joint_obst[0]*v_joint_njoint[1] - v_joint_obst[1]*v_joint_njoint[0]; //cross product for a 3D space projected on a 2D plan
        
        if (abs(angle) < 0.01 && sqrt(v_joint_obst[0]*v_joint_obst[0] + v_joint_obst[1]*v_joint_obst[1]) < sqrt(v_joint_njoint[0]*v_joint_njoint[0] + v_joint_njoint[1]*v_joint_njoint[1])){
            collision = true;}
    }

    std::vector<double> out(4);
    
    out[0] = v[0];
    out[1] = v[1];
    out[2] = v[2];
    
    if (collision){
      std::cout << "collision" << std::endl;
        out[3] = 1;}
    
    else {
        out[3] = 0;}

    return out;
};

};

#endif
#endif
