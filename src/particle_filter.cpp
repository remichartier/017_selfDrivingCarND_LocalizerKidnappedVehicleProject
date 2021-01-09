/**
 * particle_filter.cpp
 *
 * Created on: Dec 12, 2016
 * Author: Tiffany Huang
 */

/**
 * History
 * v00 : initial file
 * v01 : ParticleFilter::init() done
 * v02 : ParticleFilter::prediction() done + factorize Gaussian Noise 
 *       Distribution with applyGaussianNoise()
 * v03 : complete dataAssociation(), homogeneousTransform(),
 *			Progress on updateWeights()
 * v04 : completed updateWeights(), added filterOutOfSensorRangeLandmarks()
 * v05 : Add normalize_weights()
 * v06 : Add resample()
 * v07 : add Division by 0 Error management and debug prints
 */

#include "particle_filter.h"

#include <math.h>
#include <algorithm>
#include <iostream>
#include <iterator>
#include <numeric>
#include <random>
#include <string>
#include <vector>

#include <stdlib.h> // for exit() function
#include <limits> // for std::numeric_limits<double>::min()


#include "helper_functions.h"

using std::string;
using std::vector;

#define DEBUG01 false


void ParticleFilter::init(double x, double y, double theta, double std[]) {
  /**
   * TODO: Set the number of particles. Initialize all particles to 
   *   first position (based on estimates of x, y, theta and their uncertainties
   *   from GPS) and all weights to 1. 
   * TODO: Add random Gaussian noise to each particle.
   * NOTE: Consult particle_filter.h for more information about this method 
   *   (and others in this file).
   */
  num_particles = 100;  // TODO: Set the number of particles
  
  /*
  std::default_random_engine gen; // random engine to create a normal (Gaussian) distributions for x,y,theta
  // This line creates a normal (Gaussian) distribution for x
  //normal_distribution<double> dist_x(x, std[0]);
  // This line creates a normal (Gaussian) distribution for y
  normal_distribution<double> dist_y(y, std[1]);
  // This line creates a normal (Gaussian) distribution for theta
  normal_distribution<double> dist_theta(theta, std[2]);
  */
  for (int i = 0; i < num_particles; ++i) {
    struct Particle p;
    p.id = i;    
    // Sample from these normal distributions like this: 
    //   sample_x = dist_x(gen);
    //   where "gen" is the random engine initialized earlier.
    // p.x = dist_x(gen);
    // p.y = dist_y(gen);
    // p.theta = dist_theta(gen);
    p.x = x;
    p.y = y;
    p.theta = theta;
    applyGaussianNoise(p.x, p.y, p.theta, std);
    p.weight = 1;
    particles.push_back(p);
    
    // Also size up the Vector of weights of all particles
    weights.push_back(0.0); // may have to use this-> if not recognized ...
  }
  
  // set initialized flag
  is_initialized = true;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], 
                                double velocity, double yaw_rate) {
  /**
   * TODO: Add measurements to each particle and add random Gaussian noise.
   * NOTE: When adding noise you may find std::normal_distribution 
   *   and std::default_random_engine useful.
   *  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
   *  http://www.cplusplus.com/reference/random/default_random_engine/
   */
  
  // std_pos[] : velocity and yaw rate measurement velocities.
  if(yaw_rate == 0) {
    std::cout << "ParticleFilter::prediction() ERROR : yaw_rate=0 --> dividing by 0" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  // calculate new particle positions
  for (int i = 0; i < num_particles; ++i) {
    double d = yaw_rate * delta_t;
    double a = velocity/yaw_rate;
    double b = particles[i].theta + d;
    particles[i].x += a * (sin(b) - sin(particles[i].theta));
    particles[i].y += a * (cos(particles[i].theta) - cos(b));
    particles[i].theta += b;
    // Apply Gaussian Noise Normal Distribution to those new predicted positions.
    applyGaussianNoise(particles[i].x, particles[i].y, particles[i].theta, std_pos);
  }
}


void ParticleFilter::dataAssociation(vector<LandmarkObs> predicted,
                                     vector<LandmarkObs>& obs) {
  /**
   * TODO: Find the predicted measurement that is closest to each 
   *   observed measurement and assign the observed measurement to this 
   *   particular landmark.
   * NOTE: this method will NOT be called by the grading code. But you will 
   *   probably find it useful to implement this method and use it as a helper 
   *   during the updateWeights phase.
   */
  
  for(unsigned int o=0; o<obs.size(); ++o){  // for each observation
    // build distance vector : dist from (obs[o].x,obs[o].y) to each map_landmarks.
    vector<double> distance;
    for(unsigned int l=0; l < predicted.size(); ++l){
      distance.push_back(dist(obs[o].x, obs[o].y, predicted[l].x, predicted[l].y));
      // so here we have vector distance storing distances between observation and each landmark_list[l]
      // now need to pick the the minimum distance, retrieve the index of the landmark
      // add this index to observation[i].id...
    }
    // so here we have vector distance storing distances between observation and each landmark_list[l]
    // now need to pick the the minimum distance, retrieve the index of the landmark
    // add this index to obs[o].id...
    obs[o].id = min_index(distance);
  }
}  


void homogeneousTransform(Particle p, vector<LandmarkObs> obs, vector<LandmarkObs>& obs_map){
  /**
   * Transform observations from vehicle coordinates into maps coordinates for a particle p
   * Inputs : particule p, observations vector obs (vehicle coordinates), obs_map (copy of obs)
   * Output : observations in map coordinates relative to particule p : obs_map
   */
  for(unsigned int o=0; o<obs.size(); ++o){
      // use homogeneous equations, using x/y from the particule
      // x,y coordinates of the particule + theta heading of the particule
      // How to calculate theta ???? This is the theta from the particule ...
      // transform to map x coordinate
      obs_map[o].x = p.x + (cos(p.theta) * obs[o].x) - (sin(p.theta) * obs[o].y);
      // transform to map y coordinate
      obs_map[o].y = p.y + (sin(p.theta) * obs[o].x) + (cos(p.theta) * obs[o].y); 
    }
}

void filterOutOfSensorRangeLandmarks(const Map &map_landmarks,
                                     vector<LandmarkObs>& predicted,
                                     double sensor_range, Particle p){
  /**
   * filterOutOfSensorRangeLandmarks() filters landmarks and only keeps 
   * the ones at distance <= sensor_range of the particle filter p
   * @param map_landmarks vector of landmark coordinates and id (x,y,id)
   * @param predicted Vector of predicted landmark observations
   * @param sensor_range (sensor range in meters)
   * @param p particle x,y,id to compare distance with landmark
   *
   * Output : predicted
   */ 
  LandmarkObs obs;  
  
  // for each map landmark
  for(unsigned int l=0; l < map_landmarks.landmark_list.size(); ++l){
    double xl(map_landmarks.landmark_list[l].x_f), yl(map_landmarks.landmark_list[l].y_f);
    if(dist(p.x, p.y, xl, yl) <= sensor_range){
      obs.x = xl;
      obs.y = yl;
      obs.id = map_landmarks.landmark_list[l].id_i;
      predicted.push_back(obs);
    }
  } 
}

/**
 * normalize particle weights
 */

bool ParticleFilter::normalize_weights(double sum){
  
  if(sum==0){
  	std::cout << "ParticleFilter::normalize_weights() : ERROR division by 0" << std::endl;
    return false;
  }
  for(unsigned int i=0; i<(unsigned int)num_particles; ++i){
  	particles[i].weight /= sum;
  }
  return true;
}


void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
                                   const vector<LandmarkObs> &observations, 
                                   const Map &map_landmarks) {
  /**
   * TODO: Update the weights of each particle using a mult-variate Gaussian 
   *   distribution. You can read more about this distribution here: 
   *   https://en.wikipedia.org/wiki/Multivariate_normal_distribution
   * NOTE: The observations are given in the VEHICLE'S coordinate system. 
   *   Your particles are located according to the MAP'S coordinate system. 
   *   You will need to transform between the two systems. Keep in mind that
   *   this transformation requires both rotation AND translation (but no scaling).
   *   The following is a good resource for the theory:
   *   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
   *   and the following is a good resource for the actual equation to implement
   *   (look at equation 3.33) http://planning.cs.uiuc.edu/node99.html
   */
  // transform to map x coordinate of observations
  // but only need to do it one time, with one particule
  
  double weight_sum(0.0); // sum of weight to later normalize them before resampling
  
  for(unsigned int i=0; i<(unsigned int)num_particles; ++i){  // For each particle in the particle filter
    // need to transform each observation vehicle coordinates into map coordinates
    
    // vector<int> vect2(vect1.begin(), vect1.end()); 
    vector<LandmarkObs> obs_m(observations.begin(), observations.end()); 
    
    // homogeneousTransform(this->particles[i], observations, obs_m); 
	homogeneousTransform(particles[i], observations, obs_m); 
    // Now we have observations in map coordinates relative to particle p, in obs_m vector
    
    // Filter out landmarks out of sensor_range vs particle position
    vector<LandmarkObs> predicted;
    // filterOutOfSensorRangeLandmarks(map_landmarks, predicted, sensor_range,
    //                                this->particles[i]);
    filterOutOfSensorRangeLandmarks(map_landmarks, predicted, sensor_range,
                                    particles[i]);
    
    // We can then look for nearest landmark of this observations
    // and chose the nearest landmark for each observation, and this
    // is particle specific
    dataAssociation(predicted, obs_m);
    
    // now we have the closest landmarks to each observation on obs_m[i].id
    // corresponding to landmark index of predicted
    
    vector<double> prob(obs_m.size());
    multivariate_gaussian_probability(prob, obs_m, std_landmark,
                                      predicted);
	
    if(DEBUG01) printVectorDouble(prob, "prob"); // for debug
    
    
    // Need now to calculate particle weight ...
    // compute_particle_weight(this->particles[i].weight, prob);
    
    if(DEBUG01){
      if(particles[i].weight == 0){
        std::cout << "particles[i].weight = 0" << std::endl;
      }
    }
    
    // Also update the particle overall weight vector.
    // weights[i] = this->particles[i].weight;
    weights[i] = particles[i].weight;
    
    // to normalize weights at the end
    // weight_sum += this->particles[i].weight; 
    weight_sum += particles[i].weight; 
    
    // to work around in case weight_sum becomes 0.0 due to probs much too low : 
    if(weight_sum == 0){
         weight_sum = std::fabs(std::numeric_limits<double>::min());
    }

    
  }
  // at this point, all the particule weights have been calculated.
  // We need to normalize them to use them as probabilities in resampling step
  // so that weights are between [0,1]
  if(DEBUG01) std::cout << "weight_sum = " << weight_sum << std::endl;
  
  bool status = normalize_weights(weight_sum);
  if(!status){
    // need to print debug info to show why divide by 0
    // print prob vector
    // printVectorDouble(prob, "prob");
    std::exit(EXIT_FAILURE);
  }

}

void ParticleFilter::resample() {
  /**
   * TODO: Resample particles with replacement with probability proportional 
   *   to their weight. 
   * NOTE: You may find std::discrete_distribution helpful here.
   *   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
   */
  
  // vector<int> vect2(vect1.begin(), vect1.end()); 
  // copy particules in new vector as will draw and replace particles vector
  //vector<Particle> old_particles(this->particles.begin(),this->particles.end());
  vector<Particle> old_particles(particles.begin(),particles.end());
  
  std::default_random_engine generator;
  std::discrete_distribution<int> distribution (weights.begin(),weights.end());
  
  for(unsigned int i=0; i<(unsigned int)num_particles; ++i){
  	int number = distribution(generator);
    this->particles[i] = old_particles[number];
  }
  
  // now this->particles[i] contain new draw of num_particles from previous set
  // picked randomly according to their previous weights

}

void ParticleFilter::SetAssociations(Particle& particle, 
                                     const vector<int>& associations, 
                                     const vector<double>& sense_x, 
                                     const vector<double>& sense_y) {
  // particle: the particle to which assign each listed association, 
  //   and association's (x,y) world coordinates mapping
  // associations: The landmark id that goes along with each listed association
  // sense_x: the associations x mapping already converted to world coordinates
  // sense_y: the associations y mapping already converted to world coordinates
  particle.associations= associations;
  particle.sense_x = sense_x;
  particle.sense_y = sense_y;
}

string ParticleFilter::getAssociations(Particle best) {
  vector<int> v = best.associations;
  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<int>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}

string ParticleFilter::getSenseCoord(Particle best, string coord) {
  vector<double> v;

  if (coord == "X") {
    v = best.sense_x;
  } else {
    v = best.sense_y;
  }

  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<float>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}