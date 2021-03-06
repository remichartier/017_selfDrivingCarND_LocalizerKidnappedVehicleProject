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

#include "helper_functions.h"

using std::string;
using std::vector;


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
    return;
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

void ParticleFilter::dataAssociation(Map map_landmarks, 
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
    for(unsigned int l=0; l < map_landmarks.landmark_list.size(); ++l){
      double xl(map_landmarks.landmark_list[l].x_f), yl(map_landmarks.landmark_list[l].y_f);
      distance.push_back(dist(obs[o].x, obs[o].y, xl, yl));
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
  
  for(unsigned int i=0; i<(unsigned int)num_particles; ++i){  // For each particle in the particle filter
    // need to transform each observation vehicle coordinates into map coordinates
    
    // vector<int> vect2(vect1.begin(), vect1.end()); 
    vector<LandmarkObs> obs_m(observations.begin(), observations.end()); 
    
    homogeneousTransform(this->particles[i], observations, obs_m); 
    // Now we have observations in map coordinates relative to particle p, in obs_m vector
    
    // We can then look for nearest landmark of this observations
    // and chose the nearest landmark for each observation, and this
    // is particle specific
    
    dataAssociation(map_landmarks, obs_m);
    
    // now we have the closest landmarks to each observation on obs_m[i].id

    // Need now to calculate particle weight ...

  }  
}

void ParticleFilter::resample() {
  /**
   * TODO: Resample particles with replacement with probability proportional 
   *   to their weight. 
   * NOTE: You may find std::discrete_distribution helpful here.
   *   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
   */

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