/**
 * particle_filter.cpp
 *
 * Created on: Dec 12, 2016
 * Author: Tiffany Huang
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
using std::normal_distribution;
using std::discrete_distribution;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
  /**
   * TODO: Set the number of particles. Initialize all particles to 
   *   first position (based on estimates of x, y, theta and their uncertainties
   *   from GPS) and all weights to 1. 
   * TODO: Add random Gaussian noise to each particle.
   * NOTE: Consult particle_filter.h for more information about this method 
   *   (and others in this file).
   */
  
  //std::cout << "ParticleFilter::init()" << std::endl;

  num_particles = 9;  // TODO: Set the number of particles

  // Normal (Gaussian) distribution for x, y and theta using corresponding standard deviations
  normal_distribution<double> dist_x(x, std[0]);
  normal_distribution<double> dist_y(y, std[1]);
  normal_distribution<double> dist_theta(theta, std[2]);

  std::default_random_engine gen;
  Particle new_particle;
  
  // Initialize all weights to 1
  new_particle.weight = 1;

  for(int i = 0; i < num_particles; i++)
  {
    new_particle.id = i;
    // Sample from the normal distributions
    new_particle.x = dist_x(gen);
    new_particle.y = dist_y(gen);
    new_particle.theta = dist_theta(gen);
    
    // Print initial particle's position to the terminal.
    //std::cout << "Particle's initial position " << i << ": x=" << new_particle.x << " y=" << new_particle.y << " theta=" << new_particle.theta << std::endl;

    particles.push_back(new_particle);
  }

  weights.resize(num_particles);

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
  
  //std::cout << "ParticleFilter::prediction()" << std::endl;
  //std::cout << "delta_t=" << delta_t << " velocity=" << velocity << " yaw_rate=" << yaw_rate << " std_pos[0]=" << std_pos[0] << " std_pos[1]=" << std_pos[1] << " std_pos[2]=" << std_pos[2] << std::endl;

  std::default_random_engine gen;

  for(int i = 0; i < num_particles; i++)
  {
    // Print the original particle's position to the terminal.
    //std::cout << "Particle's original position " << i << " x=" << particles[i].x << " y=" << particles[i].y << " theta=" << particles[i].theta << std::endl;
    
    // Calculate new value of x using motion model
    double new_x = particles[i].x;
    if(yaw_rate != 0)
      new_x += (velocity / yaw_rate) * (sin(particles[i].theta + yaw_rate * delta_t) - sin(particles[i].theta));
    else
      new_x += velocity * delta_t * cos(particles[i].theta);
    // Sample from the normal (Gaussian) distribution around the new value and update the particle
    normal_distribution<double> dist_x(new_x, std_pos[0]);
    particles[i].x = dist_x(gen);
    
    // Calculate new value of y using motion model
    double new_y = particles[i].y;
    if(yaw_rate != 0)
      new_y += (velocity / yaw_rate) * (cos(particles[i].theta) - cos(particles[i].theta + yaw_rate * delta_t));
    else
      new_y += velocity * delta_t * sin(particles[i].theta);
    // Sample from the normal (Gaussian) distribution around the new value and update the particle
    normal_distribution<double> dist_y(new_y, std_pos[1]);
    particles[i].y = dist_y(gen);
    
    // Calculate new value of theta using motion model
    double new_theta = particles[i].theta + yaw_rate * delta_t;
    // Sample from the normal (Gaussian) distribution around the new value and update the particle
    normal_distribution<double> dist_theta(new_theta, std_pos[2]);
    particles[i].theta = dist_theta(gen);
    
    // Print the new particle's position to the terminal.
    //std::cout << "Particle's new position " << i << " x=" << particles[i].x << " y=" << particles[i].y << " theta=" << particles[i].theta << std::endl;
  }
}

void ParticleFilter::dataAssociation(vector<LandmarkObs> predicted, 
                                     vector<LandmarkObs>& observations) {
  /**
   * TODO: Find the predicted measurement that is closest to each 
   *   observed measurement and assign the observed measurement to this 
   *   particular landmark.
   * NOTE: this method will NOT be called by the grading code. But you will 
   *   probably find it useful to implement this method and use it as a helper 
   *   during the updateWeights phase.
   */
  
  //std::cout << "ParticleFilter::dataAssociation()" << std::endl;

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
  
  //std::cout << "ParticleFilter::updateWeights()" << std::endl;
  
  double total_weight = 0;

  for(int i = 0; i < num_particles; i++)
  {
    //std::cout << "Particle #" << i << std::endl;

    // init weight for the current particle
    particles[i].weight = 1.0;
    
    // reset associations and senor measurements for the current particle
    particles[i].associations.resize(0);
    particles[i].sense_x.resize(0);
    particles[i].sense_y.resize(0);

    for(size_t j = 0; j < observations.size(); j++)
    {
      // Convert observations from vehicle's (particle's) coordinate system into map's coordinate system
      double x_p_m = particles[i].x + observations[j].x * cos(particles[i].theta) - observations[j].y * sin(particles[i].theta);
      double y_p_m = particles[i].y + observations[j].x * sin(particles[i].theta) + observations[j].y * cos(particles[i].theta);
      //std::cout << "Transformed coordinates of the observation #" << j << " x=" << x_p_m << " y=" << y_p_m << std::endl;

      // Find the landmark which is closest to this observation
      double distance = -1;
      Map::single_landmark_s landmark;
      for(size_t k = 0; k < map_landmarks.landmark_list.size(); k++)
      {
        // Consider only the landmarks which are located in the range of the sensor from the current particle
        double particle_to_landmark_distance = dist(particles[i].x, particles[i].y, map_landmarks.landmark_list[k].x_f, map_landmarks.landmark_list[k].y_f);
        if(particle_to_landmark_distance <= sensor_range)
        {
          // Calculate distance between current observation and each landmark
          double distance_k = dist(x_p_m, y_p_m, map_landmarks.landmark_list[k].x_f, map_landmarks.landmark_list[k].y_f);
          // Find landmark at the minimal distabnce to the observation
          if(distance < 0 || distance_k < distance)
          {
            distance = distance_k;
            landmark = map_landmarks.landmark_list[k];
            //std::cout << "Using landmark #" << k << " x=" << landmark.x_f << " y=" << landmark.y_f << std::endl;
          }
        }
        //else
        //  std::cout << "Landmark " << k << "is out of range for particle" << i << std::endl;
      }

      double obs_weight = 0;
      // if the landmark was found
      if(distance >= 0)
      {
        // Calculate weight of the observation
        double weight_denom = (1 / (2 * M_PI * std_landmark[0] * std_landmark[1]));
        double dx = x_p_m - landmark.x_f;
        double sq_x = dx * dx / (2 * std_landmark[0] * std_landmark[0]);
        double dy = y_p_m - landmark.y_f;
        double sq_y = dy * dy / (2 * std_landmark[1] * std_landmark[1]);
        double exp_xy = exp(-(sq_x + sq_y));
        obs_weight = weight_denom * exp_xy;
        //std::cout << "weight_denom = " << weight_denom << " sq_x = " << sq_x << " sq_y = " << sq_y << " exp_xy = " << exp_xy << " obs_weight = " << obs_weight << std::endl;
        // fill in association to the current landmark and sensor measurements
        particles[i].associations.push_back(landmark.id_i);
        particles[i].sense_x.push_back(x_p_m);
        particles[i].sense_y.push_back(y_p_m);
      }
      //else
      //  std::cout << "No landmark found for particle " << i << ", observation " << j << std::endl;
      
      // Multiply weights of all observations to compute total weight of the particle
      particles[i].weight *= obs_weight;
      //std::cout << "Observation weight = " << obs_weight << " particle weight = " << particles[i].weight << std::endl;
    }
    
    // Calculate total weight of all particles - for later normalization of the weights
    total_weight += particles[i].weight;
  }
  
  //std::cout << "Total weight = " << total_weight << std::endl;
  
  if(total_weight > 0)
  {
    // Normalize weights of all particles
    for(int i = 0; i < num_particles; i++)
    {
      particles[i].weight /= total_weight;
      
      // Print particle's weight to the terminal.
      //std::cout << "Particle's " << i << " weight=" << particles[i].weight << std::endl;
    }
  }
}

void ParticleFilter::resample() {
  /**
   * TODO: Resample particles with replacement with probability proportional 
   *   to their weight. 
   * NOTE: You may find std::discrete_distribution helpful here.
   *   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
   */
  
  //std::cout << "ParticleFilter::resample()" << std::endl;
  
  for(int i = 0; i < num_particles; i++)
    weights[i] = particles[i].weight;

  std::default_random_engine gen;
  std::discrete_distribution<int> distribution(weights.begin(), weights.end());

  std::vector<Particle> new_particles;
  
  for(int i = 0; i < num_particles; i++)
  {
    int number = distribution(gen);
    new_particles.push_back(particles[number]);
    
    // Print selected particle's position to the terminal.
    //std::cout << "Selected particle " << number << " for position " << i << ", coordinates x=" << particles[number].x << " y=" << particles[number].y << " theta=" << particles[number].theta << std::endl;
  }

  particles = new_particles;
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
  
  //std::cout << "ParticleFilter::SetAssociations()" << std::endl;
  
  particle.associations= associations;
  particle.sense_x = sense_x;
  particle.sense_y = sense_y;
}

string ParticleFilter::getAssociations(Particle best) {
  
  //std::cout << "ParticleFilter::getAssociations()" << std::endl;
  
  vector<int> v = best.associations;
  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<int>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}

string ParticleFilter::getSenseCoord(Particle best, string coord) {
  
  //std::cout << "ParticleFilter::getSenseCoord()" << std::endl;
  
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
