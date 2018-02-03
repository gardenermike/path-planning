#include <fstream>
#include <math.h>
#include <uWS/uWS.h>
#include <chrono>
#include <iostream>
#include <thread>
#include <vector>
#include "Eigen-3.3/Eigen/Core"
#include "Eigen-3.3/Eigen/QR"
#include "json.hpp"
#include "spline.h"

using namespace std;

// for convenience
using json = nlohmann::json;

// For converting back and forth between radians and degrees.
constexpr double pi() { return M_PI; }
double deg2rad(double x) { return x * pi() / 180; }
double rad2deg(double x) { return x * 180 / pi(); }
const double meters_per_timestep = 0.43; // a little bit less than 50 mph
const double timestep_seconds = 0.05;

// Checks if the SocketIO event has JSON data.
// If there is data the JSON object in string format will be returned,
// else the empty string "" will be returned.
string hasData(string s) {
  auto found_null = s.find("null");
  auto b1 = s.find_first_of("[");
  auto b2 = s.find_first_of("}");
  if (found_null != string::npos) {
    return "";
  } else if (b1 != string::npos && b2 != string::npos) {
    return s.substr(b1, b2 - b1 + 2);
  }
  return "";
}

double distance(double x1, double y1, double x2, double y2)
{
	return sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1));
}
int ClosestWaypoint(double x, double y, const vector<double> &maps_x, const vector<double> &maps_y)
{

	double closestLen = 100000; //large number
	int closestWaypoint = 0;

	for(int i = 0; i < maps_x.size(); i++)
	{
		double map_x = maps_x[i];
		double map_y = maps_y[i];
		double dist = distance(x,y,map_x,map_y);
		if(dist < closestLen)
		{
			closestLen = dist;
			closestWaypoint = i;
		}

	}

	return closestWaypoint;

}

int NextWaypoint(double x, double y, double theta, const vector<double> &maps_x, const vector<double> &maps_y)
{

	int closestWaypoint = ClosestWaypoint(x,y,maps_x,maps_y);

	double map_x = maps_x[closestWaypoint];
	double map_y = maps_y[closestWaypoint];

	double heading = atan2((map_y-y),(map_x-x));

	double angle = fabs(theta-heading);
  angle = min(2*pi() - angle, angle);

  if(angle > pi()/4)
  {
    closestWaypoint++;
  if (closestWaypoint == maps_x.size())
  {
    closestWaypoint = 0;
  }
  }

  return closestWaypoint;
}

// Transform from Cartesian x,y coordinates to Frenet s,d coordinates
vector<double> getFrenet(double x, double y, double theta, const vector<double> &maps_x, const vector<double> &maps_y)
{
	int next_wp = NextWaypoint(x,y, theta, maps_x,maps_y);

	int prev_wp;
	prev_wp = next_wp-1;
	if(next_wp == 0)
	{
		prev_wp  = maps_x.size()-1;
	}

	double n_x = maps_x[next_wp]-maps_x[prev_wp];
	double n_y = maps_y[next_wp]-maps_y[prev_wp];
	double x_x = x - maps_x[prev_wp];
	double x_y = y - maps_y[prev_wp];

	// find the projection of x onto n
	double proj_norm = (x_x*n_x+x_y*n_y)/(n_x*n_x+n_y*n_y);
	double proj_x = proj_norm*n_x;
	double proj_y = proj_norm*n_y;

	double frenet_d = distance(x_x,x_y,proj_x,proj_y);

	//see if d value is positive or negative by comparing it to a center point

	double center_x = 1000-maps_x[prev_wp];
	double center_y = 2000-maps_y[prev_wp];
	double centerToPos = distance(center_x,center_y,x_x,x_y);
	double centerToRef = distance(center_x,center_y,proj_x,proj_y);

	if(centerToPos <= centerToRef)
	{
		frenet_d *= -1;
	}

	// calculate s value
	double frenet_s = 0;
	for(int i = 0; i < prev_wp; i++)
	{
		frenet_s += distance(maps_x[i],maps_y[i],maps_x[i+1],maps_y[i+1]);
	}

	frenet_s += distance(0,0,proj_x,proj_y);

	return {frenet_s,frenet_d};

}

// Transform from Frenet s,d coordinates to Cartesian x,y
vector<double> getXY(double s, double d, const vector<double> &maps_s, const vector<double> &maps_x, const vector<double> &maps_y)
{
	int prev_wp = -1;

	while(s > maps_s[prev_wp+1] && (prev_wp < (int)(maps_s.size()-1) ))
	{
		prev_wp++;
	}

	int wp2 = (prev_wp+1)%maps_x.size();

	double heading = atan2((maps_y[wp2]-maps_y[prev_wp]),(maps_x[wp2]-maps_x[prev_wp]));
	// the x,y,s along the segment
	double seg_s = (s-maps_s[prev_wp]);

	double seg_x = maps_x[prev_wp]+seg_s*cos(heading);
	double seg_y = maps_y[prev_wp]+seg_s*sin(heading);

	double perp_heading = heading-pi()/2;

	double x = seg_x + d*cos(perp_heading);
	double y = seg_y + d*sin(perp_heading);

	return {x,y};

}

/** 
 *  find a trajectory between points that minimizes jerk
 *  
 *  start: the start state given as a length three vectory
 *  corresponding to initial values of [s, s_dot, s_double_dot]
 *
 *  end: the desired end state. Like "start," this is a
 *  length three vector.
 *
 *  T: The duration, in seconds, over which this maneuver should occur.
 *
 *  OUTPUT 
 *  an array of length 6, each value corresponding to a coefficent in the polynomial 
 *  s(t) = a_0 + a_1 * t + a_2 * t**2 + a_3 * t**3 + a_4 * t**4 + a_5 * t**5
 */
vector<double> jerkMinimizingTrajectory(vector<double> start, vector<double> end, double T) {
  double si = start[0];
  double si_dot = start[1];
  double si_dot_dot = start[2];
  double sf = end[0];
  double sf_dot = end[1];
  double sf_dot_dot = end[2];
  
  double computed_1 = sf - (si + si_dot * T + 0.5 * si_dot_dot * T * T);
  double computed_2 = sf_dot - (si_dot + si_dot_dot * T);
  double computed_3 = sf_dot_dot - si_dot_dot;
  
  Eigen::MatrixXd A(3, 3);
  A << pow(T, 3), pow(T, 4), pow(T, 5),
       3 * pow(T, 2), 4 * pow(T, 3), 5 * pow(T, 4),
       6 * T, 12 * pow(T, 2), 20 * pow(T, 3);
  Eigen::VectorXd b(3);
  b << computed_1, computed_2, computed_3;
  Eigen::VectorXd x = A.colPivHouseholderQr().solve(b);
  
  vector<double> output = {si, si_dot, si_dot_dot / 2, x(0), x(1), x(2)};
  return output;
}

int main() {
  uWS::Hub h;

  // Load up map values for waypoint's x,y,s and d normalized normal vectors
  vector<double> map_waypoints_x;
  vector<double> map_waypoints_y;
  vector<double> map_waypoints_s;
  vector<double> map_waypoints_dx;
  vector<double> map_waypoints_dy;

  // Waypoint map to read from
  string map_file_ = "../data/highway_map.csv";
  // The max s value before wrapping around the track back to 0
  double max_s = 6945.554;

  ifstream in_map_(map_file_.c_str(), ifstream::in);

  string line;
  while (getline(in_map_, line)) {
  	istringstream iss(line);
  	double x;
  	double y;
  	float s;
  	float d_x;
  	float d_y;
  	iss >> x;
  	iss >> y;
  	iss >> s;
  	iss >> d_x;
  	iss >> d_y;
  	map_waypoints_x.push_back(x);
  	map_waypoints_y.push_back(y);
  	map_waypoints_s.push_back(s);
  	map_waypoints_dx.push_back(d_x);
  	map_waypoints_dy.push_back(d_y);
  }

  int lane = 1;
  auto lane_change_timestamp = std::chrono::system_clock::now();

  h.onMessage(
      [
        &map_waypoints_x,
        &map_waypoints_y,
        &map_waypoints_s,
        &map_waypoints_dx,
        &map_waypoints_dy,
        &lane,
        &lane_change_timestamp
      ](
        uWS::WebSocket<uWS::SERVER> ws,
        char *data,
        size_t length,
        uWS::OpCode opCode
        ) {
    // "42" at the start of the message means there's a websocket message event.
    // The 4 signifies a websocket message
    // The 2 signifies a websocket event
    //auto sdata = string(data).substr(0, length);
    //cout << sdata << endl;
    if (length && length > 2 && data[0] == '4' && data[1] == '2') {

      auto s = hasData(data);

      if (s != "") {
        auto j = json::parse(s);
        
        string event = j[0].get<string>();
        
        if (event == "telemetry") {
          // j[1] is the data JSON object
          
        	// Main car's localization Data
          double car_x = j[1]["x"];
          double car_y = j[1]["y"];
          double car_s = j[1]["s"];
          double car_d = j[1]["d"];
          double car_yaw = j[1]["yaw"];
          double car_speed = j[1]["speed"];

          // Previous path data given to the Planner
          auto previous_path_x = j[1]["previous_path_x"];
          auto previous_path_y = j[1]["previous_path_y"];
          // Previous path's end s and d values 
          double end_path_s = j[1]["end_path_s"];
          double end_path_d = j[1]["end_path_d"];

          // Sensor Fusion Data, a list of all other cars on the same side of the road.
          auto sensor_fusion = j[1]["sensor_fusion"];

          // simple state machine states
          bool straight = true;
          bool slow_down = false;
          bool turn_left = false;
          bool turn_right = false;
          double other_vehicle_velocity = car_speed;
          double other_vehicle_distance = 0;

          vector<vector<double>> vehicles_lane_left;
          vector<vector<double>> vehicles_my_lane;
          vector<vector<double>> vehicles_lane_right;

          // figure out state based on other vehicles
          // only change plans if we have not already committed to a lane change
          if ((std::chrono::system_clock::now() - lane_change_timestamp).count() > 1) {
            for (int i = 0; i < sensor_fusion.size(); i++) {
              auto other = sensor_fusion[i];
              //cout << other << endl;
              int other_index = other[0];
              double other_x = other[1];
              double other_y = other[2];
              double other_vx = other[3];
              double other_vy = other[4];
              double other_velocity = sqrt(pow(other_vx, 2) + pow(other_vy, 2)); // scalar velocity
              double other_s = other[5];
              double other_d = other[6];

              //for any vehicles on my side of the road
              if (other_d > 0) {
                int other_lane = (int)(other_d) / 4;

                if (other_lane == lane - 1) {
                  vehicles_lane_left.push_back(other);
                } else if (other_lane == lane) {
                  vehicles_my_lane.push_back(other);
                } else if (other_lane == lane + 1) {
                  vehicles_lane_right.push_back(other);
                }

              }
            }

            // are there any vehicles in my lane?
            double safe_distance = 50;
            for (int i = 0; i < vehicles_my_lane.size(); i++) {
              auto other = vehicles_my_lane[i];
              int other_index = other[0];
              double other_x = other[1];
              double other_y = other[2];
              double other_vx = other[3];
              double other_vy = other[4];
              double other_velocity = sqrt(pow(other_vx, 2) + pow(other_vy, 2)); // scalar velocity
              double other_s = other[5];
              double other_d = other[6];

              double s_distance = other_s - car_s;

              // Are we close enough to worry?
              if (s_distance > 0 && s_distance < safe_distance) {
                if (other_velocity < other_vehicle_velocity) {
                  //cout << "vehicle detected! distance: " << s_distance << "\n";
                  other_vehicle_velocity = other_velocity;
                  other_vehicle_distance = s_distance;
                  straight = false;
                  slow_down = true;
                  turn_left = false;
                  turn_right = false;
                }
              }
            }

            if (slow_down) {
              // first, try left lane
              if (lane > 0) {
                bool other_vehicle_found = false;
                for (int i = 0; i < vehicles_lane_left.size(); i++) {
                  auto other = vehicles_lane_left[i];
                  int other_index = other[0];
                  double other_x = other[1];
                  double other_y = other[2];
                  double other_vx = other[3];
                  double other_vy = other[4];
                  double other_velocity = sqrt(pow(other_vx, 2) + pow(other_vy, 2)); // scalar velocity
                  double other_s = other[5];
                  double other_d = other[6];

                  // is the car in a dangerous range?
                  if (other_s > car_s - 20 && other_s < car_s + safe_distance) {
                    // second check... is it going faster?
                    if (other_s > car_s + 20 && other_velocity > car_speed) {
                      // ok after all
                    } else {
                      other_vehicle_found = true;
                      break;
                    }
                  }
                }

                if (!other_vehicle_found) {
                  straight = false;
                  slow_down = false;
                  turn_left = true;
                  turn_right = false;
                  lane_change_timestamp = std::chrono::system_clock::now();
                }
              }

              // try right lane turn
              if (lane < 2 && !turn_left) {
                bool other_vehicle_found = false;
                for (int i = 0; i < vehicles_lane_right.size(); i++) {
                  auto other = vehicles_lane_right[i];
                  int other_index = other[0];
                  double other_x = other[1];
                  double other_y = other[2];
                  double other_vx = other[3];
                  double other_vy = other[4];
                  double other_velocity = sqrt(pow(other_vx, 2) + pow(other_vy, 2)); // scalar velocity
                  double other_s = other[5];
                  double other_d = other[6];

                  // is the car in a dangerous range?
                  if (other_s > car_s - 20 && other_s < car_s + safe_distance) {
                    // second check... is it going faster?
                    //if (!(other_s < car_s + 15 && other_velocity > car_speed)) {
                    if (other_s > car_s + 20 && other_velocity > car_speed) {
                      // ok after all
                    } else {
                      other_vehicle_found = true;
                      break;
                    }
                  }
                }

                if (!other_vehicle_found) {
                  straight = false;
                  slow_down = false;
                  turn_left = false;
                  turn_right = true;
                  lane_change_timestamp = std::chrono::system_clock::now();
                }
              }
            }
          }

          vector<double> next_x_vals;
          vector<double> next_y_vals;


          vector<double> s_d;
          double current_x = car_x;
          double current_y = car_y;
          double car_yaw_radians = deg2rad(car_yaw);

          // convert waypoints into a vehicle-origin coordinate system
          // and map out the next several waypoints in that system
          int waypoint_count = 5;
          int waypoint_distance = 20;
          bool in_lane_change = (int)car_d / 4 != lane;
          if (turn_left || turn_right || in_lane_change) {
            waypoint_distance = 60; // smoother transitions while turning
          }
          std::vector<double> vehicle_x = {0.0}; // start at the vehicle location
          std::vector<double> vehicle_y = {0.0}; // start at the vehicle location
          for (int i = 0; i < waypoint_count; ++i) {
            // find the (x, y) coordinates of the frenet coordinate i positions into the future
            vector<double> waypoint = getXY(
              car_s + waypoint_distance * (i + 1),
              2 + 4 * lane,
              map_waypoints_s,
              map_waypoints_x,
              map_waypoints_y
            );

            double map_x = waypoint[0];
            double map_y = waypoint[1];

            double delta_x = map_x - car_x;
            double delta_y = map_y - car_y;
            // add the coordinates to the spline list, each transformed to vehicle coordinates
            vehicle_x.push_back(delta_x * cos(-car_yaw_radians) - delta_y * sin(-car_yaw_radians));
            vehicle_y.push_back(delta_x * sin(-car_yaw_radians) + delta_y * cos(-car_yaw_radians));
          }

          tk::spline waypoint_spline;
          waypoint_spline.set_points(vehicle_x, vehicle_y);

          // if we have a previous trajectory, use the first point in it as the starting point
          // for the new trajectory
          if (previous_path_x.size() > 1) {
            current_x = previous_path_x[0];
            current_y = previous_path_y[0];
          }

          //cout << "!(" << current_x << ", " << current_y << ")! ";
          double current_delta_x = current_x - car_x;
          double current_delta_y = current_y - car_y;
          double current_vehicle_x = current_delta_x * cos(-car_yaw_radians) - current_delta_y * sin(-car_yaw_radians);
          double current_vehicle_y = current_delta_x * sin(-car_yaw_radians) + current_delta_y * cos(-car_yaw_radians);

          //cout << "Car yaw radians: " << car_yaw_radians << " " << car_yaw <<  "\n";
          //cout << "^(" << current_vehicle_x << ", " << current_vehicle_y << ")^\n";
          
          int steps_in_trajectory = 128;
          int remaining_points = steps_in_trajectory - previous_path_x.size();
          int existing_size = previous_path_x.size();
          double current_target_distance = car_speed * 0.02 * 0.44704; // divide by timestep and convert to meters/second

          double target_meters_per_timestep = meters_per_timestep;

          // state machine actions
          if (slow_down) {
            target_meters_per_timestep = other_vehicle_velocity * 0.02 * 0.44704;
          } else if (turn_left) {
            cout << "lane change left!" << endl;
            lane -= 1;
          } else if (turn_right) {
            cout << "lane change right!" << endl;
            lane += 1;
          }

          // We have everything we need now! Create the trajectory.

          // how far have we gone in the trajectory?
          double distance_forecast = 0;

          // add the trajectory points
          for (int i = 0; i < steps_in_trajectory; i++) {
            if (i > 0) { // don't add distance for first step in previous path, which already has unknown distance
              if (current_target_distance < target_meters_per_timestep && straight) {
                current_target_distance += (target_meters_per_timestep - current_target_distance) / 250;
              } else {
                if (slow_down) {
                  // match velocity of vehicle ahead within fraction of distance to it
                  double target_diff = 0.1;
                  double velocity_diff = other_vehicle_velocity * 0.44704 - car_speed * 0.44704;
                  double seconds_to_impact = other_vehicle_distance / velocity_diff;
                  double timesteps_to_impact = seconds_to_impact / timestep_seconds; 
                  current_target_distance -= ((velocity_diff * timestep_seconds) / timesteps_to_impact) * target_diff;
                } else if (straight) {
                  current_target_distance *= 0.999; // gently reduce velocity
                }
              }
              distance_forecast += current_target_distance;
            }
            
            double x = current_vehicle_x + distance_forecast;
            double y = waypoint_spline(x);

            // convert back to map coordinates
            double transformed_x = x * cos(car_yaw_radians) - y * sin(car_yaw_radians);
            double transformed_y = x * sin(car_yaw_radians) + y * cos(car_yaw_radians);
            transformed_x += car_x;
            transformed_y += car_y;
          
            // blend with existing path
            if (i < existing_size) {
              double new_weight = i / (double)existing_size;
              double existing_weight = 1. - new_weight;
              transformed_x = transformed_x * new_weight + (double)previous_path_x[i] * existing_weight;
              transformed_y = transformed_y * new_weight + (double)previous_path_y[i] * existing_weight;
            }
            next_x_vals.push_back(transformed_x);
            next_y_vals.push_back(transformed_y);
          }


          // log trajectory
          //for (int i = 0; i < next_x_vals.size(); ++i) {
          //  cout << "(" << next_x_vals[i] << "," << next_y_vals[i] << ") \n";
          //}
          //cout << endl;

          // now that we have a reasonable trajectory,
          // it may have unreasonable acceleration at points.
          // cut it into pieces and jerk-minimize each component
          int trajectory_components = 4;
          int points_per_component = steps_in_trajectory / trajectory_components;
          //double current_velocity = car_speed * 0.44704; // convert to meters/second
          //double current_vehicle_x = current_delta_x * cos(-car_yaw_radians) - current_delta_y * sin(-car_yaw_radians);
          //double current_vehicle_y = current_delta_x * sin(-car_yaw_radians) + current_delta_y * cos(-car_yaw_radians);
          
          for (int i = 0; i < trajectory_components; i++) {
            // start point
            int start_index = i * points_per_component;
            double start_x_val = next_x_vals[start_index];
            double start_y_val = next_y_vals[start_index];

            //convert to vehicle coordinates
            double delta_x = start_x_val - car_x;
            double delta_y = start_y_val - car_y;
            double start_vehicle_x = delta_x * cos(-car_yaw_radians) - delta_y * sin(-car_yaw_radians);
            double start_vehicle_y = delta_x * sin(-car_yaw_radians) + delta_y * cos(-car_yaw_radians);
            // convert subsequent two points to determine velocity and acceleration
            double next_delta_x = next_x_vals[start_index + 1] - car_x;
            double next_delta_y = next_y_vals[start_index + 1] - car_y;
            double next_vehicle_x = next_delta_x * cos(-car_yaw_radians) - next_delta_y * sin(-car_yaw_radians);
            double next_vehicle_y = next_delta_x * sin(-car_yaw_radians) + next_delta_y * cos(-car_yaw_radians);
            double start_vehicle_x_velocity = (next_vehicle_x - start_vehicle_x) / timestep_seconds;
            double start_vehicle_y_velocity = (next_vehicle_y - start_vehicle_y) / timestep_seconds;
            double subsequent_delta_x = next_x_vals[start_index + 2] - car_x;
            double subsequent_delta_y = next_y_vals[start_index + 2] - car_y;
            double subsequent_vehicle_x = subsequent_delta_x * cos(-car_yaw_radians) - subsequent_delta_y * sin(-car_yaw_radians);
            double subsequent_vehicle_y = subsequent_delta_x * sin(-car_yaw_radians) + subsequent_delta_y * cos(-car_yaw_radians);
            double next_vehicle_x_velocity = (subsequent_vehicle_x - next_vehicle_x) / timestep_seconds;
            double next_vehicle_y_velocity = (subsequent_vehicle_y - next_vehicle_y) / timestep_seconds;
            double start_vehicle_x_acceleration = next_vehicle_x_velocity - start_vehicle_x_velocity;
            double start_vehicle_y_acceleration = next_vehicle_y_velocity - start_vehicle_y_velocity;

            // check for boundary cases
            if (isinf(start_vehicle_x_velocity) || isnan(start_vehicle_x_velocity)) start_vehicle_x_velocity = 0;
            if (isinf(start_vehicle_y_velocity) || isnan(start_vehicle_y_velocity)) start_vehicle_y_velocity = 0;
            if (isinf(start_vehicle_x_acceleration) || isnan(start_vehicle_x_acceleration)) start_vehicle_x_acceleration = 0;
            if (isinf(start_vehicle_y_acceleration) || isnan(start_vehicle_y_acceleration)) start_vehicle_y_acceleration = 0;

            // end point
            int end_index = start_index + points_per_component - 1;
            double end_x_val = next_x_vals[end_index];
            double end_y_val = next_y_vals[end_index];

            //convert to vehicle coordinates
            delta_x = end_x_val - car_x;
            delta_y = end_y_val - car_y;
            double end_vehicle_x = delta_x * cos(-car_yaw_radians) - delta_y * sin(-car_yaw_radians);
            double end_vehicle_y = delta_x * sin(-car_yaw_radians) + delta_y * cos(-car_yaw_radians);
            // convert previous two points to determine velocity and acceleration
            double previous_delta_x = next_x_vals[end_index - 1] - car_x;
            double previous_delta_y = next_y_vals[end_index - 1] - car_y;
            double previous_vehicle_x = previous_delta_x * cos(-car_yaw_radians) - previous_delta_y * sin(-car_yaw_radians);
            double previous_vehicle_y = previous_delta_x * sin(-car_yaw_radians) + previous_delta_y * cos(-car_yaw_radians);
            double end_vehicle_x_velocity = (previous_vehicle_x - end_vehicle_x) / timestep_seconds;
            double end_vehicle_y_velocity = (previous_vehicle_y - end_vehicle_y) / timestep_seconds;
            double antecedent_delta_x = next_x_vals[end_index - 2] - car_x;
            double antecedent_delta_y = next_y_vals[end_index - 2] - car_y;
            double antecedent_vehicle_x = antecedent_delta_x * cos(-car_yaw_radians) - antecedent_delta_y * sin(-car_yaw_radians);
            double antecedent_vehicle_y = antecedent_delta_x * sin(-car_yaw_radians) + antecedent_delta_y * cos(-car_yaw_radians);
            double previous_vehicle_x_velocity = (antecedent_vehicle_x - previous_vehicle_x) / timestep_seconds;
            double previous_vehicle_y_velocity = (antecedent_vehicle_y - previous_vehicle_y) / timestep_seconds;
            double end_vehicle_x_acceleration = previous_vehicle_x_velocity - end_vehicle_x_velocity;
            double end_vehicle_y_acceleration = previous_vehicle_y_velocity - end_vehicle_y_velocity;


            if (isinf(end_vehicle_x_velocity) || isnan(end_vehicle_x_velocity)) end_vehicle_x_velocity = 0;
            if (isinf(end_vehicle_y_velocity) || isnan(end_vehicle_y_velocity)) end_vehicle_y_velocity = 0;
            if (isinf(end_vehicle_x_acceleration) || isnan(end_vehicle_x_acceleration)) end_vehicle_x_acceleration = 0;
            if (isinf(end_vehicle_y_acceleration) || isnan(end_vehicle_y_acceleration)) end_vehicle_y_acceleration = 0;

            vector<double> start_x = {start_vehicle_x, start_vehicle_x_velocity, start_vehicle_x_acceleration};
            vector<double> end_x = {end_vehicle_x, end_vehicle_x_velocity, end_vehicle_x_acceleration};
            vector<double> start_y = {start_vehicle_y, start_vehicle_y_velocity, start_vehicle_y_acceleration};
            vector<double> end_y = {end_vehicle_y, end_vehicle_y_velocity, end_vehicle_y_acceleration};
            vector<double> x_trajectory_coefficients = jerkMinimizingTrajectory(start_x, end_x, points_per_component * timestep_seconds);

            /*for (int j = 0; j < start_x.size(); j++) {
              cout << start_x[j] << " ";
            }
            cout << endl;
            for (int j = 0; j < end_x.size(); j++) {
              cout << end_x[j] << " ";
            }
            cout << endl;
            for (int j = 0; j < x_trajectory_coefficients.size(); j++) {
              cout << x_trajectory_coefficients[j] << " ";
            }
            cout << endl;*/
            vector<double> y_trajectory_coefficients = jerkMinimizingTrajectory(start_y, end_y, points_per_component * timestep_seconds);

            for (int timestep = 1; i < points_per_component - 1; i++) {
              double elapsed_time = timestep * timestep_seconds;
              double vehicle_x = 
                x_trajectory_coefficients[0] +
                x_trajectory_coefficients[1] * elapsed_time +
                x_trajectory_coefficients[2] * pow(elapsed_time, 2) +
                x_trajectory_coefficients[3] * pow(elapsed_time, 3) +
                x_trajectory_coefficients[4] * pow(elapsed_time, 4) +
                x_trajectory_coefficients[5] * pow(elapsed_time, 5);
              double vehicle_y =
                y_trajectory_coefficients[0] +
                y_trajectory_coefficients[1] * elapsed_time +
                y_trajectory_coefficients[2] * pow(elapsed_time, 2) +
                y_trajectory_coefficients[3] * pow(elapsed_time, 3) +
                y_trajectory_coefficients[4] * pow(elapsed_time, 4) +
                y_trajectory_coefficients[5] * pow(elapsed_time, 5);
              //cout << "(" << vehicle_x << ", " << vehicle_y << ") ";


              // convert back to map coordinates
              double transformed_x = vehicle_x * cos(car_yaw_radians) - vehicle_y * sin(car_yaw_radians);
              double transformed_y = vehicle_x * sin(car_yaw_radians) + vehicle_y * cos(car_yaw_radians);
              transformed_x += car_x;
              transformed_y += car_y;
              int current_index = start_index + timestep;
              //if (transformed_x != next_x_vals[current_index]) {
              //  cout << transformed_x - next_x_vals[current_index] << endl;
              //}
              next_x_vals[current_index] = transformed_x;
              next_y_vals[current_index] = transformed_y;
            }
          }
          

          // log trajectory
          /*for (int i = 0; i < next_x_vals.size(); ++i) {
            cout << "(" << next_x_vals[i] << "," << next_y_vals[i] << ") \n";
          }
          cout << endl;*/

          json msgJson;
          msgJson["next_x"] = next_x_vals;
          msgJson["next_y"] = next_y_vals;

          auto msg = "42[\"control\","+ msgJson.dump()+"]";

          //this_thread::sleep_for(chrono::milliseconds(1000));
          ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
          
        }
      } else {
        // Manual driving
        std::string msg = "42[\"manual\",{}]";
        ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
      }
    }
  });

  // We don't need this since we're not using HTTP but if it's removed the
  // program
  // doesn't compile :-(
  h.onHttpRequest([](uWS::HttpResponse *res, uWS::HttpRequest req, char *data,
                     size_t, size_t) {
    const std::string s = "<h1>Hello world!</h1>";
    if (req.getUrl().valueLength == 1) {
      res->end(s.data(), s.length());
    } else {
      // i guess this should be done more gracefully?
      res->end(nullptr, 0);
    }
  });

  h.onConnection([&h](uWS::WebSocket<uWS::SERVER> ws, uWS::HttpRequest req) {
    std::cout << "Connected!!!" << std::endl;
  });

  h.onDisconnection([&h](uWS::WebSocket<uWS::SERVER> ws, int code,
                         char *message, size_t length) {
    ws.close();
    std::cout << "Disconnected" << std::endl;
  });

  int port = 4567;
  if (h.listen(port)) {
    std::cout << "Listening to port " << port << std::endl;
  } else {
    std::cerr << "Failed to listen to port" << std::endl;
    return -1;
  }
  h.run();
}
