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

  h.onMessage(
      [
        &map_waypoints_x,
        &map_waypoints_y,
        &map_waypoints_s,
        &map_waypoints_dx,
        &map_waypoints_dy,
        &lane
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

          	json msgJson;

          	vector<double> next_x_vals;
          	vector<double> next_y_vals;


          	// TODO: define a path made up of (x,y) points that the car will visit sequentially every .02 seconds
            //
            // 1. Define a smooth path around the waypoints
            // 2. Find a trajectory close to the smooth path
            // 3. Plan a course to follow the trajectory
            // 4. Return points along the course
            //


            vector<double> s_d;
            double current_x = car_x;
            double current_y = car_y;
            double car_yaw_radians = deg2rad(car_yaw);

            int waypoint_index = NextWaypoint(current_x, current_y, car_yaw_radians, map_waypoints_x, map_waypoints_y);
					
            // convert waypoints into a vehicle-origin coordinate system
            // and map out the next several waypoints in that system
            int waypoint_count = 4;
            int waypoint_distance = 25;
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
							vehicle_x.push_back(delta_x * cos(-car_yaw_radians) - delta_y * sin(-car_yaw_radians));
							vehicle_y.push_back(delta_x * sin(-car_yaw_radians) + delta_y * cos(-car_yaw_radians));
						}

            tk::spline waypoint_spline;
            waypoint_spline.set_points(vehicle_x, vehicle_y);

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
            double transform_x_test = current_vehicle_x * cos(car_yaw_radians) - current_vehicle_y * sin(car_yaw_radians);
            double transform_y_test = current_vehicle_x * sin(car_yaw_radians) + current_vehicle_y * cos(car_yaw_radians);
            transform_x_test += car_x;
            transform_y_test += car_y;
            //cout << "%(" << transform_x_test << "," << transform_y_test << ")%\n";
            
            
            int steps_in_trajectory = 100;
            int remaining_points = steps_in_trajectory - previous_path_x.size();
            int existing_size = previous_path_x.size();
            double current_target_distance = car_speed * 0.02 * 0.44704; // divide by timestep and convert to meters/second

            double distance_forecast = 0;
            for (int i = 0; i < steps_in_trajectory; i++) {
              if (i > 0) { // don't add distance for first step in previous path, which already has unknown distance
                if (current_target_distance < meters_per_timestep) {
                  current_target_distance += (meters_per_timestep - current_target_distance) / 100;
                } else {
                  current_target_distance = meters_per_timestep;
                }
                distance_forecast += current_target_distance;
              }
              
              double x = current_vehicle_x + distance_forecast; // (meters_per_timestep * i);
              double y = waypoint_spline(x);

              // convert back to map coordinates
              double transformed_x = x * cos(car_yaw_radians) - y * sin(car_yaw_radians);
              double transformed_y = x * sin(car_yaw_radians) + y * cos(car_yaw_radians);
              transformed_x += car_x;
              transformed_y += car_y;
            
              // blend with existing path
              if (i < existing_size) {
                transformed_x = transformed_x * (i / (double)existing_size) + (double)previous_path_x[i] * (1 - (i / (double)existing_size));
                transformed_y = transformed_y * (i / (double)existing_size) + (double)previous_path_y[i] * (1 - (i / (double)existing_size));
              }
              next_x_vals.push_back(transformed_x);
              next_y_vals.push_back(transformed_y);
            }
            

            //for (int i = 0; i < next_x_vals.size(); ++i) {
            //  cout << "(" << next_x_vals[i] << "," << next_y_vals[i] << ") \n";
            //}
            //cout << endl;

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
