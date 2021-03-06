#include "boat_dynamics/boat_dynamics.h"
using namespace Eigen;
using namespace std;

namespace boat_dynamics {

    BoatDynamics::BoatDynamics() : t_prev_(0.0), t_initialized_(false), nh_private_("~")
        ///
        /// Initializes all necessary variables and calculates persistent physical characteristics for the ship.
        ///
        /// Along with initializing all of the ship's values/characteristics, this function also sets markers for the simulation.  
        /// User input is provided for entering an intial velocity vector before the simulation starts, options will be added 
        /// for inputting force, torque, and orientations without first having to hardcode them in.  All reference frames
        /// are also initialized here and are used throughout the calculation and visualizations later on.
        /// 
        /// Parameters:
        /// None
        ///
        /// Returns:
        /// None
        ///
    {
        nh_ = ros::NodeHandle();

        // Initialize persistent physical class members for calculations 
        boat_height_m_ = 2.0;
        boat_mass_kg_ = 227600;
        boat_length_m_ = 33.2232;
        boat_width_m_ = 8.504;
        grav_ = Vector3d(0., 0., -9.81);
        
        // Approximating ship as an ellipsoid
        inertia_matrix_ << (2 * boat_width_m_ * boat_width_m_), 0., 0., //comma-initialization
                        0., ((boat_width_m_ * boat_width_m_) + (boat_length_m_ * boat_length_m_)), 0.,
                        0., 0., ((boat_width_m_ * boat_width_m_) + (boat_length_m_ * boat_length_m_));
        boat_inertia_ = (1.0 / 5.0) * boat_mass_kg_ * inertia_matrix_;
        boat_inertia_inv_ = boat_inertia_.inverse();
        
        // Initialize Xformd's and orientation quaternions
        T_0_boat_ = Xformd((Vector3d() << 0.0, 0.0, boat_height_m_).finished(), Quatd::Identity());
        T_0_boatNED_ = Xformd((Vector3d() << 0.0, 0.0, boat_height_m_).finished(), Quatd::from_euler(M_PI, 0.0, 0.0));
        T_NED_0_ = Xformd((Vector3d() << 0.0, 0.0, 0.0).finished(), Quatd::from_euler(M_PI, 0.0, 0.0)).inverse();
        
        // User input for initial conditions, simplifying and adding more input options in the future
        cout << "--Please Specify Initial Conditions--\n";
        cout << "Enter an initial X velocity: ";
        cin >> X_Velocity;
        cout << "Enter an initial Y velocity: ";
        cin >> Y_Velocity;
        cout << "Enter an initial Z velocity: ";
        cin >> Z_Velocity;
        
        // Initialize Six Degree of Freedom State method for boat
        Current_State_.X = (T_0_boat_);
        Current_State_.v = Vector3d(X_Velocity, Y_Velocity, Z_Velocity);       
        Current_State_.w = Vector3d(0., 0., 0.);
        Current_State_.q = Quatd::Identity();
        
        // Initialize wrench
        u_.F = Vector3d(0.0, 0.0, 9.81*boat_mass_kg_);//-grav_ * boat_mass_kg_;
        u_.T = Vector3d(0., 0., 50000.);
        
        truth_pub_ = nh_.advertise<geometry_msgs::PoseStamped>("boat_truth_NED", 1);
        marker_pub_ = nh_.advertise<visualization_msgs::Marker>("boat_marker", 1);

        transform_.header.frame_id = "world";
        transform_.child_frame_id = "boat";
        transformNED_.header.frame_id = "world";
        transformNED_.child_frame_id = "boatNED";

        truth_.header.frame_id = "NED";

        marker_.header.frame_id = "NED";
        marker_.id = 1;
        marker_.type = visualization_msgs::Marker::MESH_RESOURCE;
        marker_.action = visualization_msgs::Marker::ADD;
        marker_.mesh_resource = "package://boat_dynamics/mesh/ship_nwu_GCS.dae";
        double boat_scale = 20.0;
        marker_.scale.x = boat_scale;
        marker_.scale.y = boat_scale;
        marker_.scale.z = boat_scale;
        marker_.color.r = 1.0;
        marker_.color.g = 1.0;
        marker_.color.b = 1.0;
        marker_.color.a = 1.0;

        ros::Time Rt = ros::Time::now();
        setMessageStates(Rt);

        timer_ = nh_.createTimer(ros::Duration(ros::Rate(20)), &BoatDynamics::onUpdate, this);
    }

    // Beginning adaptation
    modeling::State6DOF BoatDynamics::RKintegrate(modeling::State6DOF& State, const modeling::Wrench& u, const double& mass,
        const Eigen::Matrix3d& inertia, const Eigen::Matrix3d& inertia_inv, const double& dt)
        ///
        /// Implementation of 4th order Runge-Kutta integration that takes care of the bigger picture (calculating the change of state, phi in the explanation).
        ///
        /// This is done by combining the x's calculated by BoatDynamics::dynamics() with the State parameter and calculating the overall state change (phi in explanation, 
        /// dx variable) that is then returned.  The function makes four calls to the BoatDynamics::dynamics() equation, one at each step, where most of the 
        /// calculation is performed.  
        /// 
        /// Parameters:
        /// State (State6DOF struct): the current state of the object (position and orientation) in one convenient package
        /// u (const wrench): the force and torque that will be applied to the boat, currently a constant value
        /// mass (const double): the mass of the ship
        /// inertia (const Matrix3d): Eigen 3x3 double matrix for the inertia of the ship
        /// inertia_inv (const Matrix3d): Eigen 3x3 double matrix, inverse of inertia parameter
        /// dt (const double): the time change
        ///
        /// Returns:
        /// State dx (State6DOF struct): the updated state of the ship
        ///
    {
        // Runge-Kutta integration
        modeling::ErrorState6DOF k1 = dynamics(State, u, mass, inertia, inertia_inv);

        modeling::State6DOF      x2 = State + k1 * (dt / 2.0);
        modeling::ErrorState6DOF k2 = dynamics(x2, u, mass, inertia, inertia_inv);

        modeling::State6DOF      x3 = State + k2 * (dt / 2.0);
        modeling::ErrorState6DOF k3 = dynamics(x3, u, mass, inertia, inertia_inv);

        modeling::State6DOF      x4 = State + k3 * dt;
        modeling::ErrorState6DOF k4 = dynamics(x4, u, mass, inertia, inertia_inv);

        modeling::ErrorState6DOF dx = (k1 + k2 * 2.0 + k3 * 2.0 + k4) * (dt / 6.0);

        return State + dx;

    }

    modeling::ErrorState6DOF BoatDynamics::dynamics(const modeling::State6DOF& x, const modeling::Wrench& u, 
        const double& m, const Eigen::Matrix3d& J, const Eigen::Matrix3d& J_inv)
        ///
        /// Implementation of 4th order Runge-Kutta dynamics calculation to calculate dx (or f(x,u) depending on name used in explanation)
        ///
        /// The code calculates the Runge-Kutta dx matrix and stores it in an ErrorState6DOF struct so that it can be added to the initial state using
        /// the manifold addition that is implemented in the struct.
        /// 
        /// Parameters:
        /// x (State6DOF struct): the current state of the object (position and orientation) in one convenient package
        /// u (const wrench): the force and torque that will be applied to the ship, currently a constant value
        /// m (const double): the mass of the ship
        /// J (const Matrix3d): Eigen 3x3 double matrix for the inertia of the ship
        /// J_inv (const Matrix3d): Eigen 3x3 double matrix, inverse of inertia parameter
        ///
        /// Returns:
        /// dx (ErrorState6DOF struct): the state change for the ship
        ///
    {
        // Runge-Kutta Dynamics
        modeling::ErrorState6DOF dx;
        dx.p = x.q.rota(x.v);
        dx.v = u.F / m + x.q.rotp(grav_) - x.w.cross(x.v);
        dx.q = x.w;
        dx.w = J_inv * (u.T - x.w.cross(J * x.w));

        return dx;
    }

    void BoatDynamics::onUpdate(const ros::TimerEvent&)
        ///
        /// Calls the Runge-Kutta state update function and performs other necessary calculations for the simulation state
        ///
        /// Keeps track of the dt for each step, and at each step calculates the change in state of the ship before updating and sending the ROS messages
        /// 
        /// Parameters:
        /// Rt (const ros::TimerEvent): ROS timing for the event
        ///
        /// Returns:
        /// None, void function
        ///
    {
        ros::Time Rt = ros::Time::now();
        double t = Rt.toSec();

        if (!t_initialized_)
        {
            t_prev_ = t;
            t_initialized_ = true;
            return;
        }

        double dt = t - t_prev_;
        t_prev_ = t;

        Current_State_ = RKintegrate(Current_State_, u_, boat_mass_kg_, boat_inertia_, boat_inertia_inv_, dt);

        // Update and send messages
        setMessageStates(Rt);
        tbr_.sendTransform(transform_);
        tbr_.sendTransform(transformNED_);
        truth_pub_.publish(truth_);
        marker_pub_.publish(marker_);
    }

    void BoatDynamics::setMessageStates(ros::Time& rt)
        ///
        /// Sets ROS message states based of the updated ship state for each dt to be passed on to the simulation.
        ///
        /// The message state for each coordinate frame is calculated and set from the ship's updated state at each time step.
        /// All values are pulled from Current_State_ which is a State6DOF struct and so contains all necessary position and orientation info.
        /// 
        /// Parameters:
        /// rt (ros::Time): time from the ROS simulation
        ///
        /// Returns:
        /// None, void function
        ///
    {
        // Update simulation state for visuals
        transform_.header.stamp = rt;
        transform_.transform.translation.x = Current_State_.p.x();
        transform_.transform.translation.y = Current_State_.p.y();
        transform_.transform.translation.z = Current_State_.p.z();
        transform_.transform.rotation.w = Current_State_.q.w();
        transform_.transform.rotation.x = Current_State_.q.x();
        transform_.transform.rotation.y = Current_State_.q.y();
        transform_.transform.rotation.z = Current_State_.q.z();

        Xformd T_NED_boat = T_NED_0_ * Current_State_.X;

        truth_.header.stamp = rt;
        truth_.pose.position.x = T_NED_boat.t_(0);
        truth_.pose.position.y = T_NED_boat.t_(1);
        truth_.pose.position.z = T_NED_boat.t_(2);
        truth_.pose.orientation.w = T_NED_boat.q_.w();
        truth_.pose.orientation.x = T_NED_boat.q_.x();
        truth_.pose.orientation.y = T_NED_boat.q_.y();
        truth_.pose.orientation.z = T_NED_boat.q_.z();

        marker_.header.stamp = rt;
        marker_.pose = truth_.pose;

        transformNED_.header.stamp = rt;
        transformNED_.transform.translation.x = T_0_boatNED_.t_(0);
        transformNED_.transform.translation.y = T_0_boatNED_.t_(1);
        transformNED_.transform.translation.z = T_0_boatNED_.t_(2);
        transformNED_.transform.rotation.w = T_0_boatNED_.q_.w();
        transformNED_.transform.rotation.x = T_0_boatNED_.q_.x();
        transformNED_.transform.rotation.y = T_0_boatNED_.q_.y();
        transformNED_.transform.rotation.z = T_0_boatNED_.q_.z();
    }

} // end namespace boat_dynamics
