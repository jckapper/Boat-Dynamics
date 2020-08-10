#include "boat_dynamics/boat_dynamics.h"
using namespace Eigen;

namespace boat_dynamics {

    BoatDynamics::BoatDynamics() : t_prev_(0.0), t_initialized_(false), nh_private_("~")
    {
        nh_ = ros::NodeHandle();

        // Starting here: initialize persistent class members for calculations 
        boat_height_m_ = 2.0;
        boat_mass_kg = 227600;
        boat_length_m = 33.2232;
        boat_width_m = 8.504;
        Vector3d grav = Vector3d(0, 0, -9.81);
        Matrix3d inertia_matrix;//check syntax
        inertia_matrix << 2 * pow(boat_width_m, 2), 0, 0, //comma-initialization
            0, pow(boat_width_m, 2) + pow(boat_length_m, 2), 0,
            0, 0, pow(boat_width_m, 2) + pow(boat_length_m, 2);
        Matrix3d boat_inertia = (1.0 / 5) * boat_mass_kg * inertia_matrix;
        Matrix3d boat_inertia_inv = boat_inertia.inverse();

        //define wrench here
        modeling::Wrech Initial_wrench;
        Initial_wrench.F = grav*boat_mass_kg;
        Initial_wrench.T = Vector3d(0, 0, 250);

        //    boat_speed_mps_ = 0.5;
        // boat_speed_mps_ = nh_private_.param<double>("boat_speed", 0.0); // Is this speed needed

        T_0_boat_ = Xformd((Vector3d() << 0.0, 0.0, boat_height_m_).finished(), Quatd::Identity());
        T_0_boatNED_ = Xformd((Vector3d() << 0.0, 0.0, boat_height_m_).finished(), Quatd::from_euler(M_PI, 0.0, 0.0));
        T_NED_0_ = Xformd((Vector3d() << 0.0, 0.0, 0.0).finished(), Quatd::from_euler(M_PI, 0.0, 0.0)).inverse();
        modeling::State6DOF Current_State;
        Current_State.X = (T_0_boat_);
        Current_State.v = Vector3d(2.0, 0.0, 0.0);
        Current_State.w = Vector3d(0.0, 0.0, 0.0);

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
    modeling::ErrorState6DOF BoatDynamics::RKintegrate(modeling::State6DOF& State, const modeling::Wrench& u, const double& mass,
        const Eigen::Matrix3d& inertia, const Eigen::Matrix3d& inertia_inv, const double& dt)
    {
        modeling::ErrorState6DOF k1 = dynamics(State, u, mass, inertia, inertia_inv);

        modeling::State6DOF      x2 = State + k1 * (dt / 2.0);
        modeling::ErrorState6DOF k2 = dynamics(x2, u, mass, inertia, inertia_inv);

        modeling::State6DOF      x3 = State + k2 * (dt / 2.0);
        modeling::ErrorState6DOF k3 = dynamics(x3, u, mass, inertia, inertia_inv);

        modeling::State6DOF      x4 = State + k3 * dt;
        modeling::ErrorState6DOF k4 = dynamics(x4, u, mass, inertia, inertia_inv);

        modeling::ErrorState6DOF dx = (k1 + k2 * 2.0 + k3 * 2.0 + k4) * (dt / 6.0);

        return dx;

    }

    modeling::ErrorState6DOF BoatDynamics::dynamics(const modeling::State6DOF& x, const modeling::Wrench& u, 
        const double& m, const Eigen::Matrix3d& J, const Eigen::Matrix3d& J_inv)
    {
        modeling::ErrorState6DOF dx;
        dx.p = x.q.rota(x.v);
        dx.v = u.F / m + x.q.rotp(grav) - x.w.cross(x.v);
        dx.q = x.w;
        dx.w = J_inv * (u.T - x.w.cross(J * x.w));

        return dx;
    }

    void BoatDynamics::onUpdate(const ros::TimerEvent&)
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
        
        // Need to update boat state
        // T_0_boat_.t_(0) += boat_speed_mps_ * dt;
        Current_State += RKintegrate(Current_State, Initial_wrench, boat_mass_kg, boat_inertia, boat_inertia_inv, dt);

        // ++++ update orientation using from_two_unit_vectors ++++

        // update and send messages
        setMessageStates(Rt);
        tbr_.sendTransform(transform_);
        tbr_.sendTransform(transformNED_);
        truth_pub_.publish(truth_);
        marker_pub_.publish(marker_);
    }

    void BoatDynamics::setMessageStates(ros::Time& rt)
    {
        transform_.header.stamp = rt;
        transform_.transform.translation.x = Current_State.v(0)*dt;
        transform_.transform.translation.y = Current_State.v(1)*dt;
        transform_.transform.translation.z = Current_State.v(2)*dt;
        transform_.transform.rotation.w = Current_State.X.q_.w();
        transform_.transform.rotation.x = Current_State.X.q_.x();
        transform_.transform.rotation.y = Current_State.X.q_.y();
        transform_.transform.rotation.z = Current_State.X.q_.z();

        Xformd T_NED_boat = T_NED_0_ * Current_State.X;

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