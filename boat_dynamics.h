#include <ros/ros.h>
#include "geometry-utils-lib/xform.h"
#include <tf2_ros/transform_broadcaster.h>
#include <geometry_msgs/TransformStamped.h>
#include <geometry_msgs/PoseStamped.h>
#include <visualization_msgs/Marker.h>
#include <Eigen/Core>

using namespace transforms;
using namespace Eigen;

namespace modeling {
    typedef Matrix<double, 6, 1> Vector6d;
    struct ErrorState6DOF
    {
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW
            enum
        {
            SIZE = 12
        };

        Matrix<double, SIZE, 1> arr;
        Map<Vector6d> X;
        Map<Vector3d> p;
        Map<Vector3d> q;
        Map<Vector3d> v;
        Map<Vector3d> w;

        ErrorState6DOF() :
            X(arr.data()),
            p(arr.data()),
            q(arr.data() + 3),
            v(arr.data() + 6),
            w(arr.data() + 9)
        {
            arr.setZero();
        }

        ErrorState6DOF(const ErrorState6DOF& obj) :
            X(arr.data()),
            p(arr.data()),
            q(arr.data() + 3),
            v(arr.data() + 6),
            w(arr.data() + 9)
        {
            arr = obj.arr;
        }

        ErrorState6DOF& operator= (const ErrorState6DOF& obj)
        {
            arr = obj.arr;
            return *this;
        }

        ErrorState6DOF operator* (const double& s)
        {
            ErrorState6DOF out;
            out.arr = s * arr;
            return out;
        }

        ErrorState6DOF operator+ (const ErrorState6DOF& obj)
        {
            ErrorState6DOF out;
            out.arr = obj.arr + arr;
            return out;
        }
    };

    inline std::ostream& operator<< (std::ostream& os, const ErrorState6DOF& e)
    {
        os << "dt: " << e.p.transpose() << " dq: " << e.q.transpose() << " dv: " << e.v.transpose() << " dw: " << e.w.transpose();
        return os;
    }

    struct State6DOF
    {
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW
            enum
        {
            SIZE = 13
        };

        Matrix<double, SIZE, 1> arr;
        Xformd X;
        Map<Vector3d> p;
        Quatd q;
        Map<Vector3d> v;
        Map<Vector3d> w;

        State6DOF() :
            X(arr.data()),
            p(arr.data()),
            q(arr.data() + 3),
            v(arr.data() + 7),
            w(arr.data() + 10)
        {
            arr.setZero();
            q = Quatd::Identity();
        }

        State6DOF(const State6DOF& x) :
            X(arr.data()),
            p(arr.data()),
            q(arr.data() + 3),
            v(arr.data() + 7),
            w(arr.data() + 10)
        {
            arr = x.arr;
        }

        State6DOF& operator= (const State6DOF& obj)
        {

            arr = obj.arr;
            return *this;
        }

        State6DOF operator+(const ErrorState6DOF& dx) const
        {
            State6DOF xp;
            xp.p = p + dx.p;
            xp.q = q + dx.q;
            xp.v = v + dx.v;
            xp.w = w + dx.w;
            return xp;
        }

        State6DOF& operator+=(const ErrorState6DOF& dx)
        {
            p = p + dx.p;
            q = q + dx.q;
            v = v + dx.v;
            w = w + dx.w;
            return *this;
        }

        ErrorState6DOF operator-(const State6DOF& x) const
        {
            ErrorState6DOF dx;
            dx.p = p - x.p;
            dx.q = q - x.q;
            dx.v = v - x.v;
            dx.w = w - x.w;
            return dx;
        }
    };
    struct Wrench
    {
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW
            enum
        {
            SIZE = 6
        };

        Matrix<double, SIZE, 1> arr;
        Map<Vector3d> F;
        Map<Vector3d> T;

        Wrench() :
            F(arr.data()),
            T(arr.data() + 3)
        {
            arr.setZero();
        }

        Wrench(const Wrench& x) :
            F(arr.data()),
            T(arr.data() + 3)
        {
            arr = x.arr;
        }

        Wrench& operator= (const Wrench& obj)
        {
            arr = obj.arr;
            return *this;
        }

        Wrench operator+(const Wrench& dx) const
        {
            Wrench xp;
            xp.F = F + dx.F;
            xp.T = T + dx.T;
            return xp;
        }

        Wrench& operator+=(const Wrench& dx)
        {
            F = F + dx.F;
            T = T + dx.T;
            return *this;
        }

        Wrench operator-(const Wrench& x) const
        {
            Wrench dx;
            dx.F = F - x.F;
            dx.T = T - x.T;
            return dx;
        }
    };

    inline std::ostream& operator<< (std::ostream& os, const Wrench& s)
    {
        os << "F: " << s.F.transpose() << " T: " << s.T.transpose();
        return os;
    }
}

namespace boat_dynamics {

class BoatDynamics
{
public:
    BoatDynamics();

private:
    void onUpdate(const ros::TimerEvent &event);
    void setMessageStates(ros::Time &rt);
    void RKintegrate(modeling::State6DOF& State, const modeling::Wrench& u, const double& mass,const Eigen::Matrix3d& inertia, const Eigen::Matrix3d& inertia_inv, const double& dt);
    modeling::ErrorState6DOF dynamics(const modeling::State6DOF& x, const modeling::Wrench& u,const double& m, const Eigen::Matrix3d& J, const Eigen::Matrix3d& J_inv);

    ros::NodeHandle nh_;
    ros::NodeHandle nh_private_;
    ros::Timer timer_;

    double boat_height_m_;
    int boat_mass_kg;
    double boat_length_m;
    double boat_width_m;
    Vector3d grav;
    Matrix3d inertia_matrix;
    Matrix3d boat_inertia;
    Matrix3d boat_inertia_inv;
    double boat_speed_mps_;
    double t_prev_;
    bool t_initialized_;
    modeling::State6DOF Current_State;
    modeling::Wrench u;


    Xformd T_0_boat_;
    Xformd T_0_boatNED_;
    Xformd T_NED_0_;

    tf2_ros::TransformBroadcaster tbr_;
    ros::Publisher truth_pub_;
    ros::Publisher marker_pub_;

    geometry_msgs::TransformStamped transform_;
    geometry_msgs::TransformStamped transformNED_;
    geometry_msgs::PoseStamped truth_;
    visualization_msgs::Marker marker_;
};

} // end namespace boat_dynamics