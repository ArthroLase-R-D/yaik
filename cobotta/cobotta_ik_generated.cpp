#include "yaik_cpp_common.h"

using namespace yaik_cpp;

struct cobotta_ik {

// Constants for solver
static constexpr int robot_nq = 6;
static constexpr int max_n_solutions = 128;
static constexpr int n_tree_nodes = 36;
static constexpr int intermediate_solution_size = 7;
static constexpr double pose_tolerance = 1e-6;
static constexpr double pose_tolerance_degenerate = 1e-4;
static constexpr double zero_tolerance = 1e-6;
using RawIKWorksace = ::yaik_cpp::RawIkWorkspaceTemplate<intermediate_solution_size, max_n_solutions, robot_nq>;

// Robot parameters
static constexpr double a_0 = 710.0;
static constexpr double d_1 = -50.0;
static constexpr double d_2 = 590.0;
static constexpr double d_3 = 150.0;
static constexpr double d_4 = 85.0;
static constexpr double post_transform_special_symbol_23 = 75.0;
static constexpr double pre_transform_special_symbol_13 = 220.0;

// Unknown offsets from original unknown value to raw value
// Original value are the ones corresponded to robot (usually urdf/sdf)
// Raw value are the ones used in the solver
// unknown_i_raw = unknown_i_original + unknown_i_offset_original2raw
static constexpr double th_0_offset_original2raw = 0.0;
static constexpr double th_1_offset_original2raw = -1.5707963267948966;
static constexpr double th_2_offset_original2raw = 0.0;
static constexpr double th_3_offset_original2raw = 3.141592653589793;
static constexpr double th_4_offset_original2raw = 3.141592653589793;
static constexpr double th_5_offset_original2raw = -1.5707963267948966;

// The transformation between raw and original ee target
// Original value are the ones corresponded to robot (usually urdf/sdf)
// Raw value are the ones used in the solver
// ee_original = pre_transform * ee_raw * post_transform
// ee_raw = dh_forward_transform(theta_raw)
static Eigen::Matrix4d endEffectorTargetOriginalToRaw(const Eigen::Matrix4d& T_ee)
{
    const double r_11 = T_ee(0, 0);
    const double r_12 = T_ee(0, 1);
    const double r_13 = T_ee(0, 2);
    const double Px = T_ee(0, 3);
    const double r_21 = T_ee(1, 0);
    const double r_22 = T_ee(1, 1);
    const double r_23 = T_ee(1, 2);
    const double Py = T_ee(1, 3);
    const double r_31 = T_ee(2, 0);
    const double r_32 = T_ee(2, 1);
    const double r_33 = T_ee(2, 2);
    const double Pz = T_ee(2, 3);
    Eigen::Matrix4d ee_transformed;
    ee_transformed.setIdentity();
    ee_transformed(0, 0) = 1.0*r_31;
    ee_transformed(0, 1) = 1.0*r_32;
    ee_transformed(0, 2) = 1.0*r_33;
    ee_transformed(0, 3) = 1.0*Pz - 1.0*post_transform_special_symbol_23*r_33;
    ee_transformed(1, 0) = 1.0*r_11;
    ee_transformed(1, 1) = 1.0*r_12;
    ee_transformed(1, 2) = 1.0*r_13;
    ee_transformed(1, 3) = 1.0*Px - 1.0*post_transform_special_symbol_23*r_13;
    ee_transformed(2, 0) = 1.0*r_21;
    ee_transformed(2, 1) = 1.0*r_22;
    ee_transformed(2, 2) = 1.0*r_23;
    ee_transformed(2, 3) = 1.0*Py - 1.0*post_transform_special_symbol_23*r_23 - 1.0*pre_transform_special_symbol_13;
    return ee_transformed;
}

static Eigen::Matrix4d endEffectorTargetRawToOriginal(const Eigen::Matrix4d& T_ee)
{
    const double r_11 = T_ee(0, 0);
    const double r_12 = T_ee(0, 1);
    const double r_13 = T_ee(0, 2);
    const double Px = T_ee(0, 3);
    const double r_21 = T_ee(1, 0);
    const double r_22 = T_ee(1, 1);
    const double r_23 = T_ee(1, 2);
    const double Py = T_ee(1, 3);
    const double r_31 = T_ee(2, 0);
    const double r_32 = T_ee(2, 1);
    const double r_33 = T_ee(2, 2);
    const double Pz = T_ee(2, 3);
    Eigen::Matrix4d ee_transformed;
    ee_transformed.setIdentity();
    ee_transformed(0, 0) = 1.0*r_21;
    ee_transformed(0, 1) = 1.0*r_22;
    ee_transformed(0, 2) = 1.0*r_23;
    ee_transformed(0, 3) = 1.0*Py + 1.0*post_transform_special_symbol_23*r_23;
    ee_transformed(1, 0) = 1.0*r_31;
    ee_transformed(1, 1) = 1.0*r_32;
    ee_transformed(1, 2) = 1.0*r_33;
    ee_transformed(1, 3) = 1.0*Pz + 1.0*post_transform_special_symbol_23*r_33 + 1.0*pre_transform_special_symbol_13;
    ee_transformed(2, 0) = 1.0*r_11;
    ee_transformed(2, 1) = 1.0*r_12;
    ee_transformed(2, 2) = 1.0*r_13;
    ee_transformed(2, 3) = 1.0*Px + 1.0*post_transform_special_symbol_23*r_13;
    return ee_transformed;
}

///************* Below are the actual FK and IK implementations *******************
static Eigen::Matrix4d computeFK(const std::array<double, robot_nq>& theta_input_original)
{
    // Extract the variables
    const double th_0 = theta_input_original[0] + th_0_offset_original2raw;
    const double th_1 = theta_input_original[1] + th_1_offset_original2raw;
    const double th_2 = theta_input_original[2] + th_2_offset_original2raw;
    const double th_3 = theta_input_original[3] + th_3_offset_original2raw;
    const double th_4 = theta_input_original[4] + th_4_offset_original2raw;
    const double th_5 = theta_input_original[5] + th_5_offset_original2raw;
    
    // Temp variable for efficiency
    const double x0 = std::sin(th_5);
    const double x1 = std::sin(th_0);
    const double x2 = std::cos(th_3);
    const double x3 = std::sin(th_3);
    const double x4 = std::cos(th_0);
    const double x5 = std::sin(th_1);
    const double x6 = std::sin(th_2);
    const double x7 = x5*x6;
    const double x8 = std::cos(th_1);
    const double x9 = std::cos(th_2);
    const double x10 = -x4*x7 + x4*x8*x9;
    const double x11 = x1*x2 - x10*x3;
    const double x12 = std::cos(th_5);
    const double x13 = std::sin(th_4);
    const double x14 = x5*x9;
    const double x15 = x6*x8;
    const double x16 = -x14*x4 - x15*x4;
    const double x17 = std::cos(th_4);
    const double x18 = x1*x3 + x10*x2;
    const double x19 = -x13*x16 + x17*x18;
    const double x20 = -x13*x18 - x16*x17;
    const double x21 = a_0*x8;
    const double x22 = -x1*x7 + x1*x8*x9;
    const double x23 = -x2*x4 - x22*x3;
    const double x24 = -x1*x14 - x1*x15;
    const double x25 = x2*x22 - x3*x4;
    const double x26 = -x13*x24 + x17*x25;
    const double x27 = -x13*x25 - x17*x24;
    const double x28 = -x14 - x15;
    const double x29 = x28*x3;
    const double x30 = x7 - x8*x9;
    const double x31 = -x13*x30 + x17*x2*x28;
    const double x32 = -x13*x2*x28 - x17*x30;
    // End of temp variables
    Eigen::Matrix4d ee_pose_raw;
    ee_pose_raw.setIdentity();
    ee_pose_raw(0, 0) = -x0*x11 + x12*x19;
    ee_pose_raw(0, 1) = -x0*x19 - x11*x12;
    ee_pose_raw(0, 2) = x20;
    ee_pose_raw(0, 3) = -d_1*x1 + d_2*x16 + d_3*x11 + d_4*x20 + x21*x4;
    ee_pose_raw(1, 0) = -x0*x23 + x12*x26;
    ee_pose_raw(1, 1) = -x0*x26 - x12*x23;
    ee_pose_raw(1, 2) = x27;
    ee_pose_raw(1, 3) = d_1*x4 + d_2*x24 + d_3*x23 + d_4*x27 + x1*x21;
    ee_pose_raw(2, 0) = x0*x29 + x12*x31;
    ee_pose_raw(2, 1) = -x0*x31 + x12*x28*x3;
    ee_pose_raw(2, 2) = x32;
    ee_pose_raw(2, 3) = -a_0*x5 + d_2*x30 - d_3*x29 + d_4*x32;
    return endEffectorTargetRawToOriginal(ee_pose_raw);
}

static void computeTwistJacobian(const std::array<double, robot_nq>& theta_input_original, Eigen::Matrix<double, 6, robot_nq>& jacobian)
{
    // Extract the variables
    const double th_0 = theta_input_original[0] + th_0_offset_original2raw;
    const double th_1 = theta_input_original[1] + th_1_offset_original2raw;
    const double th_2 = theta_input_original[2] + th_2_offset_original2raw;
    const double th_3 = theta_input_original[3] + th_3_offset_original2raw;
    const double th_4 = theta_input_original[4] + th_4_offset_original2raw;
    const double th_5 = theta_input_original[5] + th_5_offset_original2raw;
    
    // Temp variable for efficiency
    const double x0 = std::cos(th_0);
    const double x1 = 1.0*x0;
    const double x2 = std::cos(th_2);
    const double x3 = std::sin(th_1);
    const double x4 = std::sin(th_0);
    const double x5 = 1.0*x4;
    const double x6 = x3*x5;
    const double x7 = std::sin(th_2);
    const double x8 = std::cos(th_1);
    const double x9 = x5*x8;
    const double x10 = -x2*x6 - x7*x9;
    const double x11 = std::cos(th_3);
    const double x12 = std::sin(th_3);
    const double x13 = 1.0*x2*x4*x8 - x6*x7;
    const double x14 = -x1*x11 - x12*x13;
    const double x15 = std::cos(th_4);
    const double x16 = std::sin(th_4);
    const double x17 = -x10*x15 - x16*(-x1*x12 + x11*x13);
    const double x18 = 1.0*x7;
    const double x19 = 1.0*x2;
    const double x20 = x18*x3 - x19*x8;
    const double x21 = -x18*x8 - x19*x3;
    const double x22 = x12*x21;
    const double x23 = -x11*x16*x21 - x15*x20;
    const double x24 = -x5;
    const double x25 = x1*x3;
    const double x26 = x1*x8;
    const double x27 = -x2*x25 - x26*x7;
    const double x28 = 1.0*x0*x2*x8 - x25*x7;
    const double x29 = 1.0*x11*x4 - x12*x28;
    const double x30 = -x15*x27 - x16*(x11*x28 + x12*x5);
    const double x31 = 1.0*a_0*x3;
    const double x32 = -pre_transform_special_symbol_13 + x31;
    const double x33 = -x32;
    const double x34 = d_2*x20 + pre_transform_special_symbol_13 - x31;
    const double x35 = a_0*x26 - d_1*x5;
    const double x36 = d_2*x27 + x35;
    const double x37 = d_3*x22;
    const double x38 = d_2*x20 - x32 - x37;
    const double x39 = d_3*x29 + x36;
    const double x40 = d_4*x23 + x34 - x37;
    const double x41 = d_4*x30 + x39;
    const double x42 = a_0*x9 + d_1*x1;
    const double x43 = d_2*x10 + x42;
    const double x44 = d_3*x14 + x43;
    const double x45 = d_4*x17 + x44;
    // End of temp variables
    
    jacobian.setZero();
    jacobian(0, 1) = x1;
    jacobian(0, 2) = x1;
    jacobian(0, 3) = x10;
    jacobian(0, 4) = x14;
    jacobian(0, 5) = x17;
    jacobian(1, 0) = 1.0;
    jacobian(1, 3) = x20;
    jacobian(1, 4) = -x22;
    jacobian(1, 5) = x23;
    jacobian(2, 1) = x24;
    jacobian(2, 2) = x24;
    jacobian(2, 3) = x27;
    jacobian(2, 4) = x29;
    jacobian(2, 5) = x30;
    jacobian(3, 1) = -pre_transform_special_symbol_13*x5;
    jacobian(3, 2) = -x33*x5;
    jacobian(3, 3) = -x20*x36 + x27*x34;
    jacobian(3, 4) = x22*x39 + x29*x38;
    jacobian(3, 5) = -x23*x41 + x30*x40;
    jacobian(4, 2) = x1*x35 + x42*x5;
    jacobian(4, 3) = x10*x36 - x27*x43;
    jacobian(4, 4) = x14*x39 - x29*x44;
    jacobian(4, 5) = x17*x41 - x30*x45;
    jacobian(5, 1) = -pre_transform_special_symbol_13*x1;
    jacobian(5, 2) = -x1*x33;
    jacobian(5, 3) = -x10*x34 + x20*x43;
    jacobian(5, 4) = -x14*x38 - x22*x44;
    jacobian(5, 5) = -x17*x40 + x23*x45;
    return;
}

static void computeAngularVelocityJacobian(const std::array<double, robot_nq>& theta_input_original, Eigen::Matrix<double, 3, robot_nq>& jacobian)
{
    // Extract the variables
    const double th_0 = theta_input_original[0] + th_0_offset_original2raw;
    const double th_1 = theta_input_original[1] + th_1_offset_original2raw;
    const double th_2 = theta_input_original[2] + th_2_offset_original2raw;
    const double th_3 = theta_input_original[3] + th_3_offset_original2raw;
    const double th_4 = theta_input_original[4] + th_4_offset_original2raw;
    const double th_5 = theta_input_original[5] + th_5_offset_original2raw;
    
    // Temp variable for efficiency
    const double x0 = std::cos(th_0);
    const double x1 = 1.0*x0;
    const double x2 = std::cos(th_2);
    const double x3 = std::sin(th_1);
    const double x4 = std::sin(th_0);
    const double x5 = 1.0*x4;
    const double x6 = x3*x5;
    const double x7 = std::sin(th_2);
    const double x8 = std::cos(th_1);
    const double x9 = -x2*x6 - x5*x7*x8;
    const double x10 = std::cos(th_3);
    const double x11 = std::sin(th_3);
    const double x12 = 1.0*x2*x4*x8 - x6*x7;
    const double x13 = std::cos(th_4);
    const double x14 = std::sin(th_4);
    const double x15 = 1.0*x7;
    const double x16 = 1.0*x2;
    const double x17 = x15*x3 - x16*x8;
    const double x18 = -x15*x8 - x16*x3;
    const double x19 = -x5;
    const double x20 = x1*x3;
    const double x21 = -x1*x7*x8 - x2*x20;
    const double x22 = 1.0*x0*x2*x8 - x20*x7;
    // End of temp variables
    
    jacobian.setZero();
    jacobian(0, 1) = x1;
    jacobian(0, 2) = x1;
    jacobian(0, 3) = x9;
    jacobian(0, 4) = -x1*x10 - x11*x12;
    jacobian(0, 5) = -x13*x9 - x14*(-x1*x11 + x10*x12);
    jacobian(1, 0) = 1.0;
    jacobian(1, 3) = x17;
    jacobian(1, 4) = -x11*x18;
    jacobian(1, 5) = -x10*x14*x18 - x13*x17;
    jacobian(2, 1) = x19;
    jacobian(2, 2) = x19;
    jacobian(2, 3) = x21;
    jacobian(2, 4) = 1.0*x10*x4 - x11*x22;
    jacobian(2, 5) = -x13*x21 - x14*(x10*x22 + x11*x5);
    return;
}

static void computeTransformPointJacobian(const std::array<double, robot_nq>& theta_input_original, Eigen::Vector3d& point_on_ee, Eigen::Matrix<double, 3, robot_nq>& jacobian)
{
    // Extract the variables
    const double th_0 = theta_input_original[0] + th_0_offset_original2raw;
    const double th_1 = theta_input_original[1] + th_1_offset_original2raw;
    const double th_2 = theta_input_original[2] + th_2_offset_original2raw;
    const double th_3 = theta_input_original[3] + th_3_offset_original2raw;
    const double th_4 = theta_input_original[4] + th_4_offset_original2raw;
    const double th_5 = theta_input_original[5] + th_5_offset_original2raw;
    const double p_on_ee_x = point_on_ee[0];
    const double p_on_ee_y = point_on_ee[1];
    const double p_on_ee_z = point_on_ee[2];
    
    // Temp variable for efficiency
    const double x0 = 1.0*p_on_ee_z;
    const double x1 = std::sin(th_0);
    const double x2 = 1.0*x1;
    const double x3 = p_on_ee_y*x2;
    const double x4 = std::sin(th_1);
    const double x5 = 1.0*x4;
    const double x6 = a_0*x5;
    const double x7 = -pre_transform_special_symbol_13 + x6;
    const double x8 = -x7;
    const double x9 = std::sin(th_2);
    const double x10 = x5*x9;
    const double x11 = std::cos(th_2);
    const double x12 = std::cos(th_1);
    const double x13 = 1.0*x12;
    const double x14 = x10 - x11*x13;
    const double x15 = std::cos(th_0);
    const double x16 = x11*x5;
    const double x17 = x13*x9;
    const double x18 = -x15*x16 - x15*x17;
    const double x19 = d_2*x14 + pre_transform_special_symbol_13 - x6;
    const double x20 = a_0*x13*x15 - d_1*x2;
    const double x21 = d_2*x18 + x20;
    const double x22 = std::sin(th_3);
    const double x23 = -x16 - x17;
    const double x24 = x22*x23;
    const double x25 = std::cos(th_3);
    const double x26 = -x10*x15 + 1.0*x11*x12*x15;
    const double x27 = 1.0*x1*x25 - x22*x26;
    const double x28 = d_3*x24;
    const double x29 = d_2*x14 - x28 - x7;
    const double x30 = d_3*x27 + x21;
    const double x31 = std::cos(th_4);
    const double x32 = std::sin(th_4);
    const double x33 = -x14*x31 - x23*x25*x32;
    const double x34 = -x18*x31 - x32*(x2*x22 + x25*x26);
    const double x35 = d_4*x33 + x19 - x28;
    const double x36 = d_4*x34 + x30;
    const double x37 = 1.0*p_on_ee_x;
    const double x38 = x0*x15 + x1*x37;
    const double x39 = 1.0*x15;
    const double x40 = x12*x2;
    const double x41 = a_0*x40 + d_1*x39;
    const double x42 = x2*x4;
    const double x43 = -x11*x42 - x40*x9;
    const double x44 = d_2*x43 + x41;
    const double x45 = 1.0*x1*x11*x12 - x42*x9;
    const double x46 = -x22*x45 - x25*x39;
    const double x47 = d_3*x46 + x44;
    const double x48 = -x31*x43 - x32*(-x22*x39 + x25*x45);
    const double x49 = d_4*x48 + x47;
    const double x50 = p_on_ee_y*x39;
    // End of temp variables
    
    jacobian.setZero();
    jacobian(0, 0) = x0;
    jacobian(0, 1) = -pre_transform_special_symbol_13*x2 + x3;
    jacobian(0, 2) = -x2*x8 + x3;
    jacobian(0, 3) = -p_on_ee_y*x18 + p_on_ee_z*x14 - x14*x21 + x18*x19;
    jacobian(0, 4) = -p_on_ee_y*x27 - p_on_ee_z*x24 + x22*x23*x30 + x27*x29;
    jacobian(0, 5) = -p_on_ee_y*x34 + p_on_ee_z*x33 - x33*x36 + x34*x35;
    jacobian(1, 1) = -x38;
    jacobian(1, 2) = 1.0*x1*x41 + 1.0*x15*x20 - x38;
    jacobian(1, 3) = p_on_ee_x*x18 - p_on_ee_z*x43 - x18*x44 + x21*x43;
    jacobian(1, 4) = p_on_ee_x*x27 - p_on_ee_z*x46 - x27*x47 + x30*x46;
    jacobian(1, 5) = p_on_ee_x*x34 - p_on_ee_z*x48 - x34*x49 + x36*x48;
    jacobian(2, 0) = -x37;
    jacobian(2, 1) = -pre_transform_special_symbol_13*x39 + x50;
    jacobian(2, 2) = -x39*x8 + x50;
    jacobian(2, 3) = -p_on_ee_x*x14 + p_on_ee_y*x43 + x14*x44 - x19*x43;
    jacobian(2, 4) = p_on_ee_x*x24 + p_on_ee_y*x46 - x24*x47 - x29*x46;
    jacobian(2, 5) = -p_on_ee_x*x33 + p_on_ee_y*x48 + x33*x49 - x35*x48;
    return;
}

static void computeRawIK(const Eigen::Matrix4d& T_ee, SolutionQueue<intermediate_solution_size, max_n_solutions>& solution_queue, NodeIndexWorkspace<max_n_solutions>& node_index_workspace, std::vector<std::array<double, robot_nq>>& ik_output)
{
    // Extracting the ik target symbols
    const double r_11 = T_ee(0, 0);
    const double r_12 = T_ee(0, 1);
    const double r_13 = T_ee(0, 2);
    const double Px = T_ee(0, 3);
    const double r_21 = T_ee(1, 0);
    const double r_22 = T_ee(1, 1);
    const double r_23 = T_ee(1, 2);
    const double Py = T_ee(1, 3);
    const double r_31 = T_ee(2, 0);
    const double r_32 = T_ee(2, 1);
    const double r_33 = T_ee(2, 2);
    const double Pz = T_ee(2, 3);
    const auto& ee_rotation = T_ee.block<3, 3>(0, 0);
    const auto& ee_translation = T_ee.block<3, 1>(0, 3);
    const Eigen::Vector3d inv_ee_translation = - ee_rotation.transpose() * ee_translation;
    const double inv_Px = inv_ee_translation(0);
    const double inv_Py = inv_ee_translation(1);
    const double inv_Pz = inv_ee_translation(2);
    
    solution_queue.reset();
    node_index_workspace.reset(n_tree_nodes);
    
    using RawSolution = IntermediateSolution<intermediate_solution_size>;
    auto make_raw_solution = []() -> RawSolution { return {}; };
    
    auto append_solution_to_queue = [&solution_queue](RawSolution solution_2_add) -> int {
        return solution_queue.push_solution(solution_2_add);
    };
    
    auto add_input_index_to = [&node_index_workspace](int node_idx, int solution_idx) -> void {
        if (solution_idx < 0) return;
        node_index_workspace.append_index_to_node(node_idx, solution_idx);
    };
    
    // Code for non-branch dispatcher node 0
    // Actually, there is no code
    
    // Code for general_6dof solution node 1, solved variable is th_2
    auto General6DoFNumericalReduceSolutionNode_node_1_solve_th_2_processor = [&]() -> void
    {
        const auto& this_node_input_index = node_index_workspace.get_input_indices_for_node(0);
        const bool this_input_valid = node_index_workspace.is_input_indices_valid_for_node(0);
        if (!this_input_valid)
            return;
        
        // The general 6-dof solution of root node
        Eigen::Matrix<double, 8, 8> R_l;
        R_l.setZero();
        R_l(0, 0) = d_3*r_21;
        R_l(0, 1) = d_3*r_22;
        R_l(0, 2) = d_3*r_11;
        R_l(0, 3) = d_3*r_12;
        R_l(0, 4) = Py - d_4*r_23;
        R_l(0, 5) = Px - d_4*r_13;
        R_l(1, 0) = d_3*r_11;
        R_l(1, 1) = d_3*r_12;
        R_l(1, 2) = -d_3*r_21;
        R_l(1, 3) = -d_3*r_22;
        R_l(1, 4) = Px - d_4*r_13;
        R_l(1, 5) = -Py + d_4*r_23;
        R_l(2, 6) = d_3*r_31;
        R_l(2, 7) = d_3*r_32;
        R_l(3, 0) = -r_21;
        R_l(3, 1) = -r_22;
        R_l(3, 2) = -r_11;
        R_l(3, 3) = -r_12;
        R_l(4, 0) = -r_11;
        R_l(4, 1) = -r_12;
        R_l(4, 2) = r_21;
        R_l(4, 3) = r_22;
        R_l(5, 6) = 2*Px*d_3*r_11 + 2*Py*d_3*r_21 + 2*Pz*d_3*r_31 - 2*d_3*d_4*r_11*r_13 - 2*d_3*d_4*r_21*r_23 - 2*d_3*d_4*r_31*r_33;
        R_l(5, 7) = 2*Px*d_3*r_12 + 2*Py*d_3*r_22 + 2*Pz*d_3*r_32 - 2*d_3*d_4*r_12*r_13 - 2*d_3*d_4*r_22*r_23 - 2*d_3*d_4*r_32*r_33;
        R_l(6, 0) = -Px*r_31 + Pz*r_11 - d_4*r_11*r_33 + d_4*r_13*r_31;
        R_l(6, 1) = -Px*r_32 + Pz*r_12 - d_4*r_12*r_33 + d_4*r_13*r_32;
        R_l(6, 2) = Py*r_31 - Pz*r_21 + d_4*r_21*r_33 - d_4*r_23*r_31;
        R_l(6, 3) = Py*r_32 - Pz*r_22 + d_4*r_22*r_33 - d_4*r_23*r_32;
        R_l(7, 0) = Py*r_31 - Pz*r_21 + d_4*r_21*r_33 - d_4*r_23*r_31;
        R_l(7, 1) = Py*r_32 - Pz*r_22 + d_4*r_22*r_33 - d_4*r_23*r_32;
        R_l(7, 2) = Px*r_31 - Pz*r_11 + d_4*r_11*r_33 - d_4*r_13*r_31;
        R_l(7, 3) = Px*r_32 - Pz*r_12 + d_4*r_12*r_33 - d_4*r_13*r_32;
        Eigen::Matrix<double, 8, 8> R_l_mat_inv = R_l.inverse();
        for(auto r = 0; r < R_l_mat_inv.rows(); r++) {
            for(auto c = 0; c < R_l_mat_inv.cols(); c++) {
                if (std::isnan(R_l_mat_inv(r, c)) || (!std::isfinite(R_l_mat_inv(r, c)))) return;
            }
        }
        
        const double R_l_inv_00 = R_l_mat_inv(0, 0);
        const double R_l_inv_01 = R_l_mat_inv(0, 1);
        const double R_l_inv_02 = R_l_mat_inv(0, 2);
        const double R_l_inv_03 = R_l_mat_inv(0, 3);
        const double R_l_inv_04 = R_l_mat_inv(0, 4);
        const double R_l_inv_05 = R_l_mat_inv(0, 5);
        const double R_l_inv_06 = R_l_mat_inv(0, 6);
        const double R_l_inv_07 = R_l_mat_inv(0, 7);
        const double R_l_inv_10 = R_l_mat_inv(1, 0);
        const double R_l_inv_11 = R_l_mat_inv(1, 1);
        const double R_l_inv_12 = R_l_mat_inv(1, 2);
        const double R_l_inv_13 = R_l_mat_inv(1, 3);
        const double R_l_inv_14 = R_l_mat_inv(1, 4);
        const double R_l_inv_15 = R_l_mat_inv(1, 5);
        const double R_l_inv_16 = R_l_mat_inv(1, 6);
        const double R_l_inv_17 = R_l_mat_inv(1, 7);
        const double R_l_inv_20 = R_l_mat_inv(2, 0);
        const double R_l_inv_21 = R_l_mat_inv(2, 1);
        const double R_l_inv_22 = R_l_mat_inv(2, 2);
        const double R_l_inv_23 = R_l_mat_inv(2, 3);
        const double R_l_inv_24 = R_l_mat_inv(2, 4);
        const double R_l_inv_25 = R_l_mat_inv(2, 5);
        const double R_l_inv_26 = R_l_mat_inv(2, 6);
        const double R_l_inv_27 = R_l_mat_inv(2, 7);
        const double R_l_inv_30 = R_l_mat_inv(3, 0);
        const double R_l_inv_31 = R_l_mat_inv(3, 1);
        const double R_l_inv_32 = R_l_mat_inv(3, 2);
        const double R_l_inv_33 = R_l_mat_inv(3, 3);
        const double R_l_inv_34 = R_l_mat_inv(3, 4);
        const double R_l_inv_35 = R_l_mat_inv(3, 5);
        const double R_l_inv_36 = R_l_mat_inv(3, 6);
        const double R_l_inv_37 = R_l_mat_inv(3, 7);
        const double R_l_inv_40 = R_l_mat_inv(4, 0);
        const double R_l_inv_41 = R_l_mat_inv(4, 1);
        const double R_l_inv_42 = R_l_mat_inv(4, 2);
        const double R_l_inv_43 = R_l_mat_inv(4, 3);
        const double R_l_inv_44 = R_l_mat_inv(4, 4);
        const double R_l_inv_45 = R_l_mat_inv(4, 5);
        const double R_l_inv_46 = R_l_mat_inv(4, 6);
        const double R_l_inv_47 = R_l_mat_inv(4, 7);
        const double R_l_inv_50 = R_l_mat_inv(5, 0);
        const double R_l_inv_51 = R_l_mat_inv(5, 1);
        const double R_l_inv_52 = R_l_mat_inv(5, 2);
        const double R_l_inv_53 = R_l_mat_inv(5, 3);
        const double R_l_inv_54 = R_l_mat_inv(5, 4);
        const double R_l_inv_55 = R_l_mat_inv(5, 5);
        const double R_l_inv_56 = R_l_mat_inv(5, 6);
        const double R_l_inv_57 = R_l_mat_inv(5, 7);
        const double R_l_inv_60 = R_l_mat_inv(6, 0);
        const double R_l_inv_61 = R_l_mat_inv(6, 1);
        const double R_l_inv_62 = R_l_mat_inv(6, 2);
        const double R_l_inv_63 = R_l_mat_inv(6, 3);
        const double R_l_inv_64 = R_l_mat_inv(6, 4);
        const double R_l_inv_65 = R_l_mat_inv(6, 5);
        const double R_l_inv_66 = R_l_mat_inv(6, 6);
        const double R_l_inv_67 = R_l_mat_inv(6, 7);
        const double R_l_inv_70 = R_l_mat_inv(7, 0);
        const double R_l_inv_71 = R_l_mat_inv(7, 1);
        const double R_l_inv_72 = R_l_mat_inv(7, 2);
        const double R_l_inv_73 = R_l_mat_inv(7, 3);
        const double R_l_inv_74 = R_l_mat_inv(7, 4);
        const double R_l_inv_75 = R_l_mat_inv(7, 5);
        const double R_l_inv_76 = R_l_mat_inv(7, 6);
        const double R_l_inv_77 = R_l_mat_inv(7, 7);
        
        // Temp variable for efficiency
        const double x0 = R_l_inv_64*r_31;
        const double x1 = R_l_inv_74*r_32;
        const double x2 = R_l_inv_60*r_31 + R_l_inv_70*r_32;
        const double x3 = a_0*x2;
        const double x4 = R_l_inv_62*r_31 + R_l_inv_72*r_32;
        const double x5 = d_2*x4;
        const double x6 = x0 + x1 + x3 + x5;
        const double x7 = R_l_inv_66*r_31;
        const double x8 = R_l_inv_76*r_32;
        const double x9 = x7 + x8;
        const double x10 = d_2*x9;
        const double x11 = d_1*(R_l_inv_61*r_31 + R_l_inv_71*r_32);
        const double x12 = d_4*r_33;
        const double x13 = Pz - x12;
        const double x14 = x13*x4;
        const double x15 = R_l_inv_65*r_31 + R_l_inv_75*r_32;
        const double x16 = std::pow(Px, 2);
        const double x17 = std::pow(Py, 2);
        const double x18 = std::pow(Pz, 2);
        const double x19 = std::pow(d_1, 2);
        const double x20 = -x19;
        const double x21 = std::pow(d_3, 2);
        const double x22 = std::pow(r_11, 2);
        const double x23 = x21*x22;
        const double x24 = std::pow(r_21, 2);
        const double x25 = x21*x24;
        const double x26 = std::pow(r_31, 2);
        const double x27 = x21*x26;
        const double x28 = std::pow(d_4, 2);
        const double x29 = std::pow(r_13, 2)*x28;
        const double x30 = std::pow(r_23, 2)*x28;
        const double x31 = std::pow(r_33, 2)*x28;
        const double x32 = d_4*r_13;
        const double x33 = 2*Px;
        const double x34 = x32*x33;
        const double x35 = d_4*r_23;
        const double x36 = 2*Py;
        const double x37 = x35*x36;
        const double x38 = 2*Pz;
        const double x39 = x12*x38;
        const double x40 = std::pow(a_0, 2);
        const double x41 = std::pow(d_2, 2);
        const double x42 = -x40 - x41;
        const double x43 = x16 + x17 + x18 + x20 + x23 + x25 + x27 + x29 + x30 + x31 - x34 - x37 - x39 + x42;
        const double x44 = x15*x43;
        const double x45 = x11 + x14 + x44;
        const double x46 = x10 + x45;
        const double x47 = 2*d_2;
        const double x48 = x2*x47;
        const double x49 = 2*a_0;
        const double x50 = x4*x49;
        const double x51 = x49*x9;
        const double x52 = x50 + x51;
        const double x53 = x0 + x1 - x3;
        const double x54 = -x5;
        const double x55 = -x10;
        const double x56 = x45 + x55;
        const double x57 = x54 + x56;
        const double x58 = R_l_inv_67*r_31 + R_l_inv_77*r_32;
        const double x59 = 2*R_l_inv_63;
        const double x60 = 2*R_l_inv_73;
        const double x61 = r_31*x59 + r_32*x60;
        const double x62 = 2*d_2*x58 - x61;
        const double x63 = 4*d_1;
        const double x64 = x63*x9;
        const double x65 = x64 - 4;
        const double x66 = x47*x58 + x61;
        const double x67 = -x0 - x1 + x3;
        const double x68 = -x11 - x14 - x44;
        const double x69 = Px*r_11;
        const double x70 = Py*r_21;
        const double x71 = Pz*r_31;
        const double x72 = r_11*x32;
        const double x73 = r_21*x35;
        const double x74 = r_31*x12;
        const double x75 = x69 + x70 + x71 - x72 - x73 - x74;
        const double x76 = Px*r_12;
        const double x77 = Py*r_22;
        const double x78 = Pz*r_32;
        const double x79 = r_12*x32;
        const double x80 = r_22*x35;
        const double x81 = r_32*x12;
        const double x82 = x76 + x77 + x78 - x79 - x80 - x81;
        const double x83 = R_l_inv_60*x75 + R_l_inv_70*x82;
        const double x84 = a_0*x83;
        const double x85 = R_l_inv_62*x75 + R_l_inv_72*x82;
        const double x86 = d_2*x85;
        const double x87 = R_l_inv_66*x75;
        const double x88 = R_l_inv_76*x82;
        const double x89 = x87 + x88;
        const double x90 = d_2*x89;
        const double x91 = x84 + x86 + x90;
        const double x92 = d_1*(R_l_inv_61*x75 + R_l_inv_71*x82);
        const double x93 = x13*x85;
        const double x94 = R_l_inv_65*x75 + R_l_inv_75*x82;
        const double x95 = x43*x94;
        const double x96 = R_l_inv_64*x75;
        const double x97 = R_l_inv_74*x82;
        const double x98 = -d_1 + x96 + x97;
        const double x99 = -d_3*x22 - d_3*x24 - d_3*x26 + x92 + x93 + x95 + x98;
        const double x100 = x47*x83;
        const double x101 = x49*x85;
        const double x102 = x49*x89;
        const double x103 = x101 + x102;
        const double x104 = d_3*x22;
        const double x105 = d_3*x24;
        const double x106 = d_3*x26;
        const double x107 = x104 + x105 + x106 - x92 - x93 - x95;
        const double x108 = R_l_inv_67*x75 + R_l_inv_77*x82;
        const double x109 = -2*d_2*x108;
        const double x110 = -x49;
        const double x111 = x59*x75 + x60*x82;
        const double x112 = x110 + x111;
        const double x113 = x63*x89;
        const double x114 = x108*x47;
        const double x115 = x111 + x49;
        const double x116 = -x84;
        const double x117 = x116 + x90;
        const double x118 = -x86;
        const double x119 = x107 + x98;
        const double x120 = x118 + x119;
        const double x121 = -x90;
        const double x122 = x121 + x84;
        const double x123 = x119 + x86;
        const double x124 = Px*r_21;
        const double x125 = Py*r_11;
        const double x126 = r_11*x35 - r_21*x32 + x124 - x125;
        const double x127 = R_l_inv_64*x126;
        const double x128 = Px*r_22;
        const double x129 = Py*r_12;
        const double x130 = r_12*x35 - r_22*x32 + x128 - x129;
        const double x131 = R_l_inv_74*x130;
        const double x132 = d_1*(R_l_inv_61*x126 + R_l_inv_71*x130);
        const double x133 = R_l_inv_62*x126 + R_l_inv_72*x130;
        const double x134 = x13*x133;
        const double x135 = R_l_inv_65*x126 + R_l_inv_75*x130;
        const double x136 = x135*x43;
        const double x137 = x127 + x131 + x132 + x134 + x136;
        const double x138 = R_l_inv_60*x126 + R_l_inv_70*x130;
        const double x139 = a_0*x138;
        const double x140 = d_2*x133;
        const double x141 = x139 + x140;
        const double x142 = R_l_inv_66*x126;
        const double x143 = R_l_inv_76*x130;
        const double x144 = x142 + x143;
        const double x145 = d_2*x144;
        const double x146 = a_0 + x145;
        const double x147 = -x47;
        const double x148 = x144*x49;
        const double x149 = x138*x47;
        const double x150 = x133*x49;
        const double x151 = -x149 + x150;
        const double x152 = -a_0;
        const double x153 = -x145;
        const double x154 = x152 + x153;
        const double x155 = -x139;
        const double x156 = -x140;
        const double x157 = x155 + x156;
        const double x158 = R_l_inv_67*x126 + R_l_inv_77*x130;
        const double x159 = x158*x47;
        const double x160 = 2*d_1;
        const double x161 = -x160;
        const double x162 = x126*x59 + x130*x60 + x161;
        const double x163 = -x159 + x162;
        const double x164 = x144*x63;
        const double x165 = -x164;
        const double x166 = -x159 - x162;
        const double x167 = -x127 - x131 + x132 + x134 + x136;
        const double x168 = -x148;
        const double x169 = 2*x35;
        const double x170 = 2*x32;
        const double x171 = 2*r_11;
        const double x172 = r_23*x28;
        const double x173 = r_13*x172;
        const double x174 = 2*r_31;
        const double x175 = r_33*x172;
        const double x176 = std::pow(r_21, 3)*x21 - r_21*x16 + r_21*x17 - r_21*x18 + r_21*x23 + r_21*x27 - r_21*x29 + r_21*x30 - r_21*x31 + r_21*x39 + x124*x170 - x125*x170 - x169*x69 - x169*x70 - x169*x71 + x171*x173 + x174*x175 + x36*x69 + x36*x71 - x36*x74;
        const double x177 = R_l_inv_04*x176;
        const double x178 = r_21*x21;
        const double x179 = x171*x178;
        const double x180 = x174*x178;
        const double x181 = 2*x173;
        const double x182 = 2*r_32;
        const double x183 = r_12*x179 + r_12*x181 - r_22*x16 + r_22*x17 - r_22*x18 + r_22*x23 + 3*r_22*x25 + r_22*x27 - r_22*x29 + r_22*x30 - r_22*x31 + r_22*x39 + r_32*x180 + x128*x170 - x129*x170 - x169*x76 - x169*x77 - x169*x78 + x175*x182 + x36*x76 + x36*x78 - x36*x81;
        const double x184 = R_l_inv_14*x183;
        const double x185 = r_13*r_33*x28;
        const double x186 = std::pow(r_11, 3)*x21 + r_11*x16 - r_11*x17 - r_11*x18 + r_11*x25 + r_11*x27 + r_11*x29 - r_11*x30 - r_11*x31 + r_11*x39 + r_21*x181 - x124*x169 + x125*x169 - x170*x69 - x170*x70 - x170*x71 + x174*x185 + x33*x70 + x33*x71 - x33*x74;
        const double x187 = R_l_inv_24*x186;
        const double x188 = r_31*x171*x21;
        const double x189 = r_12*x16 - r_12*x17 - r_12*x18 + 3*r_12*x23 + r_12*x25 + r_12*x27 + r_12*x29 - r_12*x30 - r_12*x31 + r_12*x39 + r_22*x179 + r_22*x181 + r_32*x188 - x128*x169 + x129*x169 - x170*x76 - x170*x77 - x170*x78 + x182*x185 + x33*x77 + x33*x78 - x33*x81;
        const double x190 = R_l_inv_34*x189;
        const double x191 = Py*x22 + Py*x24 + Py*x26 - x22*x35 - x24*x35 - x26*x35;
        const double x192 = 2*d_3;
        const double x193 = x191*x192;
        const double x194 = R_l_inv_46*x193;
        const double x195 = Px*x22 + Px*x24 + Px*x26 - x22*x32 - x24*x32 - x26*x32;
        const double x196 = x192*x195;
        const double x197 = R_l_inv_56*x196;
        const double x198 = R_l_inv_06*x176;
        const double x199 = R_l_inv_26*x186;
        const double x200 = R_l_inv_16*x183;
        const double x201 = R_l_inv_36*x189;
        const double x202 = x194 + x197 + x198 + x199 + x200 + x201;
        const double x203 = d_2*x202;
        const double x204 = d_1*x49;
        const double x205 = R_l_inv_44*x193;
        const double x206 = R_l_inv_54*x196;
        const double x207 = x177 + x184 + x187 + x190 + x203 + x204 + x205 + x206;
        const double x208 = R_l_inv_00*x176 + R_l_inv_10*x183 + R_l_inv_20*x186 + R_l_inv_30*x189 + R_l_inv_40*x193 + R_l_inv_50*x196;
        const double x209 = a_0*x208;
        const double x210 = R_l_inv_02*x176 + R_l_inv_12*x183 + R_l_inv_22*x186 + R_l_inv_32*x189 + R_l_inv_42*x193 + R_l_inv_52*x196;
        const double x211 = d_2*x210;
        const double x212 = x209 + x211;
        const double x213 = R_l_inv_01*x176 + R_l_inv_11*x183 + R_l_inv_21*x186 + R_l_inv_31*x189 + R_l_inv_41*x193 + R_l_inv_51*x196;
        const double x214 = R_l_inv_05*x176 + R_l_inv_15*x183 + R_l_inv_25*x186 + R_l_inv_35*x189 + R_l_inv_45*x193 + R_l_inv_55*x196;
        const double x215 = d_1*x213 + x13*x210 + x214*x43;
        const double x216 = d_2*x63;
        const double x217 = -x216;
        const double x218 = x202*x49;
        const double x219 = x208*x47;
        const double x220 = x210*x49;
        const double x221 = -x219 + x220;
        const double x222 = -x209;
        const double x223 = -x211;
        const double x224 = x222 + x223;
        const double x225 = -x203;
        const double x226 = -x204;
        const double x227 = x177 + x184 + x187 + x190 + x205 + x206 + x225 + x226;
        const double x228 = R_l_inv_07*x176 + R_l_inv_17*x183 + R_l_inv_27*x186 + R_l_inv_37*x189 + R_l_inv_47*x193 + R_l_inv_57*x196;
        const double x229 = x228*x47;
        const double x230 = 2*x176;
        const double x231 = 2*x183;
        const double x232 = 2*x186;
        const double x233 = 2*x189;
        const double x234 = 4*d_3;
        const double x235 = x191*x234;
        const double x236 = x195*x234;
        const double x237 = R_l_inv_03*x230 + R_l_inv_13*x231 + R_l_inv_23*x232 + R_l_inv_33*x233 + R_l_inv_43*x235 + R_l_inv_53*x236 - 2*x19 + 2*x40 - 2*x41;
        const double x238 = -x229 + x237;
        const double x239 = a_0*d_2;
        const double x240 = 8*x239;
        const double x241 = x202*x63;
        const double x242 = -x229 - x237;
        const double x243 = -d_1*x213 - x13*x210 - x214*x43;
        const double x244 = -x218;
        const double x245 = R_l_inv_40*x196;
        const double x246 = R_l_inv_50*x193;
        const double x247 = R_l_inv_00*x186;
        const double x248 = R_l_inv_20*x176;
        const double x249 = R_l_inv_10*x189;
        const double x250 = R_l_inv_30*x183;
        const double x251 = x245 - x246 + x247 - x248 + x249 - x250;
        const double x252 = a_0*x251;
        const double x253 = R_l_inv_02*x186 + R_l_inv_12*x189 - R_l_inv_22*x176 - R_l_inv_32*x183 + R_l_inv_42*x196 - R_l_inv_52*x193;
        const double x254 = d_2*x253;
        const double x255 = R_l_inv_46*x196;
        const double x256 = R_l_inv_56*x193;
        const double x257 = R_l_inv_06*x186;
        const double x258 = R_l_inv_26*x176;
        const double x259 = R_l_inv_16*x189;
        const double x260 = R_l_inv_36*x183;
        const double x261 = x255 - x256 + x257 - x258 + x259 - x260;
        const double x262 = d_2*x261;
        const double x263 = x252 + x254 + x262;
        const double x264 = R_l_inv_04*x186;
        const double x265 = R_l_inv_14*x189;
        const double x266 = R_l_inv_24*x176;
        const double x267 = R_l_inv_34*x183;
        const double x268 = R_l_inv_54*x193;
        const double x269 = R_l_inv_44*x196;
        const double x270 = R_l_inv_01*x186 + R_l_inv_11*x189 - R_l_inv_21*x176 - R_l_inv_31*x183 + R_l_inv_41*x196 - R_l_inv_51*x193;
        const double x271 = R_l_inv_05*x186 + R_l_inv_15*x189 - R_l_inv_25*x176 - R_l_inv_35*x183 + R_l_inv_45*x196 - R_l_inv_55*x193;
        const double x272 = d_1*x270 + x13*x253 + x271*x43;
        const double x273 = x19 + x264 + x265 - x266 - x267 - x268 + x269 + x272 + x42;
        const double x274 = x251*x47;
        const double x275 = x253*x49;
        const double x276 = x261*x49;
        const double x277 = x275 + x276;
        const double x278 = x20 - x264 - x265 + x266 + x267 + x268 - x269 + x40 + x41;
        const double x279 = R_l_inv_07*x186 + R_l_inv_17*x189 - R_l_inv_27*x176 - R_l_inv_37*x183 + R_l_inv_47*x196 - R_l_inv_57*x193;
        const double x280 = x279*x47;
        const double x281 = a_0*x63;
        const double x282 = R_l_inv_23*x230;
        const double x283 = R_l_inv_33*x231;
        const double x284 = R_l_inv_03*x232;
        const double x285 = R_l_inv_13*x233;
        const double x286 = R_l_inv_53*x235;
        const double x287 = R_l_inv_43*x236;
        const double x288 = -x282 - x283 + x284 + x285 - x286 + x287;
        const double x289 = x281 + x288;
        const double x290 = x261*x63;
        const double x291 = -x290;
        const double x292 = -x262;
        const double x293 = x252 + x292;
        const double x294 = x272 + x278;
        const double x295 = x254 + x294;
        const double x296 = -x252;
        const double x297 = x262 + x296;
        const double x298 = -x254;
        const double x299 = x294 + x298;
        const double x300 = 2*x12;
        const double x301 = 2*x175;
        const double x302 = r_21*x301 + std::pow(r_31, 3)*x21 - r_31*x16 - r_31*x17 + r_31*x18 + r_31*x23 + r_31*x25 - r_31*x29 - r_31*x30 + r_31*x31 + r_31*x34 + r_31*x37 + x171*x185 - x300*x69 - x300*x70 - x300*x71 + x38*x69 + x38*x70 - x38*x72 - x38*x73;
        const double x303 = R_l_inv_64*x302;
        const double x304 = 2*r_12*x185 + r_12*x188 + r_22*x180 + r_22*x301 - r_32*x16 - r_32*x17 + r_32*x18 + r_32*x23 + r_32*x25 + 3*r_32*x27 - r_32*x29 - r_32*x30 + r_32*x31 + r_32*x34 + r_32*x37 - x300*x76 - x300*x77 - x300*x78 + x38*x76 + x38*x77 - x38*x79 - x38*x80;
        const double x305 = R_l_inv_74*x304;
        const double x306 = R_l_inv_60*x302 + R_l_inv_70*x304;
        const double x307 = a_0*x306;
        const double x308 = R_l_inv_66*x302 + R_l_inv_76*x304;
        const double x309 = d_2*x308;
        const double x310 = d_1*x47;
        const double x311 = x303 + x305 + x307 + x309 + x310;
        const double x312 = R_l_inv_62*x302;
        const double x313 = R_l_inv_72*x304;
        const double x314 = x312 + x313;
        const double x315 = d_2*x314;
        const double x316 = R_l_inv_61*x302 + R_l_inv_71*x304;
        const double x317 = R_l_inv_65*x302 + R_l_inv_75*x304;
        const double x318 = x104*x38;
        const double x319 = x105*x38;
        const double x320 = x106*x38;
        const double x321 = d_1*x316 + x104*x300 + x105*x300 + x106*x300 + x13*x314 + x317*x43 - x318 - x319 - x320;
        const double x322 = x315 + x321;
        const double x323 = x306*x47;
        const double x324 = -x323;
        const double x325 = x281 + x308*x49;
        const double x326 = x314*x49 + x325;
        const double x327 = -x307;
        const double x328 = -x315;
        const double x329 = x321 + x328;
        const double x330 = x327 + x329;
        const double x331 = -x309 - x310;
        const double x332 = x303 + x305;
        const double x333 = x331 + x332;
        const double x334 = R_l_inv_67*x302 + R_l_inv_77*x304;
        const double x335 = x334*x47;
        const double x336 = -x335;
        const double x337 = 4*x239;
        const double x338 = x302*x59 + x304*x60;
        const double x339 = x337 + x338;
        const double x340 = 4*x40;
        const double x341 = 4*x19 + x308*x63 + 4*x41;
        const double x342 = -x340 + x341;
        const double x343 = -x303 - x305;
        const double x344 = x307 + x343;
        const double x345 = -2*a_0*x314 + x325;
        const double x346 = x309 + x310;
        const double x347 = x15*x49;
        const double x348 = x47*(x2 - x347);
        const double x349 = 4*d_2;
        const double x350 = x47*(-x2 - x347);
        const double x351 = 4*a_0;
        const double x352 = x351*x58;
        const double x353 = 8*R_l_inv_63;
        const double x354 = 8*R_l_inv_73;
        const double x355 = x49*x94;
        const double x356 = x47*(-x355 + x83);
        const double x357 = x47*(-x355 - x83);
        const double x358 = x108*x351;
        const double x359 = -x113;
        const double x360 = -2*a_0*x135;
        const double x361 = x138 + 1;
        const double x362 = x135*x49;
        const double x363 = x158*x351;
        const double x364 = x138 - 1;
        const double x365 = -2*a_0*x214;
        const double x366 = x160 + x208;
        const double x367 = x214*x49;
        const double x368 = x228*x351;
        const double x369 = -x241;
        const double x370 = 16*d_3;
        const double x371 = R_l_inv_43*x370;
        const double x372 = 8*R_l_inv_03;
        const double x373 = 8*R_l_inv_13;
        const double x374 = x161 + x208;
        const double x375 = x271*x49;
        const double x376 = x375 + x49;
        const double x377 = x279*x351;
        const double x378 = x47*(2*a_0*x317 - x306);
        const double x379 = x160 + x308;
        const double x380 = x47*(x306 + x317*x49);
        const double x381 = x340 + x341;
        const double x382 = x0 + x1 + x3;
        const double x383 = x139 + x156;
        const double x384 = a_0 + x153;
        const double x385 = x149 + x150;
        const double x386 = x145 + x152;
        const double x387 = x140 + x155;
        const double x388 = x177 + x184 + x187 + x190 + x205 + x206 + x215;
        const double x389 = x204 + x225;
        const double x390 = x209 + x223;
        const double x391 = x219 + x220;
        const double x392 = x211 + x222;
        const double x393 = x203 + x226;
        const double x394 = x177 + x184 + x187 + x190 + x205 + x206 + x243;
        const double x395 = -x337 + x338;
        
        Eigen::Matrix<double, 6, 9> A;
        A.setZero();
        A(0, 0) = -x46 - x6;
        A(0, 1) = x48 - x52;
        A(0, 2) = -x53 - x57;
        A(0, 3) = x62;
        A(0, 4) = x65;
        A(0, 5) = x66;
        A(0, 6) = -x5 - x56 - x67;
        A(0, 7) = x48 - x50 + x51;
        A(0, 8) = x55 + x6 + x68;
        A(1, 0) = -x91 - x99;
        A(1, 1) = x100 - x103;
        A(1, 2) = d_1 + x107 + x91 - x96 - x97;
        A(1, 3) = -x109 - x112;
        A(1, 4) = x113;
        A(1, 5) = x114 + x115;
        A(1, 6) = x117 + x120;
        A(1, 7) = x100 - x101 + x102;
        A(1, 8) = x122 + x123;
        A(2, 0) = x137 + x141 + x146;
        A(2, 1) = x147 + x148 + x151;
        A(2, 2) = x137 + x154 + x157;
        A(2, 3) = x163;
        A(2, 4) = x165;
        A(2, 5) = x166;
        A(2, 6) = x141 + x154 + x167;
        A(2, 7) = x151 + x168 + x47;
        A(2, 8) = x146 + x157 + x167;
        A(3, 0) = x207 + x212 + x215;
        A(3, 1) = x217 + x218 + x221;
        A(3, 2) = x215 + x224 + x227;
        A(3, 3) = x238;
        A(3, 4) = -x240 - x241;
        A(3, 5) = x242;
        A(3, 6) = -x207 - x224 - x243;
        A(3, 7) = x216 + x221 + x244;
        A(3, 8) = -x212 - x227 - x243;
        A(4, 0) = x263 + x273;
        A(4, 1) = -x274 + x277;
        A(4, 2) = d_1*x270 + x13*x253 - x263 + x271*x43 - x278;
        A(4, 3) = -x280 + x289;
        A(4, 4) = x291;
        A(4, 5) = -x280 + x281 - x288;
        A(4, 6) = x293 + x295;
        A(4, 7) = 2*a_0*x253 - x274 - x276;
        A(4, 8) = x297 + x299;
        A(5, 0) = x311 + x322;
        A(5, 1) = x324 + x326;
        A(5, 2) = x330 + x333;
        A(5, 3) = x336 + x339;
        A(5, 4) = -x342;
        A(5, 5) = -x335 - x339;
        A(5, 6) = x322 + x331 + x344;
        A(5, 7) = -x323 - x345;
        A(5, 8) = x330 + x343 + x346;
        
        Eigen::Matrix<double, 6, 9> B;
        B.setZero();
        B(0, 0) = x348;
        B(0, 1) = x349*(x4 + x9);
        B(0, 2) = x350;
        B(0, 3) = -x352 + x64 - 4;
        B(0, 4) = r_31*x353 + r_32*x354;
        B(0, 5) = -x352 - x65;
        B(0, 6) = x348;
        B(0, 7) = x349*(x4 - x7 - x8);
        B(0, 8) = x350;
        B(1, 0) = x356;
        B(1, 1) = x349*(x85 + x89);
        B(1, 2) = x357;
        B(1, 3) = -x358 - x359;
        B(1, 4) = x353*x75 + x354*x82;
        B(1, 5) = -x113 - x358;
        B(1, 6) = x356;
        B(1, 7) = x349*(x85 - x87 - x88);
        B(1, 8) = x357;
        B(2, 0) = x47*(-x360 - x361);
        B(2, 1) = x349*(-x133 - x144);
        B(2, 2) = x47*(x361 + x362);
        B(2, 3) = x165 + x363;
        B(2, 4) = 8*d_1 - x126*x353 - x130*x354;
        B(2, 5) = x164 + x363;
        B(2, 6) = x47*(-x360 - x364);
        B(2, 7) = x349*(-x133 + x142 + x143);
        B(2, 8) = x47*(x362 + x364);
        B(3, 0) = x47*(-x365 - x366);
        B(3, 1) = x349*(-x202 - x210);
        B(3, 2) = x47*(x366 + x367);
        B(3, 3) = x368 + x369;
        B(3, 4) = -8*R_l_inv_23*x186 - 8*R_l_inv_33*x189 - R_l_inv_53*x195*x370 - x176*x372 - x183*x373 + 8*x19 - x191*x371 + 8*x40 + 8*x41;
        B(3, 5) = x241 + x368;
        B(3, 6) = x47*(-x365 - x374);
        B(3, 7) = x349*(x194 + x197 + x198 + x199 + x200 + x201 - x210);
        B(3, 8) = x47*(x367 + x374);
        B(4, 0) = x47*(-x245 + x246 - x247 + x248 - x249 + x250 + x376);
        B(4, 1) = x349*(-x253 - x261);
        B(4, 2) = x47*(x251 + x376);
        B(4, 3) = x291 + x377;
        B(4, 4) = 8*R_l_inv_23*x176 + 8*R_l_inv_33*x183 + 16*R_l_inv_53*d_3*x191 - x186*x372 - x189*x373 - x195*x371;
        B(4, 5) = x290 + x377;
        B(4, 6) = x47*(2*a_0*x271 - x251 - x49);
        B(4, 7) = x349*(-x253 + x255 - x256 + x257 - x258 + x259 - x260);
        B(4, 8) = x47*(x110 + x251 + x375);
        B(5, 0) = x378;
        B(5, 1) = x349*(-x314 - x379);
        B(5, 2) = x380;
        B(5, 3) = 4*a_0*x334 - x381;
        B(5, 4) = -x302*x353 - x304*x354;
        B(5, 5) = x334*x351 + x381;
        B(5, 6) = x378;
        B(5, 7) = x349*(-x312 - x313 + x379);
        B(5, 8) = x380;
        
        Eigen::Matrix<double, 6, 9> C;
        C.setZero();
        C(0, 0) = -x382 - x57;
        C(0, 1) = -x48 - x52;
        C(0, 2) = -x46 - x5 - x53;
        C(0, 3) = x66;
        C(0, 4) = -x65;
        C(0, 5) = x62;
        C(0, 6) = -x46 - x54 - x67;
        C(0, 7) = -x48 - x50 + x51;
        C(0, 8) = x10 + x382 + x54 + x68;
        C(1, 0) = -x118 - x122 - x99;
        C(1, 1) = -x100 - x103;
        C(1, 2) = -x117 - x86 - x99;
        C(1, 3) = x112 + x114;
        C(1, 4) = x359;
        C(1, 5) = -x109 - x115;
        C(1, 6) = x116 + x121 + x123;
        C(1, 7) = -x100 - x101 + x102;
        C(1, 8) = x120 + x84 + x90;
        C(2, 0) = x137 + x383 + x384;
        C(2, 1) = x148 + x385 + x47;
        C(2, 2) = x137 + x386 + x387;
        C(2, 3) = x166;
        C(2, 4) = x164;
        C(2, 5) = x163;
        C(2, 6) = x167 + x383 + x386;
        C(2, 7) = x147 + x168 + x385;
        C(2, 8) = x167 + x384 + x387;
        C(3, 0) = x388 + x389 + x390;
        C(3, 1) = x216 + x218 + x391;
        C(3, 2) = x388 + x392 + x393;
        C(3, 3) = x242;
        C(3, 4) = -x240 - x369;
        C(3, 5) = x238;
        C(3, 6) = -x389 - x392 - x394;
        C(3, 7) = x217 + x244 + x391;
        C(3, 8) = -x390 - x393 - x394;
        C(4, 0) = x273 + x293 + x298;
        C(4, 1) = x274 + x277;
        C(4, 2) = x254 + x273 + x297;
        C(4, 3) = -x280 - x289;
        C(4, 4) = x290;
        C(4, 5) = -x280 - x281 - x282 - x283 + x284 + x285 - x286 + x287;
        C(4, 6) = x252 + x262 + x299;
        C(4, 7) = x274 + x275 - x276;
        C(4, 8) = x292 + x295 + x296;
        C(5, 0) = x307 + x329 + x333;
        C(5, 1) = x323 + x326;
        C(5, 2) = x322 + x327 + x332 + x346;
        C(5, 3) = -x335 - x395;
        C(5, 4) = x342;
        C(5, 5) = x336 + x395;
        C(5, 6) = x329 + x344 + x346;
        C(5, 7) = -x324 - x345;
        C(5, 8) = d_1*x316 + 2*d_3*d_4*r_33*x22 + 2*d_3*d_4*r_33*x24 + 2*d_3*d_4*r_33*x26 + x13*x314 - x311 + x317*x43 - x318 - x319 - x320 - x328;
        
        // Invoke the solver
        std::array<double, 16> solution_buffer;
        int n_solutions = yaik_cpp::general_6dof_internal::computeSolutionFromTanhalfLME(A, B, C, &solution_buffer);
        
        for(auto i = 0; i < n_solutions; i++)
        {
            auto solution_i = make_raw_solution();
            solution_i[3] = solution_buffer[i];
            int appended_idx = append_solution_to_queue(solution_i);
            add_input_index_to(2, appended_idx);
        };
    };
    // Invoke the processor
    General6DoFNumericalReduceSolutionNode_node_1_solve_th_2_processor();
    // Finish code for general_6dof solution node 0
    
    // Code for non-branch dispatcher node 2
    // Actually, there is no code
    
    // Code for explicit solution node 3, solved variable is th_3
    auto ExplicitSolutionNode_node_3_solve_th_3_processor = [&]() -> void
    {
        const auto& this_node_input_index = node_index_workspace.get_input_indices_for_node(2);
        const bool this_input_valid = node_index_workspace.is_input_indices_valid_for_node(2);
        if (!this_input_valid)
            return;
        
        // The solution of non-root node 3
        for(int i = 0; i < this_node_input_index.size(); i++)
        {
            int node_input_i_idx_in_queue = this_node_input_index[i];
            if (!solution_queue.is_solution_valid(node_input_i_idx_in_queue))
                continue;
            const auto& this_solution = solution_queue.get_solution(node_input_i_idx_in_queue);
            const double th_2 = this_solution[3];
            
            const bool condition_0 = 2*std::fabs(d_1*d_3) >= zero_tolerance || 2*std::fabs(a_0*d_3*std::cos(th_2)) >= zero_tolerance || std::fabs(std::pow(Px, 2) - 2*Px*d_4*r_13 + std::pow(Py, 2) - 2*Py*d_4*r_23 + std::pow(Pz, 2) - 2*Pz*d_4*r_33 - std::pow(a_0, 2) + 2*a_0*d_2*std::sin(th_2) - std::pow(d_1, 2) - std::pow(d_2, 2) - std::pow(d_3, 2) + std::pow(d_4, 2)*std::pow(r_13, 2) + std::pow(d_4, 2)*std::pow(r_23, 2) + std::pow(d_4, 2)*std::pow(r_33, 2)) >= zero_tolerance;
            if (condition_0)
            {
                // Temp variable for efficiency
                const double x0 = std::cos(th_2);
                const double x1 = 2*d_3;
                const double x2 = std::atan2(-a_0*x0*x1, -d_1*x1);
                const double x3 = std::pow(d_1, 2);
                const double x4 = std::pow(d_3, 2);
                const double x5 = 4*x4;
                const double x6 = std::pow(a_0, 2);
                const double x7 = 2*d_4;
                const double x8 = std::pow(d_4, 2);
                const double x9 = std::pow(Px, 2) - Px*r_13*x7 + std::pow(Py, 2) - Py*r_23*x7 + std::pow(Pz, 2) - Pz*r_33*x7 + 2*a_0*d_2*std::sin(th_2) - std::pow(d_2, 2) + std::pow(r_13, 2)*x8 + std::pow(r_23, 2)*x8 + std::pow(r_33, 2)*x8 - x3 - x4 - x6;
                const double x10 = safe_sqrt(std::pow(x0, 2)*x5*x6 + x3*x5 - std::pow(x9, 2));
                // End of temp variables
                RawSolution solution_0(this_solution);
                solution_0[4] = x2 + std::atan2(x10, x9);
                int appended_idx = append_solution_to_queue(solution_0);
                add_input_index_to(4, appended_idx);
            }
            
            const bool condition_1 = 2*std::fabs(d_1*d_3) >= zero_tolerance || 2*std::fabs(a_0*d_3*std::cos(th_2)) >= zero_tolerance || std::fabs(std::pow(Px, 2) - 2*Px*d_4*r_13 + std::pow(Py, 2) - 2*Py*d_4*r_23 + std::pow(Pz, 2) - 2*Pz*d_4*r_33 - std::pow(a_0, 2) + 2*a_0*d_2*std::sin(th_2) - std::pow(d_1, 2) - std::pow(d_2, 2) - std::pow(d_3, 2) + std::pow(d_4, 2)*std::pow(r_13, 2) + std::pow(d_4, 2)*std::pow(r_23, 2) + std::pow(d_4, 2)*std::pow(r_33, 2)) >= zero_tolerance;
            if (condition_1)
            {
                // Temp variable for efficiency
                const double x0 = std::cos(th_2);
                const double x1 = 2*d_3;
                const double x2 = std::atan2(-a_0*x0*x1, -d_1*x1);
                const double x3 = std::pow(d_1, 2);
                const double x4 = std::pow(d_3, 2);
                const double x5 = 4*x4;
                const double x6 = std::pow(a_0, 2);
                const double x7 = 2*d_4;
                const double x8 = std::pow(d_4, 2);
                const double x9 = std::pow(Px, 2) - Px*r_13*x7 + std::pow(Py, 2) - Py*r_23*x7 + std::pow(Pz, 2) - Pz*r_33*x7 + 2*a_0*d_2*std::sin(th_2) - std::pow(d_2, 2) + std::pow(r_13, 2)*x8 + std::pow(r_23, 2)*x8 + std::pow(r_33, 2)*x8 - x3 - x4 - x6;
                const double x10 = safe_sqrt(std::pow(x0, 2)*x5*x6 + x3*x5 - std::pow(x9, 2));
                // End of temp variables
                const double tmp_sol_value = x2 + std::atan2(-x10, x9);
                solution_queue.get_solution(node_input_i_idx_in_queue)[4] = tmp_sol_value;
                add_input_index_to(4, node_input_i_idx_in_queue);
            }
            else
            {
                solution_queue.solutions_validity[node_input_i_idx_in_queue] = false;
            }
        }
    };
    // Invoke the processor
    ExplicitSolutionNode_node_3_solve_th_3_processor();
    // Finish code for explicit solution node 2
    
    // Code for equation all-zero dispatcher node 4
    auto EquationAllZeroDispatcherNode_node_4_processor = [&]() -> void
    {
        const auto& this_node_input_index = node_index_workspace.get_input_indices_for_node(4);
        const bool this_input_valid = node_index_workspace.is_input_indices_valid_for_node(4);
        if (!this_input_valid)
            return;
        for(int i = 0; i < this_node_input_index.size(); i++)
        {
            int node_input_i_idx_in_queue = this_node_input_index[i];
            if (!solution_queue.is_solution_valid(node_input_i_idx_in_queue))
                continue;
            const auto& this_solution = solution_queue.get_solution(node_input_i_idx_in_queue);
            const double th_3 = this_solution[4];
            const bool checked_result = std::fabs(Px - d_4*r_13) <= 9.9999999999999995e-7 && std::fabs(Py - d_4*r_23) <= 9.9999999999999995e-7 && std::fabs(d_1 - d_3*std::cos(th_3)) <= 9.9999999999999995e-7;
            if (!checked_result)  // To non-degenerate node
                add_input_index_to(5, node_input_i_idx_in_queue);
        }
    };
    
    // Invoke the processor
    EquationAllZeroDispatcherNode_node_4_processor();
    // Finish code for equation all-zero dispatcher node 4
    
    // Code for explicit solution node 5, solved variable is th_0
    auto ExplicitSolutionNode_node_5_solve_th_0_processor = [&]() -> void
    {
        const auto& this_node_input_index = node_index_workspace.get_input_indices_for_node(5);
        const bool this_input_valid = node_index_workspace.is_input_indices_valid_for_node(5);
        if (!this_input_valid)
            return;
        
        // The solution of non-root node 5
        for(int i = 0; i < this_node_input_index.size(); i++)
        {
            int node_input_i_idx_in_queue = this_node_input_index[i];
            if (!solution_queue.is_solution_valid(node_input_i_idx_in_queue))
                continue;
            const auto& this_solution = solution_queue.get_solution(node_input_i_idx_in_queue);
            const double th_3 = this_solution[4];
            
            const bool condition_0 = std::fabs(Px - d_4*r_13) >= zero_tolerance || std::fabs(Py - d_4*r_23) >= zero_tolerance || std::fabs(d_1 - d_3*std::cos(th_3)) >= zero_tolerance;
            if (condition_0)
            {
                // Temp variable for efficiency
                const double x0 = Px - d_4*r_13;
                const double x1 = -Py + d_4*r_23;
                const double x2 = std::atan2(x0, x1);
                const double x3 = -d_1 + d_3*std::cos(th_3);
                const double x4 = safe_sqrt(std::pow(x0, 2) + std::pow(x1, 2) - std::pow(x3, 2));
                // End of temp variables
                RawSolution solution_0(this_solution);
                solution_0[0] = x2 + std::atan2(x4, x3);
                int appended_idx = append_solution_to_queue(solution_0);
                add_input_index_to(6, appended_idx);
            }
            
            const bool condition_1 = std::fabs(Px - d_4*r_13) >= zero_tolerance || std::fabs(Py - d_4*r_23) >= zero_tolerance || std::fabs(d_1 - d_3*std::cos(th_3)) >= zero_tolerance;
            if (condition_1)
            {
                // Temp variable for efficiency
                const double x0 = Px - d_4*r_13;
                const double x1 = -Py + d_4*r_23;
                const double x2 = std::atan2(x0, x1);
                const double x3 = -d_1 + d_3*std::cos(th_3);
                const double x4 = safe_sqrt(std::pow(x0, 2) + std::pow(x1, 2) - std::pow(x3, 2));
                // End of temp variables
                const double tmp_sol_value = x2 + std::atan2(-x4, x3);
                solution_queue.get_solution(node_input_i_idx_in_queue)[0] = tmp_sol_value;
                add_input_index_to(6, node_input_i_idx_in_queue);
            }
            else
            {
                solution_queue.solutions_validity[node_input_i_idx_in_queue] = false;
            }
        }
    };
    // Invoke the processor
    ExplicitSolutionNode_node_5_solve_th_0_processor();
    // Finish code for explicit solution node 5
    
    // Code for solved_variable dispatcher node 6
    auto SolvedVariableDispatcherNode_node_6_processor = [&]() -> void
    {
        const auto& this_node_input_index = node_index_workspace.get_input_indices_for_node(6);
        const bool this_input_valid = node_index_workspace.is_input_indices_valid_for_node(6);
        if (!this_input_valid)
            return;
        
        for(int i = 0; i < this_node_input_index.size(); i++)
        {
            int node_input_i_idx_in_queue = this_node_input_index[i];
            if (!solution_queue.is_solution_valid(node_input_i_idx_in_queue))
                continue;
            const auto& this_solution = solution_queue.get_solution(node_input_i_idx_in_queue);
            bool taken_by_degenerate = false;
            const double th_3 = this_solution[4];
            
            const bool degenerate_valid_0 = std::fabs(th_3) <= 9.9999999999999995e-7;
            if (degenerate_valid_0)
            {
                taken_by_degenerate = true;
                add_input_index_to(12, node_input_i_idx_in_queue);
            }
            
            const bool degenerate_valid_1 = std::fabs(th_3 - M_PI) <= 9.9999999999999995e-7;
            if (degenerate_valid_1)
            {
                taken_by_degenerate = true;
                add_input_index_to(17, node_input_i_idx_in_queue);
            }
            
            if (!taken_by_degenerate)
                add_input_index_to(7, node_input_i_idx_in_queue);
        }
    };
    
    // Invoke the processor
    SolvedVariableDispatcherNode_node_6_processor();
    // Finish code for solved_variable dispatcher node 6
    
    // Code for explicit solution node 17, solved variable is th_1th_2th_4_soa
    auto ExplicitSolutionNode_node_17_solve_th_1th_2th_4_soa_processor = [&]() -> void
    {
        const auto& this_node_input_index = node_index_workspace.get_input_indices_for_node(17);
        const bool this_input_valid = node_index_workspace.is_input_indices_valid_for_node(17);
        if (!this_input_valid)
            return;
        
        // The solution of non-root node 17
        for(int i = 0; i < this_node_input_index.size(); i++)
        {
            int node_input_i_idx_in_queue = this_node_input_index[i];
            if (!solution_queue.is_solution_valid(node_input_i_idx_in_queue))
                continue;
            const auto& this_solution = solution_queue.get_solution(node_input_i_idx_in_queue);
            const double th_0 = this_solution[0];
            
            const bool condition_0 = 1 >= zero_tolerance || std::fabs(r_33) >= zero_tolerance || std::fabs(r_13*std::cos(th_0) + r_23*std::sin(th_0)) >= zero_tolerance;
            if (condition_0)
            {
                // Temp variable for efficiency
                // End of temp variables
                const double tmp_sol_value = std::atan2(r_13*std::cos(th_0) + r_23*std::sin(th_0), r_33);
                solution_queue.get_solution(node_input_i_idx_in_queue)[2] = tmp_sol_value;
                add_input_index_to(18, node_input_i_idx_in_queue);
            }
            else
            {
                solution_queue.solutions_validity[node_input_i_idx_in_queue] = false;
            }
        }
    };
    // Invoke the processor
    ExplicitSolutionNode_node_17_solve_th_1th_2th_4_soa_processor();
    // Finish code for explicit solution node 17
    
    // Code for non-branch dispatcher node 18
    // Actually, there is no code
    
    // Code for explicit solution node 19, solved variable is th_5
    auto ExplicitSolutionNode_node_19_solve_th_5_processor = [&]() -> void
    {
        const auto& this_node_input_index = node_index_workspace.get_input_indices_for_node(18);
        const bool this_input_valid = node_index_workspace.is_input_indices_valid_for_node(18);
        if (!this_input_valid)
            return;
        
        // The solution of non-root node 19
        for(int i = 0; i < this_node_input_index.size(); i++)
        {
            int node_input_i_idx_in_queue = this_node_input_index[i];
            if (!solution_queue.is_solution_valid(node_input_i_idx_in_queue))
                continue;
            const auto& this_solution = solution_queue.get_solution(node_input_i_idx_in_queue);
            const double th_0 = this_solution[0];
            
            const bool condition_0 = 1 >= zero_tolerance || std::fabs(r_11*std::sin(th_0) - r_21*std::cos(th_0)) >= zero_tolerance || std::fabs(r_12*std::sin(th_0) - r_22*std::cos(th_0)) >= zero_tolerance;
            if (condition_0)
            {
                // Temp variable for efficiency
                const double x0 = std::sin(th_0);
                const double x1 = std::cos(th_0);
                // End of temp variables
                const double tmp_sol_value = std::atan2(r_11*x0 - r_21*x1, r_12*x0 - r_22*x1);
                solution_queue.get_solution(node_input_i_idx_in_queue)[6] = tmp_sol_value;
                add_input_index_to(20, node_input_i_idx_in_queue);
            }
            else
            {
                solution_queue.solutions_validity[node_input_i_idx_in_queue] = false;
            }
        }
    };
    // Invoke the processor
    ExplicitSolutionNode_node_19_solve_th_5_processor();
    // Finish code for explicit solution node 18
    
    // Code for solved_variable dispatcher node 20
    auto SolvedVariableDispatcherNode_node_20_processor = [&]() -> void
    {
        const auto& this_node_input_index = node_index_workspace.get_input_indices_for_node(20);
        const bool this_input_valid = node_index_workspace.is_input_indices_valid_for_node(20);
        if (!this_input_valid)
            return;
        
        for(int i = 0; i < this_node_input_index.size(); i++)
        {
            int node_input_i_idx_in_queue = this_node_input_index[i];
            if (!solution_queue.is_solution_valid(node_input_i_idx_in_queue))
                continue;
            const auto& this_solution = solution_queue.get_solution(node_input_i_idx_in_queue);
            bool taken_by_degenerate = false;
            const double th_2 = this_solution[3];
            
            const bool degenerate_valid_0 = std::fabs(th_2 - 1.0/2.0*M_PI) <= 9.9999999999999995e-7;
            if (degenerate_valid_0)
            {
                taken_by_degenerate = true;
                add_input_index_to(30, node_input_i_idx_in_queue);
            }
            
            const bool degenerate_valid_1 = std::fabs(th_2 + M_PI_2) <= 9.9999999999999995e-7;
            if (degenerate_valid_1)
            {
                taken_by_degenerate = true;
                add_input_index_to(33, node_input_i_idx_in_queue);
            }
            
            if (!taken_by_degenerate)
                add_input_index_to(21, node_input_i_idx_in_queue);
        }
    };
    
    // Invoke the processor
    SolvedVariableDispatcherNode_node_20_processor();
    // Finish code for solved_variable dispatcher node 20
    
    // Code for explicit solution node 33, solved variable is th_1
    auto ExplicitSolutionNode_node_33_solve_th_1_processor = [&]() -> void
    {
        const auto& this_node_input_index = node_index_workspace.get_input_indices_for_node(33);
        const bool this_input_valid = node_index_workspace.is_input_indices_valid_for_node(33);
        if (!this_input_valid)
            return;
        
        // The solution of non-root node 33
        for(int i = 0; i < this_node_input_index.size(); i++)
        {
            int node_input_i_idx_in_queue = this_node_input_index[i];
            if (!solution_queue.is_solution_valid(node_input_i_idx_in_queue))
                continue;
            const auto& this_solution = solution_queue.get_solution(node_input_i_idx_in_queue);
            const double th_0 = this_solution[0];
            const double th_3 = this_solution[4];
            
            const bool condition_0 = std::fabs(d_3*std::sin(th_3)) >= 9.9999999999999995e-7 || std::fabs(a_0 + d_2) >= 9.9999999999999995e-7;
            if (condition_0)
            {
                // Temp variable for efficiency
                const double x0 = a_0 + d_2;
                const double x1 = Pz - d_4*r_33;
                const double x2 = std::cos(th_0);
                const double x3 = std::sin(th_0);
                const double x4 = Px*x2 + Py*x3 - d_4*r_13*x2 - d_4*r_23*x3;
                const double x5 = d_3*std::sin(th_3);
                // End of temp variables
                const double tmp_sol_value = std::atan2(-x0*x1 - x4*x5, x0*x4 - x1*x5);
                solution_queue.get_solution(node_input_i_idx_in_queue)[1] = tmp_sol_value;
                add_input_index_to(34, node_input_i_idx_in_queue);
            }
            else
            {
                solution_queue.solutions_validity[node_input_i_idx_in_queue] = false;
            }
        }
    };
    // Invoke the processor
    ExplicitSolutionNode_node_33_solve_th_1_processor();
    // Finish code for explicit solution node 33
    
    // Code for non-branch dispatcher node 34
    // Actually, there is no code
    
    // Code for explicit solution node 35, solved variable is th_4
    auto ExplicitSolutionNode_node_35_solve_th_4_processor = [&]() -> void
    {
        const auto& this_node_input_index = node_index_workspace.get_input_indices_for_node(34);
        const bool this_input_valid = node_index_workspace.is_input_indices_valid_for_node(34);
        if (!this_input_valid)
            return;
        
        // The solution of non-root node 35
        for(int i = 0; i < this_node_input_index.size(); i++)
        {
            int node_input_i_idx_in_queue = this_node_input_index[i];
            if (!solution_queue.is_solution_valid(node_input_i_idx_in_queue))
                continue;
            const auto& this_solution = solution_queue.get_solution(node_input_i_idx_in_queue);
            const double th_0 = this_solution[0];
            const double th_1 = this_solution[1];
            const double th_3 = this_solution[4];
            
            const bool condition_0 = 1 >= zero_tolerance || std::fabs(r_13*(std::sin(th_0)*std::sin(th_3) + std::sin(th_1)*std::cos(th_0)*std::cos(th_3)) + r_23*(std::sin(th_0)*std::sin(th_1)*std::cos(th_3) - std::sin(th_3)*std::cos(th_0)) + r_33*std::cos(th_1)*std::cos(th_3)) >= zero_tolerance || std::fabs(r_13*std::cos(th_0)*std::cos(th_1) + r_23*std::sin(th_0)*std::cos(th_1) - r_33*std::sin(th_1)) >= zero_tolerance;
            if (condition_0)
            {
                // Temp variable for efficiency
                const double x0 = std::cos(th_1);
                const double x1 = std::cos(th_3);
                const double x2 = std::sin(th_0);
                const double x3 = std::sin(th_3);
                const double x4 = std::cos(th_0);
                const double x5 = std::sin(th_1);
                const double x6 = x1*x5;
                // End of temp variables
                const double tmp_sol_value = std::atan2(-r_13*(x2*x3 + x4*x6) - r_23*(x2*x6 - x3*x4) - r_33*x0*x1, -r_13*x0*x4 - r_23*x0*x2 + r_33*x5);
                solution_queue.get_solution(node_input_i_idx_in_queue)[5] = tmp_sol_value;
            }
            else
            {
                solution_queue.solutions_validity[node_input_i_idx_in_queue] = false;
            }
        }
    };
    // Invoke the processor
    ExplicitSolutionNode_node_35_solve_th_4_processor();
    // Finish code for explicit solution node 34
    
    // Code for explicit solution node 30, solved variable is th_1
    auto ExplicitSolutionNode_node_30_solve_th_1_processor = [&]() -> void
    {
        const auto& this_node_input_index = node_index_workspace.get_input_indices_for_node(30);
        const bool this_input_valid = node_index_workspace.is_input_indices_valid_for_node(30);
        if (!this_input_valid)
            return;
        
        // The solution of non-root node 30
        for(int i = 0; i < this_node_input_index.size(); i++)
        {
            int node_input_i_idx_in_queue = this_node_input_index[i];
            if (!solution_queue.is_solution_valid(node_input_i_idx_in_queue))
                continue;
            const auto& this_solution = solution_queue.get_solution(node_input_i_idx_in_queue);
            const double th_0 = this_solution[0];
            const double th_3 = this_solution[4];
            
            const bool condition_0 = std::fabs(d_3*std::sin(th_3)) >= 9.9999999999999995e-7 || std::fabs(a_0 - d_2) >= 9.9999999999999995e-7;
            if (condition_0)
            {
                // Temp variable for efficiency
                const double x0 = a_0 - d_2;
                const double x1 = Pz - d_4*r_33;
                const double x2 = std::cos(th_0);
                const double x3 = std::sin(th_0);
                const double x4 = Px*x2 + Py*x3 - d_4*r_13*x2 - d_4*r_23*x3;
                const double x5 = d_3*std::sin(th_3);
                // End of temp variables
                const double tmp_sol_value = std::atan2(-x0*x1 + x4*x5, x0*x4 + x1*x5);
                solution_queue.get_solution(node_input_i_idx_in_queue)[1] = tmp_sol_value;
                add_input_index_to(31, node_input_i_idx_in_queue);
            }
            else
            {
                solution_queue.solutions_validity[node_input_i_idx_in_queue] = false;
            }
        }
    };
    // Invoke the processor
    ExplicitSolutionNode_node_30_solve_th_1_processor();
    // Finish code for explicit solution node 30
    
    // Code for non-branch dispatcher node 31
    // Actually, there is no code
    
    // Code for explicit solution node 32, solved variable is th_4
    auto ExplicitSolutionNode_node_32_solve_th_4_processor = [&]() -> void
    {
        const auto& this_node_input_index = node_index_workspace.get_input_indices_for_node(31);
        const bool this_input_valid = node_index_workspace.is_input_indices_valid_for_node(31);
        if (!this_input_valid)
            return;
        
        // The solution of non-root node 32
        for(int i = 0; i < this_node_input_index.size(); i++)
        {
            int node_input_i_idx_in_queue = this_node_input_index[i];
            if (!solution_queue.is_solution_valid(node_input_i_idx_in_queue))
                continue;
            const auto& this_solution = solution_queue.get_solution(node_input_i_idx_in_queue);
            const double th_0 = this_solution[0];
            const double th_1 = this_solution[1];
            const double th_3 = this_solution[4];
            
            const bool condition_0 = 1 >= zero_tolerance || std::fabs(-r_13*(std::sin(th_0)*std::sin(th_3) - std::sin(th_1)*std::cos(th_0)*std::cos(th_3)) + r_23*(std::sin(th_0)*std::sin(th_1)*std::cos(th_3) + std::sin(th_3)*std::cos(th_0)) + r_33*std::cos(th_1)*std::cos(th_3)) >= zero_tolerance || std::fabs(r_13*std::cos(th_0)*std::cos(th_1) + r_23*std::sin(th_0)*std::cos(th_1) - r_33*std::sin(th_1)) >= zero_tolerance;
            if (condition_0)
            {
                // Temp variable for efficiency
                const double x0 = std::cos(th_1);
                const double x1 = std::cos(th_3);
                const double x2 = std::sin(th_3);
                const double x3 = std::cos(th_0);
                const double x4 = std::sin(th_0);
                const double x5 = std::sin(th_1);
                const double x6 = x1*x5;
                // End of temp variables
                const double tmp_sol_value = std::atan2(-r_13*(x2*x4 - x3*x6) + r_23*(x2*x3 + x4*x6) + r_33*x0*x1, r_13*x0*x3 + r_23*x0*x4 - r_33*x5);
                solution_queue.get_solution(node_input_i_idx_in_queue)[5] = tmp_sol_value;
            }
            else
            {
                solution_queue.solutions_validity[node_input_i_idx_in_queue] = false;
            }
        }
    };
    // Invoke the processor
    ExplicitSolutionNode_node_32_solve_th_4_processor();
    // Finish code for explicit solution node 31
    
    // Code for explicit solution node 21, solved variable is th_1
    auto ExplicitSolutionNode_node_21_solve_th_1_processor = [&]() -> void
    {
        const auto& this_node_input_index = node_index_workspace.get_input_indices_for_node(21);
        const bool this_input_valid = node_index_workspace.is_input_indices_valid_for_node(21);
        if (!this_input_valid)
            return;
        
        // The solution of non-root node 21
        for(int i = 0; i < this_node_input_index.size(); i++)
        {
            int node_input_i_idx_in_queue = this_node_input_index[i];
            if (!solution_queue.is_solution_valid(node_input_i_idx_in_queue))
                continue;
            const auto& this_solution = solution_queue.get_solution(node_input_i_idx_in_queue);
            const double th_0 = this_solution[0];
            const double th_1th_2th_4_soa = this_solution[2];
            const double th_2 = this_solution[3];
            
            const bool condition_0 = std::fabs(d_2*std::cos(th_2)) >= 9.9999999999999995e-7 || std::fabs(a_0 - d_2*std::sin(th_2)) >= 9.9999999999999995e-7;
            if (condition_0)
            {
                // Temp variable for efficiency
                const double x0 = Pz - d_4*r_33;
                const double x1 = a_0 - d_2*std::sin(th_2);
                const double x2 = Px*std::cos(th_0) + Py*std::sin(th_0) - d_4*std::sin(th_1th_2th_4_soa);
                const double x3 = d_2*std::cos(th_2);
                // End of temp variables
                const double tmp_sol_value = std::atan2(-x0*x1 - x2*x3, -x0*x3 + x1*x2);
                solution_queue.get_solution(node_input_i_idx_in_queue)[1] = tmp_sol_value;
                add_input_index_to(22, node_input_i_idx_in_queue);
            }
            else
            {
                solution_queue.solutions_validity[node_input_i_idx_in_queue] = false;
            }
        }
    };
    // Invoke the processor
    ExplicitSolutionNode_node_21_solve_th_1_processor();
    // Finish code for explicit solution node 21
    
    // Code for non-branch dispatcher node 22
    // Actually, there is no code
    
    // Code for explicit solution node 23, solved variable is th_4
    auto ExplicitSolutionNode_node_23_solve_th_4_processor = [&]() -> void
    {
        const auto& this_node_input_index = node_index_workspace.get_input_indices_for_node(22);
        const bool this_input_valid = node_index_workspace.is_input_indices_valid_for_node(22);
        if (!this_input_valid)
            return;
        
        // The solution of non-root node 23
        for(int i = 0; i < this_node_input_index.size(); i++)
        {
            int node_input_i_idx_in_queue = this_node_input_index[i];
            if (!solution_queue.is_solution_valid(node_input_i_idx_in_queue))
                continue;
            const auto& this_solution = solution_queue.get_solution(node_input_i_idx_in_queue);
            const double th_1 = this_solution[1];
            const double th_1th_2th_4_soa = this_solution[2];
            const double th_2 = this_solution[3];
            
            const bool condition_0 = 1 >= zero_tolerance;
            if (condition_0)
            {
                // Temp variable for efficiency
                // End of temp variables
                const double tmp_sol_value = -th_1 + th_1th_2th_4_soa - th_2;
                solution_queue.get_solution(node_input_i_idx_in_queue)[5] = tmp_sol_value;
            }
            else
            {
                solution_queue.solutions_validity[node_input_i_idx_in_queue] = false;
            }
        }
    };
    // Invoke the processor
    ExplicitSolutionNode_node_23_solve_th_4_processor();
    // Finish code for explicit solution node 22
    
    // Code for explicit solution node 12, solved variable is th_5
    auto ExplicitSolutionNode_node_12_solve_th_5_processor = [&]() -> void
    {
        const auto& this_node_input_index = node_index_workspace.get_input_indices_for_node(12);
        const bool this_input_valid = node_index_workspace.is_input_indices_valid_for_node(12);
        if (!this_input_valid)
            return;
        
        // The solution of non-root node 12
        for(int i = 0; i < this_node_input_index.size(); i++)
        {
            int node_input_i_idx_in_queue = this_node_input_index[i];
            if (!solution_queue.is_solution_valid(node_input_i_idx_in_queue))
                continue;
            const auto& this_solution = solution_queue.get_solution(node_input_i_idx_in_queue);
            const double th_0 = this_solution[0];
            
            const bool condition_0 = 1 >= zero_tolerance || std::fabs(r_11*std::sin(th_0) - r_21*std::cos(th_0)) >= zero_tolerance || std::fabs(r_12*std::sin(th_0) - r_22*std::cos(th_0)) >= zero_tolerance;
            if (condition_0)
            {
                // Temp variable for efficiency
                const double x0 = std::sin(th_0);
                const double x1 = std::cos(th_0);
                // End of temp variables
                const double tmp_sol_value = std::atan2(-r_11*x0 + r_21*x1, -r_12*x0 + r_22*x1);
                solution_queue.get_solution(node_input_i_idx_in_queue)[6] = tmp_sol_value;
                add_input_index_to(13, node_input_i_idx_in_queue);
            }
            else
            {
                solution_queue.solutions_validity[node_input_i_idx_in_queue] = false;
            }
        }
    };
    // Invoke the processor
    ExplicitSolutionNode_node_12_solve_th_5_processor();
    // Finish code for explicit solution node 12
    
    // Code for solved_variable dispatcher node 13
    auto SolvedVariableDispatcherNode_node_13_processor = [&]() -> void
    {
        const auto& this_node_input_index = node_index_workspace.get_input_indices_for_node(13);
        const bool this_input_valid = node_index_workspace.is_input_indices_valid_for_node(13);
        if (!this_input_valid)
            return;
        
        for(int i = 0; i < this_node_input_index.size(); i++)
        {
            int node_input_i_idx_in_queue = this_node_input_index[i];
            if (!solution_queue.is_solution_valid(node_input_i_idx_in_queue))
                continue;
            const auto& this_solution = solution_queue.get_solution(node_input_i_idx_in_queue);
            bool taken_by_degenerate = false;
            const double th_2 = this_solution[3];
            
            const bool degenerate_valid_0 = std::fabs(th_2 - 1.0/2.0*M_PI) <= 9.9999999999999995e-7;
            if (degenerate_valid_0)
            {
                taken_by_degenerate = true;
                add_input_index_to(24, node_input_i_idx_in_queue);
            }
            
            const bool degenerate_valid_1 = std::fabs(th_2 + M_PI_2) <= 9.9999999999999995e-7;
            if (degenerate_valid_1)
            {
                taken_by_degenerate = true;
                add_input_index_to(27, node_input_i_idx_in_queue);
            }
            
            if (!taken_by_degenerate)
                add_input_index_to(14, node_input_i_idx_in_queue);
        }
    };
    
    // Invoke the processor
    SolvedVariableDispatcherNode_node_13_processor();
    // Finish code for solved_variable dispatcher node 13
    
    // Code for explicit solution node 27, solved variable is th_1
    auto ExplicitSolutionNode_node_27_solve_th_1_processor = [&]() -> void
    {
        const auto& this_node_input_index = node_index_workspace.get_input_indices_for_node(27);
        const bool this_input_valid = node_index_workspace.is_input_indices_valid_for_node(27);
        if (!this_input_valid)
            return;
        
        // The solution of non-root node 27
        for(int i = 0; i < this_node_input_index.size(); i++)
        {
            int node_input_i_idx_in_queue = this_node_input_index[i];
            if (!solution_queue.is_solution_valid(node_input_i_idx_in_queue))
                continue;
            const auto& this_solution = solution_queue.get_solution(node_input_i_idx_in_queue);
            const double th_0 = this_solution[0];
            const double th_3 = this_solution[4];
            
            const bool condition_0 = std::fabs(d_3*std::sin(th_3)) >= 9.9999999999999995e-7 || std::fabs(a_0 + d_2) >= 9.9999999999999995e-7;
            if (condition_0)
            {
                // Temp variable for efficiency
                const double x0 = a_0 + d_2;
                const double x1 = Pz - d_4*r_33;
                const double x2 = std::cos(th_0);
                const double x3 = std::sin(th_0);
                const double x4 = Px*x2 + Py*x3 - d_4*r_13*x2 - d_4*r_23*x3;
                const double x5 = d_3*std::sin(th_3);
                // End of temp variables
                const double tmp_sol_value = std::atan2(-x0*x1 - x4*x5, x0*x4 - x1*x5);
                solution_queue.get_solution(node_input_i_idx_in_queue)[1] = tmp_sol_value;
                add_input_index_to(28, node_input_i_idx_in_queue);
            }
            else
            {
                solution_queue.solutions_validity[node_input_i_idx_in_queue] = false;
            }
        }
    };
    // Invoke the processor
    ExplicitSolutionNode_node_27_solve_th_1_processor();
    // Finish code for explicit solution node 27
    
    // Code for non-branch dispatcher node 28
    // Actually, there is no code
    
    // Code for explicit solution node 29, solved variable is th_4
    auto ExplicitSolutionNode_node_29_solve_th_4_processor = [&]() -> void
    {
        const auto& this_node_input_index = node_index_workspace.get_input_indices_for_node(28);
        const bool this_input_valid = node_index_workspace.is_input_indices_valid_for_node(28);
        if (!this_input_valid)
            return;
        
        // The solution of non-root node 29
        for(int i = 0; i < this_node_input_index.size(); i++)
        {
            int node_input_i_idx_in_queue = this_node_input_index[i];
            if (!solution_queue.is_solution_valid(node_input_i_idx_in_queue))
                continue;
            const auto& this_solution = solution_queue.get_solution(node_input_i_idx_in_queue);
            const double th_0 = this_solution[0];
            const double th_1 = this_solution[1];
            const double th_3 = this_solution[4];
            
            const bool condition_0 = 1 >= zero_tolerance || std::fabs(r_13*(std::sin(th_0)*std::sin(th_3) + std::sin(th_1)*std::cos(th_0)*std::cos(th_3)) + r_23*(std::sin(th_0)*std::sin(th_1)*std::cos(th_3) - std::sin(th_3)*std::cos(th_0)) + r_33*std::cos(th_1)*std::cos(th_3)) >= zero_tolerance || std::fabs(r_13*std::cos(th_0)*std::cos(th_1) + r_23*std::sin(th_0)*std::cos(th_1) - r_33*std::sin(th_1)) >= zero_tolerance;
            if (condition_0)
            {
                // Temp variable for efficiency
                const double x0 = std::cos(th_1);
                const double x1 = std::cos(th_3);
                const double x2 = std::sin(th_0);
                const double x3 = std::sin(th_3);
                const double x4 = std::cos(th_0);
                const double x5 = std::sin(th_1);
                const double x6 = x1*x5;
                // End of temp variables
                const double tmp_sol_value = std::atan2(-r_13*(x2*x3 + x4*x6) - r_23*(x2*x6 - x3*x4) - r_33*x0*x1, -r_13*x0*x4 - r_23*x0*x2 + r_33*x5);
                solution_queue.get_solution(node_input_i_idx_in_queue)[5] = tmp_sol_value;
            }
            else
            {
                solution_queue.solutions_validity[node_input_i_idx_in_queue] = false;
            }
        }
    };
    // Invoke the processor
    ExplicitSolutionNode_node_29_solve_th_4_processor();
    // Finish code for explicit solution node 28
    
    // Code for explicit solution node 24, solved variable is th_1
    auto ExplicitSolutionNode_node_24_solve_th_1_processor = [&]() -> void
    {
        const auto& this_node_input_index = node_index_workspace.get_input_indices_for_node(24);
        const bool this_input_valid = node_index_workspace.is_input_indices_valid_for_node(24);
        if (!this_input_valid)
            return;
        
        // The solution of non-root node 24
        for(int i = 0; i < this_node_input_index.size(); i++)
        {
            int node_input_i_idx_in_queue = this_node_input_index[i];
            if (!solution_queue.is_solution_valid(node_input_i_idx_in_queue))
                continue;
            const auto& this_solution = solution_queue.get_solution(node_input_i_idx_in_queue);
            const double th_0 = this_solution[0];
            const double th_3 = this_solution[4];
            
            const bool condition_0 = std::fabs(d_3*std::sin(th_3)) >= 9.9999999999999995e-7 || std::fabs(a_0 - d_2) >= 9.9999999999999995e-7;
            if (condition_0)
            {
                // Temp variable for efficiency
                const double x0 = a_0 - d_2;
                const double x1 = Pz - d_4*r_33;
                const double x2 = std::cos(th_0);
                const double x3 = std::sin(th_0);
                const double x4 = Px*x2 + Py*x3 - d_4*r_13*x2 - d_4*r_23*x3;
                const double x5 = d_3*std::sin(th_3);
                // End of temp variables
                const double tmp_sol_value = std::atan2(-x0*x1 + x4*x5, x0*x4 + x1*x5);
                solution_queue.get_solution(node_input_i_idx_in_queue)[1] = tmp_sol_value;
                add_input_index_to(25, node_input_i_idx_in_queue);
            }
            else
            {
                solution_queue.solutions_validity[node_input_i_idx_in_queue] = false;
            }
        }
    };
    // Invoke the processor
    ExplicitSolutionNode_node_24_solve_th_1_processor();
    // Finish code for explicit solution node 24
    
    // Code for non-branch dispatcher node 25
    // Actually, there is no code
    
    // Code for explicit solution node 26, solved variable is th_4
    auto ExplicitSolutionNode_node_26_solve_th_4_processor = [&]() -> void
    {
        const auto& this_node_input_index = node_index_workspace.get_input_indices_for_node(25);
        const bool this_input_valid = node_index_workspace.is_input_indices_valid_for_node(25);
        if (!this_input_valid)
            return;
        
        // The solution of non-root node 26
        for(int i = 0; i < this_node_input_index.size(); i++)
        {
            int node_input_i_idx_in_queue = this_node_input_index[i];
            if (!solution_queue.is_solution_valid(node_input_i_idx_in_queue))
                continue;
            const auto& this_solution = solution_queue.get_solution(node_input_i_idx_in_queue);
            const double th_0 = this_solution[0];
            const double th_1 = this_solution[1];
            const double th_3 = this_solution[4];
            
            const bool condition_0 = 1 >= zero_tolerance || std::fabs(-r_13*(std::sin(th_0)*std::sin(th_3) - std::sin(th_1)*std::cos(th_0)*std::cos(th_3)) + r_23*(std::sin(th_0)*std::sin(th_1)*std::cos(th_3) + std::sin(th_3)*std::cos(th_0)) + r_33*std::cos(th_1)*std::cos(th_3)) >= zero_tolerance || std::fabs(r_13*std::cos(th_0)*std::cos(th_1) + r_23*std::sin(th_0)*std::cos(th_1) - r_33*std::sin(th_1)) >= zero_tolerance;
            if (condition_0)
            {
                // Temp variable for efficiency
                const double x0 = std::cos(th_1);
                const double x1 = std::cos(th_3);
                const double x2 = std::sin(th_3);
                const double x3 = std::cos(th_0);
                const double x4 = std::sin(th_0);
                const double x5 = std::sin(th_1);
                const double x6 = x1*x5;
                // End of temp variables
                const double tmp_sol_value = std::atan2(-r_13*(x2*x4 - x3*x6) + r_23*(x2*x3 + x4*x6) + r_33*x0*x1, r_13*x0*x3 + r_23*x0*x4 - r_33*x5);
                solution_queue.get_solution(node_input_i_idx_in_queue)[5] = tmp_sol_value;
            }
            else
            {
                solution_queue.solutions_validity[node_input_i_idx_in_queue] = false;
            }
        }
    };
    // Invoke the processor
    ExplicitSolutionNode_node_26_solve_th_4_processor();
    // Finish code for explicit solution node 25
    
    // Code for explicit solution node 14, solved variable is th_1
    auto ExplicitSolutionNode_node_14_solve_th_1_processor = [&]() -> void
    {
        const auto& this_node_input_index = node_index_workspace.get_input_indices_for_node(14);
        const bool this_input_valid = node_index_workspace.is_input_indices_valid_for_node(14);
        if (!this_input_valid)
            return;
        
        // The solution of non-root node 14
        for(int i = 0; i < this_node_input_index.size(); i++)
        {
            int node_input_i_idx_in_queue = this_node_input_index[i];
            if (!solution_queue.is_solution_valid(node_input_i_idx_in_queue))
                continue;
            const auto& this_solution = solution_queue.get_solution(node_input_i_idx_in_queue);
            const double th_0 = this_solution[0];
            const double th_2 = this_solution[3];
            
            const bool condition_0 = std::fabs(d_2*std::cos(th_2)) >= 9.9999999999999995e-7 || std::fabs(a_0 - d_2*std::sin(th_2)) >= 9.9999999999999995e-7;
            if (condition_0)
            {
                // Temp variable for efficiency
                const double x0 = Pz - d_4*r_33;
                const double x1 = a_0 - d_2*std::sin(th_2);
                const double x2 = std::cos(th_0);
                const double x3 = std::sin(th_0);
                const double x4 = Px*x2 + Py*x3 - d_4*r_13*x2 - d_4*r_23*x3;
                const double x5 = d_2*std::cos(th_2);
                // End of temp variables
                const double tmp_sol_value = std::atan2(-x0*x1 - x4*x5, -x0*x5 + x1*x4);
                solution_queue.get_solution(node_input_i_idx_in_queue)[1] = tmp_sol_value;
                add_input_index_to(15, node_input_i_idx_in_queue);
            }
            else
            {
                solution_queue.solutions_validity[node_input_i_idx_in_queue] = false;
            }
        }
    };
    // Invoke the processor
    ExplicitSolutionNode_node_14_solve_th_1_processor();
    // Finish code for explicit solution node 14
    
    // Code for non-branch dispatcher node 15
    // Actually, there is no code
    
    // Code for explicit solution node 16, solved variable is th_4
    auto ExplicitSolutionNode_node_16_solve_th_4_processor = [&]() -> void
    {
        const auto& this_node_input_index = node_index_workspace.get_input_indices_for_node(15);
        const bool this_input_valid = node_index_workspace.is_input_indices_valid_for_node(15);
        if (!this_input_valid)
            return;
        
        // The solution of non-root node 16
        for(int i = 0; i < this_node_input_index.size(); i++)
        {
            int node_input_i_idx_in_queue = this_node_input_index[i];
            if (!solution_queue.is_solution_valid(node_input_i_idx_in_queue))
                continue;
            const auto& this_solution = solution_queue.get_solution(node_input_i_idx_in_queue);
            const double th_0 = this_solution[0];
            const double th_1 = this_solution[1];
            const double th_2 = this_solution[3];
            
            const bool condition_0 = 1 >= zero_tolerance || std::fabs(-r_13*(-std::sin(th_1)*std::sin(th_2) + std::cos(th_1)*std::cos(th_2))*std::cos(th_0) - r_23*(-std::sin(th_1)*std::sin(th_2) + std::cos(th_1)*std::cos(th_2))*std::sin(th_0) + r_33*(std::sin(th_1)*std::cos(th_2) + std::sin(th_2)*std::cos(th_1))) >= zero_tolerance || std::fabs(r_13*(std::sin(th_1)*std::cos(th_2) + std::sin(th_2)*std::cos(th_1))*std::cos(th_0) + r_23*(std::sin(th_1)*std::cos(th_2) + std::sin(th_2)*std::cos(th_1))*std::sin(th_0) + r_33*(-std::sin(th_1)*std::sin(th_2) + std::cos(th_1)*std::cos(th_2))) >= zero_tolerance;
            if (condition_0)
            {
                // Temp variable for efficiency
                const double x0 = std::sin(th_1);
                const double x1 = std::cos(th_2);
                const double x2 = std::sin(th_2);
                const double x3 = std::cos(th_1);
                const double x4 = x0*x1 + x2*x3;
                const double x5 = -x0*x2 + x1*x3;
                const double x6 = r_13*std::cos(th_0);
                const double x7 = r_23*std::sin(th_0);
                // End of temp variables
                const double tmp_sol_value = std::atan2(r_33*x4 - x5*x6 - x5*x7, r_33*x5 + x4*x6 + x4*x7);
                solution_queue.get_solution(node_input_i_idx_in_queue)[5] = tmp_sol_value;
            }
            else
            {
                solution_queue.solutions_validity[node_input_i_idx_in_queue] = false;
            }
        }
    };
    // Invoke the processor
    ExplicitSolutionNode_node_16_solve_th_4_processor();
    // Finish code for explicit solution node 15
    
    // Code for explicit solution node 7, solved variable is th_4
    auto ExplicitSolutionNode_node_7_solve_th_4_processor = [&]() -> void
    {
        const auto& this_node_input_index = node_index_workspace.get_input_indices_for_node(7);
        const bool this_input_valid = node_index_workspace.is_input_indices_valid_for_node(7);
        if (!this_input_valid)
            return;
        
        // The solution of non-root node 7
        for(int i = 0; i < this_node_input_index.size(); i++)
        {
            int node_input_i_idx_in_queue = this_node_input_index[i];
            if (!solution_queue.is_solution_valid(node_input_i_idx_in_queue))
                continue;
            const auto& this_solution = solution_queue.get_solution(node_input_i_idx_in_queue);
            const double th_0 = this_solution[0];
            const double th_3 = this_solution[4];
            
            const bool condition_0 = std::fabs((r_13*std::sin(th_0) - r_23*std::cos(th_0))/std::sin(th_3)) <= 1;
            if (condition_0)
            {
                // Temp variable for efficiency
                const double x0 = safe_asin((-r_13*std::sin(th_0) + r_23*std::cos(th_0))/std::sin(th_3));
                // End of temp variables
                RawSolution solution_0(this_solution);
                solution_0[5] = x0;
                int appended_idx = append_solution_to_queue(solution_0);
                add_input_index_to(8, appended_idx);
            }
            
            const bool condition_1 = std::fabs((r_13*std::sin(th_0) - r_23*std::cos(th_0))/std::sin(th_3)) <= 1;
            if (condition_1)
            {
                // Temp variable for efficiency
                const double x0 = safe_asin((-r_13*std::sin(th_0) + r_23*std::cos(th_0))/std::sin(th_3));
                // End of temp variables
                const double tmp_sol_value = M_PI - x0;
                solution_queue.get_solution(node_input_i_idx_in_queue)[5] = tmp_sol_value;
                add_input_index_to(8, node_input_i_idx_in_queue);
            }
            else
            {
                solution_queue.solutions_validity[node_input_i_idx_in_queue] = false;
            }
        }
    };
    // Invoke the processor
    ExplicitSolutionNode_node_7_solve_th_4_processor();
    // Finish code for explicit solution node 7
    
    // Code for equation all-zero dispatcher node 8
    auto EquationAllZeroDispatcherNode_node_8_processor = [&]() -> void
    {
        const auto& this_node_input_index = node_index_workspace.get_input_indices_for_node(8);
        const bool this_input_valid = node_index_workspace.is_input_indices_valid_for_node(8);
        if (!this_input_valid)
            return;
        for(int i = 0; i < this_node_input_index.size(); i++)
        {
            int node_input_i_idx_in_queue = this_node_input_index[i];
            if (!solution_queue.is_solution_valid(node_input_i_idx_in_queue))
                continue;
            const auto& this_solution = solution_queue.get_solution(node_input_i_idx_in_queue);
            const double th_2 = this_solution[3];
            const double th_3 = this_solution[4];
            const bool checked_result = std::fabs(d_2*std::cos(th_2) - d_3*std::sin(th_2)*std::sin(th_3)) <= 9.9999999999999995e-7 && std::fabs(-a_0 + d_2*std::sin(th_2) + d_3*std::sin(th_3)*std::cos(th_2)) <= 9.9999999999999995e-7;
            if (!checked_result)  // To non-degenerate node
                add_input_index_to(9, node_input_i_idx_in_queue);
        }
    };
    
    // Invoke the processor
    EquationAllZeroDispatcherNode_node_8_processor();
    // Finish code for equation all-zero dispatcher node 8
    
    // Code for explicit solution node 9, solved variable is th_1
    auto ExplicitSolutionNode_node_9_solve_th_1_processor = [&]() -> void
    {
        const auto& this_node_input_index = node_index_workspace.get_input_indices_for_node(9);
        const bool this_input_valid = node_index_workspace.is_input_indices_valid_for_node(9);
        if (!this_input_valid)
            return;
        
        // The solution of non-root node 9
        for(int i = 0; i < this_node_input_index.size(); i++)
        {
            int node_input_i_idx_in_queue = this_node_input_index[i];
            if (!solution_queue.is_solution_valid(node_input_i_idx_in_queue))
                continue;
            const auto& this_solution = solution_queue.get_solution(node_input_i_idx_in_queue);
            const double th_0 = this_solution[0];
            const double th_2 = this_solution[3];
            const double th_3 = this_solution[4];
            
            const bool condition_0 = std::fabs(d_2*std::cos(th_2) - d_3*std::sin(th_2)*std::sin(th_3)) >= 9.9999999999999995e-7 || std::fabs(-a_0 + d_2*std::sin(th_2) + d_3*std::sin(th_3)*std::cos(th_2)) >= 9.9999999999999995e-7;
            if (condition_0)
            {
                // Temp variable for efficiency
                const double x0 = Pz - d_4*r_33;
                const double x1 = std::sin(th_2);
                const double x2 = std::cos(th_2);
                const double x3 = std::sin(th_3);
                const double x4 = a_0 - d_2*x1 - d_3*x2*x3;
                const double x5 = -d_2*x2 + d_3*x1*x3;
                const double x6 = std::cos(th_0);
                const double x7 = std::sin(th_0);
                const double x8 = Px*x6 + Py*x7 - d_4*r_13*x6 - d_4*r_23*x7;
                // End of temp variables
                const double tmp_sol_value = std::atan2(-x0*x4 + x5*x8, x0*x5 + x4*x8);
                solution_queue.get_solution(node_input_i_idx_in_queue)[1] = tmp_sol_value;
                add_input_index_to(10, node_input_i_idx_in_queue);
            }
            else
            {
                solution_queue.solutions_validity[node_input_i_idx_in_queue] = false;
            }
        }
    };
    // Invoke the processor
    ExplicitSolutionNode_node_9_solve_th_1_processor();
    // Finish code for explicit solution node 9
    
    // Code for non-branch dispatcher node 10
    // Actually, there is no code
    
    // Code for explicit solution node 11, solved variable is th_5
    auto ExplicitSolutionNode_node_11_solve_th_5_processor = [&]() -> void
    {
        const auto& this_node_input_index = node_index_workspace.get_input_indices_for_node(10);
        const bool this_input_valid = node_index_workspace.is_input_indices_valid_for_node(10);
        if (!this_input_valid)
            return;
        
        // The solution of non-root node 11
        for(int i = 0; i < this_node_input_index.size(); i++)
        {
            int node_input_i_idx_in_queue = this_node_input_index[i];
            if (!solution_queue.is_solution_valid(node_input_i_idx_in_queue))
                continue;
            const auto& this_solution = solution_queue.get_solution(node_input_i_idx_in_queue);
            const double th_0 = this_solution[0];
            const double th_1 = this_solution[1];
            const double th_2 = this_solution[3];
            const double th_3 = this_solution[4];
            
            const bool condition_0 = 1 >= zero_tolerance || std::fabs(r_11*((std::sin(th_1)*std::sin(th_2) - std::cos(th_1)*std::cos(th_2))*std::sin(th_3)*std::cos(th_0) + std::sin(th_0)*std::cos(th_3)) + r_21*((std::sin(th_1)*std::sin(th_2) - std::cos(th_1)*std::cos(th_2))*std::sin(th_0)*std::sin(th_3) - std::cos(th_0)*std::cos(th_3)) + r_31*(std::sin(th_1)*std::cos(th_2) + std::sin(th_2)*std::cos(th_1))*std::sin(th_3)) >= zero_tolerance || std::fabs(r_12*((std::sin(th_1)*std::sin(th_2) - std::cos(th_1)*std::cos(th_2))*std::sin(th_3)*std::cos(th_0) + std::sin(th_0)*std::cos(th_3)) + r_22*((std::sin(th_1)*std::sin(th_2) - std::cos(th_1)*std::cos(th_2))*std::sin(th_0)*std::sin(th_3) - std::cos(th_0)*std::cos(th_3)) + r_32*(std::sin(th_1)*std::cos(th_2) + std::sin(th_2)*std::cos(th_1))*std::sin(th_3)) >= zero_tolerance;
            if (condition_0)
            {
                // Temp variable for efficiency
                const double x0 = std::sin(th_3);
                const double x1 = std::sin(th_1);
                const double x2 = std::cos(th_2);
                const double x3 = std::sin(th_2);
                const double x4 = std::cos(th_1);
                const double x5 = x0*(x1*x2 + x3*x4);
                const double x6 = std::sin(th_0);
                const double x7 = std::cos(th_3);
                const double x8 = std::cos(th_0);
                const double x9 = x1*x3 - x2*x4;
                const double x10 = x0*x8*x9 + x6*x7;
                const double x11 = -x0*x6*x9 + x7*x8;
                // End of temp variables
                const double tmp_sol_value = std::atan2(-r_11*x10 + r_21*x11 - r_31*x5, -r_12*x10 + r_22*x11 - r_32*x5);
                solution_queue.get_solution(node_input_i_idx_in_queue)[6] = tmp_sol_value;
            }
            else
            {
                solution_queue.solutions_validity[node_input_i_idx_in_queue] = false;
            }
        }
    };
    // Invoke the processor
    ExplicitSolutionNode_node_11_solve_th_5_processor();
    // Finish code for explicit solution node 10
    
    // Collect the output
    for(int i = 0; i < solution_queue.size(); i++)
    {
        if(!solution_queue.solutions_validity[i])
            continue;
        const auto& raw_ik_out_i = solution_queue.get_solution(i);
        std::array<double, robot_nq> new_ik_i;
        const double value_at_0 = raw_ik_out_i[0];  // th_0
        new_ik_i[0] = value_at_0;
        const double value_at_1 = raw_ik_out_i[1];  // th_1
        new_ik_i[1] = value_at_1;
        const double value_at_2 = raw_ik_out_i[3];  // th_2
        new_ik_i[2] = value_at_2;
        const double value_at_3 = raw_ik_out_i[4];  // th_3
        new_ik_i[3] = value_at_3;
        const double value_at_4 = raw_ik_out_i[5];  // th_4
        new_ik_i[4] = value_at_4;
        const double value_at_5 = raw_ik_out_i[6];  // th_5
        new_ik_i[5] = value_at_5;
        ik_output.push_back(new_ik_i);
    }
}

static void computeRawIK(const Eigen::Matrix4d& T_ee_raw, RawIKWorksace& workspace)
{
    workspace.raw_ik_out.clear();
    workspace.raw_ik_out.reserve(max_n_solutions);
    computeRawIK(T_ee_raw, workspace.solution_queue, workspace.node_index_workspace, workspace.raw_ik_out);
}

static void computeIKUnChecked(const Eigen::Matrix4d& T_ee, RawIKWorksace& workspace, std::vector<std::array<double, robot_nq>>& ik_output)
{
    const Eigen::Matrix4d& T_ee_raw = endEffectorTargetOriginalToRaw(T_ee);
    computeRawIK(T_ee_raw, workspace);
    const auto& raw_ik_out = workspace.raw_ik_out;
    ik_output.clear();
    for(int i = 0; i < raw_ik_out.size(); i++)
    {
        auto ik_out_i = raw_ik_out[i];
        ik_out_i[0] -= th_0_offset_original2raw;
        ik_out_i[1] -= th_1_offset_original2raw;
        ik_out_i[2] -= th_2_offset_original2raw;
        ik_out_i[3] -= th_3_offset_original2raw;
        ik_out_i[4] -= th_4_offset_original2raw;
        ik_out_i[5] -= th_5_offset_original2raw;
        ik_output.push_back(ik_out_i);
    }
}

static void wrapAngleToPi(std::vector<std::array<double, robot_nq>>& ik_output)
{
    for(int i = 0; i < ik_output.size(); i++)
    {
        // Wrap angular value to [-pi, pi]
        auto& solution_i = ik_output[i];
        // Revolute unknown th_0
        while(solution_i[0] > M_PI)
            solution_i[0] -= 2 * M_PI;
        while(solution_i[0] < - M_PI)
            solution_i[0] += 2 * M_PI;
        // Revolute unknown th_1
        while(solution_i[1] > M_PI)
            solution_i[1] -= 2 * M_PI;
        while(solution_i[1] < - M_PI)
            solution_i[1] += 2 * M_PI;
        // Revolute unknown th_2
        while(solution_i[2] > M_PI)
            solution_i[2] -= 2 * M_PI;
        while(solution_i[2] < - M_PI)
            solution_i[2] += 2 * M_PI;
        // Revolute unknown th_3
        while(solution_i[3] > M_PI)
            solution_i[3] -= 2 * M_PI;
        while(solution_i[3] < - M_PI)
            solution_i[3] += 2 * M_PI;
        // Revolute unknown th_4
        while(solution_i[4] > M_PI)
            solution_i[4] -= 2 * M_PI;
        while(solution_i[4] < - M_PI)
            solution_i[4] += 2 * M_PI;
        // Revolute unknown th_5
        while(solution_i[5] > M_PI)
            solution_i[5] -= 2 * M_PI;
        while(solution_i[5] < - M_PI)
            solution_i[5] += 2 * M_PI;
    }
}

static void computeIK(const Eigen::Matrix4d& T_ee, RawIKWorksace& workspace, std::vector<std::array<double, robot_nq>>& ik_output)
{
    const Eigen::Matrix4d& T_ee_raw = endEffectorTargetOriginalToRaw(T_ee);
    computeRawIK(T_ee_raw, workspace);
    const auto& raw_ik_out = workspace.raw_ik_out;
    ik_output.clear();
    for(int i = 0; i < raw_ik_out.size(); i++)
    {
        auto ik_out_i = raw_ik_out[i];
        ik_out_i[0] -= th_0_offset_original2raw;
        ik_out_i[1] -= th_1_offset_original2raw;
        ik_out_i[2] -= th_2_offset_original2raw;
        ik_out_i[3] -= th_3_offset_original2raw;
        ik_out_i[4] -= th_4_offset_original2raw;
        ik_out_i[5] -= th_5_offset_original2raw;
        const Eigen::Matrix4d ee_pose_i = computeFK(ik_out_i);
        double ee_pose_diff = (ee_pose_i - T_ee).squaredNorm();
        if (ee_pose_diff < pose_tolerance)
            ik_output.push_back(ik_out_i);
    }
    
    if (!ik_output.empty())
    {
        wrapAngleToPi(ik_output);
        removeDuplicate<robot_nq>(ik_output, zero_tolerance);
        return;
    }
    
    // Disturbing method for degenerate handling
    Eigen::Matrix4d T_ee_raw_disturbed = yaik_cpp::disturbTransform(T_ee_raw);
    Eigen::Matrix4d T_ee_disturbed = endEffectorTargetRawToOriginal(T_ee_raw_disturbed);
    computeRawIK(T_ee_raw_disturbed, workspace);
    const auto& raw_ik_out_disturb = workspace.raw_ik_out;
    for(int i = 0; i < raw_ik_out_disturb.size(); i++)
    {
        auto ik_out_i = raw_ik_out_disturb[i];
        ik_out_i[0] -= th_0_offset_original2raw;
        ik_out_i[1] -= th_1_offset_original2raw;
        ik_out_i[2] -= th_2_offset_original2raw;
        ik_out_i[3] -= th_3_offset_original2raw;
        ik_out_i[4] -= th_4_offset_original2raw;
        ik_out_i[5] -= th_5_offset_original2raw;
        Eigen::Matrix4d ee_pose_i = computeFK(ik_out_i);
        double ee_pose_diff = (ee_pose_i - T_ee_disturbed).squaredNorm();
        if (ee_pose_diff > pose_tolerance_degenerate)
            continue;
        
        // Try numerical refinement
        yaik_cpp::numericalRefinement<robot_nq>(computeFK, computeTwistJacobian, T_ee, ik_out_i);
        ee_pose_i = computeFK(ik_out_i);
        ee_pose_diff = (ee_pose_i - T_ee).squaredNorm();
        if (ee_pose_diff < pose_tolerance_degenerate)
            ik_output.push_back(ik_out_i);
    }
    
    wrapAngleToPi(ik_output);
    removeDuplicate<robot_nq>(ik_output, zero_tolerance);
}

static std::vector<std::array<double, robot_nq>> computeIK(const Eigen::Matrix4d& T_ee)
{
    std::vector<std::array<double, robot_nq>> ik_output;
    RawIKWorksace raw_ik_workspace;
    computeIK(T_ee, raw_ik_workspace, ik_output);
    return ik_output;
}

}; // struct cobotta_ik

// Code below for debug
void test_ik_solve_cobotta()
{
    std::array<double, cobotta_ik::robot_nq> theta;
    std::random_device rd;
    std::uniform_real_distribution<double> distribution;
    for(auto i = 0; i < theta.size(); i++)
        theta[i] = distribution(rd);
    const Eigen::Matrix4d ee_pose = cobotta_ik::computeFK(theta);
    auto ik_output = cobotta_ik::computeIK(ee_pose);
    for(int i = 0; i < ik_output.size(); i++)
    {
        Eigen::Matrix4d ee_pose_i = cobotta_ik::computeFK(ik_output[i]);
        double ee_pose_diff = (ee_pose_i - ee_pose).norm();
        std::cout << "For solution " << i << " Pose different with ground-truth " << ee_pose_diff << std::endl;
    }
}

int main()
{
    test_ik_solve_cobotta();
}
