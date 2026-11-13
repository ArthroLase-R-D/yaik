#include "yaik_cpp_common.h"

using namespace yaik_cpp;

struct articulating_arm_ik {

// Constants for solver
static constexpr int robot_nq = 6;
static constexpr int max_n_solutions = 128;
static constexpr int n_tree_nodes = 27;
static constexpr int intermediate_solution_size = 7;
static constexpr double pose_tolerance = 1e-6;
static constexpr double pose_tolerance_degenerate = 1e-4;
static constexpr double zero_tolerance = 1e-6;
using RawIKWorksace = ::yaik_cpp::RawIkWorkspaceTemplate<intermediate_solution_size, max_n_solutions, robot_nq>;

// Robot parameters
static constexpr double a_0 = 380.0;
static constexpr double d_1 = 43.0;
static constexpr double d_2 = 400.0;
static constexpr double d_3 = 57.0;
static constexpr double post_transform_special_symbol_03 = -57.0;
static constexpr double post_transform_special_symbol_23 = 57.0;
static constexpr double pre_transform_special_symbol_13 = 200.0;

// Unknown offsets from original unknown value to raw value
// Original value are the ones corresponded to robot (usually urdf/sdf)
// Raw value are the ones used in the solver
// unknown_i_raw = unknown_i_original + unknown_i_offset_original2raw
static constexpr double th_0_offset_original2raw = 0.0;
static constexpr double th_1_offset_original2raw = -1.5707963267948966;
static constexpr double th_2_offset_original2raw = -1.5707963267948966;
static constexpr double th_3_offset_original2raw = 3.141592653589793;
static constexpr double th_4_offset_original2raw = 3.141592653589793;
static constexpr double th_5_offset_original2raw = 1.5707963267948966;

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
    ee_transformed(0, 0) = -1.0*r_31;
    ee_transformed(0, 1) = 1.0*r_33;
    ee_transformed(0, 2) = 1.0*r_32;
    ee_transformed(0, 3) = 1.0*Pz + 1.0*post_transform_special_symbol_03*r_31 - 1.0*post_transform_special_symbol_23*r_32;
    ee_transformed(1, 0) = -1.0*r_11;
    ee_transformed(1, 1) = 1.0*r_13;
    ee_transformed(1, 2) = 1.0*r_12;
    ee_transformed(1, 3) = 1.0*Px + 1.0*post_transform_special_symbol_03*r_11 - 1.0*post_transform_special_symbol_23*r_12;
    ee_transformed(2, 0) = -1.0*r_21;
    ee_transformed(2, 1) = 1.0*r_23;
    ee_transformed(2, 2) = 1.0*r_22;
    ee_transformed(2, 3) = 1.0*Py + 1.0*post_transform_special_symbol_03*r_21 - 1.0*post_transform_special_symbol_23*r_22 - 1.0*pre_transform_special_symbol_13;
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
    ee_transformed(0, 0) = -1.0*r_21;
    ee_transformed(0, 1) = 1.0*r_23;
    ee_transformed(0, 2) = 1.0*r_22;
    ee_transformed(0, 3) = 1.0*Py + 1.0*post_transform_special_symbol_03*r_21 + 1.0*post_transform_special_symbol_23*r_23;
    ee_transformed(1, 0) = -1.0*r_31;
    ee_transformed(1, 1) = 1.0*r_33;
    ee_transformed(1, 2) = 1.0*r_32;
    ee_transformed(1, 3) = 1.0*Pz + 1.0*post_transform_special_symbol_03*r_31 + 1.0*post_transform_special_symbol_23*r_33 + 1.0*pre_transform_special_symbol_13;
    ee_transformed(2, 0) = -1.0*r_11;
    ee_transformed(2, 1) = 1.0*r_13;
    ee_transformed(2, 2) = 1.0*r_12;
    ee_transformed(2, 3) = 1.0*Px + 1.0*post_transform_special_symbol_03*r_11 + 1.0*post_transform_special_symbol_23*r_13;
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
    const double x20 = a_0*x8;
    const double x21 = -x1*x7 + x1*x8*x9;
    const double x22 = -x2*x4 - x21*x3;
    const double x23 = -x1*x14 - x1*x15;
    const double x24 = x2*x21 - x3*x4;
    const double x25 = -x13*x23 + x17*x24;
    const double x26 = -x14 - x15;
    const double x27 = x26*x3;
    const double x28 = x7 - x8*x9;
    const double x29 = -x13*x28 + x17*x2*x26;
    // End of temp variables
    Eigen::Matrix4d ee_pose_raw;
    ee_pose_raw.setIdentity();
    ee_pose_raw(0, 0) = -x0*x11 + x12*x19;
    ee_pose_raw(0, 1) = -x0*x19 - x11*x12;
    ee_pose_raw(0, 2) = -x13*x18 - x16*x17;
    ee_pose_raw(0, 3) = -d_1*x1 + d_2*x16 + d_3*x11 + x20*x4;
    ee_pose_raw(1, 0) = -x0*x22 + x12*x25;
    ee_pose_raw(1, 1) = -x0*x25 - x12*x22;
    ee_pose_raw(1, 2) = -x13*x24 - x17*x23;
    ee_pose_raw(1, 3) = d_1*x4 + d_2*x23 + d_3*x22 + x1*x20;
    ee_pose_raw(2, 0) = x0*x27 + x12*x29;
    ee_pose_raw(2, 1) = -x0*x29 + x12*x26*x3;
    ee_pose_raw(2, 2) = -x13*x2*x26 - x17*x28;
    ee_pose_raw(2, 3) = -a_0*x5 + d_2*x28 - d_3*x27;
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
    const double x37 = d_2*x20 - d_3*x22 - x32;
    const double x38 = d_3*x29 + x36;
    const double x39 = a_0*x9 + d_1*x1;
    const double x40 = d_2*x10 + x39;
    const double x41 = d_3*x14 + x40;
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
    jacobian(3, 4) = x22*x38 + x29*x37;
    jacobian(3, 5) = -x23*x38 + x30*x37;
    jacobian(4, 2) = x1*x35 + x39*x5;
    jacobian(4, 3) = x10*x36 - x27*x40;
    jacobian(4, 4) = x14*x38 - x29*x41;
    jacobian(4, 5) = x17*x38 - x30*x41;
    jacobian(5, 1) = -pre_transform_special_symbol_13*x1;
    jacobian(5, 2) = -x1*x33;
    jacobian(5, 3) = -x10*x34 + x20*x40;
    jacobian(5, 4) = -x14*x37 - x22*x41;
    jacobian(5, 5) = -x17*x37 + x23*x41;
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
    const double x28 = d_2*x14 - d_3*x24 - x7;
    const double x29 = d_3*x27 + x21;
    const double x30 = std::cos(th_4);
    const double x31 = std::sin(th_4);
    const double x32 = -x14*x30 - x23*x25*x31;
    const double x33 = -x18*x30 - x31*(x2*x22 + x25*x26);
    const double x34 = 1.0*p_on_ee_x;
    const double x35 = x0*x15 + x1*x34;
    const double x36 = 1.0*x15;
    const double x37 = x12*x2;
    const double x38 = a_0*x37 + d_1*x36;
    const double x39 = x2*x4;
    const double x40 = -x11*x39 - x37*x9;
    const double x41 = d_2*x40 + x38;
    const double x42 = 1.0*x1*x11*x12 - x39*x9;
    const double x43 = -x22*x42 - x25*x36;
    const double x44 = d_3*x43 + x41;
    const double x45 = -x30*x40 - x31*(-x22*x36 + x25*x42);
    const double x46 = p_on_ee_y*x36;
    // End of temp variables
    
    jacobian.setZero();
    jacobian(0, 0) = x0;
    jacobian(0, 1) = -pre_transform_special_symbol_13*x2 + x3;
    jacobian(0, 2) = -x2*x8 + x3;
    jacobian(0, 3) = -p_on_ee_y*x18 + p_on_ee_z*x14 - x14*x21 + x18*x19;
    jacobian(0, 4) = -p_on_ee_y*x27 - p_on_ee_z*x24 + x22*x23*x29 + x27*x28;
    jacobian(0, 5) = -p_on_ee_y*x33 + p_on_ee_z*x32 + x28*x33 - x29*x32;
    jacobian(1, 1) = -x35;
    jacobian(1, 2) = 1.0*x1*x38 + 1.0*x15*x20 - x35;
    jacobian(1, 3) = p_on_ee_x*x18 - p_on_ee_z*x40 - x18*x41 + x21*x40;
    jacobian(1, 4) = p_on_ee_x*x27 - p_on_ee_z*x43 - x27*x44 + x29*x43;
    jacobian(1, 5) = p_on_ee_x*x33 - p_on_ee_z*x45 + x29*x45 - x33*x44;
    jacobian(2, 0) = -x34;
    jacobian(2, 1) = -pre_transform_special_symbol_13*x36 + x46;
    jacobian(2, 2) = -x36*x8 + x46;
    jacobian(2, 3) = -p_on_ee_x*x14 + p_on_ee_y*x40 + x14*x41 - x19*x40;
    jacobian(2, 4) = p_on_ee_x*x24 + p_on_ee_y*x43 - x24*x44 - x28*x43;
    jacobian(2, 5) = -p_on_ee_x*x32 + p_on_ee_y*x45 - x28*x45 + x32*x44;
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
        R_l(0, 4) = Py;
        R_l(0, 5) = Px;
        R_l(1, 0) = d_3*r_11;
        R_l(1, 1) = d_3*r_12;
        R_l(1, 2) = -d_3*r_21;
        R_l(1, 3) = -d_3*r_22;
        R_l(1, 4) = Px;
        R_l(1, 5) = -Py;
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
        R_l(5, 6) = 2*Px*d_3*r_11 + 2*Py*d_3*r_21 + 2*Pz*d_3*r_31;
        R_l(5, 7) = 2*Px*d_3*r_12 + 2*Py*d_3*r_22 + 2*Pz*d_3*r_32;
        R_l(6, 0) = -Px*r_31 + Pz*r_11;
        R_l(6, 1) = -Px*r_32 + Pz*r_12;
        R_l(6, 2) = Py*r_31 - Pz*r_21;
        R_l(6, 3) = Py*r_32 - Pz*r_22;
        R_l(7, 0) = Py*r_31 - Pz*r_21;
        R_l(7, 1) = Py*r_32 - Pz*r_22;
        R_l(7, 2) = Px*r_31 - Pz*r_11;
        R_l(7, 3) = Px*r_32 - Pz*r_12;
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
        const double x0 = R_l_inv_60*r_31 + R_l_inv_70*r_32;
        const double x1 = a_0*x0;
        const double x2 = R_l_inv_62*r_31 + R_l_inv_72*r_32;
        const double x3 = d_2*x2;
        const double x4 = x1 + x3;
        const double x5 = Pz*x2;
        const double x6 = d_1*(R_l_inv_61*r_31 + R_l_inv_71*r_32);
        const double x7 = R_l_inv_66*r_31;
        const double x8 = R_l_inv_76*r_32;
        const double x9 = x7 + x8;
        const double x10 = d_2*x9;
        const double x11 = R_l_inv_65*r_31 + R_l_inv_75*r_32;
        const double x12 = std::pow(Px, 2);
        const double x13 = std::pow(Py, 2);
        const double x14 = std::pow(Pz, 2);
        const double x15 = std::pow(d_1, 2);
        const double x16 = -x15;
        const double x17 = std::pow(d_3, 2);
        const double x18 = std::pow(r_11, 2);
        const double x19 = x17*x18;
        const double x20 = std::pow(r_21, 2);
        const double x21 = x17*x20;
        const double x22 = std::pow(r_31, 2);
        const double x23 = x17*x22;
        const double x24 = std::pow(a_0, 2);
        const double x25 = std::pow(d_2, 2);
        const double x26 = -x24 - x25;
        const double x27 = x12 + x13 + x14 + x16 + x19 + x21 + x23 + x26;
        const double x28 = x11*x27;
        const double x29 = x10 + x28 + x5 + x6;
        const double x30 = R_l_inv_64*r_31;
        const double x31 = R_l_inv_74*r_32;
        const double x32 = x30 + x31;
        const double x33 = x29 + x32;
        const double x34 = 2*d_2;
        const double x35 = x0*x34;
        const double x36 = 2*a_0;
        const double x37 = x2*x36;
        const double x38 = x36*x9;
        const double x39 = x37 + x38;
        const double x40 = -x1;
        const double x41 = -x3;
        const double x42 = x40 + x41;
        const double x43 = -x10 + x28 + x5 + x6;
        const double x44 = x32 + x43;
        const double x45 = R_l_inv_67*r_31 + R_l_inv_77*r_32;
        const double x46 = 2*R_l_inv_63;
        const double x47 = 2*R_l_inv_73;
        const double x48 = r_31*x46 + r_32*x47;
        const double x49 = 2*d_2*x45 - x48;
        const double x50 = 4*d_1;
        const double x51 = x50*x9;
        const double x52 = x51 - 4;
        const double x53 = x34*x45 + x48;
        const double x54 = -x30 - x31;
        const double x55 = x43 + x54;
        const double x56 = x29 + x54;
        const double x57 = Px*r_11;
        const double x58 = Py*r_21;
        const double x59 = Pz*r_31;
        const double x60 = x57 + x58 + x59;
        const double x61 = Px*r_12;
        const double x62 = Py*r_22;
        const double x63 = Pz*r_32;
        const double x64 = x61 + x62 + x63;
        const double x65 = R_l_inv_60*x60 + R_l_inv_70*x64;
        const double x66 = a_0*x65;
        const double x67 = R_l_inv_62*x60 + R_l_inv_72*x64;
        const double x68 = d_2*x67;
        const double x69 = R_l_inv_66*x60;
        const double x70 = R_l_inv_76*x64;
        const double x71 = x69 + x70;
        const double x72 = d_2*x71;
        const double x73 = x66 + x68 + x72;
        const double x74 = R_l_inv_64*x60;
        const double x75 = R_l_inv_74*x64;
        const double x76 = -d_1 + x74 + x75;
        const double x77 = Pz*x67;
        const double x78 = d_1*(R_l_inv_61*x60 + R_l_inv_71*x64);
        const double x79 = R_l_inv_65*x60 + R_l_inv_75*x64;
        const double x80 = x27*x79;
        const double x81 = -d_3*x18 - d_3*x20 - d_3*x22 + x77 + x78 + x80;
        const double x82 = x76 + x81;
        const double x83 = x34*x65;
        const double x84 = x36*x67;
        const double x85 = x36*x71;
        const double x86 = x84 + x85;
        const double x87 = d_1 - x74 - x75;
        const double x88 = d_3*x18;
        const double x89 = d_3*x20;
        const double x90 = d_3*x22;
        const double x91 = -x77 - x78 - x80 + x88 + x89 + x90;
        const double x92 = R_l_inv_67*x60 + R_l_inv_77*x64;
        const double x93 = -2*d_2*x92;
        const double x94 = -x36;
        const double x95 = x46*x60 + x47*x64;
        const double x96 = x94 + x95;
        const double x97 = x50*x71;
        const double x98 = x34*x92;
        const double x99 = x36 + x95;
        const double x100 = x66 - x72;
        const double x101 = x100 + x68;
        const double x102 = x81 + x87;
        const double x103 = x76 + x91;
        const double x104 = Px*r_21 - Py*r_11;
        const double x105 = Px*r_22 - Py*r_12;
        const double x106 = R_l_inv_62*x104 + R_l_inv_72*x105;
        const double x107 = Pz*x106;
        const double x108 = R_l_inv_64*x104;
        const double x109 = R_l_inv_74*x105;
        const double x110 = d_1*(R_l_inv_61*x104 + R_l_inv_71*x105);
        const double x111 = R_l_inv_65*x104 + R_l_inv_75*x105;
        const double x112 = x111*x27;
        const double x113 = x107 + x108 + x109 + x110 + x112;
        const double x114 = R_l_inv_60*x104 + R_l_inv_70*x105;
        const double x115 = a_0*x114;
        const double x116 = d_2*x106;
        const double x117 = x115 + x116;
        const double x118 = R_l_inv_66*x104;
        const double x119 = R_l_inv_76*x105;
        const double x120 = x118 + x119;
        const double x121 = d_2*x120;
        const double x122 = a_0 + x121;
        const double x123 = -x34;
        const double x124 = x120*x36;
        const double x125 = x114*x34;
        const double x126 = x106*x36;
        const double x127 = -x125 + x126;
        const double x128 = -a_0;
        const double x129 = -x121;
        const double x130 = x128 + x129;
        const double x131 = -x115;
        const double x132 = -x116;
        const double x133 = x131 + x132;
        const double x134 = R_l_inv_67*x104 + R_l_inv_77*x105;
        const double x135 = x134*x34;
        const double x136 = 2*d_1;
        const double x137 = -x136;
        const double x138 = x104*x46 + x105*x47 + x137;
        const double x139 = -x135 + x138;
        const double x140 = x120*x50;
        const double x141 = -x140;
        const double x142 = -x135 - x138;
        const double x143 = x107 - x108 - x109 + x110 + x112;
        const double x144 = -x124;
        const double x145 = 2*Py;
        const double x146 = std::pow(r_21, 3)*x17 - r_21*x12 + r_21*x13 - r_21*x14 + r_21*x19 + r_21*x23 + x145*x57 + x145*x59;
        const double x147 = R_l_inv_04*x146;
        const double x148 = 2*r_21*x17;
        const double x149 = r_11*x148;
        const double x150 = r_31*r_32;
        const double x151 = r_12*x149 - r_22*x12 + r_22*x13 - r_22*x14 + r_22*x19 + 3*r_22*x21 + r_22*x23 + x145*x61 + x145*x63 + x148*x150;
        const double x152 = R_l_inv_14*x151;
        const double x153 = 2*Px;
        const double x154 = std::pow(r_11, 3)*x17 + r_11*x12 - r_11*x13 - r_11*x14 + r_11*x21 + r_11*x23 + x153*x58 + x153*x59;
        const double x155 = R_l_inv_24*x154;
        const double x156 = 2*r_11*x17;
        const double x157 = r_12*x12 - r_12*x13 - r_12*x14 + 3*r_12*x19 + r_12*x21 + r_12*x23 + r_22*x149 + x150*x156 + x153*x62 + x153*x63;
        const double x158 = R_l_inv_34*x157;
        const double x159 = x18 + x20 + x22;
        const double x160 = d_3*x159;
        const double x161 = x153*x160;
        const double x162 = R_l_inv_56*x161;
        const double x163 = x145*x160;
        const double x164 = R_l_inv_46*x163;
        const double x165 = R_l_inv_06*x146;
        const double x166 = R_l_inv_26*x154;
        const double x167 = R_l_inv_16*x151;
        const double x168 = R_l_inv_36*x157;
        const double x169 = x162 + x164 + x165 + x166 + x167 + x168;
        const double x170 = d_2*x169;
        const double x171 = d_1*x36;
        const double x172 = R_l_inv_54*x161;
        const double x173 = R_l_inv_44*x163;
        const double x174 = x147 + x152 + x155 + x158 + x170 + x171 + x172 + x173;
        const double x175 = R_l_inv_00*x146 + R_l_inv_10*x151 + R_l_inv_20*x154 + R_l_inv_30*x157 + R_l_inv_40*x163 + R_l_inv_50*x161;
        const double x176 = a_0*x175;
        const double x177 = R_l_inv_02*x146 + R_l_inv_12*x151 + R_l_inv_22*x154 + R_l_inv_32*x157 + R_l_inv_42*x163 + R_l_inv_52*x161;
        const double x178 = d_2*x177;
        const double x179 = x176 + x178;
        const double x180 = R_l_inv_01*x146 + R_l_inv_11*x151 + R_l_inv_21*x154 + R_l_inv_31*x157 + R_l_inv_41*x163 + R_l_inv_51*x161;
        const double x181 = R_l_inv_05*x146 + R_l_inv_15*x151 + R_l_inv_25*x154 + R_l_inv_35*x157 + R_l_inv_45*x163 + R_l_inv_55*x161;
        const double x182 = Pz*x177 + d_1*x180 + x181*x27;
        const double x183 = d_2*x50;
        const double x184 = -x183;
        const double x185 = x169*x36;
        const double x186 = x175*x34;
        const double x187 = x177*x36;
        const double x188 = -x186 + x187;
        const double x189 = -x176;
        const double x190 = -x178;
        const double x191 = x189 + x190;
        const double x192 = -x170;
        const double x193 = -x171;
        const double x194 = x147 + x152 + x155 + x158 + x172 + x173 + x192 + x193;
        const double x195 = R_l_inv_07*x146 + R_l_inv_17*x151 + R_l_inv_27*x154 + R_l_inv_37*x157 + R_l_inv_47*x163 + R_l_inv_57*x161;
        const double x196 = x195*x34;
        const double x197 = 2*x146;
        const double x198 = 2*x151;
        const double x199 = 2*x154;
        const double x200 = 2*x157;
        const double x201 = 4*x160;
        const double x202 = R_l_inv_53*x201;
        const double x203 = R_l_inv_43*x201;
        const double x204 = Px*x202 + Py*x203 + R_l_inv_03*x197 + R_l_inv_13*x198 + R_l_inv_23*x199 + R_l_inv_33*x200 - 2*x15 + 2*x24 - 2*x25;
        const double x205 = -x196 + x204;
        const double x206 = a_0*d_2;
        const double x207 = 8*x206;
        const double x208 = x169*x50;
        const double x209 = -x196 - x204;
        const double x210 = -Pz*x177 - d_1*x180 - x181*x27;
        const double x211 = -x185;
        const double x212 = R_l_inv_40*x161;
        const double x213 = R_l_inv_50*x163;
        const double x214 = R_l_inv_00*x154;
        const double x215 = R_l_inv_20*x146;
        const double x216 = R_l_inv_10*x157;
        const double x217 = R_l_inv_30*x151;
        const double x218 = x212 - x213 + x214 - x215 + x216 - x217;
        const double x219 = a_0*x218;
        const double x220 = R_l_inv_02*x154 + R_l_inv_12*x157 - R_l_inv_22*x146 - R_l_inv_32*x151 + R_l_inv_42*x161 - R_l_inv_52*x163;
        const double x221 = d_2*x220;
        const double x222 = R_l_inv_46*x161;
        const double x223 = R_l_inv_56*x163;
        const double x224 = R_l_inv_06*x154;
        const double x225 = R_l_inv_26*x146;
        const double x226 = R_l_inv_16*x157;
        const double x227 = R_l_inv_36*x151;
        const double x228 = x222 - x223 + x224 - x225 + x226 - x227;
        const double x229 = d_2*x228;
        const double x230 = x219 + x221 + x229;
        const double x231 = R_l_inv_04*x154;
        const double x232 = R_l_inv_14*x157;
        const double x233 = R_l_inv_24*x146;
        const double x234 = R_l_inv_34*x151;
        const double x235 = R_l_inv_54*x163;
        const double x236 = R_l_inv_44*x161;
        const double x237 = R_l_inv_01*x154 + R_l_inv_11*x157 - R_l_inv_21*x146 - R_l_inv_31*x151 + R_l_inv_41*x161 - R_l_inv_51*x163;
        const double x238 = R_l_inv_05*x154 + R_l_inv_15*x157 - R_l_inv_25*x146 - R_l_inv_35*x151 + R_l_inv_45*x161 - R_l_inv_55*x163;
        const double x239 = Pz*x220 + d_1*x237 + x238*x27;
        const double x240 = x15 + x231 + x232 - x233 - x234 - x235 + x236 + x239 + x26;
        const double x241 = x218*x34;
        const double x242 = x220*x36;
        const double x243 = x228*x36;
        const double x244 = x242 + x243;
        const double x245 = x16 - x231 - x232 + x233 + x234 + x235 - x236 + x24 + x25;
        const double x246 = R_l_inv_07*x154 + R_l_inv_17*x157 - R_l_inv_27*x146 - R_l_inv_37*x151 + R_l_inv_47*x161 - R_l_inv_57*x163;
        const double x247 = x246*x34;
        const double x248 = a_0*x50;
        const double x249 = R_l_inv_23*x197;
        const double x250 = R_l_inv_33*x198;
        const double x251 = R_l_inv_03*x199;
        const double x252 = R_l_inv_13*x200;
        const double x253 = Py*x202;
        const double x254 = Px*x203;
        const double x255 = -x249 - x250 + x251 + x252 - x253 + x254;
        const double x256 = x248 + x255;
        const double x257 = x228*x50;
        const double x258 = -x257;
        const double x259 = -x229;
        const double x260 = x219 + x259;
        const double x261 = x239 + x245;
        const double x262 = x221 + x261;
        const double x263 = -x219;
        const double x264 = x229 + x263;
        const double x265 = -x221;
        const double x266 = x261 + x265;
        const double x267 = 2*Pz;
        const double x268 = std::pow(r_31, 3)*x17 - r_31*x12 - r_31*x13 + r_31*x14 + r_31*x19 + r_31*x21 + x267*x57 + x267*x58;
        const double x269 = r_12*r_31*x156 + r_22*r_31*x148 - r_32*x12 - r_32*x13 + r_32*x14 + r_32*x19 + r_32*x21 + 3*r_32*x23 + x267*x61 + x267*x62;
        const double x270 = R_l_inv_60*x268 + R_l_inv_70*x269;
        const double x271 = a_0*x270;
        const double x272 = R_l_inv_62*x268;
        const double x273 = R_l_inv_72*x269;
        const double x274 = x272 + x273;
        const double x275 = d_2*x274;
        const double x276 = R_l_inv_66*x268 + R_l_inv_76*x269;
        const double x277 = d_2*x276;
        const double x278 = d_1*x34;
        const double x279 = x275 + x277 + x278;
        const double x280 = x271 + x279;
        const double x281 = R_l_inv_61*x268 + R_l_inv_71*x269;
        const double x282 = R_l_inv_65*x268 + R_l_inv_75*x269;
        const double x283 = x267*x88;
        const double x284 = x267*x89;
        const double x285 = x267*x90;
        const double x286 = R_l_inv_64*x268;
        const double x287 = R_l_inv_74*x269;
        const double x288 = x286 + x287;
        const double x289 = Pz*x274 + d_1*x281 + x27*x282 - x283 - x284 - x285 + x288;
        const double x290 = x270*x34;
        const double x291 = -x290;
        const double x292 = x248 + x276*x36;
        const double x293 = x274*x36 + x292;
        const double x294 = -Pz*x274 - d_1*x281 - x27*x282 + x283 + x284 + x285;
        const double x295 = -x286 - x287 + x294;
        const double x296 = R_l_inv_67*x268 + R_l_inv_77*x269;
        const double x297 = x296*x34;
        const double x298 = -x297;
        const double x299 = 4*x206;
        const double x300 = x268*x46 + x269*x47;
        const double x301 = x299 + x300;
        const double x302 = 4*x24;
        const double x303 = 4*x15 + 4*x25 + x276*x50;
        const double x304 = -x302 + x303;
        const double x305 = -x271;
        const double x306 = x288 + x294;
        const double x307 = x305 + x306;
        const double x308 = -x275 + x277 + x278;
        const double x309 = -2*a_0*x274 + x292;
        const double x310 = x271 + x306;
        const double x311 = x275 - x277 - x278;
        const double x312 = x11*x36;
        const double x313 = x34*(x0 - x312);
        const double x314 = 4*d_2;
        const double x315 = x34*(-x0 - x312);
        const double x316 = 4*a_0;
        const double x317 = x316*x45;
        const double x318 = 8*R_l_inv_63;
        const double x319 = 8*R_l_inv_73;
        const double x320 = x36*x79;
        const double x321 = x34*(-x320 + x65);
        const double x322 = x34*(-x320 - x65);
        const double x323 = x316*x92;
        const double x324 = -x97;
        const double x325 = -2*a_0*x111;
        const double x326 = x114 + 1;
        const double x327 = x111*x36;
        const double x328 = x134*x316;
        const double x329 = x114 - 1;
        const double x330 = -2*a_0*x181;
        const double x331 = x136 + x175;
        const double x332 = x181*x36;
        const double x333 = x195*x316;
        const double x334 = -x208;
        const double x335 = 16*x160;
        const double x336 = Px*x335;
        const double x337 = 8*R_l_inv_03;
        const double x338 = 8*R_l_inv_13;
        const double x339 = x137 + x175;
        const double x340 = x238*x36;
        const double x341 = x340 + x36;
        const double x342 = x246*x316;
        const double x343 = x34*(2*a_0*x282 - x270);
        const double x344 = x136 + x276;
        const double x345 = x34*(x270 + x282*x36);
        const double x346 = x302 + x303;
        const double x347 = x1 + x41;
        const double x348 = x3 + x40;
        const double x349 = -x68;
        const double x350 = x349 + x66 + x72;
        const double x351 = x115 + x132;
        const double x352 = a_0 + x129;
        const double x353 = x125 + x126;
        const double x354 = x121 + x128;
        const double x355 = x116 + x131;
        const double x356 = x147 + x152 + x155 + x158 + x172 + x173 + x182;
        const double x357 = x171 + x192;
        const double x358 = x176 + x190;
        const double x359 = x186 + x187;
        const double x360 = x178 + x189;
        const double x361 = x170 + x193;
        const double x362 = x147 + x152 + x155 + x158 + x172 + x173 + x210;
        const double x363 = x279 + x305;
        const double x364 = -x299 + x300;
        
        Eigen::Matrix<double, 6, 9> A;
        A.setZero();
        A(0, 0) = -x33 - x4;
        A(0, 1) = x35 - x39;
        A(0, 2) = -x42 - x44;
        A(0, 3) = x49;
        A(0, 4) = x52;
        A(0, 5) = x53;
        A(0, 6) = -x4 - x55;
        A(0, 7) = x35 - x37 + x38;
        A(0, 8) = -x42 - x56;
        A(1, 0) = -x73 - x82;
        A(1, 1) = x83 - x86;
        A(1, 2) = x73 + x87 + x91;
        A(1, 3) = -x93 - x96;
        A(1, 4) = x97;
        A(1, 5) = x98 + x99;
        A(1, 6) = -x101 - x102;
        A(1, 7) = x83 - x84 + x85;
        A(1, 8) = x101 + x103;
        A(2, 0) = x113 + x117 + x122;
        A(2, 1) = x123 + x124 + x127;
        A(2, 2) = x113 + x130 + x133;
        A(2, 3) = x139;
        A(2, 4) = x141;
        A(2, 5) = x142;
        A(2, 6) = x117 + x130 + x143;
        A(2, 7) = x127 + x144 + x34;
        A(2, 8) = x122 + x133 + x143;
        A(3, 0) = x174 + x179 + x182;
        A(3, 1) = x184 + x185 + x188;
        A(3, 2) = x182 + x191 + x194;
        A(3, 3) = x205;
        A(3, 4) = -x207 - x208;
        A(3, 5) = x209;
        A(3, 6) = -x174 - x191 - x210;
        A(3, 7) = x183 + x188 + x211;
        A(3, 8) = -x179 - x194 - x210;
        A(4, 0) = x230 + x240;
        A(4, 1) = -x241 + x244;
        A(4, 2) = Pz*x220 + d_1*x237 - x230 + x238*x27 - x245;
        A(4, 3) = -x247 + x256;
        A(4, 4) = x258;
        A(4, 5) = -x247 + x248 - x255;
        A(4, 6) = x260 + x262;
        A(4, 7) = 2*a_0*x220 - x241 - x243;
        A(4, 8) = x264 + x266;
        A(5, 0) = x280 + x289;
        A(5, 1) = x291 + x293;
        A(5, 2) = -x280 - x295;
        A(5, 3) = x298 + x301;
        A(5, 4) = -x304;
        A(5, 5) = -x297 - x301;
        A(5, 6) = -x307 - x308;
        A(5, 7) = -x290 - x309;
        A(5, 8) = -x310 - x311;
        
        Eigen::Matrix<double, 6, 9> B;
        B.setZero();
        B(0, 0) = x313;
        B(0, 1) = x314*(x2 + x9);
        B(0, 2) = x315;
        B(0, 3) = -x317 + x51 - 4;
        B(0, 4) = r_31*x318 + r_32*x319;
        B(0, 5) = -x317 - x52;
        B(0, 6) = x313;
        B(0, 7) = x314*(x2 - x7 - x8);
        B(0, 8) = x315;
        B(1, 0) = x321;
        B(1, 1) = x314*(x67 + x71);
        B(1, 2) = x322;
        B(1, 3) = -x323 - x324;
        B(1, 4) = x318*x60 + x319*x64;
        B(1, 5) = -x323 - x97;
        B(1, 6) = x321;
        B(1, 7) = x314*(x67 - x69 - x70);
        B(1, 8) = x322;
        B(2, 0) = x34*(-x325 - x326);
        B(2, 1) = x314*(-x106 - x120);
        B(2, 2) = x34*(x326 + x327);
        B(2, 3) = x141 + x328;
        B(2, 4) = 8*d_1 - x104*x318 - x105*x319;
        B(2, 5) = x140 + x328;
        B(2, 6) = x34*(-x325 - x329);
        B(2, 7) = x314*(-x106 + x118 + x119);
        B(2, 8) = x34*(x327 + x329);
        B(3, 0) = x34*(-x330 - x331);
        B(3, 1) = x314*(-x169 - x177);
        B(3, 2) = x34*(x331 + x332);
        B(3, 3) = x333 + x334;
        B(3, 4) = -Py*R_l_inv_43*x335 - 8*R_l_inv_23*x154 - 8*R_l_inv_33*x157 - R_l_inv_53*x336 - x146*x337 + 8*x15 - x151*x338 + 8*x24 + 8*x25;
        B(3, 5) = x208 + x333;
        B(3, 6) = x34*(-x330 - x339);
        B(3, 7) = x314*(x162 + x164 + x165 + x166 + x167 + x168 - x177);
        B(3, 8) = x34*(x332 + x339);
        B(4, 0) = x34*(-x212 + x213 - x214 + x215 - x216 + x217 + x341);
        B(4, 1) = x314*(-x220 - x228);
        B(4, 2) = x34*(x218 + x341);
        B(4, 3) = x258 + x342;
        B(4, 4) = 16*Py*R_l_inv_53*d_3*x159 + 8*R_l_inv_23*x146 + 8*R_l_inv_33*x151 - R_l_inv_43*x336 - x154*x337 - x157*x338;
        B(4, 5) = x257 + x342;
        B(4, 6) = x34*(2*a_0*x238 - x218 - x36);
        B(4, 7) = x314*(-x220 + x222 - x223 + x224 - x225 + x226 - x227);
        B(4, 8) = x34*(x218 + x340 + x94);
        B(5, 0) = x343;
        B(5, 1) = x314*(-x274 - x344);
        B(5, 2) = x345;
        B(5, 3) = 4*a_0*x296 - x346;
        B(5, 4) = -x268*x318 - x269*x319;
        B(5, 5) = x296*x316 + x346;
        B(5, 6) = x343;
        B(5, 7) = x314*(-x272 - x273 + x344);
        B(5, 8) = x345;
        
        Eigen::Matrix<double, 6, 9> C;
        C.setZero();
        C(0, 0) = -x347 - x44;
        C(0, 1) = -x35 - x39;
        C(0, 2) = -x33 - x348;
        C(0, 3) = x53;
        C(0, 4) = -x52;
        C(0, 5) = x49;
        C(0, 6) = -x347 - x56;
        C(0, 7) = -x35 - x37 + x38;
        C(0, 8) = -x348 - x55;
        C(1, 0) = -x100 - x349 - x82;
        C(1, 1) = -x83 - x86;
        C(1, 2) = x66 - x68 - x72 - x82;
        C(1, 3) = x96 + x98;
        C(1, 4) = x324;
        C(1, 5) = -x93 - x99;
        C(1, 6) = -x102 - x350;
        C(1, 7) = -x83 - x84 + x85;
        C(1, 8) = x103 + x350;
        C(2, 0) = x113 + x351 + x352;
        C(2, 1) = x124 + x34 + x353;
        C(2, 2) = x113 + x354 + x355;
        C(2, 3) = x142;
        C(2, 4) = x140;
        C(2, 5) = x139;
        C(2, 6) = x143 + x351 + x354;
        C(2, 7) = x123 + x144 + x353;
        C(2, 8) = x143 + x352 + x355;
        C(3, 0) = x356 + x357 + x358;
        C(3, 1) = x183 + x185 + x359;
        C(3, 2) = x356 + x360 + x361;
        C(3, 3) = x209;
        C(3, 4) = -x207 - x334;
        C(3, 5) = x205;
        C(3, 6) = -x357 - x360 - x362;
        C(3, 7) = x184 + x211 + x359;
        C(3, 8) = -x358 - x361 - x362;
        C(4, 0) = x240 + x260 + x265;
        C(4, 1) = x241 + x244;
        C(4, 2) = x221 + x240 + x264;
        C(4, 3) = -x247 - x256;
        C(4, 4) = x257;
        C(4, 5) = -x247 - x248 - x249 - x250 + x251 + x252 - x253 + x254;
        C(4, 6) = x219 + x229 + x266;
        C(4, 7) = x241 + x242 - x243;
        C(4, 8) = x259 + x262 + x263;
        C(5, 0) = -x295 - x363;
        C(5, 1) = x290 + x293;
        C(5, 2) = x289 + x363;
        C(5, 3) = -x297 - x364;
        C(5, 4) = x304;
        C(5, 5) = x298 + x364;
        C(5, 6) = -x307 - x311;
        C(5, 7) = -x291 - x309;
        C(5, 8) = -x308 - x310;
        
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
            
            const bool condition_0 = 2*std::fabs(d_1*d_3) >= zero_tolerance || 2*std::fabs(a_0*d_3*std::cos(th_2)) >= zero_tolerance || std::fabs(std::pow(Px, 2) + std::pow(Py, 2) + std::pow(Pz, 2) - std::pow(a_0, 2) + 2*a_0*d_2*std::sin(th_2) - std::pow(d_1, 2) - std::pow(d_2, 2) - std::pow(d_3, 2)) >= zero_tolerance;
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
                const double x7 = std::pow(Px, 2) + std::pow(Py, 2) + std::pow(Pz, 2) + 2*a_0*d_2*std::sin(th_2) - std::pow(d_2, 2) - x3 - x4 - x6;
                const double x8 = safe_sqrt(std::pow(x0, 2)*x5*x6 + x3*x5 - std::pow(x7, 2));
                // End of temp variables
                RawSolution solution_0(this_solution);
                solution_0[4] = x2 + std::atan2(x8, x7);
                int appended_idx = append_solution_to_queue(solution_0);
                add_input_index_to(4, appended_idx);
            }
            
            const bool condition_1 = 2*std::fabs(d_1*d_3) >= zero_tolerance || 2*std::fabs(a_0*d_3*std::cos(th_2)) >= zero_tolerance || std::fabs(std::pow(Px, 2) + std::pow(Py, 2) + std::pow(Pz, 2) - std::pow(a_0, 2) + 2*a_0*d_2*std::sin(th_2) - std::pow(d_1, 2) - std::pow(d_2, 2) - std::pow(d_3, 2)) >= zero_tolerance;
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
                const double x7 = std::pow(Px, 2) + std::pow(Py, 2) + std::pow(Pz, 2) + 2*a_0*d_2*std::sin(th_2) - std::pow(d_2, 2) - x3 - x4 - x6;
                const double x8 = safe_sqrt(std::pow(x0, 2)*x5*x6 + x3*x5 - std::pow(x7, 2));
                // End of temp variables
                const double tmp_sol_value = x2 + std::atan2(-x8, x7);
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
            const bool checked_result = std::fabs(Px) <= 9.9999999999999995e-7 && std::fabs(Py) <= 9.9999999999999995e-7 && std::fabs(d_1 - d_3*std::cos(th_3)) <= 9.9999999999999995e-7;
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
            
            const bool condition_0 = std::fabs(Px) >= zero_tolerance || std::fabs(Py) >= zero_tolerance || std::fabs(d_1 - d_3*std::cos(th_3)) >= zero_tolerance;
            if (condition_0)
            {
                // Temp variable for efficiency
                const double x0 = std::atan2(Px, -Py);
                const double x1 = -d_1 + d_3*std::cos(th_3);
                const double x2 = safe_sqrt(std::pow(Px, 2) + std::pow(Py, 2) - std::pow(x1, 2));
                // End of temp variables
                RawSolution solution_0(this_solution);
                solution_0[0] = x0 + std::atan2(x2, x1);
                int appended_idx = append_solution_to_queue(solution_0);
                add_input_index_to(6, appended_idx);
            }
            
            const bool condition_1 = std::fabs(Px) >= zero_tolerance || std::fabs(Py) >= zero_tolerance || std::fabs(d_1 - d_3*std::cos(th_3)) >= zero_tolerance;
            if (condition_1)
            {
                // Temp variable for efficiency
                const double x0 = std::atan2(Px, -Py);
                const double x1 = -d_1 + d_3*std::cos(th_3);
                const double x2 = safe_sqrt(std::pow(Px, 2) + std::pow(Py, 2) - std::pow(x1, 2));
                // End of temp variables
                const double tmp_sol_value = x0 + std::atan2(-x2, x1);
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
                add_input_index_to(24, node_input_i_idx_in_queue);
            }
            
            if (!taken_by_degenerate)
                add_input_index_to(21, node_input_i_idx_in_queue);
        }
    };
    
    // Invoke the processor
    SolvedVariableDispatcherNode_node_20_processor();
    // Finish code for solved_variable dispatcher node 20
    
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
            const double th_3 = this_solution[4];
            
            const bool condition_0 = std::fabs(Pz) >= zero_tolerance || std::fabs(d_3*std::sin(th_3)) >= zero_tolerance || std::fabs(a_0 - d_2) >= zero_tolerance;
            if (condition_0)
            {
                // Temp variable for efficiency
                const double x0 = -a_0 + d_2;
                const double x1 = std::sin(th_3);
                const double x2 = std::atan2(x0, d_3*x1);
                const double x3 = safe_sqrt(-std::pow(Pz, 2) + std::pow(d_3, 2)*std::pow(x1, 2) + std::pow(x0, 2));
                // End of temp variables
                RawSolution solution_0(this_solution);
                solution_0[1] = x2 + std::atan2(x3, Pz);
                int appended_idx = append_solution_to_queue(solution_0);
                add_input_index_to(25, appended_idx);
            }
            
            const bool condition_1 = std::fabs(Pz) >= zero_tolerance || std::fabs(d_3*std::sin(th_3)) >= zero_tolerance || std::fabs(a_0 - d_2) >= zero_tolerance;
            if (condition_1)
            {
                // Temp variable for efficiency
                const double x0 = -a_0 + d_2;
                const double x1 = std::sin(th_3);
                const double x2 = std::atan2(x0, d_3*x1);
                const double x3 = safe_sqrt(-std::pow(Pz, 2) + std::pow(d_3, 2)*std::pow(x1, 2) + std::pow(x0, 2));
                // End of temp variables
                const double tmp_sol_value = x2 + std::atan2(-x3, Pz);
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
    
    // Code for explicit solution node 21, solved variable is th_4
    auto ExplicitSolutionNode_node_21_solve_th_4_processor = [&]() -> void
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
            
            const bool condition_0 = std::fabs(a_0*std::cos(th_2)) >= 9.9999999999999995e-7 || std::fabs(a_0*std::sin(th_2) - d_2) >= 9.9999999999999995e-7;
            if (condition_0)
            {
                // Temp variable for efficiency
                const double x0 = std::cos(th_1th_2th_4_soa);
                const double x1 = std::sin(th_1th_2th_4_soa);
                const double x2 = Px*std::cos(th_0);
                const double x3 = Py*std::sin(th_0);
                const double x4 = Pz*x0 + x1*x2 + x1*x3;
                const double x5 = a_0*std::cos(th_2);
                const double x6 = a_0*std::sin(th_2) - d_2;
                const double x7 = Pz*x1 - x0*x2 - x0*x3;
                // End of temp variables
                const double tmp_sol_value = std::atan2(x4*x5 + x6*x7, x4*x6 - x5*x7);
                solution_queue.get_solution(node_input_i_idx_in_queue)[5] = tmp_sol_value;
                add_input_index_to(22, node_input_i_idx_in_queue);
            }
            else
            {
                solution_queue.solutions_validity[node_input_i_idx_in_queue] = false;
            }
        }
    };
    // Invoke the processor
    ExplicitSolutionNode_node_21_solve_th_4_processor();
    // Finish code for explicit solution node 21
    
    // Code for non-branch dispatcher node 22
    // Actually, there is no code
    
    // Code for explicit solution node 23, solved variable is th_1
    auto ExplicitSolutionNode_node_23_solve_th_1_processor = [&]() -> void
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
            const double th_1th_2th_4_soa = this_solution[2];
            const double th_2 = this_solution[3];
            const double th_4 = this_solution[5];
            
            const bool condition_0 = 1 >= zero_tolerance;
            if (condition_0)
            {
                // Temp variable for efficiency
                // End of temp variables
                const double tmp_sol_value = th_1th_2th_4_soa - th_2 - th_4;
                solution_queue.get_solution(node_input_i_idx_in_queue)[1] = tmp_sol_value;
            }
            else
            {
                solution_queue.solutions_validity[node_input_i_idx_in_queue] = false;
            }
        }
    };
    // Invoke the processor
    ExplicitSolutionNode_node_23_solve_th_1_processor();
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
    
    // Code for equation all-zero dispatcher node 13
    auto EquationAllZeroDispatcherNode_node_13_processor = [&]() -> void
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
            const double th_0 = this_solution[0];
            const bool checked_result = std::fabs(Pz) <= 9.9999999999999995e-7 && std::fabs(Px*std::cos(th_0) + Py*std::sin(th_0)) <= 9.9999999999999995e-7;
            if (!checked_result)  // To non-degenerate node
                add_input_index_to(14, node_input_i_idx_in_queue);
        }
    };
    
    // Invoke the processor
    EquationAllZeroDispatcherNode_node_13_processor();
    // Finish code for equation all-zero dispatcher node 13
    
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
            
            const bool condition_0 = std::fabs(Pz) >= 9.9999999999999995e-7 || std::fabs(Px*std::cos(th_0) + Py*std::sin(th_0)) >= 9.9999999999999995e-7;
            if (condition_0)
            {
                // Temp variable for efficiency
                const double x0 = -a_0 + d_2*std::sin(th_2);
                const double x1 = -Px*std::cos(th_0) - Py*std::sin(th_0);
                const double x2 = d_2*std::cos(th_2);
                // End of temp variables
                const double tmp_sol_value = std::atan2(Pz*x0 + x1*x2, -Pz*x2 + x0*x1);
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
            const double th_0 = this_solution[0];
            const bool checked_result = std::fabs(Pz) <= 9.9999999999999995e-7 && std::fabs(Px*std::cos(th_0) + Py*std::sin(th_0)) <= 9.9999999999999995e-7;
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
            
            const bool condition_0 = std::fabs(Pz) >= 9.9999999999999995e-7 || std::fabs(Px*std::cos(th_0) + Py*std::sin(th_0)) >= 9.9999999999999995e-7;
            if (condition_0)
            {
                // Temp variable for efficiency
                const double x0 = std::sin(th_2);
                const double x1 = std::cos(th_2);
                const double x2 = std::sin(th_3);
                const double x3 = -a_0 + d_2*x0 + d_3*x1*x2;
                const double x4 = -Px*std::cos(th_0) - Py*std::sin(th_0);
                const double x5 = -d_2*x1 + d_3*x0*x2;
                // End of temp variables
                const double tmp_sol_value = std::atan2(Pz*x3 - x4*x5, Pz*x5 + x3*x4);
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

}; // struct articulating_arm_ik

// Code below for debug
void test_ik_solve_articulating_arm()
{
    std::array<double, articulating_arm_ik::robot_nq> theta;
    std::random_device rd;
    std::uniform_real_distribution<double> distribution;
    for(auto i = 0; i < theta.size(); i++)
        theta[i] = distribution(rd);
    const Eigen::Matrix4d ee_pose = articulating_arm_ik::computeFK(theta);
    auto ik_output = articulating_arm_ik::computeIK(ee_pose);
    for(int i = 0; i < ik_output.size(); i++)
    {
        Eigen::Matrix4d ee_pose_i = articulating_arm_ik::computeFK(ik_output[i]);
        double ee_pose_diff = (ee_pose_i - ee_pose).norm();
        std::cout << "For solution " << i << " Pose different with ground-truth " << ee_pose_diff << std::endl;
    }
}

int main()
{
    test_ik_solve_articulating_arm();
}
