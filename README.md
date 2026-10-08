# yaik

A Python package for generating analytical inverse kinematic solvers for robotics manipulators. The generated solvers support many robots in [Mech-Mind Robotics](https://www.mech-mind.com/) and have been extensively tested in real-world deployment. The name **yaik** stands for **Y**et **A**nother **A**nalytical **I**nverse **K**inematics generator.

### Features

- Closed-form, analytical inverse kinematics for serial chain robots in python/C++/Rust
- Support robots with 7 DoF or more, such as [Kuka IIwa Arm](./gallery/kuka_iiwa), [Panda Robot](./gallery/franka_panda), or arms of the Atlas humanoid
- A new factorized polynomial solution based on *groebner reduce*
- General 6-dof solution using linear matrix equations tangent of half angles
- Analysis to avoid degenerated (or conditional degenerated) solutions

### Using the generated solver

The [gallery](./gallery) contains generated IK solvers for tens of robots. Even if your robot is not included in the library, it might include a robot with the same structure and can be used after parameter modification. The file format is

```
gallery/
    puma/
        puma_ik.yaml
        puma_ik_generated.cpp
        puma_ik_generated.py
    {robot_name}/
        {robot_name}_ik.yaml
        {robot_name}_ik_generated.cpp
        {robot_name}_ik_generated.py
```

To use the generated C++ solver, you need to include the [yaik_cpp_common.h](./yaik_cpp_common.h) file, which further depends on Eigen. Please refer to the [yaik_cpp](https://bitbucket.org/gaowei19951004/yaik_cpp/src/master/) repo for an example.

To use the generated Python solver, you need to import [python_run_import.py](./python_run_import.py), which further depends on numpy.

To use a generated Rust solver (`{robot_name}_ik_generated.rs`), copy it together with [yaik_rust_common.rs](./yaik_rust_common.rs) into your crate and declare both as sibling modules; the solver imports the common file as `super::yaik_rust_common`. The host crate needs these dependencies:

```toml
[dependencies]
nalgebra = "0.35"
rand = "0.9"
glam = "0.34"
```

The solver uses nalgebra internally; its public `compute_ik` / `compute_fk` functions accept and return any type implementing `Pose4`, which covers `nalgebra::Matrix4<f64>`, row-major `[[f64; 4]; 4]`, `glam::DMat4` and `glam::Mat4`:

```rust
mod yaik_rust_common;
mod puma_ik_generated;
use puma_ik_generated as puma;

let pose: glam::DMat4 = puma::compute_fk(&[0.3, -0.7, 1.1, 0.4, -0.9, 0.25]);
let solutions: Vec<[f64; 6]> = puma::compute_ik(&pose);
```

For allocation-free repeated solves, keep a `puma::IkWorkspace` and call `compute_ik_with_workspace`. Each generated file also contains a `#[test]` that checks IK against FK on random joint values.

### Generate a solver for your robot

To generate a new solver for a robot that is not in the gallery, you need to install python3, [pyyaml](https://github.com/yaml/pyyaml), [attrs](https://www.attrs.org/en/stable/index.html) and [sympy](https://github.com/sympy/sympy/wiki/Download-Installation) on your system. Optionally, you might install [sagemath](https://www.sagemath.org/index.html), which contains a faster groebner reduce implementation.

Then, you need to express your robot in **Modified** DH format. A documentation regarding that can be found in [here](./docs/modified_dh.md). Several examples can be found in [this file](./fk/robot_models.py). Alternatively, the modified dh model can be converted from chain model in [this file](./fk/chain_models.py), or URDF using a tool in [another library](https://bitbucket.org/gaowei19951004/multibody/src/main/examples/urdf_sdf_to_dh.cpp).

After that, you can run `ik_solve.py` to generate the solution for your robot in YAML format. This step can take a LOT of time, as we perform offline analysis to find the simplest solution and avoid degeneration. Moreover, this step might fail to find a feasible solution (although the pipeline succeeded for all the robots we tested).

Finally, you can use `cpp_codegen.py/py_codegen.py/rust_codegen.py` to generate solvers in C++, Python and Rust, similar to the ones in the gallery.

### Acknowledgement

This work was done during the author's employment at Mech-Mind Robotics. We have learned a lot from existing works, in particular [IkFast](https://github.com/rdiankov/openrave) and [IKBT](https://github.com/uw-biorobotics/IKBT).

### Contact

If you have a new challenging robot that cannot be solved, please email me.
