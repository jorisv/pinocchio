# Define Pinocchio sources and headers

SET(${PROJECT_NAME}_CORE_SOURCES
  src/utils/file-explorer.cpp
  )

SET(${PROJECT_NAME}_CORE_PUBLIC_HEADERS
  include/pinocchio/algorithm/aba-derivatives.hpp
  include/pinocchio/algorithm/aba-derivatives.hxx
  include/pinocchio/algorithm/aba.hpp
  include/pinocchio/algorithm/aba.hxx
  include/pinocchio/algorithm/admm-solver.hpp
  include/pinocchio/algorithm/admm-solver.hxx
  include/pinocchio/algorithm/broadphase-callbacks.hpp
  include/pinocchio/algorithm/broadphase-callbacks.hxx
  include/pinocchio/algorithm/broadphase.hpp
  include/pinocchio/algorithm/broadphase.hxx
  include/pinocchio/algorithm/center-of-mass-derivatives.hpp
  include/pinocchio/algorithm/center-of-mass-derivatives.hxx
  include/pinocchio/algorithm/center-of-mass.hpp
  include/pinocchio/algorithm/center-of-mass.hxx
  include/pinocchio/algorithm/centroidal-derivatives.hpp
  include/pinocchio/algorithm/centroidal-derivatives.hxx
  include/pinocchio/algorithm/centroidal.hpp
  include/pinocchio/algorithm/centroidal.hxx
  include/pinocchio/algorithm/check.hpp
  include/pinocchio/algorithm/check.hxx
  include/pinocchio/algorithm/cholesky.hpp
  include/pinocchio/algorithm/cholesky.hxx
  include/pinocchio/algorithm/collision.hpp
  include/pinocchio/algorithm/collision.hxx
  include/pinocchio/algorithm/compute-all-terms.hpp
  include/pinocchio/algorithm/compute-all-terms.hxx
  include/pinocchio/algorithm/constrained-dynamics-derivatives.hpp
  include/pinocchio/algorithm/constrained-dynamics-derivatives.hxx
  include/pinocchio/algorithm/constrained-dynamics.hpp
  include/pinocchio/algorithm/constrained-dynamics.hxx
  include/pinocchio/algorithm/constrained-problem-data.hpp
  include/pinocchio/algorithm/constraints/coulomb-friction-cone.hpp
  include/pinocchio/algorithm/constraints/fwd.hpp
  include/pinocchio/algorithm/contact-cholesky.hpp
  include/pinocchio/algorithm/contact-cholesky.hxx
  include/pinocchio/algorithm/contact-dynamics.hpp
  include/pinocchio/algorithm/contact-dynamics.hxx
  include/pinocchio/algorithm/contact.hpp
  include/pinocchio/algorithm/contact.hxx
  include/pinocchio/algorithm/contact-info.hpp
  include/pinocchio/algorithm/contact-jacobian.hpp
  include/pinocchio/algorithm/contact-jacobian.hxx
  include/pinocchio/algorithm/contact-solver-base.hpp
  include/pinocchio/algorithm/contact-solver-utils.hpp
  include/pinocchio/algorithm/context.hpp
  include/pinocchio/algorithm/copy.hpp
  include/pinocchio/algorithm/crba.hpp
  include/pinocchio/algorithm/crba.hxx
  include/pinocchio/algorithm/default-check.hpp
  include/pinocchio/algorithm/delassus.hpp
  include/pinocchio/algorithm/delassus.hxx
  include/pinocchio/algorithm/delassus-operartor-ref.hpp
  include/pinocchio/algorithm/delassus-operator-base.hpp
  include/pinocchio/algorithm/delassus-operator-dense.hpp
  include/pinocchio/algorithm/delassus-operator-rigid-body.hpp
  include/pinocchio/algorithm/delassus-operator-sparse.hpp
  include/pinocchio/algorithm/distance.hpp
  include/pinocchio/algorithm/distance.hxx
  include/pinocchio/algorithm/dynamics.hpp
  include/pinocchio/algorithm/energy.hpp
  include/pinocchio/algorithm/energy.hxx
  include/pinocchio/algorithm/frames-derivatives.hpp
  include/pinocchio/algorithm/frames-derivatives.hxx
  include/pinocchio/algorithm/frames.hpp
  include/pinocchio/algorithm/frames.hxx
  include/pinocchio/algorithm/fwd.hpp
  include/pinocchio/algorithm/geometry.hpp
  include/pinocchio/algorithm/geometry.hxx
  include/pinocchio/algorithm/impulse-dynamics-derivatives.hpp
  include/pinocchio/algorithm/impulse-dynamics-derivatives.hxx
  include/pinocchio/algorithm/impulse-dynamics.hpp
  include/pinocchio/algorithm/impulse-dynamics.hxx
  include/pinocchio/algorithm/jacobian.hpp
  include/pinocchio/algorithm/jacobian.hxx
  include/pinocchio/algorithm/joint-configuration.hpp
  include/pinocchio/algorithm/joint-configuration.hxx
  include/pinocchio/algorithm/kinematics-derivatives.hpp
  include/pinocchio/algorithm/kinematics-derivatives.hxx
  include/pinocchio/algorithm/kinematics.hpp
  include/pinocchio/algorithm/kinematics.hxx
  include/pinocchio/algorithm/model.hpp
  include/pinocchio/algorithm/model.hxx
  include/pinocchio/algorithm/parallel/aba.hpp
  include/pinocchio/algorithm/parallel/broadphase.hpp
  include/pinocchio/algorithm/parallel/geometry.hpp
  include/pinocchio/algorithm/parallel/omp.hpp
  include/pinocchio/algorithm/parallel/rnea.hpp
  include/pinocchio/algorithm/pgs-solver.hpp
  include/pinocchio/algorithm/pgs-solver.hxx
  include/pinocchio/algorithm/proximal.hpp
  include/pinocchio/algorithm/regressor.hpp
  include/pinocchio/algorithm/regressor.hxx
  include/pinocchio/algorithm/rnea-derivatives.hpp
  include/pinocchio/algorithm/rnea-derivatives.hxx
  include/pinocchio/algorithm/rnea.hpp
  include/pinocchio/algorithm/rnea.hxx
  include/pinocchio/algorithm/rnea-second-order-derivatives.hpp
  include/pinocchio/algorithm/rnea-second-order-derivatives.hxx
  include/pinocchio/algorithm/utils/force.hpp
  include/pinocchio/algorithm/utils/motion.hpp
  include/pinocchio/autodiff/casadi-algo.hpp
  include/pinocchio/autodiff/casadi.hpp
  include/pinocchio/autodiff/casadi/math/matrix.hpp
  include/pinocchio/autodiff/casadi/math/quaternion.hpp
  include/pinocchio/autodiff/casadi/math/triangular-matrix.hpp
  include/pinocchio/autodiff/casadi/spatial/se3-tpl.hpp
  include/pinocchio/autodiff/casadi/utils/static-if.hpp
  include/pinocchio/autodiff/cppad/algorithm/aba.hpp
  include/pinocchio/autodiff/cppad.hpp
  include/pinocchio/autodiff/cppad/math/eigen_plugin.hpp
  include/pinocchio/autodiff/cppad/math/quaternion.hpp
  include/pinocchio/autodiff/cppad/spatial/log.hxx
  include/pinocchio/autodiff/cppad/spatial/se3-tpl.hpp
  include/pinocchio/autodiff/cppad/utils/static-if.hpp
  include/pinocchio/codegen/code-generator-algo.hpp
  include/pinocchio/codegen/code-generator-base.hpp
  include/pinocchio/codegen/cppadcg.hpp
  include/pinocchio/container/aligned-vector.hpp
  include/pinocchio/container/boost-container-limits.hpp
  include/pinocchio/context.hpp
  include/pinocchio/core/binary-op.hpp
  include/pinocchio/core/unary-op.hpp
  include/pinocchio/deprecated-macros.hpp
  include/pinocchio/deprecated-namespaces.hpp
  include/pinocchio/deprecation.hpp
  include/pinocchio/eigen-macros.hpp
  include/pinocchio/fwd.hpp
  include/pinocchio/macros.hpp
  include/pinocchio/math/casadi.hpp
  include/pinocchio/math/comparison-operators.hpp
  include/pinocchio/math/cppadcg.hpp
  include/pinocchio/math/cppad.hpp
  include/pinocchio/math/eigenvalues.hpp
  include/pinocchio/math/fwd.hpp
  include/pinocchio/math/matrix-block.hpp
  include/pinocchio/math/matrix.hpp
  include/pinocchio/math/multiprecision.hpp
  include/pinocchio/math/multiprecision-mpfr.hpp
  include/pinocchio/math/quaternion.hpp
  include/pinocchio/math/rotation.hpp
  include/pinocchio/math/rpy.hpp
  include/pinocchio/math/rpy.hxx
  include/pinocchio/math/sign.hpp
  include/pinocchio/math/sincos.hpp
  include/pinocchio/math/taylor-expansion.hpp
  include/pinocchio/math/tensor.hpp
  include/pinocchio/math/triangular-matrix.hpp
  include/pinocchio/multibody/broadphase-manager-base.hpp
  include/pinocchio/multibody/broadphase-manager.hpp
  include/pinocchio/multibody/broadphase-manager.hxx
  include/pinocchio/multibody/data.hpp
  include/pinocchio/multibody/data.hxx
  include/pinocchio/multibody/fcl.hpp
  include/pinocchio/multibody/force-set.hpp
  include/pinocchio/multibody/frame.hpp
  include/pinocchio/multibody/fwd.hpp
  include/pinocchio/multibody/geometry.hpp
  include/pinocchio/multibody/geometry.hxx
  include/pinocchio/multibody/geometry-object-filter.hpp
  include/pinocchio/multibody/geometry-object.hpp
  include/pinocchio/multibody/geometry-object.hxx
  include/pinocchio/multibody/instance-filter.hpp
  include/pinocchio/multibody/joint/fwd.hpp
  include/pinocchio/multibody/joint/joint-base.hpp
  include/pinocchio/multibody/joint/joint-basic-visitors.hpp
  include/pinocchio/multibody/joint/joint-basic-visitors.hxx
  include/pinocchio/multibody/joint/joint-collection.hpp
  include/pinocchio/multibody/joint/joint-common-operations.hpp
  include/pinocchio/multibody/joint/joint-composite.hpp
  include/pinocchio/multibody/joint/joint-composite.hxx
  include/pinocchio/multibody/joint/joint-data-base.hpp
  include/pinocchio/multibody/joint/joint-free-flyer.hpp
  include/pinocchio/multibody/joint/joint-generic.hpp
  include/pinocchio/multibody/joint/joint-helical.hpp
  include/pinocchio/multibody/joint/joint-helical-unaligned.hpp
  include/pinocchio/multibody/joint/joint-mimic.hpp
  include/pinocchio/multibody/joint/joint-model-base.hpp
  include/pinocchio/multibody/joint/joint-planar.hpp
  include/pinocchio/multibody/joint/joint-prismatic.hpp
  include/pinocchio/multibody/joint/joint-prismatic-unaligned.hpp
  include/pinocchio/multibody/joint/joint-revolute.hpp
  include/pinocchio/multibody/joint/joint-revolute-unaligned.hpp
  include/pinocchio/multibody/joint/joint-revolute-unbounded.hpp
  include/pinocchio/multibody/joint/joint-revolute-unbounded-unaligned.hpp
  include/pinocchio/multibody/joint/joints.hpp
  include/pinocchio/multibody/joint/joint-spherical.hpp
  include/pinocchio/multibody/joint/joint-spherical-ZYX.hpp
  include/pinocchio/multibody/joint/joint-translation.hpp
  include/pinocchio/multibody/joint/joint-universal.hpp
  include/pinocchio/multibody/joint-motion-subspace-base.hpp
  include/pinocchio/multibody/joint-motion-subspace-generic.hpp
  include/pinocchio/multibody/joint-motion-subspace.hpp
  include/pinocchio/multibody/liegroup/cartesian-product.hpp
  include/pinocchio/multibody/liegroup/cartesian-product-variant.hpp
  include/pinocchio/multibody/liegroup/cartesian-product-variant.hxx
  include/pinocchio/multibody/liegroup/fwd.hpp
  include/pinocchio/multibody/liegroup/liegroup-algo.hpp
  include/pinocchio/multibody/liegroup/liegroup-algo.hxx
  include/pinocchio/multibody/liegroup/liegroup-base.hpp
  include/pinocchio/multibody/liegroup/liegroup-base.hxx
  include/pinocchio/multibody/liegroup/liegroup-collection.hpp
  include/pinocchio/multibody/liegroup/liegroup-generic.hpp
  include/pinocchio/multibody/liegroup/liegroup.hpp
  include/pinocchio/multibody/liegroup/liegroup-variant-visitors.hpp
  include/pinocchio/multibody/liegroup/liegroup-variant-visitors.hxx
  include/pinocchio/multibody/liegroup/special-euclidean.hpp
  include/pinocchio/multibody/liegroup/special-orthogonal.hpp
  include/pinocchio/multibody/liegroup/vector-space.hpp
  include/pinocchio/multibody/model.hpp
  include/pinocchio/multibody/model.hxx
  include/pinocchio/multibody/model-item.hpp
  include/pinocchio/multibody/pool/broadphase-manager.hpp
  include/pinocchio/multibody/pool/fwd.hpp
  include/pinocchio/multibody/pool/geometry.hpp
  include/pinocchio/multibody/pool/model.hpp
  include/pinocchio/multibody/tree-broadphase-manager.hpp
  include/pinocchio/multibody/tree-broadphase-manager.hxx
  include/pinocchio/multibody/visitor/fusion.hpp
  include/pinocchio/multibody/visitor.hpp
  include/pinocchio/multibody/visitor/joint-binary-visitor.hpp
  include/pinocchio/multibody/visitor/joint-unary-visitor.hpp
  include/pinocchio/parsers/python.hpp
  include/pinocchio/parsers/sample-models.hpp
  include/pinocchio/parsers/sample-models.hxx
  include/pinocchio/parsers/srdf.hpp
  include/pinocchio/parsers/srdf.hxx
  include/pinocchio/parsers/utils.hpp
  include/pinocchio/serialization/aligned-vector.hpp
  include/pinocchio/serialization/archive.hpp
  include/pinocchio/serialization/data.hpp
  include/pinocchio/serialization/eigen.hpp
  include/pinocchio/serialization/force.hpp
  include/pinocchio/serialization/frame.hpp
  include/pinocchio/serialization/fwd.hpp
  include/pinocchio/serialization/geometry.hpp
  include/pinocchio/serialization/inertia.hpp
  include/pinocchio/serialization/joints-data.hpp
  include/pinocchio/serialization/joints.hpp
  include/pinocchio/serialization/joints-model.hpp
  include/pinocchio/serialization/joints-motion.hpp
  include/pinocchio/serialization/joints-motion-subspace.hpp
  include/pinocchio/serialization/joints-transform.hpp
  include/pinocchio/serialization/model.hpp
  include/pinocchio/serialization/motion.hpp
  include/pinocchio/serialization/se3.hpp
  include/pinocchio/serialization/serializable.hpp
  include/pinocchio/serialization/spatial.hpp
  include/pinocchio/serialization/static-buffer.hpp
  include/pinocchio/serialization/symmetric3.hpp
  include/pinocchio/serialization/vector.hpp
  include/pinocchio/spatial/act-on-set.hpp
  include/pinocchio/spatial/act-on-set.hxx
  include/pinocchio/spatial/cartesian-axis.hpp
  include/pinocchio/spatial/classic-acceleration.hpp
  include/pinocchio/spatial/explog.hpp
  include/pinocchio/spatial/explog-quaternion.hpp
  include/pinocchio/spatial/force-base.hpp
  include/pinocchio/spatial/force-dense.hpp
  include/pinocchio/spatial/force.hpp
  include/pinocchio/spatial/force-ref.hpp
  include/pinocchio/spatial/force-tpl.hpp
  include/pinocchio/spatial/fwd.hpp
  include/pinocchio/spatial/inertia.hpp
  include/pinocchio/spatial/log.hpp
  include/pinocchio/spatial/log.hxx
  include/pinocchio/spatial/motion-base.hpp
  include/pinocchio/spatial/motion-dense.hpp
  include/pinocchio/spatial/motion.hpp
  include/pinocchio/spatial/motion-ref.hpp
  include/pinocchio/spatial/motion-tpl.hpp
  include/pinocchio/spatial/motion-zero.hpp
  include/pinocchio/spatial/se3-base.hpp
  include/pinocchio/spatial/se3.hpp
  include/pinocchio/spatial/se3-tpl.hpp
  include/pinocchio/spatial/skew.hpp
  include/pinocchio/spatial/spatial-axis.hpp
  include/pinocchio/spatial/symmetric3.hpp
  include/pinocchio/utils/axis-label.hpp
  include/pinocchio/utils/cast.hpp
  include/pinocchio/utils/check.hpp
  include/pinocchio/utils/eigen-fix.hpp
  include/pinocchio/utils/file-explorer.hpp
  include/pinocchio/utils/file-io.hpp
  include/pinocchio/utils/helpers.hpp
  include/pinocchio/utils/openmp.hpp
  include/pinocchio/utils/shared-ptr.hpp
  include/pinocchio/utils/static-if.hpp
  include/pinocchio/utils/string-generator.hpp
  include/pinocchio/utils/string.hpp
  include/pinocchio/utils/timer2.hpp
  include/pinocchio/utils/timer.hpp
  include/pinocchio/utils/version.hpp
)

SET(${PROJECT_NAME}_URDF_SOURCES
  src/parsers/urdf/model.cpp
  src/parsers/urdf/geometry.cpp
  src/parsers/urdf/utils.cpp
  )

SET(${PROJECT_NAME}_URDF_PUBLIC_HEADERS
  include/pinocchio/parsers/urdf.hpp
  include/pinocchio/parsers/urdf/model.hxx
  include/pinocchio/parsers/urdf/geometry.hxx
  include/pinocchio/parsers/urdf/utils.hpp
  include/pinocchio/parsers/urdf/types.hpp
  )

SET(${PROJECT_NAME}_SDF_SOURCES
  src/parsers/sdf/model.cpp
  src/parsers/sdf/geometry.cpp
  )

SET(${PROJECT_NAME}_SDF_PUBLIC_HEADERS
  include/pinocchio/parsers/sdf.hpp
  include/pinocchio/parsers/sdf/model.hxx
  include/pinocchio/parsers/sdf/geometry.hxx
  )

SET(${PROJECT_NAME}_LIBPYTHON_PUBLIC_HEADERS
  include/pinocchio/parsers/python.hpp
  )

SET(${PROJECT_NAME}_HPP_FCL_PUBLIC_HEADERS
  include/pinocchio/spatial/fcl-pinocchio-conversions.hpp
  )

SET(${PROJECT_NAME}_QHULL_SOURCES
  src/algorithm/reachable-workspace.cpp
  )

SET(${PROJECT_NAME}_QHULL_PUBLIC_HEADERS
  include/pinocchio/algorithm/reachable-workspace.hpp
  include/pinocchio/algorithm/reachable-workspace.hxx
  )

SET(_binary_headers_root ${${PROJECT_NAME}_BINARY_DIR}/include/pinocchio)
SET(${PROJECT_NAME}_CORE_GENERATED_PUBLIC_HEADERS
  ${_binary_headers_root}/config.hpp
  ${_binary_headers_root}/deprecated.hpp
  ${_binary_headers_root}/warning.hpp
  )

SET(${PROJECT_NAME}_TEMPLATE_INSTANTIATION_PUBLIC_HEADERS
  include/pinocchio/algorithm/aba-derivatives.txx
  include/pinocchio/algorithm/contact-cholesky.txx
  include/pinocchio/algorithm/contact-dynamics.txx
  include/pinocchio/algorithm/frames-derivatives.txx
  include/pinocchio/algorithm/regressor.txx
  include/pinocchio/algorithm/joint-configuration.txx
  include/pinocchio/algorithm/model.txx
  include/pinocchio/algorithm/proximal.txx
  include/pinocchio/algorithm/cholesky.txx
  include/pinocchio/algorithm/rnea.txx
  include/pinocchio/algorithm/jacobian.txx
  include/pinocchio/algorithm/geometry.txx
  include/pinocchio/algorithm/center-of-mass.txx
  include/pinocchio/algorithm/energy.txx
  include/pinocchio/algorithm/aba.txx
  include/pinocchio/algorithm/frames.txx
  include/pinocchio/algorithm/centroidal.txx
  include/pinocchio/algorithm/constrained-dynamics-derivatives.txx
  include/pinocchio/algorithm/contact-jacobian.txx
  include/pinocchio/algorithm/crba.txx
  include/pinocchio/algorithm/constrained-dynamics.txx
  include/pinocchio/algorithm/centroidal-derivatives.txx
  include/pinocchio/algorithm/rnea-derivatives.txx
  include/pinocchio/algorithm/impulse-dynamics.txx
  include/pinocchio/algorithm/compute-all-terms.txx
  include/pinocchio/algorithm/kinematics.txx
  include/pinocchio/algorithm/center-of-mass-derivatives.txx
  include/pinocchio/algorithm/kinematics-derivatives.txx
  include/pinocchio/multibody/model.txx
  include/pinocchio/multibody/data.txx
  include/pinocchio/parsers/sample-models.txx
  )

SET(${PROJECT_NAME}_TEMPLATE_INSTANTIATION_SOURCES
  src/algorithm/cholesky.cpp
  src/algorithm/aba.cpp
  src/algorithm/regressor.cpp
  src/algorithm/contact-dynamics.cpp
  src/algorithm/frames-derivatives.cpp
  src/algorithm/impulse-dynamics.cpp
  src/algorithm/model.cpp
  src/algorithm/constrained-dynamics.cpp
  src/algorithm/rnea-derivatives.cpp
  src/algorithm/compute-all-terms.cpp
  src/algorithm/jacobian.cpp
  src/algorithm/energy.cpp
  src/algorithm/centroidal-derivatives.cpp
  src/algorithm/frames.cpp
  src/algorithm/constrained-dynamics-derivatives.cpp
  src/algorithm/center-of-mass.cpp
  src/algorithm/geometry.cpp
  src/algorithm/kinematics.cpp
  src/algorithm/rnea.cpp
  src/algorithm/centroidal.cpp
  src/algorithm/aba-derivatives.cpp
  src/algorithm/crba.cpp
  src/algorithm/contact-cholesky.cpp
  src/algorithm/joint-configuration.cpp
  src/algorithm/center-of-mass-derivatives.cpp
  src/algorithm/proximal.cpp
  src/algorithm/kinematics-derivatives.cpp
  src/algorithm/contact-jacobian.cpp
  src/multibody/model.cpp
  src/multibody/data.cpp
  src/parsers/sample-models.cpp
  )

# Define Pinocchio Python binding sources and headers

SET(${PROJECT_NAME}_BINDINGS_PYTHON_PUBLIC_HEADERS
  include/pinocchio/bindings/python/spatial/explog.hpp
  include/pinocchio/bindings/python/spatial/classic-acceleration.hpp
  include/pinocchio/bindings/python/spatial/symmetric3.hpp
  include/pinocchio/bindings/python/spatial/motion.hpp
  include/pinocchio/bindings/python/spatial/force.hpp
  include/pinocchio/bindings/python/spatial/inertia.hpp
  include/pinocchio/bindings/python/spatial/se3.hpp
  include/pinocchio/bindings/python/algorithm/contact-info.hpp
  include/pinocchio/bindings/python/algorithm/contact-cholesky.hpp
  include/pinocchio/bindings/python/algorithm/proximal.hpp
  include/pinocchio/bindings/python/algorithm/algorithms.hpp
  include/pinocchio/bindings/python/algorithm/constraints/coulomb-friction-cone.hpp
  include/pinocchio/bindings/python/context.hpp
  include/pinocchio/bindings/python/pybind11-all.hpp
  include/pinocchio/bindings/python/pybind11.hpp
  include/pinocchio/bindings/python/utils/constant.hpp
  include/pinocchio/bindings/python/utils/version.hpp
  include/pinocchio/bindings/python/utils/pickle-vector.hpp
  include/pinocchio/bindings/python/utils/macros.hpp
  include/pinocchio/bindings/python/utils/std-vector.hpp
  include/pinocchio/bindings/python/utils/printable.hpp
  include/pinocchio/bindings/python/utils/dependencies.hpp
  include/pinocchio/bindings/python/utils/conversions.hpp
  include/pinocchio/bindings/python/utils/address.hpp
  include/pinocchio/bindings/python/utils/copyable.hpp
  include/pinocchio/bindings/python/utils/registration.hpp
  include/pinocchio/bindings/python/utils/pickle-map.hpp
  include/pinocchio/bindings/python/utils/list.hpp
  include/pinocchio/bindings/python/utils/std-aligned-vector.hpp
  include/pinocchio/bindings/python/utils/eigen.hpp
  include/pinocchio/bindings/python/utils/comparable.hpp
  include/pinocchio/bindings/python/utils/std-map.hpp
  include/pinocchio/bindings/python/utils/cast.hpp
  include/pinocchio/bindings/python/utils/deprecation.hpp
  include/pinocchio/bindings/python/utils/namespace.hpp
  include/pinocchio/bindings/python/context/cppadcg.hpp
  include/pinocchio/bindings/python/context/cppad.hpp
  include/pinocchio/bindings/python/context/casadi.hpp
  include/pinocchio/bindings/python/context/default.hpp
  include/pinocchio/bindings/python/context/mpfr.hpp
  include/pinocchio/bindings/python/context/generic.hpp
  include/pinocchio/bindings/python/fwd.hpp
  include/pinocchio/bindings/python/serialization/serialization.hpp
  include/pinocchio/bindings/python/serialization/serializable.hpp
  include/pinocchio/bindings/python/multibody/data.hpp
  include/pinocchio/bindings/python/multibody/frame.hpp
  include/pinocchio/bindings/python/multibody/geometry-model.hpp
  include/pinocchio/bindings/python/multibody/broadphase-manager-base.hpp
  include/pinocchio/bindings/python/multibody/model.hpp
  include/pinocchio/bindings/python/multibody/joint/joint-model.hpp
  include/pinocchio/bindings/python/multibody/joint/joint-derived.hpp
  include/pinocchio/bindings/python/multibody/joint/joints-models.hpp
  include/pinocchio/bindings/python/multibody/joint/joints-variant.hpp
  include/pinocchio/bindings/python/multibody/joint/joints-datas.hpp
  include/pinocchio/bindings/python/multibody/joint/joint-data.hpp
  include/pinocchio/bindings/python/multibody/liegroups.hpp
  include/pinocchio/bindings/python/multibody/geometry-data.hpp
  include/pinocchio/bindings/python/multibody/geometry-functors.hpp
  include/pinocchio/bindings/python/multibody/geometry-object.hpp
  include/pinocchio/bindings/python/multibody/tree-broadphase-manager.hpp
  include/pinocchio/bindings/python/parsers/urdf.hpp
  include/pinocchio/bindings/python/parsers/sample-models.hpp
  include/pinocchio/bindings/python/parsers/sdf.hpp
  include/pinocchio/bindings/python/parsers/srdf.hpp
  )

SET(${PROJECT_NAME}_BINDINGS_PYTHON_SOURCES
  bindings/python/spatial/expose-symmetric3.cpp
  bindings/python/spatial/expose-force.cpp
  bindings/python/spatial/expose-inertia.cpp
  bindings/python/spatial/expose-SE3.cpp
  bindings/python/spatial/expose-motion.cpp
  bindings/python/spatial/expose-explog.cpp
  bindings/python/spatial/expose-skew.cpp
  bindings/python/algorithm/expose-impulse-dynamics.cpp
  bindings/python/algorithm/expose-model.cpp
  bindings/python/algorithm/expose-centroidal.cpp
  bindings/python/algorithm/expose-aba.cpp
  bindings/python/algorithm/expose-algorithms.cpp
  bindings/python/algorithm/admm-solver.cpp
  bindings/python/algorithm/pgs-solver.cpp
  bindings/python/algorithm/expose-com.cpp
  bindings/python/algorithm/expose-frames.cpp
  bindings/python/algorithm/expose-energy.cpp
  bindings/python/algorithm/expose-contact-dynamics.cpp
  bindings/python/algorithm/expose-constrained-dynamics-derivatives.cpp
  bindings/python/algorithm/expose-rnea-derivatives.cpp
  bindings/python/algorithm/expose-kinematics.cpp
  bindings/python/algorithm/expose-geometry.cpp
  bindings/python/algorithm/expose-aba-derivatives.cpp
  bindings/python/algorithm/expose-contact-solvers.cpp
  bindings/python/algorithm/expose-joints.cpp
  bindings/python/algorithm/expose-constrained-dynamics.cpp
  bindings/python/algorithm/expose-rnea.cpp
  bindings/python/algorithm/expose-contact-jacobian.cpp
  bindings/python/algorithm/constraints/expose-cones.cpp
  bindings/python/algorithm/expose-cholesky.cpp
  bindings/python/algorithm/expose-regressor.cpp
  bindings/python/algorithm/expose-kinematics-derivatives.cpp
  bindings/python/algorithm/expose-cat.cpp
  bindings/python/algorithm/expose-frames-derivatives.cpp
  bindings/python/algorithm/expose-centroidal-derivatives.cpp
  bindings/python/algorithm/expose-jacobian.cpp
  bindings/python/algorithm/expose-delassus.cpp
  bindings/python/algorithm/expose-impulse-dynamics-derivatives.cpp
  bindings/python/algorithm/expose-kinematic-regressor.cpp
  bindings/python/algorithm/expose-crba.cpp
  bindings/python/module.cpp
  bindings/python/utils/version.cpp
  bindings/python/utils/dependencies.cpp
  bindings/python/utils/conversions.cpp
  bindings/python/math/expose-linalg.cpp
  bindings/python/math/expose-rpy.cpp
  bindings/python/math/expose-eigen-types.cpp
  bindings/python/serialization/serialization.cpp
  bindings/python/multibody/expose-model.cpp
  bindings/python/multibody/expose-liegroups.cpp
  bindings/python/multibody/expose-geometry.cpp
  bindings/python/multibody/expose-frame.cpp
  bindings/python/multibody/expose-data.cpp
  bindings/python/multibody/joint/expose-joints.cpp
  bindings/python/parsers/sdf/model.cpp
  bindings/python/parsers/sdf/geometry.cpp
  bindings/python/parsers/expose-parsers.cpp
  bindings/python/parsers/urdf/console-bridge.cpp
  bindings/python/parsers/urdf/model.cpp
  bindings/python/parsers/urdf/geometry.cpp
  bindings/python/parsers/srdf.cpp
  bindings/python/parsers/sample-models.cpp
  )

SET(${PROJECT_NAME}_BINDINGS_PYTHON_LIBPYTHON_SOURCES
  bindings/python/parsers/python/model.cpp
  )

SET(${PROJECT_NAME}_BINDINGS_PYTHON_LIBPYTHON_PUBLIC_HEADERS
  include/pinocchio/bindings/python/parsers/python.hpp
  )

SET(${PROJECT_NAME}_BINDINGS_PYTHON_HPP_FCL_SOURCES
  bindings/python/multibody/fcl/expose-fcl.cpp
  bindings/python/algorithm/expose-broadphase.cpp
  bindings/python/algorithm/expose-broadphase-callbacks.cpp
  )

SET(${PROJECT_NAME}_BINDINGS_PYTHON_HPP_FCL_PUBLIC_HEADERS
  include/pinocchio/bindings/python/multibody/fcl/transform.hpp
  include/pinocchio/bindings/python/multibody/broadphase-manager.hpp
  )

SET(${PROJECT_NAME}_BINDINGS_PYTHON_OPENMP_SOURCES
  bindings/python/multibody/pool/expose-pool.cpp
  bindings/python/algorithm/parallel/aba.cpp
  bindings/python/algorithm/parallel/expose-parallel.cpp
  bindings/python/algorithm/parallel/rnea.cpp
  )

SET(${PROJECT_NAME}_BINDINGS_PYTHON_OPENMP_PUBLIC_HEADERS
  include/pinocchio/bindings/python/multibody/pool/model.hpp
  )

SET(${PROJECT_NAME}_BINDINGS_PYTHON_HPP_FCL_OPENMP_SOURCES
  bindings/python/algorithm/parallel/geometry.cpp
  bindings/python/algorithm/parallel/broadphase.cpp
  )

SET(${PROJECT_NAME}_BINDINGS_PYTHON_HPP_FCL_OPENMP_PUBLIC_HEADERS
  include/pinocchio/bindings/python/multibody/pool/geometry.hpp
  include/pinocchio/bindings/python/multibody/pool/broadphase-manager.hpp
  )

SET(${PROJECT_NAME}_BINDINGS_PYTHON_QHULL_SOURCES
  bindings/python/algorithm/expose-reachable-workspace.cpp
  )

SET(${PROJECT_NAME}_BINDINGS_PYTHON_HPP_FCL_OPENMP_PUBLIC_HEADERS
  include/pinocchio/bindings/python/multibody/pool/geometry.hpp
  include/pinocchio/bindings/python/multibody/pool/broadphase-manager.hpp
  include/pinocchio/bindings/python/multibody/broadphase-manager.hpp
  )


SET(${PROJECT_NAME}_BINDINGS_PYTHON_EXTRA_MPFR_PUBLIC_HEADERS
  include/pinocchio/bindings/python/math/multiprecision/boost/number.hpp
  )

SET(${PROJECT_NAME}_BINDINGS_PYTHON_EXTRA_MPFR_SOURCES
  bindings/python/extra/mpfr/boost_number.cpp
  )
