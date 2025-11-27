#ifndef __pinocchio_multibody_visitor2_hpp__
#define __pinocchio_multibody_visitor2_hpp__

#include "pinocchio/multibody/joint/joint-base.hpp"
#include <boost/variant2.hpp>

namespace pinocchio
{

  template<typename JointModelDerived>
  struct GetDataFromModel
  {
    using type =
      typename std::remove_reference_t<std::remove_cv_t<JointModelDerived>>::JointDataDerived;
  };

  template<typename JointModelDerived, typename JointDataVariant>
  typename GetDataFromModel<JointModelDerived>::type &
  getJointData(const JointModelBase<JointModelDerived> &, JointDataVariant & data)
  {
    return boost::variant2::get<typename GetDataFromModel<JointModelDerived>::type>(data);
  }

} // namespace pinocchio

#endif // ifndef __pinocchio_multibody_visitor2_hpp__
