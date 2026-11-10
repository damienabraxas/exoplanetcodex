
#include "Lightweaver.hpp"
#include "ThreadStorage.hpp"
#include "SimdFullIterationTemplates.hpp"
#include <stdexcept>
static IterationResult bride_iteration(Context& ctx, bool lambdaIterate, ExtraParams params) {
    if(ctx.Nthreads != 1) throw std::runtime_error("Bride overlap control requires one thread");
    auto mode=LwInternal::FsMode::UpdateJ | LwInternal::FsMode::UpdateRates;
    if(lambdaIterate) mode=mode | LwInternal::FsMode::PureLambdaIteration;
    return LwInternal::bride_transition_local_iteration_matrices_impl<SimdType::Scalar>(ctx, mode, params);
}
extern "C" FsIterationFns fs_iteration_fns_provider() {
    return FsIterationFns{-1,false,true,true,true,"bride_transition_local_scalar",
        bride_iteration,formal_sol_scalar,formal_sol_full_stokes_impl,
        redistribute_prd_lines_scalar,stat_eq_impl,time_dependent_update_impl,nr_post_update_impl};
}
