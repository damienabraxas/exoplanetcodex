"""Build isolated single-thread LW iteration plugin; never modify installed LW."""
import argparse,json,shutil,subprocess,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from pipeline.bride_reference import fingerprint

PROVIDER=r'''
#include "Lightweaver.hpp"
#include "ThreadStorage.hpp"
#include "SimdFullIterationTemplates.hpp"
#include <stdexcept>
static IterationResult bride_iteration(Context& ctx, bool lambdaIterate, ExtraParams params) {
    if(ctx.Nthreads != 1) throw std::runtime_error("Bride overlap control requires one thread");
    auto mode=LwInternal::FsMode::UpdateJ | LwInternal::FsMode::UpdateRates;
    if(lambdaIterate) mode=mode | LwInternal::FsMode::PureLambdaIteration;
    return LwInternal::bride_guard_iteration_matrices_impl<SimdType::Scalar>(ctx, mode, params);
}
extern "C" FsIterationFns fs_iteration_fns_provider() {
    return FsIterationFns{-1,false,true,true,true,"bride_overlap_guard_scalar",
        bride_iteration,formal_sol_scalar,formal_sol_full_stokes_impl,
        redistribute_prd_lines_scalar,stat_eq_impl,time_dependent_update_impl,nr_post_update_impl};
}
'''

GUARD=r'''
            FsMode wavelengthMode = mode;
            bool overlappingIonContinua = false;
            for (auto* active : activeAtoms) {
                int ion = -1;
                for (auto* transition : active->trans) {
                    if (transition->type != CONTINUUM || !transition->is_active(la)) continue;
                    if (ion >= 0 && ion != transition->j) overlappingIonContinua = true;
                    ion = transition->j;
                }
            }
            if (overlappingIonContinua) wavelengthMode = wavelengthMode | FsMode::PureLambdaIteration;
'''

def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--source',type=Path,required=True);p.add_argument('--output',type=Path,required=True);p.add_argument('--method',choices=['overlap_guard','transition_local'],default='overlap_guard');a=p.parse_args()
 import lightweaver.LwCompiled as compiled
 a.output.mkdir(parents=True,exist_ok=False)
 source=a.source/'Source';target=a.output/'Source';shutil.copytree(source,target)
 original=(source/'SimdFullIterationTemplates.hpp').read_text();patched=original.replace('formal_sol_iteration_matrices_impl','bride_guard_iteration_matrices_impl')
 needle='''        for (int la = 0; la < Nspect; ++la)
        {
            f64 dJ = dispatch_intensity_core_opt_<simd>(true, false, true,
                                                        storeDepthData,
                                                        iCore, la * ctx.formalSolver.width,
                                                        mode, &params);'''
 replacement='''        for (int la = 0; la < Nspect; ++la)
        {'''+GUARD+'''
            f64 dJ = dispatch_intensity_core_opt_<simd>(true, false, true,
                                                        storeDepthData,
                                                        iCore, la * ctx.formalSolver.width,
                                                        wavelengthMode, &params);'''
 if patched.count(needle)!=1:raise ValueError('Held source patch context changed')
 if a.method=='overlap_guard':patched=patched.replace(needle,replacement)
 else:
  patched=original.replace('formal_sol_iteration_matrices_impl','bride_transition_local_iteration_matrices_impl')
  first='f64 integrand = (Uji(k) + Vji(k) * Ieff(k)) - (PsiStar(k) * atom.chi(t.i, k) * atom.U(t.j, k));'
  second='integrand = (Vij(k) * Ieff(k)) - (PsiStar(k) * atom.chi(t.j, k) * atom.U(t.i, k));'
  if patched.count(first)!=1 or patched.count(second)!=1:raise ValueError('Held transition-local patch context changed')
  patched=patched.replace(first,'f64 integrand = Uji(k) + Vji(k) * I(k) - PsiStar(k) * atom.n(t.i,k) * Vij(k) * Uji(k);').replace(second,'integrand = Vij(k) * I(k) - PsiStar(k) * atom.n(t.j,k) * Vij(k) * Uji(k);')
 (target/'SimdFullIterationTemplates.hpp').write_text(patched)
 cpp=a.output/'overlap_guard.cpp';provider=PROVIDER
 if a.method=='transition_local':provider=provider.replace('bride_guard_iteration_matrices_impl','bride_transition_local_iteration_matrices_impl').replace('bride_overlap_guard_scalar','bride_transition_local_scalar')
 cpp.write_text(provider)
 includes=[target]+sorted({x.parent for name in ['TaskScheduler.h','JasPP.hpp']for x in a.source.rglob(name)})
 library=a.output/'bride_overlap_guard.dylib'
 command=['clang++','-std=c++17','-O2','-dynamiclib','-undefined','dynamic_lookup']+[arg for inc in includes for arg in ['-I',str(inc)]]+[str(cpp),'-o',str(library)]
 with (a.output/'build.log').open('x')as stream:subprocess.run(command,stdout=stream,stderr=subprocess.STDOUT,check=True)
 if library.stat().st_size>5*1024**2:raise ValueError('Plugin exceeds original 5MiB guard')
 receipt=dict(schema='codex.bride_overlap_plugin_build/1',scientific_admission='HOLD',method=a.method,scheme='bride_'+a.method+'_scalar',generator=fingerprint(__file__),backend=fingerprint(compiled.__file__),source_git_sha=subprocess.check_output(['git','-C',str(a.source),'rev-parse','HEAD'],text=True).strip(),source_dirty=subprocess.check_output(['git','-C',str(a.source),'status','--porcelain'],text=True),source_headers={str(x.relative_to(source)):fingerprint(x)for x in source.rglob('*')if x.is_file()},original_template=fingerprint(source/'SimdFullIterationTemplates.hpp'),patched_template=fingerprint(target/'SimdFullIterationTemplates.hpp'),provider=fingerprint(cpp),plugin=fingerprint(library),command=command,scope='Isolated custom iteration plugin; one thread; overlap_guard disables MALI at distinct-ion continuum overlap, transition_local uses pairwise rate-preserving MALI. No installed source/backend/default changes. ABI and physical qualification required independently.')
 (a.output/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt['plugin']))

if __name__=='__main__':main()
