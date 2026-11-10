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
#include <cmath>
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

FULL_STORAGE=r'''
                    const bool brideFullOperator = ComputeOperator && !(mode & FsMode::PureLambdaIteration);
                    std::vector<f64> brideChi;
                    std::vector<unsigned char> brideSeen, brideEmitting;
                    std::vector<int> brideRows, brideColumns;
                    if (brideFullOperator) {
                        brideChi.assign(atom.Nlevel * Nspace, 0.0);
                        brideSeen.assign(atom.Nlevel, 0);
                        brideEmitting.assign(atom.Nlevel, 0);
                    }
'''

FULL_ACCUMULATE=r'''
                        if (brideFullOperator) {
                            for (int level : {t.i, t.j}) {
                                if (!brideSeen[level]) { brideSeen[level] = 1; brideRows.push_back(level); }
                            }
                            if (!brideEmitting[t.j]) { brideEmitting[t.j] = 1; brideColumns.push_back(t.j); }
                            for (int k = 0; k < Nspace; ++k) {
                                const f64 weightedChi = (atom.n(t.i,k) * Vij(k) - atom.n(t.j,k) * Vji(k)) * atom.wla(kr,k) * wmu;
                                brideChi[t.i * Nspace + k] += weightedChi;
                                brideChi[t.j * Nspace + k] -= weightedChi;
                            }
                        }
'''

FULL_OUTER=r'''
                    if (brideFullOperator) {
                        for (int i : brideRows) for (int j : brideColumns) {
                            if (i == j) continue;
                            for (int k = 0; k < Nspace; ++k)
                                atom.Gamma(i,j,k) -= PsiStar(k) * brideChi[i * Nspace + k] * atom.U(j,k);
                        }
                    }
'''

CAP_RESPONSE=r'''
            if constexpr (ComputeOperator) {
                if (!(mode & FsMode::PureLambdaIteration) && !toObs) {
                    if (!params || !params->contains("BrideCapLengthM"))
                        throw std::runtime_error("Cap preconditioning requires explicit BrideCapLengthM");
                    if (atmos.Ndim != 1) throw std::runtime_error("Cap preconditioning requires 1D");
                    const f64 length = params->get_as<f64>("BrideCapLengthM");
                    const f64 direction = std::abs(atmos.muz(mu));
                    if (!(length >= 0.0) || !std::isfinite(length) || !(direction > 0.0) || !(chiTot(0) >= 0.0))
                        throw std::runtime_error("Invalid fixed-cap derivative coefficients");
                    const f64 response = chiTot(0) == 0.0 ? length / direction : -std::expm1(-chiTot(0) * length / direction) / chiTot(0);
                    PsiStar(0) += response;
                    if (params->contains("BrideCapResponse")) {
                        auto destination = params->get_as<F64View2D>("BrideCapResponse");
                        if (destination.shape(0) != spect.wavelength.shape(0) || destination.shape(1) != Nrays)
                            throw std::runtime_error("Wrong cap derivative diagnostic shape");
                        destination(la,mu) = response;
                    }
                }
            }
'''

def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--source',type=Path,required=True);p.add_argument('--output',type=Path,required=True);p.add_argument('--method',choices=['overlap_guard','transition_local','full_coupling','full_coupling_cap'],default='overlap_guard');a=p.parse_args()
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
  patched=original.replace('formal_sol_iteration_matrices_impl','bride_'+a.method+'_iteration_matrices_impl')
  first='f64 integrand = (Uji(k) + Vji(k) * Ieff(k)) - (PsiStar(k) * atom.chi(t.i, k) * atom.U(t.j, k));'
  second='integrand = (Vij(k) * Ieff(k)) - (PsiStar(k) * atom.chi(t.j, k) * atom.U(t.i, k));'
  if patched.count(first)!=1 or patched.count(second)!=1:raise ValueError('Held transition-local patch context changed')
  if a.method=='transition_local':
   patched=patched.replace(first,'f64 integrand = Uji(k) + Vji(k) * I(k) - PsiStar(k) * atom.n(t.i,k) * Vij(k) * Uji(k);').replace(second,'integrand = Vij(k) * I(k) - PsiStar(k) * atom.n(t.j,k) * Vij(k) * Uji(k);')
  else:
   patched=patched.replace(first,'f64 integrand = Uji(k) + Vji(k) * Ieff(k);').replace(second,'integrand = Vij(k) * Ieff(k);')
   storage='''                    auto& atom = *activeAtoms[a];
                    if constexpr (ComputeOperator)'''
   accumulate='''                        const bool computeRates = (UpdateRates && !PrdRatesOnly) ||'''
   outer='''                                                computeRates, &atom, kr, wmu, &data);
                    }
                }
            }
            if constexpr (UpdateRates)'''
   if patched.count(storage)!=1 or patched.count(accumulate)!=2 or patched.count(outer)!=1:raise ValueError('Held full-coupling patch context changed')
   patched=patched.replace(storage,'                    auto& atom = *activeAtoms[a];\n'+FULL_STORAGE+'                    if constexpr (ComputeOperator)')
   patched=patched.replace(accumulate,FULL_ACCUMULATE+accumulate,1)
   patched=patched.replace(outer,'                                                computeRates, &atom, kr, wmu, &data);\n                    }\n'+FULL_OUTER+'                }\n            }\n            if constexpr (UpdateRates)')
   if a.method=='full_coupling_cap':
    cap_hook='''            if constexpr (UpdateRates || ComputeOperator)
            {'''
    if patched.count(cap_hook)!=1:raise ValueError('Held cap-derivative patch context changed')
    patched=patched.replace(cap_hook,CAP_RESPONSE+cap_hook)
 (target/'SimdFullIterationTemplates.hpp').write_text(patched)
 cpp=a.output/'overlap_guard.cpp';provider=PROVIDER
 if a.method!='overlap_guard':provider=provider.replace('bride_guard_iteration_matrices_impl','bride_'+a.method+'_iteration_matrices_impl').replace('bride_overlap_guard_scalar','bride_'+a.method+'_scalar')
 cpp.write_text(provider)
 includes=[target]+sorted({x.parent for name in ['TaskScheduler.h','JasPP.hpp']for x in a.source.rglob(name)})
 library=a.output/'bride_overlap_guard.dylib'
 command=['clang++','-std=c++17','-O2','-dynamiclib','-undefined','dynamic_lookup']+[arg for inc in includes for arg in ['-I',str(inc)]]+[str(cpp),'-o',str(library)]
 with (a.output/'build.log').open('x')as stream:subprocess.run(command,stdout=stream,stderr=subprocess.STDOUT,check=True)
 if library.stat().st_size>5*1024**2:raise ValueError('Plugin exceeds original 5MiB guard')
 receipt=dict(schema='codex.bride_overlap_plugin_build/1',scientific_admission='HOLD',method=a.method,scheme='bride_'+a.method+'_scalar',generator=fingerprint(__file__),backend=fingerprint(compiled.__file__),source_git_sha=subprocess.check_output(['git','-C',str(a.source),'rev-parse','HEAD'],text=True).strip(),source_dirty=subprocess.check_output(['git','-C',str(a.source),'status','--porcelain'],text=True),source_headers={str(x.relative_to(source)):fingerprint(x)for x in source.rglob('*')if x.is_file()},original_template=fingerprint(source/'SimdFullIterationTemplates.hpp'),patched_template=fingerprint(target/'SimdFullIterationTemplates.hpp'),provider=fingerprint(cpp),plugin=fingerprint(library),command=command,scope='Isolated custom iteration plugin; one thread; overlap_guard disables MALI at distinct-ion continuum overlap, transition_local uses pairwise rate-preserving MALI. No installed source/backend/default changes. ABI and physical qualification required independently.')
 (a.output/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt['plugin']))

if __name__=='__main__':main()
