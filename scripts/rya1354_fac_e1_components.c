/* Observability only: the pinned FAC library is linked unmodified.
 * Records every signed term of a native E1 reduced amplitude in its ACTUAL
 * rate context: the CI levels, mixing, angular cut and photon grid left by the
 * immediately preceding unchanged TransitionTable/RawTransition calls.
 *   A rows: native AngularZMix order (the exact list TRMultipole sums).
 *   B rows: the same AngularZMix CSF-pair loop, term by term, INCLUDING the
 *           terms the native angular cut drops (flagged, never summed into A).
 *   L rows: each level's stored basis, mixing and configuration labels.
 * TRMultipole is NOT called here (its optional pair selector reads an
 * uninitialized variable); the preceding RawTransition row is the native
 * amplitude used for exact comparison. Any unsupported branch returns -1. */
static int PRawE1Components(int argc, char *argv[], int argt[], ARRAY *variables) {
  static int seen[64][2], count=0;
  int low,up,g,i,j,m,h,nz,nmk=0,naw,p1,p2,j1,j2,kmin,kmax,swapped,ns,nh,ih1,ih2,isz,flags,nt;
  double energy,aw,s,r,prod,r0,norm,*mbk=NULL,*grid=NULL;
  LEVEL *lev1,*lev2,*la,*lb,*lv;
  STATE *st,*pb1,*pb2;
  SYMMETRY *sym;
  SHAMILTON *hams;
  ANGZ_DATUM *ad;
  ANGZ_ARY *sub;
  ANGULAR_ZMIX *ang=NULL,*az,z;
  ORBITAL *o0,*o1;
  CONFIG_GROUP *grp;
  char name[LEVEL_NAME_LEN],sname[LEVEL_NAME_LEN],nc[LEVEL_NAME_LEN];
  int vnl;
  FILE *f;
  if (argc!=4 || argt[0]!=STRING || argt[1]!=NUMBER || argt[2]!=NUMBER ||
      argt[3]!=NUMBER || count>=64) return -1;
  low=atoi(argv[1]);up=atoi(argv[2]);g=atoi(argv[3]);
  if ((g!=G_COULOMB && g!=G_BABUSHKIN) || GetTransitionGauge()!=g ||
      GetTransitionMode()!=M_FR) return -1;
  for (h=0;h<count;h++) if (seen[h][0]==low && seen[h][1]==up) return -1;
  lev1=GetLevel(low);lev2=GetLevel(up);
  if (!lev1 || !lev2 || lev1->nele!=lev2->nele || lev1->iham<0 || lev2->iham<0) return -1;
  DecodePJ(lev1->pj,&p1,&j1);DecodePJ(lev2->pj,&p2,&j2);
  if (p1==p2 || !Triangle(j1,j2,2)) return -1;
  energy=lev2->energy-lev1->energy;
  if (!isfinite(energy) || energy<=0) return -1;
  aw=FINE_STRUCTURE_CONST*energy;
  sym=GetSymmetry(lev1->pj);pb1=(STATE *)ArrayGet(&(sym->states),lev1->pb);
  sym=GetSymmetry(lev2->pj);pb2=(STATE *)ArrayGet(&(sym->states),lev2->pb);
  if (!pb1 || !pb2 || pb1->kgroup<0 || pb2->kgroup<0) return -1;
  naw=GetAWGrid(&grid);
  if (naw<1 || !grid) return -1;
  hams=GetSHamilton(&nh);
  if (!hams || lev1->iham>=nh || lev2->iham>=nh) return -1;
  kmin=abs(j1-j2);kmax=j1+j2;
  if (IsOdd(kmin)) kmin++;
  kmin=Max(kmin,2);kmax=Min(kmax,2);
  if (kmin!=2 || kmax!=2) return -1;
  f=fopen(argv[0],"a");if (!f) return -1;
  fseek(f,0,SEEK_END);
  if (ftell(f)==0) {
    fprintf(f,"#\tFAC_E1_COMPONENTS_V44\tHARTREE_EV=%.17g\tALPHA=%.17g\tANGZCUT=%.17g\tMIXCUT=%.17g\n",
            HARTREE_EV,FINE_STRUCTURE_CONST,(double)ANGZCUT,(double)MIXCUT);
  }
  swapped=lev1->iham>lev2->iham;
  fprintf(f,"C\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%.17g\t%.17g\t%d\n",
          low,up,j1,j2,p1,p2,lev1->nele,g,GetTransitionMode(),lev1->iham,lev2->iham,energy,aw,swapped);
  for (h=0;h<naw;h++) fprintf(f,"G\t%d\t%.17g\n",h,grid[h]);
  for (m=0;m<2;m++) {
    lv=m?lev2:lev1;sym=GetSymmetry(lv->pj);norm=0.0;
    for (i=0;i<lv->n_basis;i++) {
      st=(STATE *)ArrayGet(&(sym->states),lv->basis[i]);
      if (!st || st->kgroup<0) {fclose(f);return -1;}
      grp=GetGroup(st->kgroup);
      name[0]=sname[0]=nc[0]='\0';
      ConstructLevelName(name,sname,nc,&vnl,st);
      fprintf(f,"L\t%d\t%d\t%d\t%d\t%.17g\t%d\t%d\t%d\t%s\t%s\t%s\t%s\n",m,i,(int)lv->ibasis[i],
              lv->basis[i],lv->mixing[i],st->kgroup,st->kcfg,st->kstate,grp->name,nc,sname,name);
      norm+=lv->mixing[i]*lv->mixing[i];
    }
    fprintf(f,"N\t%d\t%d\t%d\t%.17g\n",m,lv->n_basis,hams[lv->iham].nbasis,norm);
  }
  nz=AngularZMix(&ang,low,up,2,2,&nmk,&mbk);
  if (nmk>0) {free(mbk);if (nz>0) free(ang);fclose(f);return -1;}
  s=0.0;
  for (i=0;i<nz;i++) {
    if (ang[i].k!=2) continue;
    r=MultipoleRadialFR(aw,-1,ang[i].k0,ang[i].k1,g);
    s+=r*ang[i].coeff;
    o0=GetOrbital(ang[i].k0);o1=GetOrbital(ang[i].k1);
    fprintf(f,"A\t%d\t%d\t%d\t%d\t%d\t%.17g\t%.17g\t%.17g\t%.17g\n",i,o0->n,o0->kappa,o1->n,o1->kappa,
            ang[i].coeff,r,r*ang[i].coeff,s);
  }
  if (nz>0) free(ang);
  la=swapped?lev2:lev1;lb=swapped?lev1:lev2;
  ns=AngularZMixStates(&ad,la->iham,lb->iham);
  if (ns<=0 || !ad || !ad->angz) {fclose(f);return -1;}
  nt=0;
  for (i=0;i<la->n_basis;i++) {
    ih1=la->ibasis[i];
    for (j=0;j<lb->n_basis;j++) {
      ih2=lb->ibasis[j];
      isz=ih1*hams[lb->iham].nbasis+ih2;
      sub=(ad->angz)[isz];
      if (!sub) continue;
      prod=la->mixing[i]*lb->mixing[j];
      az=(ANGULAR_ZMIX *)sub->az;
      for (h=0;h<sub->nz;h++) {
        if (az[h].k>kmax || az[h].k<kmin) continue;
        r0=az[h].coeff*prod;
        flags=(fabs(la->mixing[i])<ANGZCUT)|((fabs(lb->mixing[j])<ANGZCUT)<<1)|
              ((fabs(prod)<ANGZCUT)<<2)|((fabs(r0)<ANGZCUT)<<3);
        z=az[h];z.coeff=r0;
        if (swapped) AngZSwapBraKet(1,&z,j1-j2);
        r=MultipoleRadialFR(aw,-1,z.k0,z.k1,g);
        o0=GetOrbital(z.k0);o1=GetOrbital(z.k1);
        fprintf(f,"B\t%d\t%d\t%d\t%.17g\t%d\t%d\t%d\t%d\t%d\t%.17g\t%d\t%.17g\t%.17g\n",
                swapped?j:i,swapped?i:j,h,az[h].coeff,az[h].k,o0->n,o0->kappa,o1->n,o1->kappa,
                z.coeff,flags,r,r*z.coeff);
        nt++;
      }
    }
  }
  fprintf(f,"E\t%d\t%d\t%.17g\n",low,nt,s);
  fclose(f);
  seen[count][0]=low;seen[count][1]=up;count++;
  return 0;
}
