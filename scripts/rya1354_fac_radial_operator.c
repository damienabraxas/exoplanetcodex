/* Observability only: original pinned FAC radial library is not modified.
 * Each ordered orbital pair is queried once in a cold, single-gauge process.
 * Native multipole cache keys omit gauge/frequency, so reuse is rejected here.
 * NR is requested ONLY in its implemented Babushkin/length form. */
static int PRawRadialE1(int argc, char *argv[], int argt[], ARRAY *variables) {
  static int seen[256][2], count=0, context_gauge=-1;
  static double context_request=-1.0;
  int n1,k1,n2,k2,g,i,j,h,naw;
  double request,e,aw,r,nr,*grid=NULL;
  ORBITAL *a,*b;
  FILE *f;
  if (argc!=7 || argt[0]!=STRING) return -1;
  for (h=1;h<7;h++) if (argt[h]!=NUMBER) return -1;
  n1=atoi(argv[1]);k1=atoi(argv[2]);n2=atoi(argv[3]);k2=atoi(argv[4]);
  request=atof(argv[5]);g=atoi(argv[6]);
  if ((g!=1 && g!=2) || !isfinite(request) || request<0 || count>=256 ||
      GetTransitionGauge()!=g || GetTransitionMode()!=0) return -1;
  if (count && (g!=context_gauge || request!=context_request)) return -1;
  i=OrbitalIndex(n1,k1,0.0);j=OrbitalIndex(n2,k2,0.0);
  a=GetOrbitalSolved(i);b=GetOrbitalSolved(j);
  if (!a || !b || !a->wfun || !b->wfun || a->n!=n1 || a->kappa!=k1 ||
      b->n!=n2 || b->kappa!=k2 || !isfinite(a->energy) || !isfinite(b->energy) ||
      a->energy>=0 || b->energy>=0 ||
      ReducedCL(GetJFromKappa(k1),2,GetJFromKappa(k2))==0.0) return -1;
  for (h=0;h<count;h++) if (seen[h][0]==i && seen[h][1]==j) return -1;
  e=request>0?request:fabs(a->energy-b->energy);
  if (!isfinite(e) || e<=0) return -1;
  aw=FINE_STRUCTURE_CONST*e;
  SetAWGrid(1,aw,aw);
  naw=GetAWGrid(&grid);
  if (naw!=1 || !grid || grid[0]!=aw) return -1;
  r=MultipoleRadialFR(aw,-1,i,j,g);
  nr=MultipoleRadialNR(-1,i,j,2);
  if (!isfinite(r) || !isfinite(nr)) return -1;
  f=fopen(argv[0],"a");if (!f) return -1;
  fseek(f,0,SEEK_END);
  if (ftell(f)==0) {
    fprintf(f,"# FAC_RADIAL_E1_V42 HARTREE_EV=%.17g ALPHA=%.17g\n",HARTREE_EV,FINE_STRUCTURE_CONST);
    fprintf(f,"# n1 kappa1 n2 kappa2 gauge mode orbital_gap_mode grid_points requested_au used_au eps1_au eps2_au aw native_aw radial_FR radial_NR_length\n");
  }
  fprintf(f,"%d %d %d %d %d %d %d %d %.17g %.17g %.17g %.17g %.17g %.17g %.17g %.17g\n",
          n1,k1,n2,k2,g,GetTransitionMode(),request==0.0,naw,
          request,e,a->energy,b->energy,aw,grid[0],r,nr);
  fclose(f);
  seen[count][0]=i;seen[count][1]=j;count++;
  context_gauge=g;context_request=request;
  return 0;
}
