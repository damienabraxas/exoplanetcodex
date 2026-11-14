/* Observability-only SFAC command; linked to the unmodified pinned FAC library.
 * TRMultipole/OscillatorStrength and constants are the same APIs used by
 * SaveTransition0/PrintTRTable. Nonzero native status is NOT mapped to zero. */
static int PRawTransition(int argc, char *argv[], int argt[], ARRAY *variables) {
  int low, up, m, rc, k0, k1, nz=0, nmk=0, jup, jlow;
  double strength=NAN, energy=0.0, gf=NAN, rate=NAN, *mbk=NULL;
  ANGULAR_ZMIX *ang=NULL;
  FILE *f;
  LEVEL *lower, *upper;
  if (argc!=4 || argt[0]!=STRING || argt[1]!=NUMBER ||
      argt[2]!=NUMBER || argt[3]!=NUMBER) return -1;
  low=atoi(argv[1]); up=atoi(argv[2]); m=atoi(argv[3]);
  lower=GetLevel(low); upper=GetLevel(up);
  if (!lower || !upper || !m) return -1;
  jup=LevelTotalJ(up); jlow=LevelTotalJ(low);
  rc=TRMultipole(&strength,&energy,m,low,up,&k0,&k1);
  nz=AngularZMix(&ang,low,up,2*abs(m),2*abs(m),&nmk,&mbk);
  if (nz>0) free(ang);
  if (nmk>0) free(mbk);
  if (!rc) {
    gf=OscillatorStrength(m,energy,strength,&rate);
    rate=rate/(jup+1.0)*RATE_AU;
  }
  f=fopen(argv[0],"a");
  if (!f) return -1;
  fseek(f,0,SEEK_END);
  if (ftell(f)==0) {
    fprintf(f,"# FAC_RAW_V1 HARTREE_EV=%.17g RATE_AU=%.17g\n",HARTREE_EV,RATE_AU);
    fprintf(f,"# upper lower upper_twoJ lower_twoJ nele multipole gauge mode status angular_nz nmk energy_au amplitude gf A_s1\n");
  }
  fprintf(f,"%d %d %d %d %d %d %d %d %d %d %d %.17g %.17g %.17g %.17g\n",
    up,low,jup,jlow,upper->nele,m,GetTransitionGauge(),GetTransitionMode(),
    rc,nz,nmk,energy,strength,gf,rate);
  fclose(f);
  return 0;
}
