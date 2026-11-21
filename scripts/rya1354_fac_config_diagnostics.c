/* Read-only configuration state; original FAC library unchanged. */
static int PConfigDiagnostics(int argc, char *argv[], int argt[], ARRAY *variables) {
  int g, i, ng;
  FILE *f;
  CONFIG_GROUP *group;
  CONFIG *cfg;
  if (argc!=1 || argt[0]!=STRING) return -1;
  f=fopen(argv[0],"wx");
  if (!f) return -1;
  ng=GetNumGroups();
  fprintf(f,"# FAC_CONFIG_V1 mode=%d groups=%d\n",ConfigEnergyMode(),ng);
  fprintf(f,"# group config nele energy delta shift name\n");
  for (g=0;g<ng;g++) {
    group=GetGroup(g);
    for (i=0;i<group->n_cfgs;i++) {
      cfg=(CONFIG *)ArrayGet(&(group->cfg_list),i);
      fprintf(f,"%d %d %d %.17g %.17g %.17g %s\n",g,i,cfg->n_electrons,
              cfg->energy,cfg->delta,cfg->shift,group->name);
    }
  }
  fclose(f);
  return 0;
}
