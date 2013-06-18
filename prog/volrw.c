/* reweight distribution */
#define ZCOM_PICK
#define ZCOM_UTIL
#define ZCOM_AV
//#define ZCOM_HIST
#define ZCOM_ARGOPT
#include "zcom.h"
#include "avp.h"

//char *fnhis = "vol.his";
char *fnavp = "avp.dat";

/* handle input arguments */
static void doargs(int argc, char **argv)
{
  argopt_t *ao = argopt_open(0);
  argopt_add(ao, "-f", NULL, &fnavp,    "avp file");
  //argopt_add(ao, "-s", NULL, &fnhis,    "histogram file");
  argopt_addhelp(ao, "-h");
  argopt_parse(ao, argc, argv);

  argopt_dump(ao);
  argopt_close(ao);
}

/* load avp data from file */
INLINE avp_t *avp_load(const char *fn)
{
  FILE *fp;
  avp_t *avp = NULL;
  char s[1024] = "";
  int i, n;
  double rhomin, rhomax, rhodel, rho, vis, p, p0, w, w0, cs;
  /* default parameters */
  int np = 108;
  int dof = 318;
  int ensexp = 2;
  double prv0 = 0.06;

  xfopen(fp, fn, "r", return NULL);
  fgets(s, sizeof s, fp);
  if (s[0] != '#') goto ERR;
  if (4 > sscanf(s+1, "%d%lf%lf%lf%d%d%d%lf",
      &n, &rhomin, &rhodel, &vis, &np, &dof, &ensexp, &prv0)) {
    fprintf(stderr, "first line of %s broken\n", fn);
    goto ERR;
  }
  rhomax = rhomin + n * rhodel;
  avp = avp_open(rhomin, rhomax, rhodel, np, dof, ensexp, prv0, 100);
  for (i = 0; i < n; i++) {
    if (fgets(s, sizeof s, fp) == NULL) goto ERR;
    w0 = 1;
    if (sscanf(s, "%lf%lf%lf%lf%lf%lf%lf",
          &rho, &p, &w, &p0, &cs, &avp->vis[i], &w0) < 6) goto ERR;
    avp->av[i].s = w;
    avp->av[i].sx = w*p;
    avp->av0[i].s = w0;
    avp->av0[i].sx = w0*p0;
  }
  fclose(fp);
  return avp;
ERR:
  fclose(fp);
  return NULL;
}

static void output(avp_t *avp, const char *fn)
{
  FILE *fp;
  int i, j, id, m = 10, nm;
  double dvol, *fe, *pr, rho0, rho1, p;

  nm = avp->n * m;
  xnew(pr, avp->n*m + 1);
  xnew(fe, avp->n*m + 1);
  /* compute the free energy */
  fe[0] = 0;
  for (i = 0; i < avp->n; i++) {
    p = av_getave( &avp->av[i] );
    for (j = 0; j < m; j++) { /* add finer bins */
      rho0 = avp->rhomin + (i + 1.0*j/m) * avp->rhodel;
      rho1 = avp->rhomin + (i + 1.0*(j+1)/m) * avp->rhodel;
      dvol = avp->np/rho1 - avp->np/rho0;
      id = i*m + j;
      pr[id] = p;
      fe[id + 1] = fe[id] - pr[id] * dvol;
    }
  }
  pr[nm] = pr[nm-1];

  xfopen(fp, fn, "w", return);
  for (i = 0; i <= nm; i++) {
    rho0 = avp->rhomin + i * avp->rhodel/m;
    fprintf(fp, "%g %g %g\n", rho0, fe[i], pr[i]);
  }
  fclose(fp);
  free(pr);
  free(fe);
}


int main(int argc, char **argv)
{
  avp_t *avp;
  //hist_t *hs;
  //int row, ver;
  //unsigned fflags;
  //double vmin, vmax, vdel;

  doargs(argc, argv);

  die_if ((avp = avp_load(fnavp)) == NULL, "cannot load avp from %s\n", fnavp);
  //histgetinfo(fnhis, &row, &vmin, &vmax, &vdel, &ver, &fflags);
  //hs = hs_open(1, vmin, vmax, vdel);
  //die_if (hs_load(hs, fnhis, HIST_VERBOSE) != 0, "cannot load histogram from %s\n", fnhis);
  output(avp, "fe.dat");
  //hs_close(hs);
  avp_close(avp);
  return 0;
}
