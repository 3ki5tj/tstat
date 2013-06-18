/* reweight distribution */
#define ZCOM_PICK
#define ZCOM_UTIL
#define ZCOM_DEF
#define ZCOM_ARGOPT
#define ZCOM_AV
#define ZCOM_HIST
#define ZCOM_RNG
#define ZCOM_SPECFUNC
#include "zcom.h"
#include "avb.h"

char *fnhis = "epot.his";
char *fnavb = "avb.dat";

/* handle input arguments */
static void doargs(int argc, char **argv)
{
  argopt_t *ao = argopt_open(0);
  argopt_add(ao, "-f", NULL, &fnavb,    "avb file");
  argopt_add(ao, "-s", NULL, &fnhis,    "histogram file");
  argopt_addhelp(ao, "-h");
  argopt_parse(ao, argc, argv);

  argopt_dump(ao);
  argopt_close(ao);
}

/* load avb data from file */
INLINE avb_t *avb_load(const char *fn)
{
  FILE *fp;
  avb_t *avb = NULL;
  char s[1024] = "";
  int i, n;
  double emin, emax, edel, ene, vis, b, b0, w, w0, uav, cs;
  /* default parameters */
  int dof = 318;
  double bet0 = 1.0;

  if (strcmp(fn, "goavb.dat") == 0) { /* hack for 1VII simulation */
    dof = 102; bet0 = 0.5;
  }

  xfopen(fp, fn, "r", return NULL);
  fgets(s, sizeof s, fp);
  if (s[0] != '#') goto ERR;
  if (4 > sscanf(s+1, "%d%lf%lf%lf%d%lf",
      &n, &emin, &edel, &vis, &dof, &bet0)) {
    fprintf(stderr, "first line of %s broken\n", fn);
    goto ERR;
  }
  emax = emin + n * edel;
  avb = avb_open(emin, emax, edel, bet0, dof);
  for (i = 0; i < n; i++) {
    if (fgets(s, sizeof s, fp) == NULL) goto ERR;
    w0 = 1;
    if (sscanf(s, "%lf%lf%lf%lf%lf%lf%lf",
          &ene, &b, &w, &b0, &cs, &uav, &w0) < 7) goto ERR;
    avb->av[i].s = w;
    avb->av[i].sx = w*b;
    avb->av0[i].s = w0;
    avb->av0[i].sx = w0*b0;
    avb->avu[i].s = w0;
    avb->avu[i].sx = w0*uav;
  }
  fclose(fp);
  return avb;
ERR:
  fclose(fp);
  return NULL;
}

/* compute the density of energy states */
static void outpute(avb_t *avb, const char *fn)
{
  FILE *fp;
  int i, j, id, m = 10, nm;
  double *en, *bt, ene, b;

  nm = avb->n * m;
  xnew(bt, nm + 1);
  xnew(en, nm + 1);
  /* compute the free energy */
  en[0] = 0;
  for (i = 0; i < avb->n; i++) {
    b = av_getave( &avb->av[i] );
    //printf("i %d\n", i); getchar();
    for (j = 0; j < m; j++) { /* add finer bins */
      id = i * m + j;
      bt[id] = b;
      en[id + 1] = en[id] + bt[id] * avb->edel/m;
    }
  }
  bt[nm] = bt[nm-1];

  xfopen(fp, fn, "w", return);
  for (i = 0; i <= nm; i++) {
    ene = avb->emin + (i + .5) * avb->edel/m;
    fprintf(fp, "%g %g %g\n", ene, en[i], bt[i]);
  }
  fclose(fp);
  free(bt);
  free(en);
}

/* get the aggregate weight for a potential energy u
 * using incomplete gamma function */
INLINE double avb_getlnwtot(avb_t *avb, double u)
{
  int i;
  double x0, x1, bet, hdof = .5 * avb->dof;
  double lng, lng1, lng0, lnw = -1e300, lnomg = 0;

  for (i = -1; i <= avb->n; i++) {
    x0 = avb->emin + i * avb->edel - u;
    if (i < 0) x0 = -1e300;
    x1 = avb->emin + (i + 1) * avb->edel - u;
    if (i >= avb->n) x1 = 1e300;
    if (i >= 0 && i < avb->n) {
      bet = av_getave( &avb->av[i] );
    } else {
      bet = avb->bet0;
    }
    if (bet < 0) {
      fprintf(stderr, "bad bet %g, i %d/n %d\n", bet, i, avb->n);
      continue;
    }
    if (x1 > 0) {
      if (i == avb->n) {
        lng1 = lngam(hdof);
      } else {
        lng1 = lnincgam(hdof, bet * x1);
      }
      if (x0 > 0) lng0 = lnincgam(hdof, bet * x0);
      else lng0 = -1e300;
      if (lng1 < lng0) {
        printf("bad x0 %g, lng0 %g, x1 %g, lng1 %g, i %d, bet %g\n",
            x0, lng0, x1, lng1, i, bet);
      } else {
        if (i == 0) {
          lng = lnincgam(hdof, bet*x1);
        }else if (i == avb->n) {
          lng = lnincgamup(hdof, bet*x0);
        } else {
          lng = lndif(lng1, lng0);
        }
        lng += -hdof * log(bet) + bet * x0 - lnomg;
        lnw = lnadd(lnw, lng);
      }
    }
    lnomg += bet * avb->edel;
  }
  return lnw;
}

/* get the aggregate weight for a potential energy u
 * (backup routine) approximate */
INLINE double avb_getlnwtot0(avb_t *avb, double u)
{
  int i;
  double x, bet, hdof = .5 * avb->dof;
  double lng, lnw = -1e300, lnomg = 0;

  for (i = 0; i < avb->n; i++) {
    x = avb->emin + (i + .5) * avb->edel - u;
    bet = av_getave( &avb->av[i] );
    lnomg += bet * avb->edel;
    if (x > 0) {
      lng = (hdof-1) * log(x) - lnomg + log(avb->edel);
      lnw = lnadd(lnw, lng);
    }
  }
  return lnw;
}

/* write the density of potential-energy states */
static void outputu(avb_t *avb, hist_t *hs, const char *fn)
{
  FILE *fp;
  int i;
  double s, u, *lnw, *lng;

  xnew(lnw, hs->n + 1);
  xnew(lng, hs->n + 1);
  for (i = 0; i < hs->n; i++) {
    s = hs->arr[i];
    if (s <= 0) continue;
    u = hs->xmin + (i + .5) * hs->dx;
    lnw[i] = avb_getlnwtot(avb, u);
    lng[i] = log(s) - lnw[i];
  }
  //for (i = hs->n; i >= 0; i--) lng[i] -= lng[0];

  xfopen(fp, fn, "w", return);
  for (i = 0; i < hs->n; i++) {
    if (hs->arr[i] <= 0) continue;
    u = hs->xmin + (i + .5) * hs->dx;
    fprintf(fp, "%g %g %g %g\n", u, lng[i], hs->arr[i], lnw[i]);
  }
  fclose(fp);
  free(lnw);
  free(lng);
}



int main(int argc, char **argv)
{
  avb_t *avb;
  hist_t *hs;
  int row, ver;
  unsigned fflags;
  double vmin, vmax, vdel;

  doargs(argc, argv);

  die_if ((avb = avb_load(fnavb)) == NULL,
      "cannot load avb from %s\n", fnavb);
  histgetinfo(fnhis, &row, &vmin, &vmax, &vdel, &ver, &fflags);
  hs = hs_open(1, vmin, vmax, vdel);
  die_if (hs_load(hs, fnhis, HIST_VERBOSE) != 0,
      "cannot load histogram from %s\n", fnhis);
  outpute(avb, "se.dat");
  outputu(avb, hs, "su.dat");
  hs_close(hs);
  avb_close(avb);
  return 0;
}

