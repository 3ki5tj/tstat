#ifndef AVP_H__
#define AVP_H__

/* routines for sampling a flat density distribution
  Copyright (c) 2012 Cheng Zhang

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public License
  as published by the Free Software Foundation; either
  version 2.1 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this program; if not, write to the Free Software
  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301 USA
*/

typedef struct {
  double vmin, vmax;
  int np; /* number of particles */
  int dof; /* degrees of freedom */
  double rhomin, rhomax, rhodel;
  double prv0; /* default pressure */
  double wlimit;
  int n;
  int ensexp; /* dV / V^ensexp */
  double *vis; /* # of visits */
  av_t *av, *av0;
} avp_t;

INLINE avp_t *avp_open(double rhomin, double rhomax, double rhodel,
    int np, int dof, int ensexp, double prv0, double wlimit)
{
  avp_t *avp;
  int i;

  die_if(rhomin >= rhomax, "rhomin %g >= rhomax %g\n", rhomin, rhomax);
  xnew(avp, 1);
  avp->np = np;
  avp->rhomin = rhomin;
  avp->rhomax = rhomax;
  avp->n = (int)((avp->rhomax - avp->rhomin)/rhodel + .5);
  avp->rhodel = (avp->rhomax - avp->rhomin)/avp->n;
  avp->vmin = np/rhomax;
  avp->vmax = np/rhomin;
  avp->dof = dof;
  avp->ensexp = ensexp;
  avp->prv0 = prv0;
  avp->wlimit = wlimit;
  xnew(avp->vis, avp->n + 1);
  xnew(avp->av,  avp->n + 1);
  xnew(avp->av0, avp->n + 1);
  for (i = 0; i <= avp->n; i++) {
    avp->vis[i] = 0;
    av_clear( &(avp->av[i]) );
    av_clear( &(avp->av0[i]) );
  }
  return avp;
}

INLINE void avp_close(avp_t *avp)
{
  free(avp->vis);
  free(avp->av);
  free(avp->av0);
  free(avp);
}

/* return 0 if successful */
INLINE int avp_add(avp_t *avp, double vol, double pr, double ngam)
{
  int i, j;
  double va, rho;

  if (vol < avp->vmin || vol >= avp->vmax) return 1;
  rho = avp->np / vol;
  i = (int)((rho - avp->rhomin) / avp->rhodel);
  die_if(i < 0 || i >= avp->n, "index %d out of range, n %d, vol %g, rho %g\n", i, avp->n, vol, rho);
  /* if the ensemble weight is 1/V^a, then
   * each point should be weighted as V^a */
  for (va = 1.0, j = 0; j < avp->ensexp; j++) va *= vol;
  avp->vis[i] += 1.0;
  av_gaddw( &(avp->av[i]), pr, va, ngam);
  av_addw( &(avp->av0[i]), pr, va);
  return 0;
}

#define AVP_GETPRV(avp, i) (avp->vis[i] < avp->wlimit ? avp->prv0 : av_getave( &(avp->av[i]) ))

/* get the average pressure */
INLINE double avp_getprv(avp_t *avp, double vol)
{
  int i;
  double rho;

  if (vol < avp->vmin || vol >= avp->vmax) return avp->prv0;
  rho = avp->np / vol;
  i = (int)( (rho - avp->rhomin) / avp->rhodel );
  die_if(i < 0 || i >= avp->n, "index %d out of range, n %d, vol %g, rho %g\n", i, avp->n, vol, rho);
  return AVP_GETPRV(avp, i);
}

/* get \int_vol1^vol2 p dV */
INLINE double avp_getdW(avp_t *avp, double vol1, double vol2)
{
  int i1, i2, i, sgn = 1, np = avp->np;
  double dW = 0, rho1, rho2;

  if (vol1 > vol2) { sgn = -1; dblswap(vol1, vol2); }
  if (vol2 < avp->vmin || vol1 >= avp->vmax) return (vol2 - vol1) * sgn * avp->prv0;

  /* switch to rho */
  rho1 = np / vol2;
  rho2 = np / vol1;

  if (rho1 < avp->rhomin) { dW += (np/rho1 - np/avp->rhomin) * avp->prv0; rho1 = avp->rhomin; }
  if (rho2 >= avp->rhomax) { dW += (np/avp->vmax - np/rho2) * avp->prv0; rho2 = avp->rhomax - 1e-8; }

  i1 = (int)( (rho1 - avp->rhomin) / avp->rhodel );
  i2 = (int)( (rho2 - avp->rhomin) / avp->rhodel );
  die_if(i1 < 0 || i1 >= avp->n || i2 < 0 || i2 >= avp->n,
    "i1 %d, i2 %d out of range, n %d, rho1 %g, rho2 %g\n", i1, i2, avp->n, rho1, rho2);

  if (i1 == i2) { /* same bin */
    dW += (vol2 - vol1) * AVP_GETPRV(avp, i1);
  } else { /* crosses multiple bins */
    dW += (vol2 - np / (avp->rhomin + (i1+1)*avp->rhodel)) * AVP_GETPRV(avp, i1);
    dW += (np / (avp->rhomin + i2*avp->rhodel) - vol1) * AVP_GETPRV(avp, i2);
    for (i = i1 + 1; i < i2; i++)
      dW += AVP_GETPRV(avp, i) * (np/(avp->rhomin + i*avp->rhodel)
         - np/(avp->rhomin + (i+1)*avp->rhodel));
  }
  return sgn * dW;
}

/* write the output */
INLINE int avp_write(avp_t *avp, const char *fn)
{
  int i;
  FILE *fp;
  double prv, w, prv0, cnt0, tot, scal, rho, w0;

  xfopen(fp, fn, "w", return -1);
  /* count the total */
  for (tot = 0, i = 0; i < avp->n; i++)
    tot += avp->vis[i];
  scal = 1.0/tot;
  fprintf(fp, "# %d %g %g %g %d %d %d %g\n", avp->n,
    avp->rhomin, avp->rhodel, tot, avp->np, avp->dof, avp->ensexp, avp->prv0);

  for (i = 0; i < avp->n; i++) {
    prv = av_getave( &(avp->av[i]) );
    w = avp->av[i].s;
    prv0 = av_getave( &(avp->av0[i]) );
    w0 = avp->av0[i].s;
    cnt0 = avp->vis[i];
    rho = avp->rhomin + i * avp->rhodel;
    fprintf(fp, "%g %g %g %g %g %g %g\n",
      rho, prv, w, prv0, cnt0*scal/avp->rhodel, cnt0, w0);
  }
  fclose(fp);
  return 0;
}

#endif

