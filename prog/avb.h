#ifndef AVB_H__
#define AVB_H__

/* routines for sampling a flat energy distribution
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
  double emin, emax, edel;
  double bet0; /* default temperature */
  int n, dof;
  av_t *av, *av0, *avu;
} avb_t;

INLINE avb_t *avb_open(double emin, double emax, double edel,
    double bet0, int dof)
{
  avb_t *avb;
  int i;

  xnew(avb, 1);
  avb->emin = emin;
  avb->emax = emax;
  avb->n = (int)((emax - emin)/edel + .5);
  avb->bet0 = bet0;
  avb->dof = dof;
  avb->edel = (emax - emin)/avb->n;
  xnew(avb->av,  avb->n + 1);
  xnew(avb->av0, avb->n + 1);
  xnew(avb->avu, avb->n + 1);
  for (i = 0; i <= avb->n; i++) {
    av_clear( &(avb->av[i]) );
    av_clear( &(avb->av0[i]) );
    av_clear( &(avb->avu[i]) );
  }
  return avb;
}

INLINE void avb_close(avb_t *avb)
{
  free(avb->av);
  free(avb->av0);
  free(avb->avu);
  free(avb);
}

/* return 0 if successful
   for total-energy entropic sampling */
INLINE int avb_add(avb_t *avb, double etot, double ekin, double ngam)
{
  int i;

  if (etot < avb->emin || etot >= avb->emax) return 1;
  i = (int)((etot - avb->emin) / avb->edel);
  die_if(i < 0 || i >= avb->n, "index %d out of range, n %d, etot %g\n", i, avb->n, etot);
  av_gadd( &(avb->av[i]), (.5*avb->dof - 1)/ekin, ngam);
  av_add( &(avb->av0[i]), (.5*avb->dof - 1)/ekin);
  av_add( &(avb->avu[i]), etot - ekin);
  return 0;
}

/* return 0 if successful
   for potential-energy entropic sampling */
INLINE int avb_addbet(avb_t *avb, double etot, double bet, double f2, double ngam)
{
  int i;

  if (etot < avb->emin || etot >= avb->emax) return 1;
  i = (int)((etot - avb->emin) / avb->edel);
  die_if(i < 0 || i >= avb->n, "index %d out of range, n %d, etot %g\n", i, avb->n, etot);
  av_gadd( &(avb->av[i]), bet, ngam);
  av_add( &(avb->av0[i]), bet);
  av_add( &(avb->avu[i]), f2);
  return 0;
}


/* return 0 if successful (ratio method)
   for the formula bet = < Lap U > / < F.F > */
INLINE int avb_addbetrat(avb_t *avb, double etot, double bet, double f2, double ngam)
{
  int i;

  if (etot < avb->emin || etot >= avb->emax) return 1;
  i = (int)((etot - avb->emin) / avb->edel);
  die_if(i < 0 || i >= avb->n, "index %d out of range, n %d, etot %g\n", i, avb->n, etot);
  av_gaddw( &(avb->av[i]), bet, f2, ngam);
  av_addw( &(avb->av0[i]), bet, f2);
  av_add( &(avb->avu[i]), f2);
  return 0;
}


#define AVB_GETBET(avb, i) ( avb->avu[i].s < 100.0 ? avb->bet0 : av_getave( &(avb->av[i])) )

/* get the local temperature */
INLINE double avb_getbet(avb_t *avb, double etot)
{
  int i;

  if (etot < avb->emin || etot >= avb->emax) return avb->bet0;
  i = (int)( (etot - avb->emin) / avb->edel );
  die_if(i < 0 || i >= avb->n, "index %d out of range, n %d, etot %g\n", i, avb->n, etot);
  return AVB_GETBET(avb, i);
}

/* get the entropy difference */ 
INLINE double avb_getdS(avb_t *avb, double e1, double e2)
{
  int i1, i2, i, sgn = 1;
  double dS = 0;
  
  if (e1 > e2) { sgn = -1; dblswap(e1, e2); }
  if (e2 < avb->emin || e1 >= avb->emax) return (e2 - e1) * sgn * avb->bet0;

  if (e1 < avb->emin) { dS += (avb->emin - e1) * avb->bet0; e1 = avb->emin; }
  if (e2 >= avb->emax) { dS += (e2 - avb->emax) * avb->bet0; e2 = avb->emax - 1e-8; }
  
  i1 = (int)( (e1 - avb->emin) / avb->edel );
  i2 = (int)( (e2 - avb->emin) / avb->edel );
  die_if(i1 < 0 || i1 >= avb->n || i2 < 0 || i2 >= avb->n,
    "i1 %d, i2 %d out of range, n %d, e1 %g, e2 %g\n", i1, i2, avb->n, e1, e2);

  if (i1 == i2) {
    dS += (e2 - e1) * AVB_GETBET(avb, i1);
  } else {
    dS += ((avb->emin + (i1+1)*avb->edel) - e1) * AVB_GETBET(avb, i1);
    dS += (e2 - (avb->emin + i2*avb->edel)) * AVB_GETBET(avb, i2);
    for (i = i1 + 1; i < i2; i++)
      dS += AVB_GETBET(avb, i) * avb->edel;
  }
  return sgn * dS;
}



/* write the output */
INLINE int avb_write(avb_t *avb, const char *fn)
{
  int i;
  FILE *fp;
  double bet, cnt, bet0, cnt0, uav, tot, scal;

  xfopen(fp, fn, "w", return -1);
  /* count the total */
  for (tot = 0, i = 0; i < avb->n; i++) {
    tot += avb->avu[i].s;
  }
  scal = 1.0/(tot*avb->edel);
  fprintf(fp, "# %d %g %g %g %d %g\n", avb->n, avb->emin, avb->edel,
      tot, avb->dof, avb->bet0);
  
  for (i = 0; i < avb->n; i++) {
    bet = av_getave( &(avb->av[i]) );
    cnt = avb->av[i].s;
    bet0 = av_getave( &(avb->av0[i]) );
    cnt0 = avb->avu[i].s;
    uav = av_getave( &(avb->avu[i]) );
    fprintf(fp, "%g %g %g %g %g %g %g\n",
      avb->emin + i * avb->edel, bet, cnt, bet0, cnt0*scal, uav, cnt0);
  }
  fclose(fp);
  return 0;
}

/* change the kinetic energy by a Monte-Carlo move (exact) */
INLINE int avb_mcvrescale(avb_t *avb, real *v, int nd,
  real amp, real ep, real *ekin, real *tkin)
{
  int i;
  real ek1 = *ekin, s;
  double logek1, logek2, ek2, r, dS, etot1, etot2;

  logek1 = log(ek1);
  logek2 = logek1 + amp*(2.f*rnd0() - 1);
  ek2 = exp(logek2);
  etot1 = ek1 + ep;
  etot2 = ek2 + ep;
  /* causes troubles in LJ with the following limitation */
/*
  if ((etot1 > etotmin && etot1 < etotmax)
   && (etot2 > etotmax || etot2 < etotmin))
    return 0;
*/
  dS = avb_getdS(avb, etot1, etot2);
  r = dS - .5*avb->dof*(logek2 - logek1);
  
  if (r <= 0 || rnd0() < exp(-r)) {
    s = (real) sqrt(ek2/ek1);
    for (i = 0; i < nd; i++)
      v[i] *= s;
    *ekin = ek2;
    if (tkin) *tkin *= s*s;
    return 1;
  } else { /* do nothing otherwise */
    return 0;
  }
}

/* Exact Andersen theromstat */
INLINE int avb_mcandersen(avb_t *avb, real *v, int n, int d,
  real ep, real *ekin, real *tkin)
{
  int i, j;
  real ek1 = *ekin, ek2, eki1, eki2, vi[3];
  double r, sqtp, dS, bet1, bet2, etot1, etot2, lnrose;

  i = (int)(rnd0() * n);
  etot1 = ek1 + ep;
  bet1 = avb_getbet(avb, etot1);
  sqtp = 1.0/sqrt(bet1);
  for (eki1 = eki2 = 0, j = 0; j < d; j++) {
    eki1 += .5 * v[i*d + j] * v[i*d + j];
    vi[j] = (real) (sqtp * grand0());
    eki2 += .5 * vi[j] * vi[j];
  }
  ek2 = ek1 + eki2 - eki1;
  etot2 = ek2 + ep;
  
  bet2 = avb_getbet(avb, etot2);
  lnrose = .5*d*log(bet1/bet2);
  dS = avb_getdS(avb, etot1, etot2);
  r = bet2 * eki1 - bet1 * eki2 + dS + lnrose;

  if (r <= 0 || rnd0() < exp(-r)) {
    for (j = 0; j < d; j++)
      v[i*d + j] = vi[j];
    *ekin = ek2;
    *tkin = 2.*ek2/avb->dof; 
    return 1;
  } else {
    return 0;
  }
}

#endif

