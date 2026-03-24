#include "Interpolation/chebyshev_grid.hh"
#include <stdexcept>

namespace Interpolation
{
namespace Chebyshev
{
StandardGrid::StandardGrid(size_t p)
{
   _p = p;

   _tj.resize(_p + 1, 0.);
   _betaj.resize(_p + 1, 0.);
   for (size_t i = 0; i <= _p; i++) {
      _tj[i]    = cos(i * M_PI / (static_cast<double>(_p)));
      _betaj[i] = (i % 2 == 0) ? +1 : -1;
      if (i == 0 || i == p) _betaj[i] *= 0.5;
   }

   _Dij.resize(_p + 1, vector_d(_p + 1, 0.));
   _Dij[0][0]   = (2 * _p * _p + 1) / 6.0;
   _Dij[_p][_p] = -_Dij[0][0];
   for (size_t j = 1; j < _p; j++) {
      _Dij[j][j] = -_tj[j] * 0.5 / (1.0 - _tj[j] * _tj[j]);
   }
   for (size_t j = 0; j <= _p; j++) {
      for (size_t k = 0; k <= _p; k++) {
         if (k == j) continue;
         _Dij[j][k] = -(_betaj[k] / _betaj[j]) / (_tj[j] - _tj[k]);
      }
   }
   /// line
}

double StandardGrid::interpolate(double t, const vector_d &fj, size_t start, size_t end) const
{
   if (t < -1 || t > 1) {
      throw std::domain_error("StandardGrid::interpolate t must be in [-1, 1]");
   }
   if (end - start != _p) {
      throw std::domain_error("StandardGrid::interpolate end-start should be = to p");
   }

   double den = 0.;
   for (size_t j = 0; j <= _p; j++) {
      // if (t == _tj[j]) return fj[j + start];
      if (std::abs(t - _tj[j]) < 1.0e-15) return fj[j + start];

      den += _betaj[j] / (t - _tj[j]);
   }
   double res = 0.;
   for (size_t i = 0; i <= _p; i++) {
      res += poli_weight(t, i, den) * fj[i + start];
   }
   return res;
   /*
      double res = 0.;
      for (size_t i = 0; i <= _p; i++) {
         res += poli_weight(t, i) * fj[i + start];
      }
      return res;*/
}

double StandardGrid::poli_weight(double t, size_t j, double den) const
{
   if (std::abs(t - _tj[j]) < 1.0e-15) return 1.;

   double res = 0.;
   res        = _betaj[j] / (t - _tj[j]) / den;

   return res;
}

double StandardGrid::poli_weight(double t, size_t j) const
{
   if (std::abs(t - _tj[j]) < 1.0e-15) return 1.;

   double den = 0.;
   for (size_t j = 0; j <= _p; j++) {
      if (std::abs(t - _tj[j]) < 1.0e-15) return 0.;

      den += _betaj[j] / (t - _tj[j]);
   }

   double res = 0.;
   res        = _betaj[j] / (t - _tj[j]) / den;

   return res;
}

vector_d StandardGrid::discretize(const std::function<double(double)> &fnc) const
{
   vector_d fj(_p + 1, 0.);
   for (size_t i = 0; i <= _p; i++) {
      fj[i] = fnc(_tj[i]);
   }
   return fj;
}

} // namespace Chebyshev
} // namespace Interpolation