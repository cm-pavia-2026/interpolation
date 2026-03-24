#include "Interpolation/chebyshev_grid.hh"

namespace Interpolation
{
namespace Chebyshev
{
StandardGrid::StandardGrid(size_t p)
{
   _p = p;

   _tj.resize(_p + 1, 0.);
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
}

} // namespace Chebyshev
} // namespace Interpolation