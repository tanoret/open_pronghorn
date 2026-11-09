# Orifice of the Newton anchor (stage 3, orifice_implicit.i) at the design power: the
# cold-frit drag multiplier chi_cf = exp(lambda_j) of each of the 12 axial slices.
# Written by extract_orifice.py from ../03_ludewig_anchor/anchor_out.csv; s_power = 1.

[Functions]
  [chi_cf]
    type = PiecewiseConstant
    axis = y
    direction = left
    x = '0 0.06078333333 0.1215666667 0.18235 0.2431333333 0.3039166667 0.3647 0.4254833333 0.4862666667 0.54705 0.6078333333 0.6686166667'
    y = '0.14798053 2.6937365 5.1310677 7.1514394 9.0632585 10.828177 12.480756 14.035824 15.678109 17.234944 18.896257 20.630455'
  []
[]
