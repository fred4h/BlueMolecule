%chk=nh3.chk
#p hf/sto-3g opt

Ammonia Z-matrix with dummy atom and variables

0 1
X
N1  1  1.0
H1  N1 nh  1  hnx
H2  N1 nh  1  hnx  H1  d1
H3  N1 nh  1  hnx  H1  -d1
Variables:
nh   1.0120
hnx  112.15
Constants:
d1   120.0

B 2 3 F

