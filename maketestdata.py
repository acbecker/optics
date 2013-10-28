import sys
import random

ncenters = int(sys.argv[1])
nepochs  = int(sys.argv[2])
scatter  = float(sys.argv[3]) # in arcsec
scatter /= 3600.              # in degrees

centers = []
for i in range(ncenters):
    ra  = 360 * random.random() - 180   # -180 to 180
    dec = 180 * random.random() - 90    # -90  to 90
    centers.append((ra, dec))

testfile = 'test_optics.dat'
buf = open(testfile, 'w')

for i in range(nepochs):
    for j in range(ncenters):
        ra  = centers[j][0]
        dec = centers[j][1]

        rao  = ra  + 2 * scatter * random.random() - scatter
        deco = dec + 2 * scatter * random.random() - scatter

        buf.write('%f %f\n' % (rao, deco))

buf.close()
        
