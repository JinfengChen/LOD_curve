echo "trait distance"
python TraitDist.py --input ../input/test.trait > ../input/test.trait.dist

echo "Convert chromosome position to single genome and genenate midpoint of each chromosome"
python LOD_Genome.py --input ../input/MPR.cross.uniq.QTL.mr.table > ../input/MSU7.Chr.midpoint

echo "LOD_curve"
python LOD_curve.py --input ../input/test.trait.dist

