mkdir sofia
cp -R avx2 sofia/avx2
cp -R ref sofia/ref
cp -R scripts sofia/scripts
cp ../scripts/signature_size.py sofia/scripts
cp README.md sofia/
cp test_compatibility.sh sofia/
make -C sofia/avx2 clean > /dev/null
make -C sofia/ref/c clean > /dev/null
find sofia -name '.git' -delete
COPYFILE_DISABLE=1 tar -czv sofia
rm -r sofia
