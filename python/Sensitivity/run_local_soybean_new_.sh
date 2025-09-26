cd local_soybean_new_
for f in *.job; do
    sbatch "$f"
done