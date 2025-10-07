cd local_soybean_new_free
for f in *.job; do
    sbatch "$f"
done