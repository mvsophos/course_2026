mkdir -p ordered
i=1
for f in $(ls *.png | sort -t '=' -k4 -n); do
    cp "$f" "$(printf "ordered/frame_%04d.png" $i)"
    i=$((i+1))
done
