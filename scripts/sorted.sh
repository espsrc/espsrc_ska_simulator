#!/bin/bash
echo "Sorting files..."
if [[ -n "$1" ]]; then
    BASE_DIR="$1"
else
    BASE_DIR="."
fi
cd "$BASE_DIR" || exit 1
for file in *; do
    if [[ -f "$file" && "$file" =~ ^([0-9]{8}_[0-9]{4}) ]]; then
        folder="${BASH_REMATCH[1]}"
        mkdir -p "$folder"
        rm -rf *.tmp
        mv "$file" "$folder/"
        echo "Movido: $file → $folder/"
    fi
    if [[ -d "$file" && "$file" =~ ^([0-9]{8}_[0-9]{4}).*\.MS$ ]]; then
        folder="${BASH_REMATCH[1]}"
        mkdir -p "$folder"
        mv "$file" "$folder/"
        echo "Movido: $file → $folder/"
    fi
done



