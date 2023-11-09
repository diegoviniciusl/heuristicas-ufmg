#!/bin/bash

echo -e "File\t\tTotal distance"

for tsp_file in files/*.tsp; do
    result=$(./tsp_nearest_insertion < $tsp_file | grep "Total distance of the tour")

    distance=$(echo $result | cut -d ' ' -f 6)

    echo -e "$tsp_file\t$distance"
done
