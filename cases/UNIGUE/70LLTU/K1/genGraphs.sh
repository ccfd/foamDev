time=$(foamListTimes -latestTime)
mkdir images
gnuplot -e "time='$time'" graphs