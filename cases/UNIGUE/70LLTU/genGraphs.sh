timeR=$(foamListTimes -case R1 -latestTime)
timeK=$(foamListTimes -case K1 -latestTime)
mkdir images
gnuplot -e "timeR='$timeR'; timeK='$timeK'" graphs