# A POSIX variable
# Reset in case getopts has been used previously in the shell.
OPTIND=1         

# Initialize our own variables:

# for the -f 
custom_output_dir=false
custom_output_path=""

# for the -v option
verbose=1

# for the -n  and -l options
n_times_from_back=0
not_all_times=false

is_dry=0

show_help () {
  printf "Tool for copying and storing results from openFOAM cases.\n"
  printf "Options:\n"
  printf "\t -l : process latest time only\n"
  printf "\t -s : silent\n"
  printf "\t -f : select different output directory\n"
  printf "\t -h : display help\n"
  printf "\t -c : destination to openFoam case\n"
  printf "\t -n : process N times counting from the last time\n"
  printf "\n"
  printf "Done by Wojciech Sadowski (wojciech1sadowski@gmail.com)\n"
}

ver_msg(){
  if (( $verbose > $1 )); then printf "$2"; fi
}

error_msg()
{
  printf "\nFatal error:\n "$1""
}

dry_run()
{
  if [[ "$is_dry" -eq "0" ]]
  then 
    $1
  else
    ver_msg -1 "$1 - dry running"
  fi
}

check_case(){
  if  [ ! -d "./constant" ] | [ ! -d "./system" ]
  then
    error_msg "Folder is not openFOAM case.\nSpecify -c option if you want to run outside of openFOAM case.\n"
    exit 0;
  fi

  if [ ! -d "./postProcessing" ]
  then
    error_msg
    printf "Case does not contain postProcessing folder.\n"
    exit 0
  fi

}

make_results(){
  local name_date=$1

  if [ -d "./results" ]
  then
    ver_msg 1 "Results directory exists.\n"
  else
    ver_msg 1 "Making directory folder... "
    mkdir results
    ver_msg 1 "Done.\n"
  fi

  ver_msg 1 "Making date-time directory... "
  if [ -d "results/$name_date" ]
  then
    error_msg "Saved results with exact same date|time exists!\nRefusing to overwrite!\n"
    exit 0
  fi
  ver_msg 1 "Done.\n"

  cur_results="./results/$name_date"
  mkdir $cur_results
}

get_timesteps(){
  if [ "$not_all_times" = true ]
  then
    timesteps=$(ls $work_dir | sort -V | tail -n$n_times_from_back)
  else
    timesteps=$(ls $work_dir | sort -V)
  fi
}

while getopts "h?dlsv:f:c:n:" opt; do
    case "$opt" in
    h|\?)
        show_help
        exit 0
        ;;
    d)  is_dry=1
        ;;
    l)  not_all_times=true
        n_times_from_back=1
        ;;
    s)  verbose=0
        ;;
    v)  verbose=$OPTARG
        ;;
    f)  output_folder=$OPTARG
        ;;
    c)  cd $OPTARG
        ;;
    #r)  res_list=$OPTARG
    #    chosen_res=1
    #    ;;
    n)  not_all_times=true
        n_times_from_back=$OPTARG
        ;;
    esac
done


check_case

ver_msg 0 "Executing results script ";

if [[ "$is_dry" -eq "1" ]]
then 
  ver_msg -1 " -- Dry run -- "
fi

ver_msg -1 "\n"

#get case name and current date
case_name=${PWD##*/}
cur_date=$(date '+%d.%b.%Y-%H:%M:%S')

ver_msg 0 "Case name:.. $case_name \n"
ver_msg 0 "Date:....... $cur_date\n"

make_results $cur_date
ver_msg 0 "Results in:. $cur_results\n"

# get postProcessing direcotry in a variablw
post_dirs=$(ls postProcessing)

ver_msg 2 "Directories available in postProcessing:\n"
for p_dir in $post_dirs;
do
  ver_msg 2 " - $p_dir \n"
done

for p_dir in $post_dirs;
do
  ver_msg 0 "Processing directory $p_dir... \n"
  mkdir $cur_results/$p_dir

  work_dir=postProcessing/$p_dir
  get_timesteps
  for timestep in $timesteps
  do
    ver_msg 1 " - T = $timestep\n"
    ver_msg 2 "   : cp -r $work_dir/$timestep $cur_results/$p_dir\n"
    cp -r $work_dir/$timestep $cur_results/$p_dir
  done
  ver_msg 0 "Done\n"
done



