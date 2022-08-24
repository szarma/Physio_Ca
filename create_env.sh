## A script to install or update islets in a conda environment
## Examples:
## create_env.sh physio1          ## will create an environment physio1, install islets inside and kreate a jupyter kernel with the same name
## create_env.sh physio1 --live   ## as above, but it will a developer (live) version of the package (add --silent if you wish less output)
## create_env.sh physio1 --update ## if the venv exists, and the islets is installed (in static, non-dev, mode) updates the islets

venv_name=$1
conda_path=`which conda`
main_path=$(dirname $(dirname $conda_path))
env_path="${main_path}/envs/${venv_name}"
echo "venv_name:   ${venv_name}"
echo "conda path:  ${conda_path}"
echo "main path:   ${main_path}"
echo "env path:    ${env_path}"



ARGS=$(getopt -a --options slu --long "silent,live,update" -- "$@")
eval set -- "$ARGS"
silent="false"
live="false"
update="false"
while true; do
  case "$1" in
    -s|--silent)
      silent="true"
      shift;;
    -l|--live)
      live="true"
      shift;;
    -u|--update)
      update="true"
      shift;;
    --)
      break;;
     *)
      printf "Unknown option %s\n" "$1"
      exit 1;;
  esac
done

#echo "silent: $silent"
#echo "live: $live"
#echo "update: $update"

if [ $update == true ]
then
  version=$(<VERSION)
  sudo -E $env_path/bin/python setup.py bdist_wheel -d ./wheels
  sudo -E $env_path/bin/pip install "./wheels/islets-${version}-py3-none-any.whl"
else
  if [ $silent != true ]; then echo "============== Attempting to set up the environment '${venv_name}'"; fi
  #sudo -E -s ${conda_path} create -n $venv_name python==3.8 pip ipykernel setuptools -y
  sudo -E -s ${conda_path} create -n $venv_name python==3.8 pip ipykernel setuptools python-javabridge -c conda-forge -y

  #if [ $silent != true ]; then echo "============== Installing python-javabridge"; fi
  #sudo -E ${conda_path} install -n $venv_name python-javabridge -y -c conda-forge

  if [ $silent != true ]; then echo "============== Installing python-bioformats"; fi
  sudo $env_path/bin/pip install python-bioformats==4.0.0

  if [ $silent != true ]; then echo "============== Installing islets"; fi

  if [ $live == true ];
    then
      echo "!!!!!!!!! Installing developer (live) version in ${venv_name}"
      sudo -E $env_path/bin/pip install -e .
    else
      version=$(<VERSION)
      sudo -E $env_path/bin/python setup.py bdist_wheel -d ./wheels
      sudo -E $env_path/bin/pip install "./wheels/islets-${version}-py3-none-any.whl"
  fi

  if [ $silent != true ]; then echo "============== Attempting to set up jupyter kernel for '${venv_name}'"; fi
  sudo $env_path/bin/python -m ipykernel install --name $venv_name --display-name "Physio [${venv_name}]"
fi