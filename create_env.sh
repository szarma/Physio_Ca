## this will create developer (live) version of Physio_Ca
venv_name=$1
conda_path=`which conda`
main_path=$(dirname $(dirname $conda_path))
env_path="${main_path}/envs/${venv_name}"
echo "venv_name:   ${venv_name}"
echo "conda path:  ${conda_path}"
echo "main path:   ${main_path}"
echo "env path:    ${env_path}"



ARGS=$(getopt -a --options sl --long "silent,live" -- "$@")
eval set -- "$ARGS"
silent="false"
live="false"
while true; do
  case "$1" in
    -s|--silent)
      silent="true"
      shift;;
    -l|--live)
      live="true"
      shift;;
    --)
      break;;
     *)
      printf "Unknown option %s\n" "$1"
      exit 1;;
  esac
done

#if [ $live == true ]
#  then echo "live"
#  else echo "not live"
#fi

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