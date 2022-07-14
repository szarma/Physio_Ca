### intended to be run as the following:
### sudo -E -s ./install_in_env.sh
if ((id -u != 0)); then
  echo "Please run as a root with flags '-s -E'"
  exit
fi

path_to_env=$1

env_name=${path_to_env##*/}

echo "Environment name: ${env_name}"

n=${#env_name}

echo "Environment name length: ${n}"

if [ ${n} == 0 ]; then
    echo "Please enter the path without the last slash (/) and try again."
    exit
fi

echo 'Script for installing islets package.'
echo "Building wheel and installing it into ${path_to_env}..."
version=$(<VERSION)
path_to_pip="${path_to_env}/bin/pip"
path_to_python="${path_to_env}/bin/python"
# sudo -E $path_to_python setup.py bdist_wheel -d ./wheels/
# sudo -E $path_to_pip install "./wheels/islets-${version}-py3-none-any.whl"
echo "Successfully installed in ${path_to_env}."


sudo -E $path_to_python -m ipykernel install --name $env_name --display-name "${env_name}"