### intended to be run as the following:
### sudo -E -s ./install_in_env.sh
if ((id -u != 0)); then
  echo "Please run as a root with flags '-s -E'"
  exit
fi

env=$1

echo 'Script for installing islets package.'
echo "Building wheel and installing it into ${env}..."

version=$(<VERSION)

path_to_pip="${env}/bin/pip"
path_to_python="${env}/bin/python"
sudo -E $path_to_python setup.py bdist_wheel -d ./wheels/
sudo -E $path_to_pip install "./wheels/islets-${version}-py3-none-any.whl"
echo "Successfully installed in ${env}."