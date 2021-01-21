### intended to be run as the following:
### sudo -E -s ./install_package.sh


echo 'Script for installing islets package.'
echo 'Building wheel and installing it into current environment...'

if ((id -u != 0)); then
  echo "Please run as a root with flags '-s -E'"
  exit
fi

path_to_pip=/opt/tljh/user/envs/physio/bin/pip
path_to_python=/opt/tljh/user/envs/physio/bin/python

version=$(<VERSION)

sudo -E $path_to_python setup.py bdist_wheel -d ./wheels/
sudo -E $path_to_pip install "./wheels/islets-${version}-py3-none-any.whl"

echo 'Finished script for installing package.'
