echo 'Script for installing islets package.'
echo 'Building wheel and installing it into current environment...'

if ((id -u != 0)); then
  echo "Please run as a root"
  exit
fi

path_to_pip=/opt/tljh/user/envs/physio/bin/pip
version=$(<VERSION)

sudo -E python setup.py bdist_wheel -d ./wheels/
sudo $path_to_pip install "./wheels/islets-$version-py3-none-any.whl"

echo 'Finished script for installing package.'
