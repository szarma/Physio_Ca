echo 'Script for installing islets package.'
echo 'Building wheel and installing it into current environment...'

version=$(<VERSION)
python setup.py bdist_wheel -d ./wheels/
pip install "./wheels/islets-$version-py3-none-any.whl"

echo 'Finished script for installing package.'