## this will create developer (live) version of Physio_Ca
venv_name=$1
echo "============== Attempting to set up ${venv_name} env"
sudo /opt/tljh/user/bin/conda create -n $venv_name python==3.8 pip ipykernel setuptools -y 
echo "============== Installing python-javabridge"
sudo /opt/tljh/user/bin/conda install -n $venv_name python-javabridge -c conda-forge -y
echo "============== Installing python-bioformats"
sudo /opt/tljh/user/envs/$venv_name/bin/pip install python-bioformats
echo "============== Installing islets"
sudo -E /opt/tljh/user/envs/${venv_name}/bin/pip install -e .
sudo /opt/tljh/user/envs/$venv_name/bin/python -m ipykernel install --name $venv_name --display-name "${venv_name}"