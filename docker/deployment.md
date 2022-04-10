## Docker Deployment (still in development and to be changed in near future)
JupyterHub can be deployed as a docker image giving users individual sandboxes to work with and easily upgradable containers to developers.

### Prerequisites
To install the CTN-JupyterHub make sure to have:
1. administrative rights
2. internet connection
3. Installed dependencies:
   - docker engine version 19.03.0 or higher
   - docker-compose

### Preparation
To install CTN make sure to clone the repository with its content to one of your folders on the server.
The essential parts for the docker installation beside the notebooks and the islets module are:
- Compose file letting docker know what services to have and how to connect: ./docker-compose.yml
- Environment file defining specific necessary settings: ./.env
- Folder containing the dockerfiles and helpers for image creation: ./docker/

Before installation you may need to edit the environment file!
Inside the file there are different parameters to be set:
- SERVER_DOCKER_PATH - defines where in the data to look for the dockerfiles of the server definition. Usually there are no changes needed!
- HUB_DOCKER_PATH - defines where in the data to look for the dockerfiles of the hub definition. Usually there are no changes needed!
- HOST - IP-address to be used for referring with traefik. On test systems leave it as is is. **For productive systems change to correct URL!**
- COMPOSE_PROJECT_NAME

### Installation
Installation of the images is rather easy.
Once navigated into the cloned directory we need to build the images using the following:
```
sudo docker-compose build
```
This will create 3 different images:
- traefik (Webserver to dispatch requests from the global network and between containers)
- jupyterhub (The hub itself for user management)
- jupyterhub_server (The image for a single person calculating server)

To start up the system we just need to type:
```
sudo docker-compose up
# optional: add -d as parameter to run it as a daemon process
```
After this the server should be ready to work.
