# openbread Docker

This Docker offers a suitable environment to run openbread.

## Requirements

Make sure [docker](https://www.docker.com/) is installed and running before entering the commands below.


## Building and Runnig the Docker

First, build the container with `build.bat` or `. build`.
Then start the container with `run.bat` or `. run`.
```bash
. ./build.sh
. ./run.sh
```


## Additional information

If you are running your Docker container in a Unix-based environment, you might get permission errors on the `.sh` scripts.
This is because permissions are inherited from the host environment. You can fix this by running in the `docker` folder:
```bash
chmod +x utils/*.sh

```
