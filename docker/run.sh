#!/bin/sh
docker run --rm -it \
        -v $(pwd)/work:/home/user/work \
        -v $(pwd)/..:/openbread \
        -v $(pwd)/..:/geek \
        geek_docker
