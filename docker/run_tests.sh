#!/bin/sh
chmod -R u+X .
docker run --rm \
        -v $(pwd)/..:/geek \
        geek_docker_ci 	\
        bash -c "cd /geek && py.test"
