docker run --rm -it ^
        -v %CD%\work:/home/user/work ^
        -v %CD%/..:/openbread ^
        -v %CD%/..:/geek ^
        geek_docker %*
