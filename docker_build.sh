#!/usr/bin/env bash
set -e


DIR=$(cd $(dirname $0) && pwd)

# Registries and tags
STORAGE_URL="gitlab.ispras.ru:4567"
IMAGE_URL="gitlab.ispras.ru:4567/bioinformatics/images/hobotnica"
TOOL_NAME="hobotnica"

# hobotnica version
HOBOTNICA_PUBLIC_VERSION="7.0.1"

# Necessary tools and help text
TOOLS=(docker gcloud)
HELP="$(basename "$0") [-h|--help] [-t|--tools] -- script to build the hobotnica image and push to ISPRAS image store

where:
    -h|--help Show help text
    -v|--version of hobotnica to use (default: $HOBOTNICA_PUBLIC_VERSION)
    -t|--tools Show tools needed to run script
    "

function main(){
    for t in "${TOOLS[@]}"; do which $t >/dev/null || ok=no; done
        if [[ $ok == no ]]; then
            echo "Missing one of the following tools: "
            for t in "${TOOLS[@]}"; do echo "$t"; done
            exit 1
        fi

    while [[ $# -gt 0 ]]
    do
    key="$1"
    case $key in
        -v|--version)
        HOBOTNICA_PUBLIC_VERSION="$2"
        shift
        shift
        ;;
        -h|--help)
        echo "$HELP"
        exit 0
        ;;
        -t|--tools)
        for t in "${TOOLS[@]}"; do echo $t; done
        exit 0
        ;;
        *)
        shift
        ;;
    esac
    done

    IMAGE_TAG="$HOBOTNICA_PUBLIC_VERSION"

    echo "building and pushing Image - $IMAGE_URL:$IMAGE_TAG"
    docker build -t "$TOOL_NAME:$IMAGE_TAG" \
        --build-arg HOBOTNICA_PUBLIC_VERSION="$HOBOTNICA_PUBLIC_VERSION" "$DIR"
    docker login "$STORAGE_URL"
    docker tag "$TOOL_NAME:$IMAGE_TAG" "$IMAGE_URL:$IMAGE_TAG"
    docker push "$IMAGE_URL:$IMAGE_TAG"

    echo "$IMAGE_URL:$IMAGE_TAG" >> "$DIR/docker_versions.tsv"
    echo "done"
}

main "$@"
