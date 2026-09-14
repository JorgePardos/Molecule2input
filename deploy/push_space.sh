#!/bin/sh
# Publish the committed state of the browser interface to a Hugging Face Space.
#
#   deploy/push_space.sh https://huggingface.co/spaces/<user>/<space>
#
# The Space receives a snapshot, not this repository's history: only what the
# image needs (Dockerfile, pinned requirements, README with the Space header,
# licence, the m2i package), taken from the last commit. Two reasons. Spaces
# refuse pushes whose history holds binary files outside Git LFS -- and the
# example picture is one. And the tests, fixtures and examples have no business
# on a public server.
#
# git asks for your Hugging Face username and, as the password, an access token
# with write permission (huggingface.co/settings/tokens).
set -eu

remote="${1:?usage: deploy/push_space.sh https://huggingface.co/spaces/<user>/<space>}"
root=$(git rev-parse --show-toplevel)
revision=$(git -C "$root" rev-parse --short HEAD)

if [ -n "$(git -C "$root" status --porcelain -- Dockerfile requirements-web.txt README.md LICENSE m2i)" ]; then
    echo "error: commit your changes first; the Space gets the last commit, not the working tree." >&2
    exit 1
fi

snapshot=$(mktemp -d)
trap 'rm -rf "$snapshot"' EXIT

git -C "$root" archive HEAD Dockerfile .dockerignore requirements-web.txt README.md LICENSE m2i \
    | tar -x -C "$snapshot"

cd "$snapshot"
git init -q
git checkout -q -b main
git add -A
git commit -q -m "m2i $revision"
git push --force "$remote" main
echo "Published m2i $revision to $remote"
