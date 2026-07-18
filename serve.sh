#!/bin/bash
# Serve this Jekyll site locally, matching GitHub Pages' toolchain.
#
# Why Docker: Fedora ships Ruby 4.0, but the github-pages gem set only
# supports Ruby 3.x -- its precompiled nokogiri declares
# required_ruby_version < 3.5.dev. On Ruby 4.0 bundler can't install the
# locked gems, silently re-resolves, and backtracks all the way to
# github-pages 8 (2013), which then fails compiling C extensions that
# predate Ruby 2.4. Pinning the interpreter is the durable fix.
#
# Gems persist in the 'picup-bundle' docker volume, so only the first
# run is slow.
#
#   --security-opt label=disable  : Fedora SELinux would otherwise deny
#                                   the bind mount. Non-destructive --
#                                   does not relabel your files.
#   --user $(id -u)               : so _site/ isn't created root-owned.
#
# Site will be at: http://localhost:4000/picup-workshop/
#                  (the /picup-workshop/ suffix comes from `baseurl` in
#                   _config.yml -- plain localhost:4000 will 404)
set -e

cd "$(dirname "$0")"

# One-time: make the gem volume writable by your uid.
docker run --rm --security-opt label=disable \
  -v picup-bundle:/bundle ruby:3.3 \
  chown -R "$(id -u):$(id -g)" /bundle

exec docker run --rm -it \
  --security-opt label=disable \
  --user "$(id -u):$(id -g)" \
  -e HOME=/tmp \
  -e BUNDLE_PATH=/bundle \
  -v "$(pwd)":/srv/jekyll \
  -v picup-bundle:/bundle \
  -p 4000:4000 \
  -w /srv/jekyll \
  ruby:3.3 \
  bash -c "bundle install && bundle exec jekyll serve --host 0.0.0.0 --force_polling"
