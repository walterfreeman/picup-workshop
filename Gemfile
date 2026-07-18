source 'https://rubygems.org'

# Plain Jekyll, not the `github-pages` gem.
#
# GitHub Pages ignores this Gemfile entirely -- it builds with its own
# server-side toolchain. The github-pages gem only ever existed to
# approximate that locally, and it drags in ~60 pinned gems including a
# precompiled nokogiri that refuses to install on Ruby 3.5+. This site
# uses no plugins and no remote theme (_layouts/_includes/_sass are all
# local), so it needs none of that.
gem 'jekyll', '~> 4.4'
gem 'webrick'   # no longer in Ruby's stdlib; `jekyll serve` needs it
