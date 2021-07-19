FROM ruby:2.5

RUN apt-get update && \
  apt-get install -y \
    libfox-1.6-dev \
    libxrandr-dev \
    pkg-config

RUN gem install --verbose fxruby

WORKDIR /usr/local/app
COPY inchi_v1.4.rb ./inchi_v1.4.rb
COPY periodic_table.rb ./periodic_table.rb

CMD ["ruby", "inchi_v1.4.rb"]