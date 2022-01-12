FROM jupyter/minimal-notebook:latest

# Switch to root user in order to install Ruby.
USER root
RUN apt-get update
RUN apt-get --yes install ruby-full
RUN gem install rubocop
# Switch back to original notebook user.
USER $NB_UID