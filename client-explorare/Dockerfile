FROM node:9.2.1
MAINTAINER Marcin Piechota <piechota.marcin@gmail.com>
# use changes to package.json to force Docker not to use the cache
# when we change our application's nodejs dependencies:
ADD package.json /tmp/package.json
RUN cd /tmp && yarn install
RUN mkdir -p /usr/src/app && cp -a /tmp/node_modules /usr/src/app
# From here we load our application's code in, therefore the previous docker
# "layer" thats been cached will be used if possible
WORKDIR /usr/src/app
ADD . /usr/src/app
RUN npm run build
RUN rm -rf ./build
RUN rm -rf ./test
RUN rm -rf ./src
ENV PORT=8081
EXPOSE 8081
CMD [ "npm", "run", "serve" ]
