FROM    httpd:2.4

# Set port
RUN sed -i 's/80/8085/' /usr/local/apache2/conf/httpd.conf
#RUN sed -i 's/LoadModule mpm_event_module/#LoadModule mpm_event_module/' /usr/local/apache2/conf/httpd.conf
#RUN sed -i 's/#LoadModule mpm_prefork_module/LoadModule mpm_prefork_module/' /usr/local/apache2/conf/httpd.conf
#RUN echo "LoadModule mpm_event_module modules/mod_mpm_event.so" >> /usr/local/apache2/conf/httpd.conf


# Copy html
ADD . / /usr/local/apache2/htdocs/
#ADD httpd.conf /usr/local/apache2/conf/httpd.conf
